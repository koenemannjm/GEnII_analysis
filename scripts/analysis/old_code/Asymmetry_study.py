import ROOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
import pandas as pd
from scipy.interpolate import NearestNDInterpolator
from matplotlib.backends.backend_pdf import PdfPages

data_dirpath = "/w/halla-scshelf2102/sbs/koeneman/GEnII_analysis/outfiles/HUNTERrootfiles/"
model_dirpath = "/w/halla-scshelf2102/sbs/koeneman/GEnII_analysis/outfiles/csvfiles/"
model_filename = "he3model.csv"

df = pd.read_csv(model_dirpath + model_filename)

F1_QE_interpolated = NearestNDInterpolator(df[['Q2','X']], df['F1_QE'])
F2_QE_interpolated = NearestNDInterpolator(df[['Q2','X']], df['F2_QE'])
G1_QE_interpolated = NearestNDInterpolator(df[['Q2','X']], df['G1_QE'])
G2_QE_interpolated = NearestNDInterpolator(df[['Q2','X']], df['G2_QE'])

# Inelastic
F1_Inelastic_interpolated = NearestNDInterpolator(df[['Q2','X']], df['F1_Inelastic'])
F2_Inelastic_interpolated = NearestNDInterpolator(df[['Q2','X']], df['F2_Inelastic'])
G1_Inelastic_interpolated = NearestNDInterpolator(df[['Q2','X']], df['G1_Inelastic'])
G2_Inelastic_interpolated = NearestNDInterpolator(df[['Q2','X']], df['G2_Inelastic'])

# Inelastic
# F1_Inelastic_interpolated = NearestNDInterpolator(df[['Q2','X']], df['F1_IpQE'])
# F2_Inelastic_interpolated = NearestNDInterpolator(df[['Q2','X']], df['F2_IpQE'])
# G1_Inelastic_interpolated = NearestNDInterpolator(df[['Q2','X']], df['G1_IpQE'])
# G2_Inelastic_interpolated = NearestNDInterpolator(df[['Q2','X']], df['G2_IpQE'])

# Consants
Mp = 0.938  # Proton mass in GeV/c^2
Mn = 0.939  # Proton mass in GeV/c^2
alpha = 1.0 / 137.0 # E&M fine structure constant
B_He3 = 7.718 * 10 ** (-3) # Helium 3 binding energy
M3He = Mn + 2 * Mp - B_He3 # Helium 3 mass
electronEnergyInterval=11-.5 # G4 E'Max-E'Min in generator
hbarc=3.89e-28 # factor of hbar*c
cm2nb=1e33 # cm^2 to nb converter
MN = (Mp + Mn) / 2

def mott_cros_sec(E, th):
    """
    This Function computes the Mott Cross-section

    d^2sigma/dOmega = alpha^2 * cos^2(theta/2) / ( 4 * E^2 * sin^4(theta/2) )

    Input:
    E -> Initial electron energy (beam energy)
    theta -> the electron scattering angle in scattering plane

    output:
    Mott Cross-section value

    """
    sin_th_sqr = np.sin(th/2)**2
    numer = alpha * np.cos(th/2) #numerator of mott cross-section inside (...)^2
    denom = 2 * E * sin_th_sqr #denominator of mott cross-section (...)^2
    mott = (numer/denom)**2 #value of mott cross-section
    return mott

def unpol_cros_sec(E, th, Ep, arg):
    """
    This Function computes the Unpolarized Cross-section for Nucleon/Nucleus scattering

    d^2/dOmega dE_prime = [MOTT] * ( F2 / nu + (2/M)*tan^2(theta/2)F1 )

    Input:
    E : Initial electron energy
    theta : Angle between k and k' in scattering plane
    Ep : Final electron energy


    """
    mott_xsec = mott_cros_sec(E, th)
    massn = MN

    q2 = 2 * E * Ep * (1 - np.cos(th))
    nu = E - Ep
    x = q2 / (2 * massn * nu)

    # Computing F1 and F2
    if arg == 'QE':
      F1 = F1_QE_interpolated(q2, x)
      F2 = F2_QE_interpolated(q2, x)
    elif arg == 'IN':
      F1 = F1_Inelastic_interpolated(q2, x)
      F2 = F2_Inelastic_interpolated(q2, x)
    else:
      raise ValueError("Invalid argument. Expected 'QE' or 'IN'.")

    tan_th_sqr = np.tan(th/2) **2

    unpol = mott_xsec * ( F2/nu + 2 * tan_th_sqr * F1 / massn ) #value of unpolarized cross-section

    return unpol

def pol_cros_sec(E, th, Ep, phi_star, theta_star, arg):

    massn = MN

    nu = E - Ep #Energy transfered from electron to Nucleon in Nucleus
    q2 = 2 * E * Ep * (1 - np.cos(th))
    x = q2 / (2 * massn * nu)
    EovM = E / massn
    tau = q2 / (4 * massn**2)

    Cosbeta = (np.cos(phi_star)*np.sin(theta_star)*np.sin(th)*(EovM - 2 * tau) + np.cos(theta_star)*(EovM - (EovM - 2 * tau)*np.cos(th)))/(2*np.sqrt(tau*(1 + tau)))
    CosTheta = (np.cos(phi_star)*np.sin(theta_star)*np.sin(th)*EovM + np.cos(theta_star)*(EovM*np.cos(th) - (EovM - 2 * tau)))/(2*np.sqrt(tau*(1 + tau)))

    Gamma = -(4 * (alpha **2) * Ep) / (q2 * E * massn * nu)
    A = E * Cosbeta + Ep * CosTheta
    B = (2 * E * Ep / nu) * (CosTheta - Cosbeta)

    # Computing g1 and g2
    if arg == 'QE':
      g1 = G1_QE_interpolated(q2, x)
      g2 = G2_QE_interpolated(q2, x)
    elif arg == 'IN':
      g1 = G1_Inelastic_interpolated(q2, x)
      g2 = G2_Inelastic_interpolated(q2, x)
    else:
      raise ValueError("Invalid argument. Expected 'QE' or 'IN'.")

    pol = Gamma * (A * g1 + B * g2)

    return pol

def sgn(x):
    return (x > 0) - (x < 0)

if len(sys.argv)<1:
    print("Missing kinematic arguments: kin2, kin3, kin4a, kin4b")
    sys.exit()

kine_list = ["kin2","kin3","kin4a","kin4b"]
# kine_list = [sys.argv[1]]

kine_counter = 1

for kine in kine_list:
    print(f"Processing {kine} ... | {kine_counter}/{len(kine_list)}")
    kine_counter += 1
    if kine == "kin2":
        E = 4.291 #GeV
        target = "He3"
        datapass = "pass2"
        expname = "GEN2"
        Pkin = -1
        coin_time_mean = 128.542
        coin_time_sigma = 1.3707
        T_magnet_angle_deg = 59.34
        T_magnet_pitch_deg = 1.840
        T_magnet_angle_rad = np.deg2rad(T_magnet_angle_deg)
        T_magnet_pitch_rad = np.deg2rad(T_magnet_pitch_deg)
        eprime_min = 0.5
        eprime_max = 3.5
        etheta_min_deg = 25.0
        etheta_max_deg = 36.0
        etheta_min_rad = np.deg2rad(etheta_min_deg)
        etheta_max_rad = np.deg2rad(etheta_max_deg)
        ephi_min_rad = -0.4
        ephi_max_rad = 0.4
        dy_min_anti = -3.0
        dy_max_anti = 2.5
    elif kine == "kin3":
        E = 6.373 #GeV
        target = "He3"
        datapass = "pass2"
        expname = "GEN3"
        Pkin = 1
        coin_time_mean = 119.849
        coin_time_sigma = 1.71674
        T_magnet_angle_deg = 70.94
        T_magnet_pitch_deg = 0.510
        T_magnet_angle_rad = np.deg2rad(T_magnet_angle_deg)
        T_magnet_pitch_rad = np.deg2rad(T_magnet_pitch_deg)
        eprime_min = 0.5
        eprime_max = 3.5
        etheta_min_deg = 31.0
        etheta_max_deg = 43.0
        etheta_min_rad = np.deg2rad(etheta_min_deg)
        etheta_max_rad = np.deg2rad(etheta_max_deg)
        ephi_min_rad = -0.4
        ephi_max_rad = 0.4
        dy_min_anti = -1.3
        dy_max_anti = 1.0
    elif kine == "kin4a":
        E = 8.448 #GeV
        target = "He3"
        datapass = "pass2"
        expname = "GEN4"
        Pkin = 1
        coin_time_mean = 121.231
        coin_time_sigma = 1.74333
        T_magnet_angle_deg = 74.07
        T_magnet_pitch_deg = 0.550
        T_magnet_angle_rad = np.deg2rad(T_magnet_angle_deg)
        T_magnet_pitch_rad = np.deg2rad(T_magnet_pitch_deg)
        eprime_min = 0.5
        eprime_max = 4.0
        etheta_min_deg = 30
        etheta_max_deg = 41.0
        etheta_min_rad = np.deg2rad(etheta_min_deg)
        etheta_max_rad = np.deg2rad(etheta_max_deg)
        ephi_min_rad = -0.4
        ephi_max_rad = 0.4
        dy_min_anti = -1.1
        dy_max_anti = 1.0
    elif kine == "kin4b":
        E = 8.448 #GeV
        target = "He3"
        datapass = "pass2"
        expname = "GEN4b"
        Pkin = 1
        coin_time_mean = 185.098
        coin_time_sigma = 1.76543
        T_magnet_angle_deg = 74.07
        T_magnet_pitch_deg = 0.550
        T_magnet_angle_rad = np.deg2rad(T_magnet_angle_deg)
        T_magnet_pitch_rad = np.deg2rad(T_magnet_pitch_deg)
        T_magnet_angle_rad = np.deg2rad(T_magnet_angle_deg)
        T_magnet_pitch_rad = np.deg2rad(T_magnet_pitch_deg)
        eprime_min = 0.5
        eprime_max = 4.0
        etheta_min_deg = 30.0
        etheta_max_deg = 41.0
        etheta_min_rad = np.deg2rad(etheta_min_deg)
        etheta_max_rad = np.deg2rad(etheta_max_deg)
        ephi_min_rad = -0.4
        ephi_max_rad = 0.4
        dy_min_anti = -1.1
        dy_max_anti = 1.0

    numPoints = 10000
    ebeam = np.ones(numPoints) * E
    eprime = np.linspace(eprime_min,eprime_max,numPoints)
    etheta = np.linspace(etheta_min_rad,etheta_max_rad,numPoints)
    ephi = np.linspace(ephi_min_rad,ephi_max_rad,numPoints)
    
    theta_star = []
    phi_star = []
    for i in range(numPoints):
        ethetai = etheta[i]
        eprimei = eprime[i]
        ephii = ephi[i]
        
        ebeami = ebeam[i]
        ke = ebeami * np.array([[0, 0, 1]])
        keprime = eprimei * np.array([[np.cos(ephii)*np.sin(ethetai),
                                       np.sin(ephii)*np.sin(ethetai),
                                       np.cos(ethetai)]])
        target_spin = np.array([[np.sin(T_magnet_angle_rad)*np.cos(T_magnet_pitch_rad),
                                 np.sin(T_magnet_angle_rad)*np.sin(T_magnet_pitch_rad),
                                 np.cos(T_magnet_angle_rad)]])
        q = ke - keprime
        
        z_star = q / np.linalg.norm(q)
        y_star = np.cross(ke,keprime) / (np.linalg.norm(ke) * np.linalg.norm(keprime))
        x_star = np.cross(y_star,z_star)
        
        target_spin_x_star = np.dot(target_spin,x_star.T).item()
        target_spin_y_star = np.dot(target_spin,y_star.T).item()
        target_spin_z_star = np.dot(target_spin,z_star.T).item()
        
        target_spin_star = np.array([[target_spin_x_star,target_spin_y_star,target_spin_z_star]])
        
        theta_stari = np.arccos(target_spin_star.item(2)/np.linalg.norm(target_spin_star))
        phi_stari = sgn(target_spin_star.item(1)) * np.arccos(target_spin_star.item(0)/np.sqrt(target_spin_star.item(0)**2 + target_spin_star.item(1)**2))
            
        theta_star.append(theta_stari)
        phi_star.append(phi_stari)
            
    nu_model = ebeam - eprime
    Q2_model = 2 * eprime * ebeam * (1 - np.cos(etheta))
    X_model = Q2_model / (2 * MN * nu_model)
    W2_model = MN**2 - Q2_model + 2 * MN * nu_model
    
    theta_star_model = np.array(theta_star)
    phi_star_model = np.array(phi_star)
    
    process = 'IN'
    xz_unpol_cross_sec = unpol_cros_sec(ebeam, etheta, eprime, process)
    xz_pol_cross_sec = pol_cros_sec(ebeam, etheta, eprime, phi_star, theta_star, process)
    xz_mott_cross_sec = mott_cros_sec(ebeam, etheta)
    
    Asymmetry_model = []
    for i in range(len(xz_pol_cross_sec)):
        if xz_unpol_cross_sec[i] == 0:
            Asymmetry_model.append(0)
            continue
        Asymmetry_model.append(xz_pol_cross_sec[i]/(2*xz_unpol_cross_sec[i]))

    Asymmetry_model = np.array(Asymmetry_model)
        
    data_filename = f"Final_data_{expname}_sbs100p_nucleon_np_model2.root"
    data = data_dirpath + data_filename
    tree = "Tout"
    
    print(f"Reading in data for {kine}")
    raw_data = ROOT.RDataFrame(tree, data)
    branch_variables = ["dx","dy","W2","ePS","eSH","eHCAL","coin_time","helicity","IHWP","runnum","trP","trigbits","Q2","nu"]
    raw_data_numpy = raw_data.AsNumpy(branch_variables)
    
    print("Extracting variables")
    
    print("Loading dx...")
    dx = raw_data_numpy["dx"]
    print("Loading dy...")
    dy = raw_data_numpy["dy"]
    print("Loading W2...")
    W2 = raw_data_numpy["W2"]
    print("Loading Q2...")
    Q2 = raw_data_numpy["Q2"]
    print("Loading nu...")
    nu = raw_data_numpy["nu"]
    print("Loading ePS...")
    ePS = raw_data_numpy["ePS"]
    print("Loading eSH...")
    eSH = raw_data_numpy["eSH"]
    print("Loading eHCAL...")
    eHCAL = raw_data_numpy["eHCAL"]
    print("Loading coin_time...")
    coin_time = raw_data_numpy["coin_time"]
    print("Loading helicity...")
    helicity = raw_data_numpy["helicity"]
    print("Loading IHWP...")
    IHWP = raw_data_numpy["IHWP"]
    print("Loading runnum...")
    runnum = raw_data_numpy["runnum"]
    print("loading trP...")
    trP = raw_data_numpy["trP"]
    print("loading trigbits...")
    trigbits = raw_data_numpy["trigbits"]
    eop = (ePS + eSH)/trP

    del raw_data, raw_data_numpy
    
    print("Cutting IHWP != 1 or -1 and helicity != 1 or -1")
    
    label1 = ["QE","Anti"]
    value = [(np.abs(dy)<0.5)&(W2>0.0)&(W2<2.0)&(np.abs(dx)<0.5),(np.abs(dy - (dy_max_anti+dy_min_anti)/2)>np.abs(dy_max_anti-dy_min_anti)/2)&(W2>-4.0)&(W2<12)]
    with PdfPages(f"fig_asymmetry_study_{kine}.pdf") as pdf:
        for j in range(2):
            
            cut = (np.abs(IHWP) == 1) & (np.abs(helicity) == 1) & value[j] & (ePS>0.2) & (eHCAL>0.025) & (np.abs(coin_time - coin_time_mean)<coin_time_sigma) & (eop > 0.7) & (trigbits == 4)
            helicity_cut = helicity[cut]
            IHWP_cut = IHWP[cut]
            raw_helicity_cut = IHWP_cut*Pkin*helicity_cut
            
            dx_cut = dx[cut]
            dy_cut = dy[cut]

            bins_edges_dx = np.histogram_bin_edges(dx_cut, bins='auto')
            bins_edges_dy = np.histogram_bin_edges(dy_cut, bins='auto')

            Counts_dydx, edges_dy, edges_dx = np.histogram2d(dy_cut,dx_cut,bins=(bins_edges_dy, bins_edges_dx))

            hel = (raw_helicity_cut == 1)
            dx_cut_hel = dx_cut[hel]
            dy_cut_hel = dy_cut[hel]

            Counts_dydx_up, edges_dydx_junk, junk =  np.histogram2d(dy_cut_hel,dx_cut_hel,bins=(bins_edges_dy, bins_edges_dx))
            errCounts_dydx_up = np.sqrt(Counts_dydx_up)


            hel = (raw_helicity_cut == -1)
            dx_cut_hel = dx_cut[hel]
            dy_cut_hel = dy_cut[hel]

            Counts_dydx_dwn, edges_dydx_junk, junk =  np.histogram2d(dy_cut_hel,dx_cut_hel,bins=(bins_edges_dy, bins_edges_dx))
            errCounts_dydx_dwn = np.sqrt(Counts_dydx_dwn)

            del dx_cut, dy_cut, edges_dydx_junk, dx_cut_hel, dy_cut_hel, hel, junk

            Counts_dydx_up[Counts_dydx_up == 0] = np.nan
            Counts_dydx_dwn[Counts_dydx_dwn == 0] = np.nan

            Asymmetry_dydx = (Counts_dydx_up - Counts_dydx_dwn) / (Counts_dydx_up + Counts_dydx_dwn)
            # Asymmetry_dydx = np.abs(Asymmetry_dydx)
            # Asymmetry_dydx[Asymmetry_dydx < 0] *= -1
            errAsymmetry_dydx = 2 * np.sqrt(np.power(Counts_dydx_up*errCounts_dydx_dwn,2) + np.power(Counts_dydx_dwn*errCounts_dydx_up,2)) / np.power(Counts_dydx_up + Counts_dydx_dwn, 2)
            
            W2_cut = W2[cut]
            Q2_cut = Q2[cut]
            nu_cut = nu[cut]
            xbjk_cut = Q2_cut / (2 * MN * nu_cut)
            del nu_cut

            bins_edges_W2 = np.histogram_bin_edges(W2_cut, bins='auto')
            bins_edges_Q2 = np.histogram_bin_edges(Q2_cut, bins='auto')
            # bins_edges_xbjk = np.histogram_bin_edges(xbjk_cut, bins-'auto')

            Counts_W2Q2, edges_W2, edges_Q2 = np.histogram2d(W2_cut,Q2_cut,bins=(bins_edges_W2, bins_edges_Q2))

            hel = (raw_helicity_cut == 1)
            W2_cut_hel = W2_cut[hel]
            Q2_cut_hel = Q2_cut[hel]

            Counts_W2Q2_up, edges_W2Q2_junk, junk =  np.histogram2d(W2_cut_hel,Q2_cut_hel,bins=(bins_edges_W2, bins_edges_Q2))
            errCounts_W2Q2_up = np.sqrt(Counts_W2Q2_up)

            hel = (raw_helicity_cut == -1)
            W2_cut_hel = W2_cut[hel]
            Q2_cut_hel = Q2_cut[hel]

            Counts_W2Q2_dwn, edges_W2Q2_junk, junk =  np.histogram2d(W2_cut_hel,Q2_cut_hel,bins=(bins_edges_W2, bins_edges_Q2))
            errCounts_W2Q2_dwn = np.sqrt(Counts_W2Q2_up)

            del W2_cut, Q2_cut, edges_W2Q2_junk, W2_cut_hel, Q2_cut_hel, hel, junk

            Counts_W2Q2_up[Counts_W2Q2_up == 0] = np.nan
            Counts_W2Q2_dwn[Counts_W2Q2_dwn == 0] = np.nan

            Asymmetry_W2Q2 = (Counts_W2Q2_up - Counts_W2Q2_dwn) / (Counts_W2Q2_up + Counts_W2Q2_dwn)
            # Asymmetry_W2Q2 = np.abs(Asymmetry_W2Q2)
            # Asymmetry_W2Q2[Asymmetry_W2Q2 < 0] *= -1
            errAsymmetry_W2Q2 = 2 * np.sqrt(np.power(Counts_W2Q2_up*errCounts_W2Q2_dwn,2) + np.power(Counts_W2Q2_dwn*errCounts_W2Q2_up,2)) / np.power(Counts_W2Q2_up + Counts_W2Q2_dwn, 2)

            # -------- GENERATING PLOTS --------
            # Asymmetry_min = np.min(Asymmetry_dydx[(np.abs(Asymmetry_dydx)!=1)])
            # Asymmetry_max = np.max(Asymmetry_dydx[(np.abs(Asymmetry_dydx)!=1)])
            Asymmetry_min = -1.0
            Asymmetry_max = 1.0
            # dy - dx histograms and plots
            print(f"Making dx-dy plots {label1[j]}")
            fig, axes = plt.subplots(2, 3, figsize=(22,12), sharex=False, sharey=False)
            vmin = min(Counts_dydx.min(), Counts_dydx_up.min(), Counts_dydx_dwn.min())
            vmax = max(Counts_dydx.max(), Counts_dydx_up.max(), Counts_dydx_dwn.max())
            number_levels = 100
            levels = np.linspace(vmin,vmax, number_levels+1)
            norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)
            
            for ax, Counts, title in zip(axes[0], [Counts_dydx, Counts_dydx_up, Counts_dydx_dwn], ["All", "rawhel = +1", "rawhel = -1"]):
                mesh = ax.pcolormesh(edges_dy, edges_dx, Counts.T, cmap='rainbow', norm=norm)
                ax.set_title(f"{title} Events kine = {expname} ({label1[j]})")
                ax.set_xlabel("dy")
                ax.set_ylabel("dx")
                ax.minorticks_on()
                ax.grid(True)


            cbar = fig.colorbar(mesh, ax=axes[0].ravel().tolist(), shrink=0.8)
            cbar.set_label("Counts")

            centers_dx = 0.5*(edges_dx[:-1] + edges_dx[1:])
            centers_dy = 0.5*(edges_dy[:-1] + edges_dy[1:])

            weights_dydx = 1/np.power(errAsymmetry_dydx,2)
            
            axes[1,0].errorbar(centers_dx, np.nansum((Asymmetry_dydx*weights_dydx).T, axis=1)/np.nansum(weights_dydx.T,axis=1), yerr= np.sqrt(1/np.nansum(weights_dydx.T,axis=1)), fmt='o', markersize=3, capsize=3)
            axes[1,0].set_xlabel("dx")
            axes[1,0].set_ylabel("Asymmetry")
            axes[1,0].set_title(f"Average Asymm. for dx kine = {expname} ({label1[j]})")
            axes[1,0].minorticks_on()
            axes[1,0].set_ylim(-0.12,0.04)
            axes[1,0].grid(True)

            axes[1,1].errorbar(centers_dy, np.nansum((Asymmetry_dydx*weights_dydx).T, axis=0)/np.nansum(weights_dydx.T,axis=0), yerr= np.sqrt(1/np.nansum(weights_dydx.T,axis=0)), fmt='o', markersize=3, capsize=3)
            axes[1,1].set_xlabel("dy")
            axes[1,1].set_ylabel("Asymmetry")
            axes[1,1].set_title(f"Average Asymm. for dy kine = {expname} ({label1[j]})")
            axes[1,1].set_ylim(-0.12,0.04)
            axes[1,1].minorticks_on()
            axes[1,1].grid(True)
            
            zmin, zmax = Asymmetry_min, Asymmetry_max
            number_levels = 100
            levels = np.linspace(zmin,zmax, number_levels+1)
            norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)
            mesh = axes[1,2].pcolormesh(edges_dy, edges_dx, Asymmetry_dydx.T, cmap='rainbow', norm=norm)

            cbar = fig.colorbar(mesh, ax=axes[1,2], shrink=0.8)
            cbar.set_label("Asymmetry")
            axes[1,2].set_xlabel("dy")
            axes[1,2].set_ylabel("dx")
            axes[1,2].set_title(f"Asymmetry heatmap kine = {expname} ({label1[j]})")
            axes[1,2].minorticks_on()
            axes[1,2].grid(True)
            
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            del fig, ax, axes

            # W2 - Q2 histograms and plots
            print(f"Making W2-Q2 plots {label1[j]}")
            Asymmetry_min = np.min(Asymmetry_W2Q2[(np.abs(Asymmetry_W2Q2)!=1 )])
            Asymmetry_max = np.max(Asymmetry_W2Q2[(np.abs(Asymmetry_W2Q2)!=1 )])
            Asymmetry_min = -1.0
            Asymmetry_max = 1.0
            fig, axes = plt.subplots(2, 3, figsize=(22, 12), sharex=False, sharey=False)

            # ---- FIRST ROW: Counts Heatmaps ----
            vmin = min(Counts_W2Q2.min(), Counts_W2Q2_up.min(), Counts_W2Q2_dwn.min())
            vmax = max(Counts_W2Q2.max(), Counts_W2Q2_up.max(), Counts_W2Q2_dwn.max())
            number_levels = 100
            levels = np.linspace(vmin,vmax, number_levels+1)
            norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)
            
            for ax, Counts, title in zip(axes[0], 
                                         [Counts_W2Q2, Counts_W2Q2_up, Counts_W2Q2_dwn],
                                         ["All", "rawhel = +1", "rawhel = -1"]):
                mesh = ax.pcolormesh(edges_W2, edges_Q2, Counts.T, cmap='rainbow',norm=norm)
                ax.set_title(f"{title} Events kine = {expname} ({label1[j]})")
                ax.set_xlabel("W2")
                ax.set_ylabel("Q2")
                ax.minorticks_on()
                ax.grid(True)
                
            # Share colorbar for first row
            cbar = fig.colorbar(mesh, ax=axes[0].ravel().tolist(), shrink=0.8)
            cbar.set_label("Counts")
            
            # ---- SECOND ROW: Asymmetry Plots ----
            centers_W2 = 0.5 * (edges_W2[:-1] + edges_W2[1:])
            centers_Q2 = 0.5 * (edges_Q2[:-1] + edges_Q2[1:])
            weights_W2Q2 = 1 / np.power(errAsymmetry_W2Q2, 2)
            
            # Q2-projected asymmetry plot
            axes[1, 0].errorbar(
                centers_Q2,
                np.nansum((Asymmetry_W2Q2 * weights_W2Q2).T, axis=1) / np.nansum(weights_W2Q2.T, axis=1),
                yerr=np.sqrt(1 / np.nansum(weights_W2Q2.T, axis=1)),
                fmt='o', markersize=3, capsize=3
            )
            axes[1, 0].set_xlabel("Q2")
            axes[1, 0].set_ylabel("Asymmetry")
            axes[1, 0].set_title(f"Average Asymm. for Q2 kine = {expname} ({label1[j]})")
            axes[1, 0].set_ylim(-0.12,0.04)
            axes[1, 0].minorticks_on()
            axes[1, 0].grid(True)
            
            # W2-projected asymmetry plot
            axes[1, 1].errorbar(
                centers_W2,
                np.nansum((Asymmetry_W2Q2 * weights_W2Q2).T, axis=0) / np.nansum(weights_W2Q2.T, axis=0),
                yerr=np.sqrt(1 / np.nansum(weights_W2Q2.T, axis=0)),
                fmt='o', markersize=3, capsize=3
            )
            axes[1, 1].set_xlabel("W2")
            axes[1, 1].set_ylabel("Asymmetry")
            axes[1, 1].set_title(f"Average Asymm. for W2 kine = {expname} ({label1[j]})")
            axes[1, 1].minorticks_on()
            axes[1, 1].set_ylim(-0.12,0.04)
            axes[1, 1].grid(True)
                
            # Asymmetry heatmap
            zmin, zmax = Asymmetry_min, Asymmetry_max
            number_levels = 100
            levels = np.linspace(zmin,zmax, number_levels+1)
            norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)
            mesh = axes[1, 2].pcolormesh(edges_W2, edges_Q2, Asymmetry_W2Q2.T, cmap='rainbow', norm=norm)
            cbar = fig.colorbar(mesh, ax=axes[1, 2], shrink=0.8)
            cbar.set_label("Asymmetry")
            axes[1, 2].set_xlabel("W2")
            axes[1, 2].set_ylabel("Q2")
            axes[1, 2].set_title(f"Asymmetry heatmap kine = {expname} ({label1[j]})")
            axes[1, 2].minorticks_on()
            axes[1, 2].grid(True)
                
            # ---- Save and Close ----
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            
            del fig, ax
                             
                
    del dx, dy, W2, Q2, nu, ePS, eSH, eHCAL, coin_time, helicity, IHWP, runnum, trP, trigbits, eop
    response = input(f"Good plots for {kine}? [y,n]: ")
    if str(response) == "y":
        continue
    elif str(response) == "n":
        sys.exit()
   
            

            

            

            

            

            
