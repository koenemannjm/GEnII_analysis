import ROOT
import numpy as np
import matplotlib.pyplot as plt
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
        
    data_filename = f"Final_data_{expname}_sbs100p_nucleon_np_model2.root"
        
    data = data_dirpath + data_filename
        
    tree = "Tout"
        
    print(f"Reading in data for {kine}")
    raw_data = ROOT.RDataFrame(tree, data)
    branch_variables = ["dx","dy","W2","ePS","eSH","eHCAL","coin_time","helicity","IHWP","runnum","grinch_clus_size","grinch_clus_trackindex","trP","trigbits"]
    raw_data_numpy = raw_data.AsNumpy(branch_variables)
    
    print("Extracting variables")
    
    print("Loading dx...")
    dx = raw_data_numpy["dx"]
    print("Loading dy...")
    dy = raw_data_numpy["dy"]
    print("Loading W2...")
    W2_data = raw_data_numpy["W2"]
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
    print("Loading grinch_clus_size...")
    grinch_clus_size = raw_data_numpy["grinch_clus_size"]
    print("Loading grinch_clus_trackindex...")
    grinch_clus_trackindex = raw_data_numpy["grinch_clus_trackindex"]
    print("loading trP...")
    trP = raw_data_numpy["trP"]
    print("loading trigbits...")
    trigbits = raw_data_numpy["trigbits"]
    eop = (ePS + eSH)/trP
    
    print("Cutting IHWP != 1 or -1 and helicity != 1 or -1")
    
    label1 = ["QE","Anti"]
    value = [(np.abs(dy)<0.5)&(W2_data>0.0)&(W2_data<2.0)&(np.abs(dx)<0.5),(np.abs(dy - (dy_max_anti+dy_min_anti)/2)>np.abs(dy_max_anti-dy_min_anti)/2)&(W2_data>-4.0)&(W2_data<12)]
    with PdfPages(f"fig_asymmetry_{kine}.pdf") as pdf:
        for j in range(2):
            
            cut = (np.abs(IHWP) == 1) & (np.abs(helicity) == 1) & value[j] & (ePS>0.2) & (eHCAL>0.025) & (np.abs(coin_time - coin_time_mean)<coin_time_sigma) & (eop > 0.7) & (trigbits == 4) # & (grinch_clus_size>=3) & (grinch_clus_trackindex >=2)
            dx_cut = dx[cut]
            dy_cut = dy[cut]
            W2_cut = W2_data[cut]
            ePS_cut = ePS[cut]
            eSH_cut = eSH[cut]
            eHCAL_cut = eHCAL[cut]
            coin_time_cut = coin_time[cut]
            helicity_cut = helicity[cut]
            IHWP_cut = IHWP[cut]
            runnum_cut = runnum[cut]
            
            Nevents = len(dy_cut)
            # Nbins = int(np.sqrt(Nevents))
            bins_dy = np.histogram_bin_edges(dy_cut, bins='fd')
            bins_W2 = np.histogram_bin_edges(W2_cut, bins='fd')
            CountsW2dy, xedges_W2, yedges_dy = np.histogram2d(W2_cut, dy_cut, bins=[bins_W2,bins_dy])
            # CountsW2dy, xedges_W2, yedges_dy = np.histogram2d(W2_cut, dy_cut, bins=Nbins)
            maskedCountsW2dy = np.ma.masked_where(CountsW2dy == 0, CountsW2dy)
            
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
            
            nu = ebeam - eprime
            Q2 = 2 * eprime * ebeam * (1 - np.cos(etheta))
            X = Q2 / (2 * MN * nu)
            W2 = MN**2 - Q2 + 2 * MN * nu
            
            theta_star = np.array(theta_star)
            phi_star = np.array(phi_star)
            
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

            true_helicity_cut = IHWP_cut*Pkin*helicity_cut
            
            up_hel_cut = (true_helicity_cut == 1)
            W2_up = W2_cut[up_hel_cut]
            dx_up = dx_cut[up_hel_cut]
            dwn_hel_cut = (true_helicity_cut == -1)
            W2_dwn = W2_cut[dwn_hel_cut]
            dx_dwn = dx_cut[dwn_hel_cut]
            
            Nevents_avg = len(W2_cut)
            bins_W2_cut = np.histogram_bin_edges(W2_cut, bins='auto')
            bins_dx_cut = np.histogram_bin_edges(dx_cut, bins='auto')
            # nbins = int(np.sqrt(np.sqrt(Nevents_avg)))
            Counts_up, xedges_up = np.histogram(W2_up, bins=bins_W2_cut)
            Counts_dx_up, dxedges_up = np.histogram(dx_up,bins=bins_dx_cut)
            errCounts_up = np.sqrt(Counts_up)
            errCounts_dx_up = np.sqrt(Counts_dx_up)
            Counts_dwn, xedges_dwn = np.histogram(W2_dwn, bins=bins_W2_cut)
            Counts_dx_dwn, dxedges_dwn = np.histogram(dx_dwn, bins=bins_dx_cut)
            errCounts_dwn = np.sqrt(Counts_dwn)
            errCounts_dx_dwn = np.sqrt(Counts_dx_dwn)
                
            W2_centers_up = 0.5*(xedges_up[1:] + xedges_up[:-1])
            W2_centers_dwn = 0.5*(xedges_dwn[1:] + xedges_dwn[:-1])
            W2_centers = 0.5*(W2_centers_up + W2_centers_dwn)

            dx_centers_up = 0.5*(dxedges_up[1:] + dxedges_up[:-1])
            dx_centers_dwn = 0.5*(dxedges_dwn[1:] + dxedges_dwn[:-1])
            dx_centers = 0.5*(dx_centers_up + dx_centers_dwn)
        
            noZeros = (Counts_dwn != 0) & (Counts_up !=0)
            Counts_up = Counts_up[noZeros]
            Counts_dwn = Counts_dwn[noZeros]
            errCounts_up = errCounts_up[noZeros]
            errCounts_dwn = errCounts_dwn[noZeros]
            W2_centers = W2_centers[noZeros]
            
            Asymmetry = (Counts_up - Counts_dwn) / (Counts_up + Counts_dwn)
            errAsymmetry = 2/((Counts_up + Counts_dwn)**2)*np.sqrt((Counts_up*errCounts_dwn)**2 + (Counts_dwn*errCounts_up)**2)
            
            with np.errstate(divide='ignore', invalid='ignore'):
                relerrAsymmetry = np.where(Asymmetry != 0, errAsymmetry / Asymmetry, np.inf)
            
            smallerror = (np.abs(relerrAsymmetry) < 1.0)
            
            Asymmetry = Asymmetry[smallerror]
            W2_centers = W2_centers[smallerror]
            errAsymmetry = errAsymmetry[smallerror]
            
            print(f"Making plot {label1[j]}...")
            fig, ax = plt.subplots(figsize=(14, 8))
            axa = ax
            axb = axa.twinx()
            
            # Plot the colormesh on axb, but push it behind axa
            c = axa.pcolormesh(xedges_W2, yedges_dy, maskedCountsW2dy.T, cmap='rainbow')
            #fig.tight_layout()
            cbar = fig.colorbar(c, ax=axb, pad=0.12)
            cbar.set_label('Counts')
            
            # Plot error bars and model on top
            axb.errorbar(W2_centers, Asymmetry, yerr=errAsymmetry,
                         fmt='ko', ecolor='k', capsize=3, markersize=4, label="data", zorder=20)
            axb.plot(W2, Asymmetry_model, "m--", label="XZ model", zorder=19)
            
            # Labels and formatting
            axa.set_xlabel("W2")
            axb.set_ylabel("Asymmetry")
            axa.set_ylabel("dy")
            axa.set_title(f"Raw Asymmetry vs W2 kin = {expname} ({label1[j]})")
            axa.minorticks_on()
            axb.minorticks_on()
            axa.grid(True)
            axa.set_xlim(min(W2_centers) - 0.1, max(W2_centers) + 0.5)
            # axb.set_ylim(-1.0, 1.0)
            axb.set_ylim(-0.25, 0.25)
            axb.legend()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

            noZeros = (Counts_dx_dwn != 0) & (Counts_dx_up !=0)
            Counts_up = Counts_dx_up[noZeros]
            Counts_dwn = Counts_dx_dwn[noZeros]
            errCounts_up = errCounts_dx_up[noZeros]
            errCounts_dwn = errCounts_dx_dwn[noZeros]
            dx_centers = dx_centers[noZeros]
            
            Asymmetry = (Counts_up - Counts_dwn) / (Counts_up + Counts_dwn)
            errAsymmetry = 2/((Counts_up + Counts_dwn)**2)*np.sqrt((Counts_up*errCounts_dwn)**2 + (Counts_dwn*errCounts_up)**2)

            with np.errstate(divide='ignore', invalid='ignore'):
                relerrAsymmetry = np.where(Asymmetry != 0, errAsymmetry / Asymmetry, np.inf)
            
            smallerror = (np.abs(relerrAsymmetry) < 1.0)
            
            Asymmetry = Asymmetry[smallerror]
            dx_centers = dx_centers[smallerror]
            errAsymmetry = errAsymmetry[smallerror]
            
            fig, ax = plt.subplots(figsize=(14, 8))
            ax.errorbar(dx_centers, Asymmetry,
                         yerr=errAsymmetry,
                         fmt='o', capsize=3, markersize=4,label=f"data")
            ax.set_xlabel("dx")
            ax.set_ylabel("Asymmetry")
            ax.set_title(f"Raw Asymmetry vs dx kin = {expname} ({label1[j]})")
            ax.minorticks_on()
            ax.grid(True)
            ax.legend()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
            
            del fig
            del ax
            del axa
            del axb
            del dx_cut 
            del dy_cut
            del W2_cut
            del ePS_cut
            del eSH_cut
            del eHCAL_cut
            del coin_time_cut
            del helicity_cut
            del IHWP_cut
            del runnum_cut
            del cut
        del raw_data
        del raw_data_numpy
        del dx
        del dy
        del W2_data
        del ePS
        del eSH
        del eHCAL
        del coin_time
        del helicity
        del IHWP
        del runnum
        del grinch_clus_size
        del grinch_clus_trackindex
        del trP
        del eop
        del value
        del trigbits

print("Done! <3")
