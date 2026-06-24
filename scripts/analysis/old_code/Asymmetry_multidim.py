import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import resource
import gc
from matplotlib.backends.backend_pdf import PdfPages

data_dirpath = "/w/halla-scshelf2102/sbs/koeneman/GEnII_analysis/outfiles/HUNTERrootfiles/"

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
    print("loading trP...")
    trP = raw_data_numpy["trP"]
    print("loading trigbits...")
    trigbits = raw_data_numpy["trigbits"]
    eop = (ePS + eSH)/trP

    del raw_data_numpy, raw_data
    gc.collect()
    
    print("Cutting IHWP != 1 or -1 and helicity != 1 or -1")
    
    label1 = ["QE"]
    #label1 = ["QE","Anti"]
    value = [(np.abs(dy)<0.5)&(W2_data>0.0)&(W2_data<2.0)&(np.abs(dx)<0.5),(np.abs(dy - (dy_max_anti+dy_min_anti)/2)>np.abs(dy_max_anti-dy_min_anti)/2)&(W2_data>-4.0)&(W2_data<12)]
    
    for j in range(len(label1)):
        with PdfPages(f"fig_asymmetry_multidim_{kine}_{label1[j]}.pdf") as pdf:
            print(f"Starting {label1[j]} analysis...")

            cut = (np.abs(IHWP) == 1) & (np.abs(helicity) == 1) & value[j] & (ePS>0.2) & (eHCAL>0.025) & (np.abs(coin_time - coin_time_mean)<coin_time_sigma) & (eop > 0.7) & (trigbits == 4) # & (grinch_clus_size>=3) & (grinch_clus_trackindex >=2)
            dx_cut = dx[cut]
            dy_cut = dy[cut]
            W2_cut = W2_data[cut]
            eHCAL_cut = eHCAL[cut]
            coin_time_cut = coin_time[cut]
            helicity_cut = helicity[cut]
            IHWP_cut = IHWP[cut]
            runnum_cut = runnum[cut]
            eop_cut = eop[cut]

            del cut
            gc.collect()

            raw_helicity_cut = IHWP_cut*Pkin*helicity_cut

            unique_runnum = np.unique(runnum_cut)
            for run_number in unique_runnum:
                print(f"Run = {run_number} | ({np.where(unique_runnum == run_number)[0][0]}/{len(unique_runnum)})")
                singel_runum = (runnum_cut == run_number)
                
                dx_run = dx_cut[singel_runum]
                dy_run = dy_cut[singel_runum]
                W2_run = W2_cut[singel_runum]
                raw_helicity_run = raw_helicity_cut[singel_runum]
                coin_time_run = coin_time_cut[singel_runum]
                eop_run = eop_cut[singel_runum]
        
                bins_dx = np.histogram_bin_edges(dx_run, bins='auto')
                bins_dy = np.histogram_bin_edges(dy_run, bins='auto')
                bins_coin_time = np.histogram_bin_edges(coin_time_run, bins='auto')
                bins_W2 = np.histogram_bin_edges(W2_run, bins='auto')
                bins_eop = np.histogram_bin_edges(eop_run, bins='auto')
            
                data = np.vstack((dx_run, dy_run, coin_time_run, W2_run, eop_run)).T
                Counts, edges = np.histogramdd(data, bins=(bins_dx, bins_dy, bins_coin_time, bins_W2, bins_eop))
                del data
                gc.collect

                Counts_dx = np.sum(Counts, axis=(1,2,3,4))
                edges_dx = edges[0]
                dx_centers = 0.5*(edges_dx[:-1] + edges_dx[1:])
                Counts_dy = np.sum(Counts, axis=(0,2,3,4))
                edges_dy = edges[1]
                dy_centers = 0.5*(edges_dy[:-1] + edges_dy[1:])
                Counts_coin_time = np.sum(Counts, axis=(0,1,3,4))
                edges_coin_time = edges[2]
                coin_time_centers = 0.5*(edges_coin_time[:-1] + edges_coin_time[1:])
                Counts_W2 = np.sum(Counts, axis=(0,1,2,4))
                edges_W2 = edges[3]
                W2_centers = 0.5*(edges_W2[:-1] + edges_W2[1:])
                Counts_eop = np.sum(Counts, axis=(0,1,2,3))
                edges_eop = edges[4]
                eop_centers = 0.5*(edges_eop[:-1] + edges_eop[1:])

                # Up helicity events
                hel = (raw_helicity_run == 1)
                dx_run_hel = dx_run[hel]
                dy_run_hel = dy_run[hel]
                W2_run_hel = W2_run[hel]
                coin_time_run_hel = coin_time_run[hel]
                eop_run_hel = eop_run[hel]
            
                data = np.vstack((dx_run_hel, dy_run_hel, coin_time_run_hel, W2_run_hel, eop_run_hel)).T
                Counts_up, edges_junk = np.histogramdd(data, bins=(bins_dx, bins_dy, bins_coin_time, bins_W2, bins_eop))
                del data
                gc.collect
                
                Counts_dx_up = np.sum(Counts_up, axis=(1,2,3,4))
                Counts_dy_up = np.sum(Counts_up, axis=(0,2,3,4))
                Counts_coin_time_up = np.sum(Counts_up, axis=(0,1,3,4))
                Counts_W2_up = np.sum(Counts_up, axis=(0,1,2,4))
                Counts_eop_up = np.sum(Counts_up, axis=(0,1,2,3))
                
                # Down helicity events
                hel = (raw_helicity_run == -1)
                dx_run_hel = dx_run[hel]
                dy_run_hel = dy_run[hel]
                W2_run_hel = W2_run[hel]
                coin_time_run_hel = coin_time_run[hel]
                eop_run_hel = eop_run[hel]
            
                data = np.vstack((dx_run_hel, dy_run_hel, coin_time_run_hel, W2_run_hel, eop_run_hel)).T
                Counts_dwn, edges_junk = np.histogramdd(data, bins=(bins_dx, bins_dy, bins_coin_time, bins_W2, bins_eop))
                del data
                gc.collect

                Counts_dx_dwn = np.sum(Counts_dwn, axis=(1,2,3,4))
                Counts_dy_dwn = np.sum(Counts_dwn, axis=(0,2,3,4))
                Counts_coin_time_dwn = np.sum(Counts_dwn, axis=(0,1,3,4))
                Counts_W2_dwn = np.sum(Counts_dwn, axis=(0,1,2,4))
                Counts_eop_dwn = np.sum(Counts_dwn, axis=(0,1,2,3))

                del edges_junk
                gc.collect()
                
                # noZeros = (Counts_up != 0) & (Counts_dwn != 0)
                # Asymmetry = (Counts_up - Counts_dwn) / (Counts_up + Counts_down)

                variables = ["dx","dy","coint_time","W2", "eop"]
                Counts = [Counts_dx, Counts_dy, Counts_coin_time, Counts_W2, Counts_eop]
                Counts_up = [Counts_dx_up, Counts_dy_up, Counts_coin_time_up, Counts_W2_up, Counts_eop_up]
                Counts_dwn = [Counts_dx_dwn, Counts_dy_dwn, Counts_coin_time_dwn, Counts_W2_dwn, Counts_eop_dwn]
                centers = [dx_centers, dy_centers, coin_time_centers, W2_centers, eop_centers]
                edges = [edges_dx, edges_dy, edges_coin_time, edges_W2, edges_eop]

                fig, ax = plt.subplots(3, 2, figsize=(14,10))
                ax[0,0].bar(centers[0], Counts[0], width = np.diff(edges[0]), align='center', edgecolor='k',label='all')
                ax[0,0].bar(centers[0], Counts_up[0], width = np.diff(edges[0]), align='center', edgecolor='k',label='rawhel=1')
                ax[0,0].bar(centers[0], Counts_dwn[0], width = np.diff(edges[0]), align='center', edgecolor='k',label='rawhel=-1')
                ax[0,0].set_xlabel(variables[0])
                ax[0,0].set_ylabel('Counts')
                ax[0,0].legend()
                ax[0,0].minorticks_on()
                ax[0,0].grid(True)

                ax[0,1].bar(centers[1], Counts[1], width = np.diff(edges[1]), align='center', edgecolor='k',label='all')
                ax[0,1].bar(centers[1], Counts_up[1], width = np.diff(edges[1]), align='center', edgecolor='k',label='rawhel=1')
                ax[0,1].bar(centers[1], Counts_dwn[1], width = np.diff(edges[1]), align='center', edgecolor='k',label='rawhel=-1')
                ax[0,1].set_xlabel(variables[1])
                ax[0,1].set_ylabel('Counts')
                ax[0,1].legend()
                ax[0,1].minorticks_on()
                ax[0,1].grid(True)

                ax[1,0].bar(centers[2], Counts[2], width = np.diff(edges[2]), align='center', edgecolor='k',label='all')
                ax[1,0].bar(centers[2], Counts_up[2], width = np.diff(edges[2]), align='center', edgecolor='k',label='rawhel=1')
                ax[1,0].bar(centers[2], Counts_dwn[2], width = np.diff(edges[2]), align='center', edgecolor='k',label='rawhel=-1')
                ax[1,0].set_xlabel(variables[2])
                ax[1,0].set_ylabel('Counts')
                ax[1,0].legend()
                ax[1,0].minorticks_on()
                ax[1,0].grid(True)

                ax[1,1].bar(centers[3], Counts[3], width = np.diff(edges[3]), align='center', edgecolor='k',label='all')
                ax[1,1].bar(centers[3], Counts_up[3], width = np.diff(edges[3]), align='center', edgecolor='k',label='rawhel=1')
                ax[1,1].bar(centers[3], Counts_dwn[3], width = np.diff(edges[3]), align='center', edgecolor='k',label='rawhel=-1')
                ax[1,1].set_xlabel(variables[3])
                ax[1,1].set_ylabel('Counts')
                ax[1,1].legend()
                ax[1,1].minorticks_on()
                ax[1,1].grid(True)

                ax[2,0].bar(centers[4], Counts[4], width = np.diff(edges[4]), align='center', edgecolor='k',label='all')
                ax[2,0].bar(centers[4], Counts_up[4], width = np.diff(edges[4]), align='center', edgecolor='k',label='rawhel=1')
                ax[2,0].bar(centers[4], Counts_dwn[4], width = np.diff(edges[4]), align='center', edgecolor='k',label='rawhel=-1')
                ax[2,0].set_xlabel(variables[4])
                ax[2,0].set_ylabel('Counts')
                ax[2,0].legend()
                ax[2,0].minorticks_on()
                ax[2,0].grid(True)

                fig.suptitle(f"Helicity Dependence of Runum = {run_number} kine = {expname} ({label1[j]})")
                fig.delaxes(ax[2, 1])
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)

                #usage = resource.getrusage(resource.RUSAGE_SELF)
                #print(f"Memory used: {usage.ru_maxrss / 1024:.2f} MB")

                del fig, ax, Counts_up, Counts_dwn, Counts, edges, centers
                gc.collect()
                

            del dx_cut, dy_cut, W2_cut, eHCAL_cut, coin_time_cut, helicity_cut, IHWP_cut, runnum_cut, eop_cut, pdf
            gc.collect()

    del dx, dy, ePS, eSH, eHCAL, coin_time, helicity, IHWP, runnum, eop, trigbits, trP
    gc.collect()
    
            
