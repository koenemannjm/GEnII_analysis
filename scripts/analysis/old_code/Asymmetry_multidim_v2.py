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

            bin_edges_dx = np.histogram_bin_edges(dx_cut, bins='auto')
            bin_edges_dy = np.histogram_bin_edges(dy_cut, bins='auto')
            bin_edges_coin_time = np.histogram_bin_edges(coin_time_cut, bins='auto')
            bin_edges_W2 = np.histogram_bin_edges(W2_cut, bins='auto')
            bin_edges_eop = np.histogram_bin_edges(eop_cut, bins='auto')

            Counts_dx = []
            CountsUp_dx = []
            CountsDwn_dx = []
            Asymmetry_dx = []
            Counts_dy = []
            CountsUp_dy = []
            CountsDwn_dy = []
            Asymmetry_dy = []
            Counts_coin_time = []
            CountsUp_coin_time = []
            CountsDwn_coin_time = []
            Asymmetry_coin_time = []
            Counts_W2 = []
            CountsUp_W2 = []
            CountsDwn_W2 = []
            Asymmetry_W2 = []
            Counts_eop = []
            CountsUp_eop = []
            CountsDwn_eop = []
            Asymmetry_eop = []

            for 

            
