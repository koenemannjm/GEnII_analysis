import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
from scipy.interpolate import NearestNDInterpolator

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

sim_dirpath = "/w/halla-scshelf2102/sbs/koeneman/GEnII_analysis/outfiles/rootfiles/"
raw_dirpath = "/w/halla-scshelf2102/sbs/koeneman/GEnII_analysis/outfiles/HUNTERrootfiles/"

kine = sys.argv[1]

if kine == "kin2":
    E = 4.291 #GeV
    target = "He3"
    datapass = "pass2"
    expname = "GEN2"
    expname_sim ="GEN2"
    Pkin = -1
    coin_time_mean = 129.977
    coin_time_sigma = 1.81583
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
    expname_sim ="GEN3"
    Pkin = 1
    coin_time_mean = 120.533
    coin_time_sigma = 1.81322
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
    expname_sim ="GEN4"
    Pkin = 1
    coin_time_mean = 121.692
    coin_time_sigma = 1.94162
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
    expname_sim ="GEN4"
    Pkin = 1
    coin_time_mean = 185.835
    coin_time_sigma = 2.23273
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


raw_filename = f"Final_data_{expname}_sbs100p_nucleon_np_model2.root"
sim_filename = f"IN_sim_{expname_sim}_sbs100p_nucleon_np.root"

raw_path = raw_dirpath + raw_filename
sim_path = sim_dirpath + sim_filename

tree="Tout"

raw_data = ROOT.RDataFrame(tree,raw_path)
sim_data = ROOT.RDataFrame(tree,sim_path)

raw_branch_variables = ["dx","dy","W2","trP","ePS","eSH","eHCAL","coin_time","helicity","IHWP","runnum","grinch_clus_size","grinch_clus_trackindex"]
sim_branch_variables = ["dx","dy","MC.mc_sigma","MC.mc_sigmaold","e.kine.W2","e.kine.Q2","bb.ps.e","sbs.hcal.e","bb.etot_over_p"]

#raw_data_numpy = raw_data.AsNumpy(raw_branch_variables)
sim_data_numpy = sim_data.AsNumpy(sim_branch_variables)

etot_over_p_first = np.array([vec[0] if len(vec) > 0 else np.nan for vec in sim_data_numpy["bb.etot_over_p"]])

sim_data_numpy["bb.etot_over_p"] = etot_over_p_first


raw_dx = raw_data_numpy["dx"]
raw_dy = raw_data_numpy["dy"]
raw_W2 = raw_data_numpy["W2"]
raw_ePS = raw_data_numpy["ePS"]
raw_eHCAL = raw_data_numpy["eHCAL"]
raw_coin_time = raw_data_numpy["coin_time"]
raw_trP = raw_data_numpy["trP"]
raw_eSH = raw_data_numpy["eSH"]
raw_eop = (raw_ePS + raw_eSH)/raw_trP

sim_dx = sim_data_numpy["dx"]
sim_dy = sim_data_numpy["dy"]
sim_W2 = sim_data_numpy["e.kine.W2"]
sim_eSP = sim_data_numpy["bb.ps.e"]
sim_eHCAL = sim_data_numpy["sbs.hcal.e"]
sim_eop = sim_data_numpy["bb.etot_over_p"]
sim_XZsigma = sim_data_numpy["MC.mc_sigma"]
sim_CBsigma = sim_data_numpy["MC.mc_sigmaold"]

cut = (np.abs(sim_dx)<4)&(sim_dy<2)&(sim_dy>-10)

sim_dx = sim_dx[cut]
sim_dy = sim_dy[cut]
sim_XZsigma = sim_XZsigma[cut]
sim_CBsigma = sim_CBsigma[cut]
sim_eop = sim_eop[cut]
sim_eHCAL = sim_eHCAL[cut]
sim_W2 = sim_W2[cut]
sim_eSP = sim_eSP[cut]

sim_Nevents = len(sim_dx)
nbins = int(np.sqrt(sim_Nevents))
plt.hist(sim_dx, bins=nbins, weights=1.47e5*sim_XZsigma,label="XZ",alpha=0.5)
plt.hist(sim_dx, bins=nbins, weights=sim_CBsigma,label="CB",alpha=0.5)
plt.xlabel("dx")
plt.ylabel("d2sigma/dEprimedOmega")
plt.legend()
plt.minorticks_on()
plt.savefig(f"fig_xsec_XZ_CB_{expname_sim}.pdf", dpi=300, bbox_inches='tight')
plt.close()

anti_cut = (np.abs(sim_dy) > 1.0) & (sim_W2 < 2.0) & (sim_eop > 0.6)

sim_dx_anti = sim_dx[anti_cut]
sim_dy_anti = sim_dy[anti_cut]
sim_XZsigma_anti = sim_XZsigma[anti_cut]
sim_CBsigma_anti = sim_CBsigma[anti_cut]

sim_Nevents_anti = len(sim_dx_anti)
nbins_anti = int(np.sqrt(sim_Nevents_anti))
plt.hist(sim_dx_anti, bins=nbins_anti, weights=1.47e5*sim_XZsigma_anti,label="XZ",alpha=0.5)
plt.hist(sim_dx_anti, bins=nbins_anti, weights=sim_CBsigma_anti,label="CB",alpha=0.5)
plt.xlabel("dx")
plt.ylabel("d2sigma/dEprimedOmega")
plt.legend()
plt.minorticks_on()
plt.savefig(f"fig_xsec_XZ_CB_{expname_sim}_anti.pdf", dpi=300, bbox_inches='tight')
plt.close()
qe_cut = (np.abs(sim_dy) < 1.0) & (sim_W2 < 2.0) & (sim_eop > 0.6)

sim_dx_qe = sim_dx[qe_cut]
sim_dy_qe = sim_dy[qe_cut]
sim_XZsigma_qe = sim_XZsigma[qe_cut]
sim_CBsigma_qe = sim_CBsigma[qe_cut]

sim_Nevents_qe = len(sim_dx_qe)
nbins_qe = int(np.sqrt(sim_Nevents_qe))
plt.hist(sim_dx_qe, bins=nbins_qe, weights=1.47e5*sim_XZsigma_qe,label="XZ",alpha=0.5)
plt.hist(sim_dx_qe, bins=nbins_qe, weights=sim_CBsigma_qe,label="CB",alpha=0.5)
plt.xlabel("dx")
plt.ylabel("d2sigma/dEprimedOmega")
plt.legend()
plt.minorticks_on()
plt.savefig(f"fig_xsec_XZ_CB_{expname_sim}_qe.pdf", dpi=300, bbox_inches='tight')
plt.close()

