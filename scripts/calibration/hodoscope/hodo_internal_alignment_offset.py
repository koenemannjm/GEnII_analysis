import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

# --- file --- #
kin = "GEN3"
kin_pass = "pass2"

path = "../../outfiles/rootfiles/"
filename = f"QE_{kin_pass}_try8_v4_data_{kin}_sbs100p_nucleon_np.root"
#filename = f"hodo_timing_{kin_pass}_try8_v4_data_{kin}_He3.root"
offsetname = f"offset_hodo_internal_alignment_{kin_pass}_{kin}.txt"
filepath = path + filename

# --- Constants --- #
c = 0.299792458 # units -> m/ns

# --- Extracting Timing Branches --- #

# -- Branch Names -- #
# detector_time_variable = "bb.hodotdc.clus.tmean"
# detector_blkid_variable = "bb.hodotdc.clus.id"

detector_time_variable = "bb.hodotdc.clus.bar.tdc.meantime"
detector_blkid_variable = "bb.hodotdc.clus.bar.tdc.id"

with uproot.open(filepath) as ROOTFile:
    tree = ROOTFile["Tout"]

    Nevents = tree.num_entries

    blktime = tree[detector_time_variable].arrays(library="ak")[detector_time_variable]
    blkid = tree[detector_blkid_variable].arrays(library="ak")[detector_blkid_variable]

mask = ak.num(blktime) > 0

blktime = blktime[mask]
blkid = blkid[mask]

blktime_main_clus = ak.firsts(blktime).to_numpy()
blkid_main_clus = ak.firsts(blkid).to_numpy()
Nevents_mask = len(blktime_main_clus)

# --- Extracting Offsets --- #

offset_blkid, offset_blktime = np.loadtxt(offsetname, delimiter="\t", unpack=True)

# --- Applying Correction --- #
max_id = int(offset_blkid.max()) + 1
lookup_table = np.full(max_id, np.nan)
lookup_table[offset_blkid.astype(int)] = offset_blktime
offsets = lookup_table[blkid_main_clus.astype(int)]
print("Applying Correction")
blktime_main_clus_cor = blktime_main_clus - offsets

print("Creating plots")
savefig_name = f"plots_hodo_internal_alignment_{kin_pass}_{kin}.pdf"
with PdfPages(savefig_name) as pdf:
    
    plt.hist2d(blkid_main_clus, blktime_main_clus, bins=(200, 90), range=[[-0.01,89.01],[-10,15]], cmap='rainbow')
    plt.xlabel("Block ID")
    plt.ylabel("Block Time")
    plt.title("Main Cluster (Uncorrected)")
    pdf.savefig()  # Save this figure to the PDF
    plt.close()

    plt.hist2d(blkid_main_clus, blktime_main_clus_cor, bins=(200, 90), range=[[-0.01,89.01],[-10,15]], cmap='rainbow')
    plt.xlabel("Block ID")
    plt.ylabel("Block Time")
    plt.title("Main Cluster (Corrected)")
    pdf.savefig()
    plt.close()

    plt.hist(blktime_main_clus, bins=200, range=[-20,20], alpha=0.5, label="Uncorrected")
    plt.hist(blktime_main_clus_cor, bins=200, range=[-20,20], alpha=0.5, label="Corrected")
    plt.xlabel("Block Time")
    plt.ylabel("Counts")
    plt.legend()
    pdf.savefig()
    plt.close()

print("All Done :)")
