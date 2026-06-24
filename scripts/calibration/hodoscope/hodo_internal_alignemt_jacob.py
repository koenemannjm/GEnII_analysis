import uproot
import numpy as np
import awkward as ak

# --- file --- #
kin = "GEN2"
kin_pass = "pass2"

path = "../../outfiles/rootfiles/"
# filename = f"QE_{kin_pass}_try8_v4_data_{kin}_sbs100p_nucleon_np.root"
filename = f"hodo_timing_{kin_pass}_try8_v4_data_{kin}_He3.root"
filepath = path + filename

# --- Constants --- #
c = 0.299792458 # units -> m/ns
tdiff_max = 10 # units -> ns

# --- Extracting Timing Branches --- #

# -- Branch Names -- #
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

blktime_firsts = ak.firsts(blktime)
Nevents_mask = len(blktime_firsts)

blktime_flat = ak.flatten(blktime).to_numpy()
blkid_flat = ak.flatten(blkid).to_numpy()

timeoffset = []
offsetid = []

bar_index = np.arange(0,90,1)

for bar in bar_index:
    cut = (blkid_flat == bar)

    blktime_bar = blktime_flat[cut]
    totalblktime_bar = np.sum(blktime_bar)

    timeoffset_bar = totalblktime_bar/Nevents_mask

    timeoffset.append(timeoffset_bar)
    offsetid.append(bar)

timeoffset = np.array(timeoffset)
offsetid = np.array(offsetid)

for i, offset in enumerate(timeoffset):
    
    print(f"{i} | {offset:.8f}")
    
save_data_offset = np.column_stack((np.arange(len(timeoffset)), timeoffset))

offset_filename = f"offset_hodo_internal_alignment_jacob_{kin_pass}_{kin}.txt"

np.savetxt(offset_filename, save_data_offset, fmt=["%d","%.8f"], delimiter="\t")



