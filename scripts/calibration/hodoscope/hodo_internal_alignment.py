import ROOT
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy.linalg import svd

#########################################
#########################################
##
##  Created by: Jacob Koenemann
##  email: bxy3zr@virginia.edu
##  
##  Last updated: 08/11/2025 
##
##  Puprose: This script is meant to do
##  internal alignment of all block/bar
##  based dectors that do timing. This
##  script is based on Hunter Presley's
##  python jupyter notebook.
##
#########################################
#########################################

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

def AccumulateInternalAlignment(M, b, n_events, event_counts, blktime, blkid, nblk, n_blocks):
    """
    Inputs:
    
    M -> 
    b ->
    event_counts ->
    blktime ->
    blkid ->
    nblk ->
    n_blocks ->

    Ouputs:
    
    
    """
    print(f"Number of events = {n_events}")

    for ev in range(n_events):
        if ev % 10000 == 0:
            print(f"    Processing event {ev}/{n_events} {ev/n_events*100:.4f}%", end='\r', flush=True)

        times = blktime[ev][:int(nblk[ev])]
        ids = blkid[ev][:int(nblk[ev])]
        
        for i in range(int(nblk[ev])):

            id_i = int(ids[i])
            time_i = times[i]
            
            if i == 0:
                tmean0=time_i
            if abs(time_i-tmean0)>tdiff_max:
                continue
            if not (0 <= id_i < n_blocks) or np.isnan(time_i):
                continue
            
            for j in range(i+1, int(nblk[ev])):
                id_j = int(ids[j])
                time_j = times[j]
                if not (0 <= id_j < n_blocks) or np.isnan(time_j):
                    continue
                
                dt = time_i - time_j
                if abs(dt)>tdiff_max:
                    continue
                #if ev % 100000 == 0:
                #    print(f"nblk= {nblk[ev]}, length:{len(ids)}, dt:{dt}")
                event_counts[id_i] += 1
                event_counts[id_j] += 1
                
                b[id_i] += dt
                b[id_j] -= dt
                
                M[id_i, id_i] += 1
                M[id_j, id_j] += 1
                M[id_i, id_j] -= 1
                M[id_j, id_i] -= 1

    return M, b, event_counts
def SolveInternalAlignment(M, b, event_counts, refID, min_events=100):
    n_blocks = M.shape[0]
    for i in range(n_blocks):
        if event_counts[i] < min_events:
            b[i] = 0
            M[i, :] = 0
            M[:, i] = 0
            M[i, i] = 1

    U, S, Vt = svd(M)
    S_inv = np.array([1/s if s > 1e-10 else 0 for s in S])
    M_inv = Vt.T @ np.diag(S_inv) @ U.T
    delta = M_inv @ b

    corr = -delta[refID]
    delta += corr
    
    return delta

def ChunkwiseInternalAlignment(blktime, blkid, nblk, n_events, n_blocks, max_chunks=1, refID=45):
    M = np.zeros((n_blocks, n_blocks))
    b = np.zeros(n_blocks)
    event_counts = np.zeros(n_blocks)

    i = 0
    while True:
        if max_chunks and i >= max_chunks:
            break

        M, b, event_counts = AccumulateInternalAlignment(M, b, n_events, event_counts, blktime, blkid, nblk, n_blocks)

        i += 1

    return SolveInternalAlignment(M, b, event_counts, refID)
    

# --- Extracting Timing Branches --- #

# -- Branch Names -- #
detector_time_variable = "bb.hodotdc.clus.bar.tdc.meantime"
detector_blkid_variable = "bb.hodotdc.clus.bar.tdc.id"
detector_nblk_variable = "Ndata.bb.hodotdc.clus.bar.tdc.meantime"

with uproot.open(filepath) as ROOTFile:
    tree = ROOTFile["Tout"]

    Nevents = tree.num_entries

    blktime = tree[detector_time_variable].arrays(library="ak")[detector_time_variable]
    blkid = tree[detector_blkid_variable].arrays(library="ak")[detector_blkid_variable]
    nblk = tree[detector_nblk_variable].arrays(library="ak")[detector_nblk_variable]


mask = ak.num(blktime) > 1

blktime = blktime[mask]
blkid = blkid[mask]
nblk = nblk[mask]

blktime_firsts = ak.firsts(blktime)
Nevents_mask = len(blktime_firsts)

# samp_event_nums = 1000000
# samp_events = np.linspace(0, Nevents_mask - 1, samp_event_nums, dtype=int)

# blktime_samp = blktime[samp_events]
# blkid_samp = blkid[samp_events]
# nblk_samp = nblk[samp_events]

# --- Running the Alignment Script --- #

# offset_hodo = ChunkwiseInternalAlignment(blktime_samp, blkid_samp, nblk_samp, len(blktime_samp), 90)
offset_hodo = ChunkwiseInternalAlignment(blktime, blkid, nblk, Nevents_mask, 90)


for i, offset in enumerate(offset_hodo):
    
    print(f"{i} | {offset:.8f}")
    
save_data_offset = np.column_stack((np.arange(len(offset_hodo)), offset_hodo))

offset_filename = f"offset_hodo_internal_alignment_{kin_pass}_{kin}.txt"

np.savetxt(offset_filename, save_data_offset, fmt=["%d","%.8f"], delimiter="\t")

