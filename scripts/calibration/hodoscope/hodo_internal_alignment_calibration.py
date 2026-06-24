import ROOT
import numpy as np
from scipy.optimize import curve_fit
from numpy.linalg import svd
import os

#########################################
#########################################
##
##  Created by: Jacob Koenemann
##  email: bxy3zr@virginia.edu
##  
##  Last updated: 10/08/2025 
##
##  Puprose: This script is meant to do
##  internal alignment of all block/bar
##  based dectors that do timing. This
##  script is based on Hunter Presley's
##  python jupyter notebook.
##
#########################################
#########################################

# --- Complimentary Functions --- #
def AccumulateInternalAlignment(filepaths, M, b, event_counts, nbars):
    """
    Inputs:
    
    filepaths ->
    M -> 
    b ->
    event_counts ->
    nbars ->
    PMT ->

    Ouputs:
    
    
    """

    ### Unpacking filepaths ###
    filepath, vscintcalibration, ToTcalibration = filepaths
    
    ### Constants ###
    speed_of_light = 0.299792458 # meters/ns
    tdiff_max = 10 # units -> ns
    
    hodo_bar_height = 0.025 # meters
    hodo_bar_width = 0.60 # meters
    z_hodo = 1.854454 # meters

    BBdist = 1.63 # meters

    hodscope_height = nbars*hodo_bar_height

    etof0 = (BBdist + 3.0)/speed_of_light

    ### Forming Chain ###
    chain = ROOT.TChain("Tout")
    chain.Add(filepath)
    print(f"Initialized TChain: {filepath} with TTree = Tout")

    ### Creating TTreeReader ###
    reader = ROOT.TTreeReader(chain)

    ### Variables needed for Calibration ###
    nhits = ROOT.TTreeReaderValue("int")(reader, "Ndata.bb.hodotdc.clus.bar.tdc.id")
    yfp = ROOT.TTreeReaderArray("double")(reader, "bb.tr.y")
    phfp = ROOT.TTreeReaderArray("double")(reader, "bb.tr.ph")
    barid = ROOT.TTreeReaderArray("double")(reader, "bb.hodotdc.clus.bar.tdc.id")
    tleft = ROOT.TTreeReaderArray("double")(reader, "bb.hodotdc.clus.bar.tdc.tleft")
    tright = ROOT.TTreeReaderArray("double")(reader, "bb.hodotdc.clus.bar.tdc.tright")
    trigtime = ROOT.TTreeReaderValue("double")(reader, "bb.hodotdc.trigtime")
    totleft = ROOT.TTreeReaderArray("double")(reader, "bb.hodotdc.clus.bar.tdc.totleft")
    totright = ROOT.TTreeReaderArray("double")(reader, "bb.hodotdc.clus.bar.tdc.totright")
    etof = ROOT.TTreeReaderArray("double")(reader, "bb.hodotdc.clus.bar.tdc.etof")
    runnum = ROOT.TTreeReaderValue("double")(reader, "g.runnum")
    gevtime = ROOT.TTreeReaderValue("double")(reader, "g.evtime")

    ### Computing TBBtrig_t0 ###
    print("Computing BB trigtime offset")
    """
    while reader.Next():

        trig = trigtime.Get()[0]
        # print(trig)
        evt = int(gevtime.Get()[0])
        # print(evt)
        
        
        if (trig < 500.0):
            
            sumreftime += trig - 4.0 * (6 % (2 + 6 % (evt)))
            sumN += 1.0

    TBBtrig_t0 = sumreftime / sumN
    """
    TBBtrig_t0 = 320.38437858338546
    print(f"TBBtrig_t0 = {TBBtrig_t0}")

    ### Unpacking ToT calibration ###
    barid_element, wleft, vleft, wright, vright = np.loadtxt(
        f"./outfiles/{ToTcalibration}",
        unpack=True,
        dtype=float
    )

    ### Unpacking Vscint calibration ###
    barid_element, vscint = np.loadtxt(
        f"./outfiles/{vscintcalibration}",
        unpack=True,
        dtype=float
    )
    
    ### The event loop ###
    print("Starting Event Loop ...")
    Nevents = chain.GetEntries()
    print(f"total events = {Nevents}")
    ievent = 0
    while reader.Next():

        ievent += 1
        
        trig = trigtime.Get()[0]
        evt = int(gevtime.Get()[0])
        
        y = yfp[0] + z_hodo * phfp[0]

        dleft = min(hodo_bar_width, max(0.0, hodo_bar_width/2 - y))
        dright = min(hodo_bar_width, max(0.0, hodo_bar_width/2 + y))

        if not (trig<400):
            continue

        tref = trig - 4.0 * (6 % (2 + 6 % (evt))) - TBBtrig_t0

        if not (nhits.Get()[0]>1):
            continue
        
        for i in range(nhits.Get()[0]):
            if not ((totleft[i]>0.0) &
                    (totright[i]>0.0)
                    ):
                continue

            id_i = int(barid[i])

            timeL_i = tleft[i] + tref - wleft[id_i] * totleft[i] - dleft / vscint[id_i] - (etof[i] - etof0)
            timeR_i = tright[i] + tref - wright[id_i] * totright[i] - dright / vscint[id_i] - (etof[i] - etof0)

            time_i = 0.5 * (timeL_i + timeR_i)

            # --- LEFT PMT --- #
            
            
            if i == 0:
                tmean0=time_i
            if abs(time_i-tmean0)>tdiff_max:
                continue
            if not (0 <= id_i < nbars) or np.isnan(time_i):
                continue
            
            for j in range(i+1, nhits.Get()[0]):
                if not ((totleft[j]>0.0)
                        ):
                    continue
                id_j = int(barid[j])

                timeL_j = tleft[j] + tref - wleft[id_j] * totleft[j] - dleft / vscint[id_j] - (etof[j] - etof0)
                timeR_j = tright[j] + tref - wright[id_j] * totright[j] - dright / vscint[id_j] - (etof[j] - etof0)

                time_j = 0.5 * (timeL_j + timeR_j)

                if not (0 <= id_j < nbars) or np.isnan(time_j):
                    continue
                
                dt = time_i - time_j
                if abs(dt)>tdiff_max:
                    continue
                
                event_counts[id_i] += 1
                event_counts[id_j] += 1
                
                b[id_i] += dt
                b[id_j] -= dt
                
                M[id_i, id_i] += 1
                M[id_j, id_j] += 1
                M[id_i, id_j] -= 1
                M[id_j, id_i] -= 1

            
        if ievent%100000 == 0:
            print(f"Percent complete: {ievent/Nevents*100.0:.2f}%")
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

def ChunkwiseInternalAlignment(filepaths, nbars, refID=45):
    M = np.zeros((nbars, nbars))
    b = np.zeros(nbars)
    event_counts = np.zeros(nbars)

    M, b, event_counts = AccumulateInternalAlignment(filepaths, M, b, event_counts, nbars)

    offset = SolveInternalAlignment(M, b, event_counts, refID)
    return offset
    

# --- Internal Alignment Calibration  --- #
def InternalAlignmentCalibration(filepath,
                                 nbars,
                                 vscintcalibration="default_hodo_vscint_calibration.txt",
                                 ToTcalibration="default_hodo_ToT_calibration.txt"
                                 ):

    filepath_list = [filepath, vscintcalibration, ToTcalibration]
    
    offset_hodo = ChunkwiseInternalAlignment(filepath_list, nbars)
    
    save_data_offset = np.column_stack((np.arange(len(offset_hodo)), offset_hodo))
    offset_filename = f"./outfiles/hodo_internal_alignment_calibration.txt"
    np.savetxt(offset_filename, save_data_offset, fmt=["%-10d", "%-10.6f"], delimiter="\t")

    return offset_hodo


# --- Testing Function --- #
def main() -> None:
    print(f"Testing {os.path.basename(__file__)} ...")
    
    filename = "hodo_timing_pass2_data_GEN3_He3.root"
    dirpath = "./../../../outfiles/rootfiles/"
    filepath = dirpath + filename
    
    nbars = 90

    vscintcalibration = "fitresults_hodo_vscint_calibration.txt"
    ToTcalibration = "fitresults_hodo_ToT_calibration.txt"

    offset_hodo = InternalAlignmentCalibration(filepath,
                                               nbars,
                                               vscintcalibration=vscintcalibration,
                                               ToTcalibration=ToTcalibration
                                               )
    
    def_save_data_offset = np.column_stack((np.arange(len(offset_hodo)), np.zeros(nbars)))
    def_offset_filename = f"./outfiles/default_hodo_internal_alignment_calibration.txt"
    np.savetxt(def_offset_filename, def_save_data_offset, fmt=["%-10d", "%-10.6f"], delimiter="\t")

if __name__=="__main__":
    main()




