import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os 

def VscintCalibration(filepath: str,
                      nbars: int,
                      ToTcalibration: str ="default_hodo_ToT_calibration.txt",
                      internaloffset: str ="default_hodo_internal_alignment_calibration.txt"
                      ):

    ### Constants ###
    speed_of_light = 0.299792458 # meters/ns
    
    hodo_bar_height = 0.025 # meters
    hodo_bar_width = 0.60 # meters
    z_hodo = 1.854454 # meters
    # z_hodo = 1.63 # meters

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
    nhits = ROOT.TTreeReaderValue["int"](reader, "Ndata.bb.hodotdc.clus.bar.tdc.id")
    yfp = ROOT.TTreeReaderArray("double")(reader, "bb.tr.y")
    xfp = ROOT.TTreeReaderArray("double")(reader, "bb.tr.x")
    thfp = ROOT.TTreeReaderArray("double")(reader, "bb.tr.th")
    phfp = ROOT.TTreeReaderArray("double")(reader, "bb.tr.ph")
    barid = ROOT.TTreeReaderArray["double"](reader, "bb.hodotdc.clus.bar.tdc.id")
    tleft = ROOT.TTreeReaderArray["double"](reader, "bb.hodotdc.clus.bar.tdc.tleft")
    tright = ROOT.TTreeReaderArray["double"](reader, "bb.hodotdc.clus.bar.tdc.tright")
    trigtime = ROOT.TTreeReaderValue["double"](reader, "bb.hodotdc.trigtime")
    totleft = ROOT.TTreeReaderArray["double"](reader, "bb.hodotdc.clus.bar.tdc.totleft")
    totright = ROOT.TTreeReaderArray["double"](reader, "bb.hodotdc.clus.bar.tdc.totright")
    etof = ROOT.TTreeReaderArray("double")(reader, "bb.hodotdc.clus.bar.tdc.etof")
    runnum = ROOT.TTreeReaderValue["double"](reader, "g.runnum")
    gevtime = ROOT.TTreeReaderValue["double"](reader, "g.evtime")

    sumreftime = 0
    sumN = 0
    ### Computing TBBtrig_t0 ###
    print("Computing BB trigtime offset")
    """
    while reader.Next():
        
        if not trigtime.IsValid() or not gevtime.IsValid():
            continue

        trig = trigtime.Get()
        evt = gevtime.Get()

        if (trig<400):
            
            sumreftime += trig - 4.0 * (6 % (2 + 6 % (int(evt))))
            sumN += 1

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

    ### Unpacking internal alignment calibration ###
    barid_element, offset = np.loadtxt(
        f"./outfiles/{internaloffset}",
        unpack=True,
        dtype=float
    )

    # vscint = np.zeros(nbars)
    sum_d_times_t = np.zeros(nbars)
    sum_d_sqr = np.zeros(nbars)
    sum_t = np.zeros(nbars)
    sum_d = np.zeros(nbars)
    counts = np.zeros(nbars)

    ### Forming histograms ###
    hist2d_list = []
    hodoymin, hodoymax = -hodo_bar_width/2 , hodo_bar_width/2
    tdiff_min, tdiff_max = -20, 20
    nbins = 100
    for bar in range(nbars):
        histname = f"tdiff_vs_hodoy_bar{bar}"
        histtitle = f"tdiff vs hodoy bar={bar}"
        hist2d = ROOT.TH2F(histname,
                           histtitle,
                           nbins,
                           hodoymin,
                           hodoymax,
                           nbins,
                           tdiff_min,
                           tdiff_max
                           )
        hist2d_list.append(hist2d)
    
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
        x = xfp[0] + z_hodo * thfp[0]
        
        dleft = min(hodo_bar_width, max(0.0, hodo_bar_width/2.0 - y))
        dright = min(hodo_bar_width, max(0.0, hodo_bar_width/2.0 + y))

        if (trig<400):
            
            tref = trig - 4.0 * (6 % (2 + 6 % (evt))) - TBBtrig_t0

            for i in range(nhits.Get()[0]):

                barid_i = int(barid[i])
                barx = 0.5 * hodscope_height - hodo_bar_height * (0.5 + barid_i)
                diffx = x - barx
                
                if ((totleft[i]>0.0) &
                    (totright[i]>0.0) &
                    (abs(dleft - dright)<0.3) &
                    (abs(diffx)<0.5*hodo_bar_height)
                    ):

                    tleft_corr = tleft[i] + tref - wleft[barid_i] * totleft[i] - offset[barid_i] - (etof[i] - etof0)
                    tright_corr = tright[i] + tref - wright[barid_i] * totright[i] - offset[barid_i] - (etof[i] - etof0)
                    
                    sum_d_times_t[barid_i] += (dleft - dright) * (tleft_corr - tright_corr)
                    sum_d_sqr[barid_i] += (dleft - dright)**2
                    sum_t[barid_i] += (tleft_corr - tright_corr)
                    sum_d[barid_i] += (dleft - dright)
                    counts[barid_i] += 1.0

                    hist2d_list[barid_i].Fill(dleft-dright,tleft_corr-tright_corr)
                    
        if ievent%100000 == 0:
            print(f"Percent complete: {ievent/Nevents*100.0:.2f}%")

    print("Percent complete: 100.00%")
    print("Event Loop Finished!")
    
    inv_vscint = []
    tdiff_offset = []
    for i in range(nbars):
        if (counts[i]>0):
            
            avg_t = sum_t[i] / counts[i]
            avg_d = sum_d[i] / counts[i]
            avg_d_times_t = sum_d_times_t[i] / counts[i]
            avg_d_sqr = sum_d_sqr[i] / counts[i]

            if (avg_d_sqr - (avg_d)**2) != 0.0:
                inv_vscinti = (avg_d_times_t - avg_t * avg_d) / (avg_d_sqr - (avg_d)**2)
                someintercept = avg_t - inv_vscinti * avg_d

                if inv_vscinti < 0.0:
                    inv_vscinti = 1.0/0.16
            else :
                someintercept = 0.0
                inv_vscinti = 1.0/0.16

        else:
            inv_vscinti = 1.0/0.16
            someintercept = 0.0

        inv_vscint.append(inv_vscinti)
        tdiff_offset.append(someintercept)
    
    inv_vscint = np.array(inv_vscint)
    tdiff_offset = np.array(tdiff_offset)
    vscint = 1 / inv_vscint

    ### Plotting vscint for each bar ###
    pdf_filename = "./figures/plots_hodo_vscint_calibration.pdf"
    print(f"Creating PDF figure :{pdf_filename}")

    c = ROOT.TCanvas("c", "Vscint Correction", 800, 600)
    c.Print(f"{pdf_filename}[")
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    ROOT.gStyle.SetNumberContours(50)

    fitfuncs = []
    for bar in range(nbars):
        inv_vscinti = inv_vscint[bar]
        tdiff_offseti = tdiff_offset[bar]
        
        hist2d_list[bar].SetTitle(f"Vscint bar={bar}; dL - dR; tL - tR")
        hist2d_list[bar].Draw("COLZ")
        fitfunc = ROOT.TF1(f"vscint{bar}",
                           f"{tdiff_offseti} + {inv_vscinti}*x",
                           -hodo_bar_width/2,
                           hodo_bar_width/2
                           )
        fitfunc.SetLineColor(ROOT.kRed)
        fitfunc.SetLineWidth(2)
        fitfunc.Draw("SAME")
        fitfuncs.append(fitfunc)

        c.Print(pdf_filename)

    c.Print(f"{pdf_filename}]")

    ### Creating calibration output File ###
    txtout = "./outfiles/fitresults_hodo_vscint_calibration.txt"

    barid_element = np.arange(len(vscint))
    
    txttable = np.column_stack((barid_element, vscint))
    np.savetxt(txtout, txttable, fmt=["%-10d", "%-10.6f"])
    
    return vscint
# -------------------------------------------------------- #
### --- Testing the Function --- ###

def main() -> None:
    print(f"Testing {os.path.basename(__file__)} ...")

    filename = "hodo_timing_pass2_data_GEN3_He3.root"
    dirpath = "../../../outfiles/rootfiles/"
    filepath = dirpath + filename

    ToTcalibration = "fitresults_hodo_ToT_calibration.txt"
    internaloffset = "hodo_internal_alignment_calibration.txt"

    nbars = 90
    
    vscint = VscintCalibration(filepath,
                               nbars,
                               ToTcalibration=ToTcalibration,
                               internaloffset=internaloffset
                               )

if __name__ == "__main__":
    main()

    
    
    
    
    
    
    
    

    

