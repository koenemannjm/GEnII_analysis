import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

def ToTCalibration(filepath: str,
                   nbars: int,
                   vscintcalibration: str ="default_hodo_vscint_calibration.txt",
                   internaloffset: str ="default_hodo_internal_alignment_calibration.txt"
                   ):
    
    """
    This function corrects for the linear dependence of leading edge time vs time over threshold
    for the left and right pmts of each bar. How this is done is by a Chi^2 reduction model

    Chi^2 = sum(yi - (a*xi + b))

    where
    yi = leading edge time of left or right pmt of specific bar
    xi = time over threshold of left or right pmt of specific bar

    Grad_a(Chi^2) = 0
    Grad_b(Chi^2) = 0

    a = (avg(x*y) - avg(x) * avg(y)) / (avg(x^2) - avg(x)^2)
    b = avg(y) - a * avg(x)

    where avg(s) = sum(s) / N

    The loop adds to the sum of each average shown above (ex. sum(xy) += xi * yi)

    Inputs:
    'ToTCalibration(filepath, nbars):'

    filepath (string) -> Path to ROOT file
    nbars (integer) -> Number of bars in detector

    Outputs:
    ' return slope_left, intercept_left, slope_right, intercept_right '

    slope_left (Array | len = nbars) -> Slope of linear fit of left PMT per bar
    intercept_left (Array | len = nbars) -> Intercept of linear fit of left PMT per bar
    slope_right (Array | len = nbars) -> Slope of linear fit of right PMT per bar
    intercept_right Array | len = nbars) -> Intercept of linear fit of right PMT per bar

    """

    ### Constants ###
    speed_of_light = 0.299792458 # meters/ns
    
    hodo_bar_height = 0.025 # meters
    hodo_bar_width = 0.60 # meters
    z_hodo = 1.854454 # meters
    # z_hodo = 1.63 # meters

    BBdist = 1.63

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

    sumreftime = 0.0
    sumN = 0.0

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

    ### Unpacking Vscint calibration ###
    barid_element, vscint = np.loadtxt(
        f"./outfiles/{vscintcalibration}",
        unpack=True,
        dtype=float
    )

    ### Unpacking internal alignment calibration ###
    barid_element, offset = np.loadtxt(
        f"./outfiles/{internaloffset}",
        unpack=True,
        dtype=float
    )
    
    ### Initializing all summation terms for Chi^2 Optimizaiton for all bars ###
    ## Left PMT ##
    sum_tleft = np.zeros(nbars)
    sum_totleft = np.zeros(nbars)
    sum_totleft_sq = np.zeros(nbars)
    sum_tleft_totleft = np.zeros(nbars)
    sum_left_elements = np.zeros(nbars)

    ## Right PMT ##
    sum_tright = np.zeros(nbars)
    sum_totright = np.zeros(nbars)
    sum_totright_sq = np.zeros(nbars)
    sum_tright_totright = np.zeros(nbars)
    sum_right_elements = np.zeros(nbars)

    ### Forming histograms ###
    totmin, totmax = 0.0, 50.0
    tmin, tmax = -40, 40
    nbins= 100

    histleft2d_list = []
    graphleft_list = []
    for bar in range(nbars):
        histname = f"tleft_vs_totleft_bar{bar}"
        histtitle = f"tleft vs totleft bar={bar}"
        histleft2d = ROOT.TH2F(histname,
                               histtitle,
                               nbins,
                               totmin,
                               totmax,
                               nbins,
                               tmin,
                               tmax
                               )
        histleft2d_list.append(histleft2d)

        graph = ROOT.TGraph()
        graph.SetName(f"graphleft_bar{bar}")
        graph.SetTitle(f"tleft vs totleft bar={bar};ToT;Tleft")
        graphleft_list.append(graph)

    histright2d_list = []
    graphright_list = []
    for bar in range(nbars):
        histname = f"tright_vs_totright{bar}"
        histtitle = f"tright vs totright bar={bar}"
        histright2d = ROOT.TH2F(histname,
                                histtitle,
                                nbins,
                                totmin,
                                totmax,
                                nbins,
                                tmin,
                                tmax
                                )
        histright2d_list.append(histright2d)

        graph = ROOT.TGraph()
        graph.SetName(f"graphright_bar{bar}")
        graph.SetTitle(f"tright vs totright bar={bar};ToT;Tright")
        graphright_list.append(graph)

    ### The event loop ###
    print("Starting Event Loop ...")
    Nevents = chain.GetEntries()
    print(f"total events = {Nevents}")
    ievent = 0
    while reader.Next():
        
        trig = trigtime.Get()[0]
        evt = int(gevtime.Get()[0])
        
        y = yfp[0] + z_hodo * phfp[0]

        dleft = min(hodo_bar_width, max(0.0, hodo_bar_width/2 - y))
        dright = min(hodo_bar_width, max(0.0, hodo_bar_width/2 + y))

        if (trig<400):

            tref = trig - 4.0 * (6 % (2 + 6 % (evt))) - TBBtrig_t0

            for i in range(nhits.Get()[0]):
                if ((totleft[i]>0.0) &
                    (totright[i]>0.0) &
                    (abs(dleft-dright)<0.2)
                    ):

                    barid_i = int(barid[i])
                    # print(barid_i)

                    tleft_corr = tleft[i] + tref - dleft / vscint[barid_i] - offset[barid_i] - (etof[i] - etof0)
                    sum_tleft[barid_i] += tleft_corr
                    sum_totleft[barid_i] += totleft[i]
                    sum_totleft_sq[barid_i] += totleft[i] * totleft[i]
                    sum_tleft_totleft[barid_i] += tleft_corr * totleft[i]
                    sum_left_elements[barid_i] += 1

                    histleft2d_list[barid_i].Fill(totleft[i], tleft_corr)
                    # graphleft_list[barid_i].AddPoint(totleft[i], tleft_corr)

                    tright_corr = tright[i] + tref - dright / vscint[barid_i] - offset[barid_i] - (etof[i] - etof0)
                    sum_tright[barid_i] += tright_corr
                    sum_totright[barid_i] += totright[i]
                    sum_totright_sq[barid_i] += totright[i] * totright[i]
                    sum_tright_totright[barid_i] += tright_corr * totright[i]
                    sum_right_elements[barid_i] += 1

                    histright2d_list[barid_i].Fill(totright[i], tright_corr)
                    # graphright_list[barid_i].AddPoint(totright[i], tright_corr)

        ievent += 1

        if ievent%100000 == 0:
            print(f"Percent complete: {100.0 * ievent / Nevents:.2f}%")

    print("Percent complete: 100.00%")
    print("Event Loop Finished!")

    intercept_left = []
    slope_left = []
    intercept_right = []
    slope_right = []

    for i in range(nbars):
        """
        gL = graphleft_list[i]
        gR = graphright_list[i]
        
        fit_func_left = ROOT.TF1("fit_func_left", "[0] + [1]*x", totmin, totmax)

        fit_result_left = gL.Fit(fit_func_left, "SQ")
        intercept_lefti = fit_func_left.GetParameter(0)
        slope_lefti = fit_func_left.GetParameter(1)
        intercept_left.append(intercept_lefti)
        slope_left.append(slope_lefti)

        fit_func_right = ROOT.TF1("fit_func_right", "[0] + [1]*x", totmin, totmax)

        fit_result_righti = gR.Fit(fit_func_right, "SQ")
        intercept_righti = fit_func_right.GetParameter(0)
        slope_righti = fit_func_right.GetParameter(1)
        intercept_right.append(intercept_righti)
        slope_right.append(slope_righti)
        """
        
        if (sum_left_elements[i]>0 and sum_right_elements[i]>0):
            Nleft = sum_left_elements[i]
            average_tleft = sum_tleft[i] / Nleft
            average_totleft = sum_totleft[i] / Nleft
            average_totleft_sq = sum_totleft_sq[i] / Nleft
            average_tleft_totleft = sum_tleft_totleft[i] / Nleft
            if abs(average_totleft_sq - average_totleft * average_totleft)>0:
                slope_lefti = ((average_tleft_totleft - average_totleft * average_tleft)/
                               (average_totleft_sq - average_totleft * average_totleft))
                intercept_lefti = ((average_tleft * average_totleft_sq - average_totleft * average_tleft_totleft)/
                                   (average_totleft_sq - average_totleft * average_totleft))
            else:
                slope_lefti = 0.0
                intercept_lefti = average_tleft

            slope_left.append(slope_lefti)
            intercept_left.append(intercept_lefti)

            Nright = sum_right_elements[i]
            average_tright= sum_tright[i] / Nright
            average_totright = sum_totright[i] / Nright
            average_totright_sq = sum_totright_sq[i] / Nright
            average_tright_totright = sum_tright_totright[i] / Nright
            if abs(average_totright_sq - average_totright * average_totright)>0:
                slope_righti = ((average_tright_totright - average_totright * average_tright)/
                                (average_totright_sq - average_totright * average_totright))
                intercept_righti = ((average_tright * average_totright_sq - average_totright * average_tright_totright)/
                                   (average_totright_sq - average_totright * average_totright))
                
            else:
                slope_righti = 0.0
                intercept_righti = average_tright 

            slope_right.append(slope_righti)
            intercept_right.append(intercept_righti)

        else:
            slope_left.append(0.0)
            intercept_left.append(0.0)
            slope_right.append(0.0)
            intercept_right.append(0.0)
        

    slope_left = np.array(slope_left)
    intercept_left = np.array(intercept_left)
    slope_right = np.array(slope_right)
    intercept_right = np.array(intercept_right)

    ### Generating PDF output ###
    pdf_name = "plots_hodo_ToT_calibration.pdf"
    pdf_path = "./figures/"
    pdf_filename = pdf_path + pdf_name
    print(f"Creating PDF figure :{pdf_filename}")

    c = ROOT.TCanvas("c", "Time-Walk Correction", 1700, 2200)
    c.Print(f"{pdf_filename}[")
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    ROOT.gStyle.SetNumberContours(50)

    j = 1
    fitfuncs = []
    for bar in range(nbars):
        slopeL, interceptL = slope_left[bar], intercept_left[bar]
        slopeR, interceptR = slope_right[bar], intercept_right[bar]

        if j == 1:
            c.Clear()
            c.SetCanvasSize(1700, 2200)
            c.Divide(2,6)

        # Left
        c.cd(j)
        # ROOT.gPad.SetMargin(0.12, 0.15, 0.12, 0.08)
        histleft2d_list[bar].SetTitle(f"Time Walk left bar={bar}; ToT; Corrected Time")
        histleft2d_list[bar].Draw("COLZ")
        fitfuncL = ROOT.TF1(f"twL_bar{bar}", f"{interceptL} + {slopeL}*x", 0.0, 50.0)
        fitfuncL.SetLineColor(ROOT.kRed)
        fitfuncL.SetLineWidth(2)
        fitfuncL.Draw("SAME")
        fitfuncs.append(fitfuncL)

        # Right
        c.cd(j + 1)
        # ROOT.gPad.SetMargin(0.12, 0.15, 0.12, 0.08)
        histright2d_list[bar].SetTitle(f"Time Walk right bar={bar}; ToT; Corrected Time")
        histright2d_list[bar].Draw("COLZ")
        fitfuncR = ROOT.TF1(f"twR_bar{bar}", f"{interceptR} + {slopeR}*x", 0.0, 50.0)
        fitfuncR.SetLineColor(ROOT.kRed)
        fitfuncR.SetLineWidth(2)
        fitfuncR.Draw("SAME")
        fitfuncs.append(fitfuncR)

        j += 2

        if j > 12:
            c.Print(pdf_filename)
            j = 1

    if j !=1:
        c.Print(pdf_filename)

    c.Print(f"{pdf_filename}]")

    ### Creating calibration output File ###
    txtout = "./outfiles/fitresults_hodo_ToT_calibration.txt"

    print(f"Creating calibration table: {txtout}")

    barid_element = np.array(range(nbars))
    txttable = np.column_stack((barid_element,
                                slope_left,
                                intercept_left,
                                slope_right,
                                intercept_right)
                               )
    np.savetxt(txtout, txttable, fmt=["%-10d","%-10.6f","%-10.6f","%-10.6f","%-10.6f"])
            
    return slope_left, intercept_left, slope_right, intercept_right
# -------------------------------------------------------- #
### --- Testing the Function --- ###

def main() -> None:
    print(f"Testing {os.path.basename(__file__)} ...")
    
    filename = "hodo_timing_pass2_data_GEN3_He3.root"
    dirpath = "./../../../outfiles/rootfiles/"
    filepath = dirpath + filename

    vscintcalibration ="fitresults_hodo_vscint_calibration.txt"
    internaloffset = "hodo_internal_alignment_calibration.txt"

    nbars = 90

    wleft, vleft, wright, vright = ToTCalibration(filepath,
                                                  nbars,
                                                  vscintcalibration=vscintcalibration,
                                                  internaloffset=internaloffset
                                                  )

if __name__=="__main__":
    main()
    


    

    
