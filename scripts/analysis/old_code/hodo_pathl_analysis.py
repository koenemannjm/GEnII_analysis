import ROOT
import numpy as np
import matplotlib.pyplot as plt


#- contstants -#
c = 3*1e17 #m/ns

path  = "../../outfiles/rootfiles/"
filename  = "QE_pass2_try8_v4_data_GEN3_sbs100p_nucleon_np.root"

filepath = path + filename
tree = "Tout"

raw_data = ROOT.RDataFrame(tree, filepath)
branch_variables = ["bb.hodotdc.clus.tmean",
                    "bb.hodotdc.clus.bar.tdc.id",
                    "bb.etot_over_p",
                    "bb.tr.pathl",
                    "g.runnum"]
raw_data_numpy = raw_data.AsNumpy(branch_variables)
del raw_data

bb_hodo_tmean = raw_data_numpy[branch_variables[0]]
bb_hodo_bar_id = raw_data_numpy[branch_variables[1]]
bb_pathl = raw_data_numpy[branch_variables[2]]
bb_eop = raw_data_numpy[branch_variables[3]]
rununum = raw_data_numpy[branch_variables[4]]
del raw_data_numpy

# plt.hist2d(hodo_barid0, hodo_tmean0, bins=(100,100))
# plt.show()
