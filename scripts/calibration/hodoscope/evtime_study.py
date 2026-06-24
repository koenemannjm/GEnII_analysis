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

path = "../../../outfiles/rootfiles/"
# filename = f"QE_{kin_pass}_try8_v4_data_{kin}_sbs100p_nucleon_np.root"
filename = f"hodo_timing_{kin_pass}_try8_v4_data_{kin}_He3.root"
filepath = path + filename

# --- Constants --- #
c = 0.299792458 # units -> m/ns
tdiff_max = 10 # units -> ns

branch_name = "g.evtime"

with uproot.open(filepath) as ROOTFILE:

    tree = ROOTFILE["Tout"]

    Nevents = tree.num_entries

    gevtime = tree[branch_name].arrays(library="np")[branch_name]

sorted_gevtime = np.sort(gevtime)
dgevtime = np.diff(sorted_gevtime)
dgevtime = dgevtime[dgevtime > 0]

candidates = np.linspace(1.0, 20.0, 190001)
tol = 0.1

def score(T):
    r = np.mod(dgevtime, T)
    r = np.minimum(r, T - r)
    return np.mean(r < tol)

scores = np.array([score(T) for T in candidates])
T_hat = candidates[np.argmax(scores)]

print("Estimated clock period (ns):", T_hat)


# plt.hist(g_q, bins=int(np.sqrt(Nevents)), range=(-0.5,5))
# plt.xlabel("g.evtime[i] - g.evtime[j]")
# plt.ylabel("Counts")
# plt.savefig("plot_devtime.pdf")
