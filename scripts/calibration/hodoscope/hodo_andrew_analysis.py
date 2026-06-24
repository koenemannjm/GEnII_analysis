import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

def gaussian(x, A, mean, sigma, B):
    return A * np.exp(-0.5 * ((x - mean)/sigma)**2) + B

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

detector_time_variable = "bb.hodotdc.clus.bar.tdc.meantime" #vector
detector_blkid_variable = "bb.hodotdc.clus.bar.tdc.id" #vector

gevtime_variable = "g.evtime" #scalar
bb_sh_time_variable = "bb.sh.atimeblk" #scalar
hcal_time_variable = "sbs.hcal.atimeblk" # scalar
gem_trigtime_variable = "bb.gem.trigtime" #scalar
W2_variabale = "e.kine.W2" #scalar
ps_e_variable = "bb.ps.e" #scalar
bb_vz_variable= "bb.tr.vz" #vector
dx_variable = "dx" #scalar
dy_variable = "dy" #scalar

print("Reading in Root Variables")
with uproot.open(filepath) as ROOTFile:
    tree = ROOTFile["Tout"]

    # Nevents = tree.num_entries

    blktime = tree[detector_time_variable].arrays(library="ak")[detector_time_variable]
    blkid = tree[detector_blkid_variable].arrays(library="ak")[detector_blkid_variable]
    gevtime = tree[gevtime_variable].arrays(library="ak")[gevtime_variable]
    bb_sh_time = tree[bb_sh_time_variable].arrays(library="ak")[bb_sh_time_variable]
    hcal_time = tree[hcal_time_variable].arrays(library="ak")[hcal_time_variable]
    gem_trigtime = tree[gem_trigtime_variable].arrays(library="ak")[gem_trigtime_variable]
    W2 = tree[W2_variabale].arrays(library="ak")[W2_variabale]
    ps_e = tree[ps_e_variable].arrays(library="ak")[ps_e_variable]
    bb_vz = tree[bb_vz_variable].arrays(library="ak")[bb_vz_variable]
    dy = tree[dy_variable].arrays(library="ak")[dy_variable]
    

Nevents = len(gevtime)
    
mask = ak.num(blktime) > 0

blktime = blktime[mask]
blkid = blkid[mask]
gevtime = gevtime[mask].to_numpy()
bb_sh_time = bb_sh_time[mask].to_numpy()
hcal_time = hcal_time[mask].to_numpy()
gem_trigtime = gem_trigtime[mask].to_numpy()
W2 = W2[mask].to_numpy()
ps_e = ps_e[mask].to_numpy()
bb_vz = bb_vz[mask]
dy = dy[mask].to_numpy()

blktime_main_clus = ak.firsts(blktime).to_numpy()
bb_vz = ak.firsts(bb_vz).to_numpy()

blkid_main_clus = ak.firsts(blkid).to_list()
blkid_main_clus = np.array(blkid_main_clus)
hidden_none = ak.is_none(blkid_main_clus)
index_hidden_none = ak.where(hidden_none)[0].to_numpy()
print(index_hidden_none)

Nevents_mask = len(blktime_main_clus)
print(f"Events Lost Using Hodoscope: {Nevents_mask}/{Nevents} | {Nevents_mask/Nevents*100:.3f}%")

# savefig_name = f"plots_andrew_analysis_{kin_pass}_{kin}.pdf"


# doing QE Cuts Here:
savefig_name = f"plots_andrew_analysis_{kin_pass}_{kin}_QE.pdf"

mask = (W2 < 2.2) & (W2>-0.2) & (ps_e>0.2) & (np.abs(dy)<1.0) & (np.abs(bb_vz)<0.27)

blktime_main_clus = blktime_main_clus[mask]
blkid_main_clus = blkid_main_clus[mask]
blktime = blktime[mask]
blkid = blkid[mask]
gevtime = gevtime[mask]
bb_sh_time = bb_sh_time[mask]
hcal_time = hcal_time[mask]
gem_trigtime = gem_trigtime[mask]

hidden_none = ak.is_none(blkid_main_clus)
index_hidden_none = ak.where(hidden_none)[0].to_numpy()
print(index_hidden_none)

Nevents_QE = len(dx)

print(f"QE Events: {Nevents_QE}")



# --- Extracting Offsets --- #

offset_blkid, offset_blktime = np.loadtxt(offsetname, delimiter="\t", unpack=True)

# --- Applying Correction --- #
max_id = int(max(offset_blkid.max(), np.max(blkid_main_clus))) + 1
lookup_table = np.full(max_id, np.nan)
lookup_table[offset_blkid.astype(int)] = offset_blktime
offsets = lookup_table[blkid_main_clus.astype(int)]
print("Applying Correction")
blktime_main_clus_cor = blktime_main_clus - offsets


# --- Constructing Andrew's variable --- #

mod_evtime = np.mod(gevtime, 6)

bb_hodo_coin = blktime_main_clus + gem_trigtime - bb_sh_time
bb_hodo_coin_cor = blktime_main_clus_cor + gem_trigtime - bb_sh_time
bb_hodo_coin_stack = blktime_main_clus + gem_trigtime - bb_sh_time - 4.0*np.mod(mod_evtime + 2,6)
bb_hodo_coin_cor_stack = blktime_main_clus_cor + gem_trigtime - bb_sh_time - 4.0*np.mod(mod_evtime + 2,6)

hcal_hodo_coin = blktime_main_clus + gem_trigtime - hcal_time
hcal_hodo_coin_cor = blktime_main_clus_cor + gem_trigtime - hcal_time
hcal_hodo_coin_stack = blktime_main_clus + gem_trigtime - hcal_time - 4.0*np.mod(mod_evtime + 2,6)
hcal_hodo_coin_cor_stack = blktime_main_clus_cor + gem_trigtime - hcal_time - 4.0*np.mod(mod_evtime + 2,6)

print("Creating plots")

with PdfPages(savefig_name) as pdf:
    # -- 2D Histograms -- #
    """
    mask = ak.num(blktime) == 4
    blktime = blktime[mask]
    blkid = blkid[mask]

    event_index = 0
    fig, ax = plt.subplots(figsize=(14,2))
    
    for col in range(90):

        rect = Rectangle((col,0), 1, 0.2,
                         linewidth=1, edgecolor='black', facecolor='white')
        ax.add_patch(rect)

        if col in blkid[event_index]:
            rect.set_facecolor('red')
    ax.set_xlim(0,90)
    ax.set_ylim(0, 0.3)
    ax.set_aspect('equal')
    ax.invert_yaxis()  # Optional: so bar 0 is at bottom
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Detector Hits for Event")
    pdf.savefig(fig,dpi=300)
    plt.close(fig)
    
    """
    
    # - Uncorrected Hodoscope vs bar ID - #
    h = plt.hist2d(blkid_main_clus, blktime_main_clus, bins=(90, 200), range=[[-0.01,89.01],[-10,15]], cmap='rainbow')
    counts, xedges, yedges, img = h
    counts_masked = np.ma.masked_where(counts ==0, counts)
    plt.clf()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = plt.pcolormesh(X, Y, counts_masked.T, cmap='rainbow')
    plt.colorbar(mesh, label="Counts")
    plt.xlabel("hodo.id[0]")
    plt.ylabel("hodo.meantime[0] (ns)")
    plt.title("Hodoscope (Uncorrected) {kin}")
    pdf.savefig(dpi=300)  # Save this figure to the PDF
    plt.close()

    # - Corrected Hodoscope vs bar ID - #
    h = plt.hist2d(blkid_main_clus, blktime_main_clus_cor, bins=(90, 200), range=[[-0.01,89.01],[-10,15]], cmap='rainbow')
    counts, xedges, yedges, img = h
    counts_masked = np.ma.masked_where(counts ==0, counts)
    plt.clf()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = plt.pcolormesh(X, Y, counts_masked.T, cmap='rainbow')
    plt.colorbar(mesh, label="Counts")
    plt.xlabel("hodo.id[0]")
    plt.ylabel("hodo.meantime[0] - offset (ns)")
    plt.title("Hodoscope (Corrected w/ IA) {kin}")
    pdf.savefig(dpi=300)
    plt.close()

    # - Uncorrected Andrew Variable vs mod(evtime, 6) - #
    h = plt.hist2d(mod_evtime, bb_hodo_coin, bins=(6, 200), range=[[-0.01,5.01],[-20,20]], cmap='rainbow')
    counts, xedges, yedges, img = h
    counts_masked = np.ma.masked_where(counts ==0, counts)
    plt.clf()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = plt.pcolormesh(X, Y, counts_masked.T, cmap='rainbow')
    plt.colorbar(mesh, label="Counts")
    plt.xlabel("mod(g.evtime,6)")
    plt.ylabel("hodo.meantime[0] + bb.gem.trigtime - bb.sh.atimeblk (ns)")
    plt.title("Hodo-Shower Coin. vs mod(evtime,6) {kin}")
    pdf.savefig(dpi=300)
    plt.close()

    # - Corrected Andrew Variable vs mod(evtime, 6) - #
    h = plt.hist2d(mod_evtime, bb_hodo_coin_cor, bins=(6, 200), range=[[-0.01,5.01],[-20,20]], cmap='rainbow')
    counts, xedges, yedges, img = h
    counts_masked = np.ma.masked_where(counts ==0, counts)
    plt.clf()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = plt.pcolormesh(X, Y, counts_masked.T, cmap='rainbow')
    plt.colorbar(mesh, label="Counts")
    plt.xlabel("mod(g.evtime,6)")
    plt.ylabel("hodo.meantime[0] - offset + bb.gem.trigtime - bb.sh.atimeblk (ns)")
    plt.title("Hodo-Shower Coin. vs mod(evtime,6) (Corrected w/ IA) {kin}")
    pdf.savefig(dpi=300)
    plt.close()

    # - Uncorrected Andrew Variable vs bar ID - #
    h = plt.hist2d(blkid_main_clus, bb_hodo_coin, bins=(90, 200), range=[[-0.01,89.01],[-20,20]], cmap='rainbow')
    counts, xedges, yedges, img = h
    counts_masked = np.ma.masked_where(counts ==0, counts)
    plt.clf()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = plt.pcolormesh(X, Y, counts_masked.T, cmap='rainbow')
    plt.colorbar(mesh, label="Counts")
    plt.xlabel("hodo.id[0]")
    plt.ylabel("hodo.meantime[0] + bb.gem.trigtime - bb.sh.atimeblk (ns)")
    plt.title("Hodo-Shower Coin. vs bar id{kin}")
    pdf.savefig(dpi=300)
    plt.close()

    # - Corrected Andrew Variable vs vs bar ID - #
    h = plt.hist2d(blkid_main_clus, bb_hodo_coin_cor, bins=(90, 200), range=[[-0.01,89.01],[-20,20]], cmap='rainbow')
    counts, xedges, yedges, img = h
    counts_masked = np.ma.masked_where(counts ==0, counts)
    plt.clf()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = plt.pcolormesh(X, Y, counts_masked.T, cmap='rainbow')
    plt.colorbar(mesh, label="Counts")
    plt.xlabel("hodo.id[0]")
    plt.ylabel("hodo.meantime[0] - offset + bb.gem.trigtime - bb.sh.atimeblk (ns)")
    plt.title("Hodo-Shower Coin. vs bar id (Corrected w/ IA) {kin}")
    pdf.savefig(dpi=300)
    plt.close()

    # - Uncorrected Andrew Variable vs bar ID - #
    h = plt.hist2d(blkid_main_clus, bb_hodo_coin_stack, bins=(90, 200), range=[[-0.01,89.01],[-20,20]], cmap='rainbow')
    counts, xedges, yedges, img = h
    counts_masked = np.ma.masked_where(counts ==0, counts)
    plt.clf()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = plt.pcolormesh(X, Y, counts_masked.T, cmap='rainbow')
    plt.colorbar(mesh, label="Counts")

    x_centers = 0.5*(xedges[:-1] + xedges[1:])
    y_centers = 0.5*(yedges[:-1] + yedges[1:])

    maxy = []
    for i in range(len(x_centers)):
        hy = counts[i, :]
        if np.all(hy == 0):
            maxy.append(np.nan)
        else:
            j = np.argmax(hy)
            maxy.append(y_centers[j])
    maxy = np.array(maxy)
    plt.plot(x_centers, maxy, 'r.', label='bar peak')
    plt.xlabel("hodo.id[0]")
    plt.ylabel("hodo.meantime[0] + bb.gem.trigtime - bb.sh.atimeblk - 4((g.evtime%6+2))%6 (ns)")
    plt.title("Hodo-Shower Coin. vs bar id{kin}")
    pdf.savefig(dpi=300)
    plt.close()

    # - Corrected Andrew Variable vs vs bar ID - #
    h = plt.hist2d(blkid_main_clus, bb_hodo_coin_cor_stack, bins=(90, 200), range=[[-0.01,89.01],[-20,20]], cmap='rainbow')
    counts, xedges, yedges, img = h
    counts_masked = np.ma.masked_where(counts ==0, counts)
    plt.clf()
    X, Y = np.meshgrid(xedges, yedges)
    mesh = plt.pcolormesh(X, Y, counts_masked.T, cmap='rainbow')
    plt.colorbar(mesh, label="Counts")

    x_centers = 0.5*(xedges[:-1] + xedges[1:])
    y_centers = 0.5*(yedges[:-1] + yedges[1:])

    maxy = []
    for i in range(len(x_centers)):
        hy = counts[i, :]
        if np.all(hy == 0):
            maxy.append(np.nan)
        else:
            j = np.argmax(hy)
            maxy.append(y_centers[j])
    maxy = np.array(maxy)
    plt.plot(x_centers, maxy, 'r.', label='bar peak')
    plt.xlabel("hodo.id[0]")
    plt.ylabel("hodo.meantime[0] - offset + bb.gem.trigtime - bb.sh.atimeblk - 4((g.evtime%6+2))%6 (ns)")
    plt.title("Hodo-Shower Coin. vs bar id (Corrected w/ IA) {kin}")
    pdf.savefig(dpi=300)
    plt.close()

    # -- 1D Histograms -- #

    # - Corrected and Uncorrected Hodoscope Time - #
    plt.hist(blktime_main_clus, bins=200, range=[-20,20], alpha=0.5, label="Uncorrected")
    plt.hist(blktime_main_clus_cor, bins=200, range=[-20,20], alpha=0.5, label="Corrected")
    plt.xlabel("Hodo Time")
    plt.ylabel("Counts")
    plt.title("Hodo")
    plt.legend()
    pdf.savefig(dpi=300)
    plt.close()

    # - Uncorrected Andrew variable Time - #
    plt.hist(bb_hodo_coin, bins=200, range=[-26,20], alpha=0.5, label="Uncorrected")
    plt.hist(bb_hodo_coin_stack, bins=200, range=[-26,20], alpha=0.5, label="Corrected")
    plt.xlabel("time (ns)")
    plt.ylabel("Counts")
    plt.title("Hodo - BBCAL")
    plt.legend()
    pdf.savefig(dpi=300)
    plt.close()

    # - Corrected Andrew variable Time - #
    plt.hist(bb_hodo_coin_cor, bins=200, range=[-20,20], alpha=0.5, label="Uncorrected")
    plt.hist(bb_hodo_coin_cor_stack, bins=200, range=[-20,20], alpha=0.5, label="Corrected")
    plt.xlabel("time (ns)")
    plt.ylabel("Counts")
    plt.title("Hodo - BBCAL")
    plt.legend()
    pdf.savefig(dpi=300)
    plt.close()

    # - Corrected and Uncorrected Hodoscope Time - #
    plt.hist(bb_hodo_coin, bins=200, range=[-20,20], alpha=0.5, label="Uncorrected")
    plt.hist(bb_hodo_coin_cor, bins=200, range=[-20,20], alpha=0.5, label="Corrected")
    plt.xlabel("time (ns)")
    plt.ylabel("Counts")
    plt.title("Hodo - BBCAL")
    plt.legend()
    pdf.savefig(dpi=300)
    plt.close()

    # - Corrected and Uncorrected Hodoscope Time - #
    plt.hist(bb_hodo_coin_stack, bins=200, range=[-20,-2], alpha=0.5, label="Uncorrected")
    counts, bin_edges, _ = plt.hist(bb_hodo_coin_cor_stack, bins=200, range=[-20,-2], alpha=0.5, label="Corrected")
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    peak_idx = np.argmax(counts)
    
    A_guess = counts.max()
    mean_guess = bin_centers[peak_idx]
    sigma_guess = 2
    B_guess = counts[0]
    p0 = [A_guess, mean_guess, sigma_guess, B_guess]

    popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=p0)
    perr = np.sqrt(np.diagonal(pcov))
    plt.plot(bin_centers, gaussian(bin_centers, *popt), color='red', lw=2, label=f'Fit: mean={popt[1]:.2f}±{perr[1]:.2f}, sigma={popt[2]:.2f}±{perr[2]:.2f}')
    plt.xlabel("time (ns)")
    plt.ylabel("Counts")
    plt.title("Hodo - BBCAL")
    plt.legend()
    pdf.savefig(dpi=300)
    plt.close()

    # - Uncorrected Andrew variable Time - #
    plt.hist(hcal_hodo_coin, bins=200, range=[-180,-100], alpha=0.5, label="Uncorrected")
    plt.hist(hcal_hodo_coin_stack, bins=200, range=[-180,-100], alpha=0.5, label="Corrected")
    plt.xlabel("time (ns)")
    plt.ylabel("Counts")
    plt.title("Hodo - HCAL")
    plt.legend()
    pdf.savefig(dpi=300)
    plt.close()

    # - Corrected Andrew variable Time - #
    plt.hist(hcal_hodo_coin_cor, bins=200, range=[-180,-100], alpha=0.5, label="Uncorrected")
    plt.hist(hcal_hodo_coin_cor_stack, bins=200, range=[-180,-100], alpha=0.5, label="Corrected")
    plt.xlabel("time (ns)")
    plt.ylabel("Counts")
    plt.title("Hodo - HCAL")
    plt.legend()
    pdf.savefig(dpi=300)
    plt.close()

    # - Corrected and Uncorrected Hodoscope Time - #
    plt.hist(hcal_hodo_coin, bins=200, range=[-180,-100], alpha=0.5, label="Uncorrected")
    plt.hist(hcal_hodo_coin_cor, bins=200, range=[-180,-100], alpha=0.5, label="Corrected")
    plt.xlabel("time (ns)")
    plt.ylabel("Counts")
    plt.title("Hodo - HCAL")
    plt.legend()
    pdf.savefig(dpi=300)
    plt.close()

    # - Corrected and Uncorrected Hodoscope Time - #
    counts, bin_edges, _ = plt.hist(hcal_hodo_coin_stack, bins=200, range=[-180,-100], alpha=0.5, label="Uncorrected")
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    peak_idx = np.argmax(counts)
    
    A_guess = counts.max()
    mean_guess = bin_centers[peak_idx]
    sigma_guess = 2
    B_guess = counts[0]
    p0 = [A_guess, mean_guess, sigma_guess, B_guess]

    popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=p0)
    perr = np.sqrt(np.diagonal(pcov))
    plt.plot(bin_centers, gaussian(bin_centers, *popt), color='blue', lw=2, label=f'Fit UC')
    print("Fit: Uncorrected Hodo-Hcal")
    for i in range(len(popt)):
        print(f"{popt[i]:.2f}±{perr[i]:.2f}")
    counts, bin_edges, _ = plt.hist(hcal_hodo_coin_cor_stack, bins=200, range=[-180,-100], alpha=0.5, label="Corrected")
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    peak_idx = np.argmax(counts)
    
    A_guess = counts.max()
    mean_guess = bin_centers[peak_idx]
    sigma_guess = 2
    B_guess = counts[0]
    p0 = [A_guess, mean_guess, sigma_guess, B_guess]

    popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=p0)
    perr = np.sqrt(np.diagonal(pcov))
    plt.plot(bin_centers, gaussian(bin_centers, *popt), color='red', lw=2, label=f'Fit C')
    print("Fit: Corrected Hodo-Hcal")
    for i in range(len(popt)):
        print(f"{popt[i]:.2f}±{perr[i]:.2f}")
    plt.xlabel("time (ns)")
    plt.ylabel("Counts")
    plt.title("Hodo - HCAL")
    plt.legend()
    pdf.savefig(dpi=300)
    plt.close()

print("All Done :)")
