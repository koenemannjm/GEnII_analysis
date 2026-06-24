//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
//   Created by Jacob Koenemann & ChatGPT
//   contact: bxy3zr@virginia.edu
//                                                                     
//   Last Modified July 28, 2025   
//
//      
//   The purpose of this script is to trim down on the         
//   raw files for GEnII experiment and produce a
//   single root file that can then be analyzed.
//                                                    
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////     
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TRegexp.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TChain.h"

#include <set>
#include <map>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <regex>
#include <filesystem>

// Need for reading .cfg file
#include <iostream>
#include <cstdlib>  // for std::exit and EXIT_FAILURE
#include <fstream>

// Function for finding root files with specified run numbers
namespace fs = std::filesystem;

std::vector<std::string> getMatchingFiles(const std::string& directory, const std::vector<int>& selectedIDs) {
    std::vector<std::string> matchingFiles;

    // Loop over selected IDs and create regex patterns for each
    for (int id : selectedIDs) {
        std::string pattern = ".*_" + std::to_string(id) + "_.*\\.root";  // _####_ pattern

        // Regular expression for matching
        std::regex regexPattern(pattern);

        // Iterate through files in the directory
        for (const auto& entry : fs::directory_iterator(directory)) {
            std::string filename = entry.path().filename().string();

            // If filename matches the pattern, add it to the list
            if (std::regex_match(filename, regexPattern)) {
                matchingFiles.push_back(filename);
            }
        }
    }

    return matchingFiles;
}

void raw_try8_v4_trimming(std::string config_file) {


    // Reading in config file
    std::string config_path = "../../config/" + config_file;
    std::ifstream cfg_file(config_path);
    if (!cfg_file) {
       std::cerr << "Error: Config file does not exist at " << config_path << std::endl;
       std::exit(EXIT_FAILURE);
    }
    
    std::ifstream cfg(config_path);
    std::string line;
    std::map<std::string, std::string> settings;

    while (std::getline(cfg, line)) {
        if (line.empty() || line[0] == '#') continue;
	size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;
        std::string key = line.substr(0, eq_pos);
        std::string value = line.substr(eq_pos + 1);
        settings[key] = value;
    }

    // Pulling in Experiment naming
    std::string config = settings["config"];
    std::string exp_name = settings["exp_name"];
    std::string pass = settings["pass"];
    std::string target = settings["target"];

    // Pulling in Experiment kinematic parameters
    double ebeam = std::stod(settings["ebeam"]);
    double HCAL_angle_deg = std::stod(settings["hcal_angle"]);
    double HCAL_angle = HCAL_angle_deg * TMath::DegToRad();
    double HCAL_distance = std::stod(settings["hcal_distance"]);

    // Pulling selected run numbers
    std::istringstream iss(settings["selected_numbers"]);
    std::vector<int> selected_numbers;
    int val;
    while (iss >> val) selected_numbers.push_back(val);
    
    // Some constant(s)
    double Mp = 0.938272088; // PDG
    double Mn = 0.939565420; // PDG
    double MN = 0.5*(Mp + Mn);

    TVector3 HCAL_vector(-HCAL_distance * TMath::Sin(HCAL_angle), 0.0, HCAL_distance * TMath::Cos(HCAL_angle));
    TVector3 HCAL_unitvector_z(-TMath::Sin(HCAL_angle), 0.0, TMath::Cos(HCAL_angle));
    TVector3 HCAL_unitvector_x(0, -1.0, 0);
    TVector3 HCAL_unitvector_y = HCAL_unitvector_z.Cross(HCAL_unitvector_x);
   
   // Directory containing the files
    std::string input_dir = "/v/lustre24/expphy/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/" + pass + "/TEST/try8/" + exp_name + "/rootfiles/";  // cache directory for Raw data
    // /v/lustre24/expphy/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/pass2/TEST/try7/GEN2/rootfiles

    TSystemDirectory dir("input", input_dir.c_str());
    TList *files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "No files found in directory: " << input_dir << std::endl;
        return;
    }

    std::vector<std::string>matchingFiles = getMatchingFiles(input_dir, selected_numbers);

    //std::string output_dir = "/lustre24/expphy/volatile/halla/sbs/koeneman/data/raw/" + pass + "/test/try8/" + config + "_" + target + "/";
    std::string output_dir = "../../outfiles/rootfiles/data/";
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    std::string output_filename = output_dir + "He3" + exp_name + pass + ".root";

    // Create an output ROOT file
    TFile *outputFile = new TFile(output_filename.c_str(), "RECREATE");
    TTree *outputTree = new TTree("Tout", "Merged data from selected ROOT files");

    const int kMaxArraySize = 100;

    // Defining the output branch names
    // ----------------------------------------------------------------
    
    // vector bb.tr branches
    int bb_tr_n;
    outputTree->Branch("bb.tr.n", &bb_tr_n, "bb.tr.n/I");
    double bb_tr_p[kMaxArraySize];
    double bb_tr_px[kMaxArraySize];
    double bb_tr_py[kMaxArraySize];
    double bb_tr_pz[kMaxArraySize];
    double bb_tr_vx[kMaxArraySize];
    double bb_tr_vy[kMaxArraySize];
    double bb_tr_vz[kMaxArraySize];
    double bb_etot_over_p[kMaxArraySize];
    double bb_gem_track_nhits[kMaxArraySize];
    double bb_tr_pathl[kMaxArraySize];
    outputTree->Branch("bb.tr.p", bb_tr_p, "bb.tr.p[bb.tr.n]/D");
    outputTree->Branch("bb.tr.px", bb_tr_px, "bb.tr.px[bb.tr.n]/D");
    outputTree->Branch("bb.tr.py", bb_tr_py, "bb.tr.py[bb.tr.n]/D");
    outputTree->Branch("bb.tr.pz", bb_tr_pz, "bb.tr.pz[bb.tr.n]/D");
    outputTree->Branch("bb.tr.vx", bb_tr_vx, "bb.tr.vx[bb.tr.n]/D");
    outputTree->Branch("bb.tr.vy", bb_tr_vy, "bb.tr.vy[bb.tr.n]/D");
    outputTree->Branch("bb.tr.vz", bb_tr_vz, "bb.tr.vz[bb.tr.n]/D");
    outputTree->Branch("bb.tr.pathl",bb_tr_pathl,"bb.tr.pathl[bb.tr.n]/D");
    outputTree->Branch("bb.etot_over_p", bb_etot_over_p, "bb.etot_over_p[bb.tr.n]/D");
    outputTree->Branch("bb.gem.track.nhits", bb_gem_track_nhits, "bb.gem.track.nhits[bb.tr.n]/D");

    // vector other bb. branches
    int Ndata_bb_hodotdc_clus_bar_tdc_meantime;
    int Ndata_bb_hodotdc_clus_tmean;
    outputTree->Branch("Ndata.bb.hodotdc.clus.bar.tdc.meantime", &Ndata_bb_hodotdc_clus_bar_tdc_meantime, "Ndata.bb.hodotdc.clus.bar.tdc.meantime/I");
    outputTree->Branch("Ndata.bb.hodotdc.clus.tmean", &Ndata_bb_hodotdc_clus_tmean, "Ndata.bb.hodotdc.clus.tmean/I");
    double bb_tdctrig_tdcelemID;
    double bb_tdctrig_tdc;
    double bb_hodotdc_clus_bar_tdc_meantime[kMaxArraySize];
    double bb_hodotdc_clus_bar_tdc_id[kMaxArraySize];
    double bb_hodotdc_clus_tmean[kMaxArraySize];
    outputTree->Branch("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID, "bb.tdctrig.tdcelemID/D");
    outputTree->Branch("bb.tdctrig.tdc", &bb_tdctrig_tdc, "bb.tdctrig.tdc/D");
    outputTree->Branch("bb.hodotdc.clus.bar.tdc.meantime", bb_hodotdc_clus_bar_tdc_meantime, "bb.hodotdc.clus.bar.tdc.meantime[Ndata.bb.hodotdc.clus.bar.tdc.meantime]/D");
    outputTree->Branch("bb.hodotdc.clus.bar.tdc.id", bb_hodotdc_clus_bar_tdc_id, "bb.hodotdc.clus.bar.tdc.id[Ndata.bb.hodotdc.clus.bar.tdc.meantime]/D");
    outputTree->Branch("bb.hodotdc.clus.tmean", bb_hodotdc_clus_tmean, "bb.hodotdc.clus.tmean[Ndata.bb.hodotdc.clus.tmean]/D");

    // vector sbs.hcal branches
    int  sbs_hcal_nclus;
    outputTree->Branch("sbs.hcal.nclus", &sbs_hcal_nclus, "sbs.hcal.nclus/I");
    double sbs_hcal_clus_blk_tdctime[kMaxArraySize];
    double sbs_hcal_clus_blk_tdctime_tw[kMaxArraySize];
    double sbs_hcal_clus_blk_id[kMaxArraySize];
    double sbs_hcal_clus_x[kMaxArraySize];
    double sbs_hcal_clus_y[kMaxArraySize];
    double sbs_hcal_clus_e[kMaxArraySize];
    outputTree->Branch("sbs.hcal.clus_blk.tdctime", sbs_hcal_clus_blk_tdctime, "sbs.hcal.clus_blk.tdctime[sbs.hcal.nclus]/D");
    outputTree->Branch("sbs.hcal.clus_blk.tdctime_tw", sbs_hcal_clus_blk_tdctime_tw, "sbs.hcal.clus_blk.tdctime_tw[sbs.hcal.nclus]/D");
    outputTree->Branch("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id, "sbs.hcal.clus_blk.id[sbs.hcal.nclus]/D");
    outputTree->Branch("sbs.hcal.clus.x", sbs_hcal_clus_x, "sbs.hcal.clus.x[sbs.hcal.nclus]/D");
    outputTree->Branch("sbs.hcal.clus.y", sbs_hcal_clus_y, "sbs.hcal.clus.y[sbs.hcal.nclus]/D");
    outputTree->Branch("sbs.hcal.clus.e", sbs_hcal_clus_e, "sbs.hcal.clus.e[sbs.hcal.nclus]/D");

    // scalar sbs. branches
    double sbs_hcal_e;
    double sbs_hcal_x;
    double sbs_hcal_y;
    double sbs_hcal_rowblk;
    double sbs_hcal_colblk;
    double sbs_hcal_idblk;
    double sbs_hcal_tdctimeblk;
    double sbs_tdctrig_rftime;
    double sbs_hcal_atimeblk;
    outputTree->Branch("sbs.hcal.e", &sbs_hcal_e, "sbs.hcal.e/D");
    outputTree->Branch("sbs.hcal.x", &sbs_hcal_x, "sbs.hcal.x/D");
    outputTree->Branch("sbs.hcal.y", &sbs_hcal_y, "sbs.hcal.y/D");
    outputTree->Branch("sbs.hcal.rowblk", &sbs_hcal_rowblk, "sbs.hcal.rowblk/D");
    outputTree->Branch("sbs.hcal.colblk", &sbs_hcal_colblk, "sbs.hcal.colblk/D");
    outputTree->Branch("sbs.hcal.idblk", &sbs_hcal_idblk, "sbs.hcal.idblk/D");
    outputTree->Branch("sbs.hcal.tdctimeblk", &sbs_hcal_tdctimeblk, "sbs.hcal.tdctimeblk/D");
    outputTree->Branch("sbs.hcal.atimeblk", &sbs_hcal_atimeblk, "sbs.hcal.atimeblk/D");
    outputTree->Branch("sbs.tdctrig.rftime", &sbs_tdctrig_rftime, "sbs.tdctrig.rftime/D");
    outputTree->Branch("sbs.hcal.atimeblk", &sbs_hcal_atimeblk, "sbs.hcal.atimeblk/D");

    // scalar bb. branches
    double bb_grinch_tdc_clus_trackindex;
    double bb_grinch_tdc_clus_size;
    double bb_sh_e;
    double bb_ps_e;
    double bb_tdctrig_rftime;
    double bb_ps_atimeblk;
    double bb_sh_atimeblk;
    double bb_gem_trigtime;
    outputTree->Branch("bb.grinch_tdc.clus.trackindex", &bb_grinch_tdc_clus_trackindex, "bb.grinch_tdc.clus.trackindex/D");
    outputTree->Branch("bb.grinch_tdc.clus.size", &bb_grinch_tdc_clus_size, "bb.grind_tdc.clus.size/D");
    outputTree->Branch("bb.sh.e", &bb_sh_e, "bb.sh.e/D");
    outputTree->Branch("bb.sh.atimeblk", &bb_sh_atimeblk, "bb.sh.atimeblk/D");
    outputTree->Branch("bb.ps.e", &bb_ps_e, "bb.ps.e/D");
    outputTree->Branch("bb.tdctrig.rftime", &bb_tdctrig_rftime, "bb.tdctrig.rftime/D");
    outputTree->Branch("bb.ps.atimeblk", &bb_ps_atimeblk, "bb.ps.atimeblk/D");
    outputTree->Branch("bb.sh.atimeblk", &bb_sh_atimeblk, "bb.sh.atimeblk/D");
    outputTree->Branch("bb.gem.trigtime", &bb_gem_trigtime, "bb.gem.trigtime/D");

    // computed variables e.kine branches
    double e_kine_W2;
    double e_kine_Q2;
    double e_beam;
    double e_ppara_mag;
    double e_pperp_mag;
    double e_missing_mass2;
    outputTree->Branch("e.kine.W2", &e_kine_W2, "e.kine.W2/D");
    outputTree->Branch("e.kine.Q2", &e_kine_Q2, "e.kine.Q2/D");
    outputTree->Branch("e.beam", &e_beam, "e.beam/D");
    outputTree->Branch("e.ppara.mag", &e_ppara_mag, "e.ppara.mag/D");
    outputTree->Branch("e.pperp.mag", &e_pperp_mag, "e.pperp.mag/D");
    outputTree->Branch("e.missing_mass2", &e_missing_mass2, "e.missing_mass2/D");

    // other variables
    double g_trigbits;
    double HALLA_p;
    double hac_bcm_average;
    double g_evtime;
    double scalhel_hel;
    double IHWP;
    double g_runnum;
    outputTree->Branch("g.trigbits", &g_trigbits, "g.trigbits/D");
    outputTree->Branch("HALLA_p", &HALLA_p, "HALLA_p/D");
    outputTree->Branch("hac_bcm_average", &hac_bcm_average, "hac_bcm_average/D");
    outputTree->Branch("g.evtime", &g_evtime, "g.evtime/D");
    outputTree->Branch("scalhel.hel", &scalhel_hel, "scalhel.hel/D");
    outputTree->Branch("g.runnum", &g_runnum, "g.runnum/D");
    outputTree->Branch("IHWP", &IHWP, "IHWP/D");

    // dx and dy variables
    double dx;
    double dy;
    double sbs_hcal_x_exp;
    double sbs_hcal_y_exp;
    double dx_clus[kMaxArraySize];
    double dy_clus[kMaxArraySize];
    outputTree->Branch("sbs.hcal.x_exp", &sbs_hcal_x_exp, "sbs.hcal.x_exp/D");
    outputTree->Branch("sbs.hcal.y_exp", &sbs_hcal_y_exp, "sbs.hcal.y_exp/D");
    outputTree->Branch("dx", &dx, "dx/D");
    outputTree->Branch("dy", &dy, "dy/D");
    outputTree->Branch("dx_clus", dx_clus, "dx_clus[sbs.hcal.nclus]/D");
    outputTree->Branch("dy_clus", dy_clus, "dy_clus[sbs.hcal.nclus]/D");
    
    // ----------------------------------------------------------------

    // Creating TChain to have all "T"'s from every file
    TChain *C = new TChain("T");

    Long64_t eventTracker = 0;
    std::cout << "Merging ROOT files into TChain" << std::endl;
    int totalFiles = matchingFiles.size();
    const std::regex fourDigits{R"(_(\d{4})_)"};
    for (int j =0; j < totalFiles; ++j) {
        std::string filename = matchingFiles[j];
        std::string filePath = input_dir + filename;
	TFile *file = TFile::Open(filePath.c_str(), "READ");
        if (!file || file->IsZombie()) {
	    std::cerr << "Error opening file: " << filePath << std::endl;
	    std::exit(EXIT_FAILURE);
	}
	
	TTree* inputTree = (TTree*)file->Get("T");
	if (!inputTree) {
	    std::cerr << "Tree 'T' not found in " << filePath << std::endl;
	    file->Close();
	    delete file;
	    continue;
	}
	
	double percent_files = (1.0 * (j + 1)) / totalFiles * 100.0;
	std::cout << "Progress: " << std::fixed << std::setprecision(2) << percent_files << "% complete \r";
	std::cout.flush();
	
        C->Add(filePath.c_str());
	
	file->Close();
	delete file;
    }
    std::cout<< "\nCalculating total events\n" << "...\n" << std::endl;

    Long64_t totalEvents;
    totalEvents = C->GetEntries();
    
    std::cout<< "\nTotal number events = " << totalEvents << "\n" << std::endl;     

    C->SetBranchStatus("*",0);

    C->SetBranchStatus("bb.tr.p", 1);
    C->SetBranchStatus("bb.tr.px", 1);
    C->SetBranchStatus("bb.tr.py", 1);
    C->SetBranchStatus("bb.tr.pz", 1);
    C->SetBranchStatus("bb.tr.vx", 1);
    C->SetBranchStatus("bb.tr.vy", 1);
    C->SetBranchStatus("bb.tr.vz", 1);
    C->SetBranchStatus("bb.etot_over_p", 1);
    C->SetBranchStatus("Ndata.bb.hodotdc.clus.bar.tdc.meantime", 1);
    C->SetBranchStatus("bb.tdctrig.tdcelemID", 1);
    C->SetBranchStatus("bb.tdctrig.tdc", 1);
    C->SetBranchStatus("bb.hodotdc.clus.bar.tdc.meantime", 1);
    C->SetBranchStatus("bb.hodotdc.clus.bar.tdc.id", 1);
    C->SetBranchStatus("bb.gem.track.nhits", 1);
    C->SetBranchStatus("sbs.hcal.clus_blk.tdctime", 1);
    C->SetBranchStatus("sbs.hcal.clus_blk.tdctime_tw", 1);
    C->SetBranchStatus("sbs.hcal.clus_blk.id", 1);
    C->SetBranchStatus("sbs.hcal.nclus", 1);
    C->SetBranchStatus("sbs.hcal.clus_blk.tdctime_tw", 1);
    C->SetBranchStatus("sbs.hcal.clus_blk.tdctime_tw", 1);
    C->SetBranchStatus("sbs.hcal.e", 1);
    C->SetBranchStatus("sbs.hcal.x", 1);
    C->SetBranchStatus("sbs.hcal.y", 1);
    C->SetBranchStatus("sbs.hcal.rowblk", 1);
    C->SetBranchStatus("sbs.hcal.colblk", 1);
    C->SetBranchStatus("sbs.hcal.idblk", 1);
    C->SetBranchStatus("sbs.hcal.tdctimeblk", 1);
    C->SetBranchStatus("sbs.hcal.atimeblk", 1);
    C->SetBranchStatus("sbs.tdctrig.rftime", 1);
    C->SetBranchStatus("bb.grinch_tdc.clus.trackindex", 1);
    C->SetBranchStatus("bb.grinch_tdc.clus.size", 1);
    C->SetBranchStatus("bb.sh.e", 1);
    C->SetBranchStatus("bb.sh.atimeblk", 1);
    C->SetBranchStatus("bb.ps.e", 1);
    C->SetBranchStatus("bb.tr.n", 1);
    C->SetBranchStatus("bb.tdctrig.rftime", 1);
    C->SetBranchStatus("e.kine.W2", 1);
    C->SetBranchStatus("e.kine.Q2", 1);
    C->SetBranchStatus("g.trigbits", 1);
    C->SetBranchStatus("HALLA_p", 1);
    C->SetBranchStatus("hac_bcm_average", 1);
    C->SetBranchStatus("g.evtime", 1);
    C->SetBranchStatus("bb.tr.pathl", 1);
    C->SetBranchStatus("sbs.hcal.atimeblk", 1);
    C->SetBranchStatus("bb.ps.atimeblk", 1);
    C->SetBranchStatus("bb.sh.atimeblk", 1);
    C->SetBranchStatus("scalhel.hel", 1);
    C->SetBranchStatus("sbs.hcal.clus.x", 1);
    C->SetBranchStatus("sbs.hcal.clus.y", 1);
    C->SetBranchStatus("sbs.hcal.clus.e", 1);
    C->SetBranchStatus("bb.gem.trigtime", 1);
    C->SetBranchStatus("Ndata.bb.hodotdc.clus.tmean", 1);
    C->SetBranchStatus("bb.hodotdc.clus.tmean", 1);
    C->SetBranchStatus("IGL1I00OD16_16", 1);
    C->SetBranchStatus("g.runnum", 1);

    double bb_tr_p_in[kMaxArraySize];
    double bb_tr_px_in[kMaxArraySize];
    double bb_tr_py_in[kMaxArraySize];
    double bb_tr_pz_in[kMaxArraySize];
    double bb_tr_vx_in[kMaxArraySize];
    double bb_tr_vy_in[kMaxArraySize];
    double bb_tr_vz_in[kMaxArraySize];
    double bb_etot_over_p_in[kMaxArraySize];
    double bb_tr_pathl_in[kMaxArraySize];
    C->SetBranchAddress("bb.tr.p", bb_tr_p_in);
    C->SetBranchAddress("bb.tr.px", bb_tr_px_in);
    C->SetBranchAddress("bb.tr.py", bb_tr_py_in);
    C->SetBranchAddress("bb.tr.pz", bb_tr_pz_in);
    C->SetBranchAddress("bb.tr.vx", bb_tr_vx_in);
    C->SetBranchAddress("bb.tr.vy", bb_tr_vy_in);
    C->SetBranchAddress("bb.tr.vz", bb_tr_vz_in);
    C->SetBranchAddress("bb.etot_over_p", bb_etot_over_p_in);
    C->SetBranchAddress("bb.tr.pathl", bb_tr_pathl_in);

    int Ndata_bb_hodotdc_clus_bar_tdc_meantime_in;
    int Ndata_bb_hodotdc_clus_tmean_in;
    C->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.meantime", &Ndata_bb_hodotdc_clus_bar_tdc_meantime_in);
    C->SetBranchAddress("Ndata.bb.hodotdc.clus.tmean", &Ndata_bb_hodotdc_clus_tmean_in);
    double bb_tdctrig_tdcelemID_in;
    double bb_tdctrig_tdc_in;
    double bb_hodotdc_clus_bar_tdc_meantime_in[kMaxArraySize];
    double bb_hodotdc_clus_bar_tdc_id_in[kMaxArraySize];
    double bb_gem_track_nhits_in[kMaxArraySize];
    double bb_hodotdc_clus_tmean_in[kMaxArraySize];
    C->SetBranchAddress("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID_in);
    C->SetBranchAddress("bb.tdctrig.tdc", &bb_tdctrig_tdc_in);
    C->SetBranchAddress("bb.hodotdc.clus.bar.tdc.meantime", bb_hodotdc_clus_bar_tdc_meantime_in);
    C->SetBranchAddress("bb.hodotdc.clus.bar.tdc.id", bb_hodotdc_clus_bar_tdc_id_in);
    C->SetBranchAddress("bb.gem.track.nhits", bb_gem_track_nhits_in);
    C->SetBranchAddress("bb.hodotdc.clus.tmean", bb_hodotdc_clus_tmean_in);

    double sbs_hcal_clus_blk_tdctime_in[kMaxArraySize];
    double sbs_hcal_clus_blk_tdctime_tw_in[kMaxArraySize];
    double sbs_hcal_clus_blk_id_in[kMaxArraySize];
    double sbs_hcal_clus_x_in[kMaxArraySize];
    double sbs_hcal_clus_y_in[kMaxArraySize];
    double sbs_hcal_clus_e_in[kMaxArraySize];
    C->SetBranchAddress("sbs.hcal.clus_blk.tdctime", sbs_hcal_clus_blk_tdctime_in);
    C->SetBranchAddress("sbs.hcal.clus_blk.tdctime_tw", sbs_hcal_clus_blk_tdctime_tw_in);
    C->SetBranchAddress("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id_in);
    C->SetBranchAddress("sbs.hcal.clus.x", sbs_hcal_clus_x_in);
    C->SetBranchAddress("sbs.hcal.clus.y", sbs_hcal_clus_y_in);
    C->SetBranchAddress("sbs.hcal.clus.e", sbs_hcal_clus_e_in);

    // scalar sbs. branches
    double sbs_hcal_nclus_in;
    C->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus_in);
    double sbs_hcal_e_in;
    double sbs_hcal_x_in;
    double sbs_hcal_y_in;
    double sbs_hcal_rowblk_in;
    double sbs_hcal_colblk_in;
    double sbs_hcal_idblk_in;
    double sbs_hcal_tdctimeblk_in;
    double sbs_tdctrig_rftime_in;
    double sbs_hcal_atimeblk_in;
    C->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e_in);
    C->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x_in);
    C->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y_in);
    C->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk_in);
    C->SetBranchAddress("sbs.hcal.colblk", &sbs_hcal_colblk_in);
    C->SetBranchAddress("sbs.hcal.idblk", &sbs_hcal_idblk_in);
    C->SetBranchAddress("sbs.hcal.tdctimeblk", &sbs_hcal_tdctimeblk_in);
    C->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk_in);
    C->SetBranchAddress("sbs.tdctrig.rftime", &sbs_tdctrig_rftime_in);
    C->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk_in);

    // scalar bb. branches
    double bb_grinch_tdc_clus_trackindex_in;
    double bb_grinch_tdc_clus_size_in;
    double bb_sh_e_in;
    double bb_ps_e_in;
    double bb_tr_n_in;
    double bb_tdctrig_rftime_in;
    double bb_sh_atimeblk_in;
    double bb_ps_atimeblk_in;
    double bb_gem_trigtime_in;
    C->SetBranchAddress("bb.grinch_tdc.clus.trackindex", &bb_grinch_tdc_clus_trackindex_in);
    C->SetBranchAddress("bb.grinch_tdc.clus.size", &bb_grinch_tdc_clus_size_in);
    C->SetBranchAddress("bb.sh.e", &bb_sh_e_in);
    C->SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk_in);
    C->SetBranchAddress("bb.ps.e", &bb_ps_e_in);
    C->SetBranchAddress("bb.tr.n", &bb_tr_n_in);
    C->SetBranchAddress("bb.tdctrig.rftime", &bb_tdctrig_rftime_in);
    C->SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk_in);
    C->SetBranchAddress("bb.ps.atimeblk", &bb_ps_atimeblk_in);
    C->SetBranchAddress("bb.gem.trigtime", &bb_gem_trigtime_in);

    // computed variables e.kine branches
    double e_kine_W2_in;
    double e_kine_Q2_in;
    C->SetBranchAddress("e.kine.W2", &e_kine_W2_in);
    C->SetBranchAddress("e.kine.Q2", &e_kine_Q2_in);

    // other variables
    double g_trigbits_in;
    double HALLA_p_in;
    double hac_bcm_average_in;
    double g_evtime_in;
    double scalhel_hel_in;
    double IHWP_in;
    double g_runnum_in;
    C->SetBranchAddress("g.trigbits", &g_trigbits_in);
    C->SetBranchAddress("HALLA_p", &HALLA_p_in);
    C->SetBranchAddress("hac_bcm_average", &hac_bcm_average_in);
    C->SetBranchAddress("g.evtime", &g_evtime_in);
    C->SetBranchAddress("scalhel.hel", &scalhel_hel_in);
    C->SetBranchAddress("IGL1I00OD16_16", &IHWP_in);
    C->SetBranchAddress("g.runnum", &g_runnum_in);

    // Input branches
    // Long64_t nEntries = C->GetEntries();
    for (Long64_t i = 0; i < totalEvents; ++i) {
	// reading in the tree
	C->GetEntry(i);
	int tNum = C->GetTreeNumber();
	if (bb_tr_n_in>0 && sbs_hcal_e_in > 0.025 && abs(bb_tr_vz_in[0])<=0.33){

	    // vector bb.tr branches being added
	    bb_tr_n = bb_tr_n_in;
	    for (int k = 0; k < bb_tr_n; ++k){
	        bb_tr_p[k] = bb_tr_p_in[k];
	        bb_tr_px[k] = bb_tr_px_in[k];
	        bb_tr_py[k] = bb_tr_py_in[k];
	        bb_tr_pz[k] = bb_tr_pz_in[k];
	        bb_tr_vx[k] = bb_tr_vx_in[k];
	        bb_tr_vy[k] = bb_tr_vy_in[k];
	        bb_tr_vz[k] = bb_tr_vz_in[k];
	        bb_etot_over_p[k] = bb_etot_over_p_in[k];
	        bb_gem_track_nhits[k] = bb_gem_track_nhits_in[k];
	        bb_tr_pathl[k] = bb_tr_pathl_in[k];
	    }  
	    
	    bb_tdctrig_tdcelemID = bb_tdctrig_tdcelemID_in;
	    bb_tdctrig_tdc = bb_tdctrig_tdc_in;
            
	    Ndata_bb_hodotdc_clus_bar_tdc_meantime = Ndata_bb_hodotdc_clus_bar_tdc_meantime_in;
	    for (int k = 0; k < Ndata_bb_hodotdc_clus_bar_tdc_meantime; ++k){
	        bb_hodotdc_clus_bar_tdc_meantime[k] = bb_hodotdc_clus_bar_tdc_meantime_in[k];
	        bb_hodotdc_clus_bar_tdc_id[k] = bb_hodotdc_clus_bar_tdc_id_in[k];
	    }

	    Ndata_bb_hodotdc_clus_tmean = Ndata_bb_hodotdc_clus_tmean_in;
	    for (int k = 0; k < Ndata_bb_hodotdc_clus_tmean; ++k){
	        bb_hodotdc_clus_tmean[k] = bb_hodotdc_clus_tmean_in[k];
	    }
	    
	    int int_sbs_hcal_nclus_in = (int)sbs_hcal_nclus_in;
	    sbs_hcal_nclus = int_sbs_hcal_nclus_in;
	    for (int k = 0; k < sbs_hcal_nclus; ++k){
	        sbs_hcal_clus_blk_tdctime[k] = sbs_hcal_clus_blk_tdctime_in[k];
	        sbs_hcal_clus_blk_tdctime_tw[k] = sbs_hcal_clus_blk_tdctime_tw_in[k];
	        sbs_hcal_clus_blk_id[k] = sbs_hcal_clus_blk_id_in[k];
	    }
	    
	    sbs_hcal_e = sbs_hcal_e_in;
	    sbs_hcal_x = sbs_hcal_x_in;
	    sbs_hcal_y = sbs_hcal_y_in;
	    sbs_hcal_rowblk = sbs_hcal_rowblk_in;
	    sbs_hcal_colblk = sbs_hcal_colblk_in;
	    sbs_hcal_idblk = sbs_hcal_idblk_in;
	    sbs_hcal_tdctimeblk = sbs_hcal_tdctimeblk_in;
	    sbs_hcal_atimeblk = sbs_hcal_atimeblk_in;
	    sbs_tdctrig_rftime = sbs_tdctrig_rftime_in;
	    sbs_hcal_atimeblk = sbs_hcal_atimeblk_in;
	    
	    bb_grinch_tdc_clus_trackindex = bb_grinch_tdc_clus_trackindex_in;
	    bb_grinch_tdc_clus_size = bb_grinch_tdc_clus_size_in;
	    bb_sh_e = bb_sh_e_in;
	    bb_sh_atimeblk = bb_sh_atimeblk_in;
	    bb_ps_e = bb_ps_e_in;
	    bb_tr_n = bb_tr_n_in;
	    bb_tdctrig_rftime = bb_tdctrig_rftime_in;
	    bb_sh_atimeblk = bb_sh_atimeblk_in;
	    bb_ps_atimeblk = bb_ps_atimeblk_in;
	    bb_gem_trigtime = bb_gem_trigtime_in;
	    
	    e_kine_W2 = e_kine_W2_in;
	    e_kine_Q2 = e_kine_Q2_in;
	    e_beam = ebeam;
	    
	    g_trigbits = g_trigbits_in;
	    HALLA_p = HALLA_p_in / 1000; // Converting MeV -> GeV
	    hac_bcm_average = hac_bcm_average_in;
	    g_evtime = g_evtime_in;
	    scalhel_hel = scalhel_hel_in;
	    IHWP = IHWP_in;
	    g_runnum = g_runnum_in;
	    
	    // ---- Beginning of dx and dy Calculation ----
	    
	    double keprime_x, keprime_y, keprime_z, target_x, target_y, target_z;
	    keprime_x = bb_tr_px[0];
	    keprime_y = bb_tr_py[0];
	    keprime_z = bb_tr_pz[0];
	    
	    target_x = bb_tr_vx[0];
	    target_y = bb_tr_vy[0];
	    target_z = bb_tr_vz[0];
	    
	    // Defining momentum 3-vector
	    TVector3 keprime_vec(keprime_x,keprime_y,keprime_z);
	    TVector3 target_vec(target_x,target_y,target_z);
	    
	    double keprime_mag = keprime_vec.Mag();
	    double etheta = keprime_vec.Theta();
	    double ephi = keprime_vec.Phi();
	    double eprime = keprime_vec.Mag();
	    // double ebeam = ebeam_epics/1000;
	    
	    TLorentzVector keprime(keprime_x,keprime_y,keprime_z,eprime);
	    TLorentzVector ke(0.0,0.0,ebeam,ebeam);
	    TLorentzVector P(0.0,0.0,0.0,MN);
	    TLorentzVector PHe3(0.0,0.0,0.0,(2*Mp + Mn));
	    TLorentzVector q = ke - keprime;
	    TVector3 q_vec = q.Vect();
	    TVector3 q_unitvec = q_vec.Unit();
	    TLorentzVector Pprime = q + P;
	    
	    TVector3 Pprime_vec = Pprime.Vect();
	    TVector3 Pprime_unitvec = Pprime_vec.Unit();
	    
	    double w = (HCAL_vector - target_vec).Dot(HCAL_unitvector_z) / (q_unitvec.Dot(HCAL_unitvector_z));
	    TVector3 w_vec = target_vec + w*q_unitvec;
	    TVector3 D_vec = w_vec - HCAL_vector;
	    
	    sbs_hcal_x_exp = D_vec.Dot(HCAL_unitvector_x);
	    sbs_hcal_y_exp = D_vec.Dot(HCAL_unitvector_y);
	    
	    dx = sbs_hcal_x - sbs_hcal_x_exp;
	    dy = sbs_hcal_y - sbs_hcal_y_exp;
	    
	    // Computing dx_clus and dy_clus
	    for (int l = 0; l < sbs_hcal_nclus; ++l){
	        sbs_hcal_clus_x[l] = sbs_hcal_clus_x_in[l];
	        sbs_hcal_clus_y[l] = sbs_hcal_clus_y_in[l];
	        sbs_hcal_clus_e[l] = sbs_hcal_clus_e_in[l];
	        dx_clus[l] = sbs_hcal_clus_x_in[l] - sbs_hcal_x_exp;
	        dy_clus[l] = sbs_hcal_clus_y_in[l] - sbs_hcal_y_exp;
	    }
		
	    TVector3 HCAL_detec_vec = sbs_hcal_x*HCAL_unitvector_x + sbs_hcal_y*HCAL_unitvector_y;
	    TVector3 Pprime_detec_unitvec = (HCAL_detec_vec + HCAL_vector - target_vec).Unit();
	    TVector3 Pprime_pseudo_vec = Pprime_detec_unitvec*(Pprime_vec.Mag());
	    double Pprime_pseudo_E = sqrt(pow(Pprime_vec.Mag(),2) + MN*MN);
	    TLorentzVector Pprime_pseudo(Pprime_pseudo_vec.X(), Pprime_pseudo_vec.Y(), Pprime_pseudo_vec.Z(), Pprime_pseudo_E);
	    // TVector3 pperp_vec = Pprime_detec_vec - q_vec*(Pprime_detec_vec.Dot(q_vec)/q_vec.Mag2());
	    e_ppara_mag = q_unitvec.Dot(q_vec - Pprime_pseudo_vec);
	    e_pperp_mag = (q_vec - Pprime_pseudo_vec - q_unitvec*e_ppara_mag).Mag();
	    e_missing_mass2 = (PHe3 + q - Pprime_pseudo).M2();
	    
	    // filling the output tree
	    outputTree->Fill();
        }
	++eventTracker;
	double percent_calc = 100.0 * static_cast<double>(eventTracker) / static_cast<double>(totalEvents);
	std::cout << "Progress: " << std::fixed << std::setprecision(2) << percent_calc << "% complete\r";
	std::cout.flush();
	
    }

    std::cout << "Progress: 100% complete\n" << std::endl;                  
    std::cout << "All done <3\n" << std::endl; 
    delete C;
    // Write and close the output file
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();

    std::cout << "Merged ROOT file created: " << output_filename << std::endl;
}
