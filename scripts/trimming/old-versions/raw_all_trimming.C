//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
//   Created by Jacob Koenemann & ChatGPT
//   contact: bxy3zr@virginia.edu
//                                                                     
//   Last Modified April 4, 2025   
//
//      
//   The purpose of this script is to trim down and the         
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

void raw_all_trimming(std::string config_file) {


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

    // Some constant(s)
    double MN = 0.9385;

    TVector3 HCAL_vector(-HCAL_distance*TMath::Sin(HCAL_angle), 0.0, HCAL_distance*TMath::Cos(HCAL_angle));
    TVector3 HCAL_unitvector_z(TMath::Sin(HCAL_angle), 0.0, TMath::Cos(HCAL_angle));
    TVector3 HCAL_unitvector_y(TMath::Sin(TMath::Pi() - HCAL_angle), 0.0, TMath::Cos(TMath::Pi() - HCAL_angle));
    TVector3 HCAL_unitvector_x = HCAL_unitvector_y.Cross(HCAL_unitvector_z);
   
   // Directory containing the files
    std::string input_dir = "/lustre24/expphy/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/" + pass + "/TEST/try7/" + exp_name + "/" + target  +"/rootfiles/";  // cache directory for Raw data
    // /v/lustre24/expphy/volatile/halla/sbs/sbs-gen/GEN_REPLAYS/pass2/TEST/try7/GEN2/rootfiles

    TSystemDirectory dir("input", input_dir.c_str());
    TList *files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "No files found in directory: " << input_dir << std::endl;
        return;
    }
    std::vector<std::string>matchingFiles;

    TIter next(files);
    TSystemFile* f;
    while ((f = (TSystemFile*)next())) {
        TString fname = f->GetName();
        if (!f->IsDirectory() && !fname.BeginsWith(".")) {
	  matchingFiles.push_back(fname.Data());
        }
    }

    std::string output_dir = "/lustre24/expphy/volatile/halla/sbs/koeneman/data/raw/" + pass + "/test/try7/" + config + "_" + target + "/";
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    std::string output_filename = output_dir + "QE_test_data_" + exp_name + "_sbs100p_nucleon_np" + ".root";

    // Create an output ROOT file
    TFile *outputFile = new TFile(output_filename.c_str(), "RECREATE");
    TTree *outputTree = new TTree("Tout", "Merged data from selected ROOT files");

    const int kMaxArraySize = 100;

    // Defining the output branch names
    // ----------------------------------------------------------------
    
    // vector bb.tr branches
    int bb_tr_n;
    outputTree->Branch("bb.tr.n", &bb_tr_n, "bb.tr.n/I");
    double bb_tr_p[kMaxArraySize], bb_tr_px[kMaxArraySize], bb_tr_py[kMaxArraySize], bb_tr_pz[kMaxArraySize], bb_tr_vx[kMaxArraySize], bb_tr_vy[kMaxArraySize], bb_tr_vz[kMaxArraySize], bb_etot_over_p[kMaxArraySize], bb_gem_track_nhits[kMaxArraySize];
    outputTree->Branch("bb.tr.p", bb_tr_p, "bb.tr.p[bb.tr.n]/D");
    outputTree->Branch("bb.tr.px", bb_tr_px, "bb.tr.px[bb.tr.n]/D");
    outputTree->Branch("bb.tr.py", bb_tr_py, "bb.tr.py[bb.tr.n]/D");
    outputTree->Branch("bb.tr.pz", bb_tr_pz, "bb.tr.pz[bb.tr.n]/D");
    outputTree->Branch("bb.tr.vx", bb_tr_vx, "bb.tr.vx[bb.tr.n]/D");
    outputTree->Branch("bb.tr.vy", bb_tr_vy, "bb.tr.vy[bb.tr.n]/D");
    outputTree->Branch("bb.tr.vz", bb_tr_vz, "bb.tr.vz[bb.tr.n]/D");
    outputTree->Branch("bb.etot_over_p", bb_etot_over_p, "bb.etot_over_p[bb.tr.n]/D");
    outputTree->Branch("bb.gem.track.nhits", bb_gem_track_nhits, "bb.gem.track.nhits[bb.tr.n]/D");

    // vector other bb. branches
    int Ndata_bb_hodotdc_clus_bar_tdc_meantime;
    outputTree->Branch("Ndata.bb.hodotdc.clus.bar.tdc.meantime", &Ndata_bb_hodotdc_clus_bar_tdc_meantime, "Ndata.bb.hodotdc.clus.bar.tdc.meantime/I");
    double bb_tdctrig_tdcelemID, bb_tdctrig_tdc, bb_hodotdc_clus_bar_tdc_meantime[kMaxArraySize], bb_hodotdc_clus_bar_tdc_id[kMaxArraySize];
    outputTree->Branch("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID, "bb.tdctrig.tdcelemID/D");
    outputTree->Branch("bb.tdctrig.tdc", &bb_tdctrig_tdc, "bb.tdctrig.tdc/D");
    outputTree->Branch("bb.hodotdc.clus.bar.tdc.meantime", bb_hodotdc_clus_bar_tdc_meantime, "bb.hodotdc.clus.bar.tdc.meantime[Ndata.bb.hodotdc.clus.bar.tdc.meantime]/D");
    outputTree->Branch("bb.hodotdc.clus.bar.tdc.id", bb_hodotdc_clus_bar_tdc_id, "bb.hodotdc.clus.bar.tdc.id[Ndata.bb.hodotdc.clus.bar.tdc.meantime]/D");

    // vector sbs.hcal branches
    int  sbs_hcal_nclus;
    outputTree->Branch("sbs.hcal.nclus", &sbs_hcal_nclus, "sbs.hcal.nclus/I");
    double sbs_hcal_clus_blk_tdctime[kMaxArraySize], sbs_hcal_clus_blk_tdctime_tw[kMaxArraySize], sbs_hcal_clus_blk_id[kMaxArraySize];
    outputTree->Branch("sbs.hcal.clus_blk.tdctime", sbs_hcal_clus_blk_tdctime, "sbs.hcal.clus_blk.tdctime[sbs.hcal.nclus]/D");
    outputTree->Branch("sbs.hcal.clus_blk.tdctime_tw", sbs_hcal_clus_blk_tdctime_tw, "sbs.hcal.clus_blk.tdctime_tw[sbs.hcal.nclus]/D");
    outputTree->Branch("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id, "sbs.hcal.clus_blk.id[sbs.hcal.nclus]/D");

    // scalar sbs. branches
    double sbs_hcal_e, sbs_hcal_x, sbs_hcal_y, sbs_hcal_rowblk, sbs_hcal_colblk, sbs_hcal_idblk, sbs_hcal_tdctimeblk, sbs_hcal_atimeblk, sbs_tdctrig_rftime;
    outputTree->Branch("sbs.hcal.e", &sbs_hcal_e, "sbs.hcal.e/D");
    outputTree->Branch("sbs.hcal.x", &sbs_hcal_x, "sbs.hcal.x/D");
    outputTree->Branch("sbs.hcal.y", &sbs_hcal_y, "sbs.hcal.y/D");
    outputTree->Branch("sbs.hcal.rowblk", &sbs_hcal_rowblk, "sbs.hcal.rowblk/D");
    outputTree->Branch("sbs.hcal.colblk", &sbs_hcal_colblk, "sbs.hcal.colblk/D");
    outputTree->Branch("sbs.hcal.idblk", &sbs_hcal_idblk, "sbs.hcal.idblk/D");
    outputTree->Branch("sbs.hcal.tdctimeblk", &sbs_hcal_tdctimeblk, "sbs.hcal.tdctimeblk/D");
    outputTree->Branch("sbs.hcal.atimeblk", &sbs_hcal_atimeblk, "sbs.hcal.atimeblk/D");
    outputTree->Branch("sbs.tdctrig.rftime", &sbs_tdctrig_rftime, "sbs.tdctrig.rftime/D");

    // scalar bb. branches
    double bb_grinch_tdc_clus_trackindex, bb_grinch_tdc_clus_size, bb_sh_e, bb_sh_atimeblk, bb_ps_e,  bb_tdctrig_rftime;
    outputTree->Branch("bb.grinch_tdc.clus.trackindex", &bb_grinch_tdc_clus_trackindex, "bb.grinch_tdc.clus.trackindex/D");
    outputTree->Branch("bb.grinch_tdc.clus.size", &bb_grinch_tdc_clus_size, "bb.grind_tdc.clus.size/D");
    outputTree->Branch("bb.sh.e", &bb_sh_e, "bb.sh.e/D");
    outputTree->Branch("bb.sh.atimeblk", &bb_sh_atimeblk, "bb.sh.atimeblk/D");
    outputTree->Branch("bb.ps.e", &bb_ps_e, "bb.ps.e/D");
    outputTree->Branch("bb.tdctrig.rftime", &bb_tdctrig_rftime, "bb.tdctrig.rftime/D");

    // computed variables e.kine branches
    double e_kine_W2, e_kine_Q2;
    outputTree->Branch("e.kine.W2", &e_kine_W2, "e.kine.W2/D");
    outputTree->Branch("e.kine.Q2", &e_kine_Q2, "e.kine.Q2/D");

    // other variables
    double g_trigbits, HALLA_p, hac_bcm_average, g_evtime;
    outputTree->Branch("g.trigbits", &g_trigbits, "g.trigbits/D");
    outputTree->Branch("HALLA_p", &HALLA_p, "HALLA_p/D");
    outputTree->Branch("hac_bcm_average", &hac_bcm_average, "hac_bcm_average/D");
    outputTree->Branch("g.evtime", &g_evtime, "g.evtime/D");

    // dx and dy variables
    double dx, dy, sbs_hcal_x_exp, sbs_hcal_y_exp;
    outputTree->Branch("sbs.hcal.x_exp", &sbs_hcal_x_exp, "sbs.hcal.x_exp/D");
    outputTree->Branch("sbs.hcal.y_exp", &sbs_hcal_y_exp, "sbs.hcal.y_exp/D");
    outputTree->Branch("dx", &dx, "dx/D");
    outputTree->Branch("dy", &dy, "dy/D");
    // ----------------------------------------------------------------

    std::cout<< "\nCalculating total events\n" << "...\n" << std::endl;
    // Computing total number of events for all files
    Long64_t totalEvents = 0;
    for (const auto& filename : matchingFiles) {
        std::string filePath = input_dir + filename;
        TFile* inFile = TFile::Open(filePath.c_str(),"READ");
	if (!inFile || inFile->IsZombie()) continue;
	TTree* inTree = (TTree*)inFile->Get("T");
	if (!inTree) {
	    std::cerr << "Tree 'T' not found in file: " << filename << std::endl;
	    inFile->Close();
	    continue;
	}
	totalEvents += inTree->GetEntries();
	inFile->Close();
	delete inFile;
    }
    std::cout<< "\nTotal number events = " << totalEvents << "\n" << std::endl;

    int totalFiles = matchingFiles.size();
    int eventTracker = 0;
    for (size_t j = 0; j < totalFiles; ++j) {

        const std::string& filename = matchingFiles[j];
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

	inputTree->SetBranchStatus("*",0);

	inputTree->SetBranchStatus("bb.tr.p", 1);
	inputTree->SetBranchStatus("bb.tr.px", 1);
	inputTree->SetBranchStatus("bb.tr.py", 1);
	inputTree->SetBranchStatus("bb.tr.pz", 1);
	inputTree->SetBranchStatus("bb.tr.vx", 1);
	inputTree->SetBranchStatus("bb.tr.vy", 1);
	inputTree->SetBranchStatus("bb.tr.vz", 1);
	inputTree->SetBranchStatus("bb.etot_over_p", 1);
	inputTree->SetBranchStatus("Ndata.bb.hodotdc.clus.bar.tdc.meantime", 1);
	inputTree->SetBranchStatus("bb.tdctrig.tdcelemID", 1);
	inputTree->SetBranchStatus("bb.tdctrig.tdc", 1);
	inputTree->SetBranchStatus("bb.hodotdc.clus.bar.tdc.meantime", 1);
	inputTree->SetBranchStatus("bb.hodotdc.clus.bar.tdc.id", 1);
	inputTree->SetBranchStatus("bb.gem.track.nhits", 1);
	inputTree->SetBranchStatus("sbs.hcal.clus_blk.tdctime", 1);
	inputTree->SetBranchStatus("sbs.hcal.clus_blk.tdctime_tw", 1);
	inputTree->SetBranchStatus("sbs.hcal.clus_blk.id", 1);
	inputTree->SetBranchStatus("sbs.hcal.nclus", 1);
	inputTree->SetBranchStatus("sbs.hcal.clus_blk.tdctime_tw", 1);
	inputTree->SetBranchStatus("sbs.hcal.clus_blk.tdctime_tw", 1);
	inputTree->SetBranchStatus("sbs.hcal.e", 1);
        inputTree->SetBranchStatus("sbs.hcal.x", 1);
        inputTree->SetBranchStatus("sbs.hcal.y", 1);
        inputTree->SetBranchStatus("sbs.hcal.rowblk", 1);
        inputTree->SetBranchStatus("sbs.hcal.colblk", 1);
        inputTree->SetBranchStatus("sbs.hcal.idblk", 1);
        inputTree->SetBranchStatus("sbs.hcal.tdctimeblk", 1);
        inputTree->SetBranchStatus("sbs.hcal.atimeblk", 1);
        inputTree->SetBranchStatus("sbs.tdctrig.rftime", 1);
        inputTree->SetBranchStatus("bb.grinch_tdc.clus.trackindex", 1);
        inputTree->SetBranchStatus("bb.grinch_tdc.clus.size", 1);
        inputTree->SetBranchStatus("bb.sh.e", 1);
        inputTree->SetBranchStatus("bb.sh.atimeblk", 1);
        inputTree->SetBranchStatus("bb.ps.e", 1);
        inputTree->SetBranchStatus("bb.tr.n", 1);
        inputTree->SetBranchStatus("bb.tdctrig.rftime", 1);
        inputTree->SetBranchStatus("e.kine.W2", 1);
        inputTree->SetBranchStatus("e.kine.Q2", 1);
        inputTree->SetBranchStatus("g.trigbits", 1);
	inputTree->SetBranchStatus("HALLA_p", 1);
	inputTree->SetBranchStatus("hac_bcm_average", 1);
	inputTree->SetBranchStatus("g.evtime", 1);

	double bb_tr_p_in[kMaxArraySize], bb_tr_px_in[kMaxArraySize], bb_tr_py_in[kMaxArraySize], bb_tr_pz_in[kMaxArraySize], bb_tr_vx_in[kMaxArraySize], bb_tr_vy_in[kMaxArraySize], bb_tr_vz_in[kMaxArraySize], bb_etot_over_p_in[kMaxArraySize];
        inputTree->SetBranchAddress("bb.tr.p", bb_tr_p_in);
        inputTree->SetBranchAddress("bb.tr.px", bb_tr_px_in);
        inputTree->SetBranchAddress("bb.tr.py", bb_tr_py_in);
        inputTree->SetBranchAddress("bb.tr.pz", bb_tr_pz_in);
        inputTree->SetBranchAddress("bb.tr.vx", bb_tr_vx_in);
        inputTree->SetBranchAddress("bb.tr.vy", bb_tr_vy_in);
        inputTree->SetBranchAddress("bb.tr.vz", bb_tr_vz_in);
        inputTree->SetBranchAddress("bb.etot_over_p", bb_etot_over_p_in);

	int Ndata_bb_hodotdc_clus_bar_tdc_meantime_in;
        inputTree->SetBranchAddress("Ndata.bb.hodotdc.clus.bar.tdc.meantime", &Ndata_bb_hodotdc_clus_bar_tdc_meantime_in);
	double bb_tdctrig_tdcelemID_in, bb_tdctrig_tdc_in, bb_hodotdc_clus_bar_tdc_meantime_in[kMaxArraySize], bb_hodotdc_clus_bar_tdc_id_in[kMaxArraySize], bb_gem_track_nhits_in[kMaxArraySize];
        inputTree->SetBranchAddress("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID_in);
	inputTree->SetBranchAddress("bb.tdctrig.tdc", &bb_tdctrig_tdc_in);
	inputTree->SetBranchAddress("bb.hodotdc.clus.bar.tdc.meantime", bb_hodotdc_clus_bar_tdc_meantime_in);
	inputTree->SetBranchAddress("bb.hodotdc.clus.bar.tdc.id", bb_hodotdc_clus_bar_tdc_id_in);
	inputTree->SetBranchAddress("bb.gem.track.nhits", bb_gem_track_nhits_in);

	double sbs_hcal_clus_blk_tdctime_in[kMaxArraySize], sbs_hcal_clus_blk_tdctime_tw_in[kMaxArraySize], sbs_hcal_clus_blk_id_in[kMaxArraySize];
        inputTree->SetBranchAddress("sbs.hcal.clus_blk.tdctime", sbs_hcal_clus_blk_tdctime_in);
        inputTree->SetBranchAddress("sbs.hcal.clus_blk.tdctime_tw", sbs_hcal_clus_blk_tdctime_tw_in);
        inputTree->SetBranchAddress("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id_in);

	// scalar sbs. branches
	double sbs_hcal_nclus_in;
	inputTree->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus_in);
        double sbs_hcal_e_in, sbs_hcal_x_in, sbs_hcal_y_in, sbs_hcal_rowblk_in, sbs_hcal_colblk_in, sbs_hcal_idblk_in, sbs_hcal_tdctimeblk_in, sbs_hcal_atimeblk_in,  sbs_tdctrig_rftime_in;
        inputTree->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e_in);
        inputTree->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x_in);
        inputTree->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y_in);
        inputTree->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk_in);
        inputTree->SetBranchAddress("sbs.hcal.colblk", &sbs_hcal_colblk_in);
        inputTree->SetBranchAddress("sbs.hcal.idblk", &sbs_hcal_idblk_in);
        inputTree->SetBranchAddress("sbs.hcal.tdctimeblk", &sbs_hcal_tdctimeblk_in);
        inputTree->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk_in);
        inputTree->SetBranchAddress("sbs.tdctrig.rftime", &sbs_tdctrig_rftime_in);

	// scalar bb. branches
        double bb_grinch_tdc_clus_trackindex_in, bb_grinch_tdc_clus_size_in, bb_sh_e_in, bb_sh_atimeblk_in, bb_ps_e_in, bb_tr_n_in, bb_tdctrig_rftime_in;
        inputTree->SetBranchAddress("bb.grinch_tdc.clus.trackindex", &bb_grinch_tdc_clus_trackindex_in);
        inputTree->SetBranchAddress("bb.grinch_tdc.clus.size", &bb_grinch_tdc_clus_size_in);
        inputTree->SetBranchAddress("bb.sh.e", &bb_sh_e_in);
        inputTree->SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk_in);
        inputTree->SetBranchAddress("bb.ps.e", &bb_ps_e_in);
        inputTree->SetBranchAddress("bb.tr.n", &bb_tr_n_in);
        inputTree->SetBranchAddress("bb.tdctrig.rftime", &bb_tdctrig_rftime_in);

        // computed variables e.kine branches
        double e_kine_W2_in, e_kine_Q2_in;
        inputTree->SetBranchAddress("e.kine.W2", &e_kine_W2_in);
        inputTree->SetBranchAddress("e.kine.Q2", &e_kine_Q2_in);

        // other variables
        double g_trigbits_in, HALLA_p_in, hac_bcm_average_in, g_evtime_in;
        inputTree->SetBranchAddress("g.trigbits", &g_trigbits_in);
	inputTree->SetBranchAddress("HALLA_p", &HALLA_p_in);
	inputTree->SetBranchAddress("hac_bcm_average", &hac_bcm_average);
	inputTree->SetBranchAddress("g.evtime", &g_evtime_in);

	// Input branches
        Long64_t nEntries = inputTree->GetEntries();
	for (Long64_t i = 0; i < nEntries; ++i) {
	    // reading in the tree
	    inputTree->GetEntry(i);
	   
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
	    }
	    
            
	    bb_tdctrig_tdcelemID = bb_tdctrig_tdcelemID_in;
	    bb_tdctrig_tdc = bb_tdctrig_tdc_in;

	    Ndata_bb_hodotdc_clus_bar_tdc_meantime = Ndata_bb_hodotdc_clus_bar_tdc_meantime_in;
	    for (int k = 0; k < Ndata_bb_hodotdc_clus_bar_tdc_meantime; ++k){
	        bb_hodotdc_clus_bar_tdc_meantime[k] = bb_hodotdc_clus_bar_tdc_meantime_in[k];
	        bb_hodotdc_clus_bar_tdc_id[k] = bb_hodotdc_clus_bar_tdc_id_in[k];
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

	    bb_grinch_tdc_clus_trackindex = bb_grinch_tdc_clus_trackindex_in;
	    bb_grinch_tdc_clus_size = bb_grinch_tdc_clus_size_in;
	    bb_sh_e = bb_sh_e_in;
	    bb_sh_atimeblk = bb_sh_atimeblk_in;
	    bb_ps_e = bb_ps_e_in;
	    bb_tr_n = bb_tr_n_in;
	    bb_tdctrig_rftime = bb_tdctrig_rftime_in;

	    e_kine_W2 = e_kine_W2_in;
	    e_kine_Q2 = e_kine_Q2_in;

	    g_trigbits = g_trigbits_in;
	    HALLA_p = HALLA_p_in / 1000;
	    hac_bcm_average = hac_bcm_average_in;
	    g_evtime = g_evtime_in;

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
	    TLorentzVector q = ke - keprime;
	    TLorentzVector Pprime = q + P;

	    double eprime_el = ebeam / (1 + ebeam/MN * (1 - TMath::Cos(etheta)));
	    TVector3 keprime_el_vec(eprime_el*TMath::Cos(ephi)*TMath::Sin(etheta),eprime_el*TMath::Sin(ephi)*TMath::Sin(etheta),eprime_el*TMath::Cos(etheta));
	    TLorentzVector keprime_el(keprime_el_vec.X(),keprime_el_vec.Y(),keprime_el_vec.Z(),eprime_el);

	    TLorentzVector Pprime_el = ke - keprime_el + P;


	    TVector3 Pprime_vec = Pprime.Vect();
	    TVector3 Pprime_unitvec = Pprime_vec.Unit();

	    double w = (HCAL_vector - target_vec).Dot(HCAL_unitvector_z) / (Pprime_unitvec.Dot(HCAL_unitvector_z));
	    TVector3 w_vec = target_vec + w*Pprime_unitvec;
	    TVector3 D_vec = w_vec - HCAL_vector;

	    sbs_hcal_x_exp = D_vec.Dot(HCAL_unitvector_x);
	    sbs_hcal_y_exp = D_vec.Dot(HCAL_unitvector_y);

	    dx = sbs_hcal_x - sbs_hcal_x_exp;
	    dy = sbs_hcal_y - sbs_hcal_y_exp;
	    
	    // filling the output tree
	    outputTree->Fill();
	     ++eventTracker;
	    double percent_calc = 100.0 * eventTracker / totalEvents;
	    std::cout << "Progress: " << std::fixed << std::setprecision(2) << percent_calc << "% complete\r" << std::flush;
	}
	
	

	file->Close();
	delete file;
	
    }

    
    // Write and close the output file
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();

    std::cout << "Merged ROOT file created: " << output_filename << std::endl;
}
