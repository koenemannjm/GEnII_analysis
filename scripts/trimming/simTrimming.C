//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
//   Created by Jacob Koenemann & ChatGPT
//   contact: bxy3zr@virginia.edu
//                                                                     
//   Last Modified May 12, 2025   
//
//      
//   The purpose of this script is to trim down on the         
//   sim files for GEnII experiment and produce a
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

void sim_trimming(std::string config_file) {


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
    std::string proc = settings["proc"];
    std::string target = settings["target"];
    std::string base_dir = settings["dir"];

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
    std::string input_dir = "/v" + base_dir + exp_name + "/"; 

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
        if (!f->IsDirectory() && fname.BeginsWith("replayed") && fname.EndsWith(".root")) {
	  matchingFiles.push_back(fname.Data());
        }
    }

    std::string output_dir = "/v/lustre24/expphy/volatile/halla/sbs/koeneman/data/sim/" + proc + "/"  + config + "_" + proc + "_" + target + "/";
    gSystem->mkdir(output_dir.c_str(), kTRUE);
    std::string output_filename = output_dir + proc + "_sim_" + exp_name + "_sbs100p_nucleon_np" + ".root";

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

    // vector sbs.hcal branches
    int  sbs_hcal_nclus;
    outputTree->Branch("sbs.hcal.nclus", &sbs_hcal_nclus, "sbs.hcal.nclus/I");
    double sbs_hcal_clus_blk_id[kMaxArraySize];
    outputTree->Branch("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id, "sbs.hcal.clus_blk.id[sbs.hcal.nclus]/D");

    // scalar sbs. branches
    double sbs_hcal_e, sbs_hcal_x, sbs_hcal_y, sbs_hcal_rowblk, sbs_hcal_colblk, sbs_hcal_idblk;
    outputTree->Branch("sbs.hcal.e", &sbs_hcal_e, "sbs.hcal.e/D");
    outputTree->Branch("sbs.hcal.x", &sbs_hcal_x, "sbs.hcal.x/D");
    outputTree->Branch("sbs.hcal.y", &sbs_hcal_y, "sbs.hcal.y/D");
    outputTree->Branch("sbs.hcal.rowblk", &sbs_hcal_rowblk, "sbs.hcal.rowblk/D");
    outputTree->Branch("sbs.hcal.colblk", &sbs_hcal_colblk, "sbs.hcal.colblk/D");
    outputTree->Branch("sbs.hcal.idblk", &sbs_hcal_idblk, "sbs.hcal.idblk/D");

    // scalar bb. branches
    double bb_sh_e, bb_ps_e;
    outputTree->Branch("bb.sh.e", &bb_sh_e, "bb.sh.e/D");
    outputTree->Branch("bb.ps.e", &bb_ps_e, "bb.ps.e/D");

    // computed variables e.kine branches
    double e_kine_W2, e_kine_Q2;
    outputTree->Branch("e.kine.W2", &e_kine_W2, "e.kine.W2/D");
    outputTree->Branch("e.kine.Q2", &e_kine_Q2, "e.kine.Q2/D");

    // computed variables MC. branches
    double MC_mc_sigmaold, MC_mc_sigma, MC_mc_sigmaPol, MC_mc_omega, MC_mc_weight, MC_mc_fnucl, MC_mc_THETA, MC_mc_BETA;
    outputTree->Branch("MC.mc_sigmaold", &MC_mc_sigmaold, "MC.mc_sigmaold/D");
    outputTree->Branch("MC.mc_sigma", &MC_mc_sigma, "MC.mc_sigma/D");
    outputTree->Branch("MC.mc_sigmaPol", &MC_mc_sigmaPol, "MC.mc_sigmaPol/D");
    outputTree->Branch("MC.mc_omega", &MC_mc_omega, "MC.mc_omega/D");
    outputTree->Branch("MC.mc_weight", &MC_mc_weight, "MC.mc_weight/D");
    outputTree->Branch("MC.mc_fnucl", &MC_mc_fnucl, "MC.mc_fnucl/D");
    outputTree->Branch("MC.mc_THETA", &MC_mc_THETA, "MC.mc_THETA/D");
    outputTree->Branch("MC.mc_BETA", &MC_mc_BETA, "MC.mc_BETA/D");

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
	inputTree->SetBranchStatus("bb.gem.track.nhits", 1);
	inputTree->SetBranchStatus("sbs.hcal.clus_blk.id", 1);
	inputTree->SetBranchStatus("sbs.hcal.nclus", 1);
	inputTree->SetBranchStatus("sbs.hcal.e", 1);
        inputTree->SetBranchStatus("sbs.hcal.x", 1);
        inputTree->SetBranchStatus("sbs.hcal.y", 1);
        inputTree->SetBranchStatus("sbs.hcal.rowblk", 1);
        inputTree->SetBranchStatus("sbs.hcal.colblk", 1);
        inputTree->SetBranchStatus("sbs.hcal.idblk", 1);
        inputTree->SetBranchStatus("bb.sh.e", 1);
        inputTree->SetBranchStatus("bb.ps.e", 1);
        inputTree->SetBranchStatus("bb.tr.n", 1);
	inputTree->SetBranchStatus("MC.mc_sigmaold", 1);
	inputTree->SetBranchStatus("MC.mc_sigma", 1);
	inputTree->SetBranchStatus("MC.mc_sigmaPol", 1);
	inputTree->SetBranchStatus("MC.mc_omega", 1);
	inputTree->SetBranchStatus("MC.mc_fnucl", 1);
	inputTree->SetBranchStatus("MC.mc_THETA", 1);
	inputTree->SetBranchStatus("MC.mc_BETA", 1);

	double bb_tr_p_in[kMaxArraySize], bb_tr_px_in[kMaxArraySize], bb_tr_py_in[kMaxArraySize], bb_tr_pz_in[kMaxArraySize], bb_tr_vx_in[kMaxArraySize], bb_tr_vy_in[kMaxArraySize], bb_tr_vz_in[kMaxArraySize], bb_etot_over_p_in[kMaxArraySize];
        inputTree->SetBranchAddress("bb.tr.p", bb_tr_p_in);
        inputTree->SetBranchAddress("bb.tr.px", bb_tr_px_in);
        inputTree->SetBranchAddress("bb.tr.py", bb_tr_py_in);
        inputTree->SetBranchAddress("bb.tr.pz", bb_tr_pz_in);
        inputTree->SetBranchAddress("bb.tr.vx", bb_tr_vx_in);
        inputTree->SetBranchAddress("bb.tr.vy", bb_tr_vy_in);
        inputTree->SetBranchAddress("bb.tr.vz", bb_tr_vz_in);
        inputTree->SetBranchAddress("bb.etot_over_p", bb_etot_over_p_in);

	double bb_gem_track_nhits_in[kMaxArraySize];
	inputTree->SetBranchAddress("bb.gem.track.nhits", bb_gem_track_nhits_in);

	double sbs_hcal_clus_blk_id_in[kMaxArraySize];
        inputTree->SetBranchAddress("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id_in);

	// scalar sbs. branches
	double sbs_hcal_nclus_in;
	inputTree->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus_in);
        double sbs_hcal_e_in, sbs_hcal_x_in, sbs_hcal_y_in, sbs_hcal_rowblk_in, sbs_hcal_colblk_in, sbs_hcal_idblk_in;
        inputTree->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e_in);
        inputTree->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x_in);
        inputTree->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y_in);
        inputTree->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk_in);
        inputTree->SetBranchAddress("sbs.hcal.colblk", &sbs_hcal_colblk_in);
        inputTree->SetBranchAddress("sbs.hcal.idblk", &sbs_hcal_idblk_in);

	// scalar bb. branches
        double bb_sh_e_in, bb_ps_e_in, bb_tr_n_in;
        inputTree->SetBranchAddress("bb.sh.e", &bb_sh_e_in);
        inputTree->SetBranchAddress("bb.ps.e", &bb_ps_e_in);
        inputTree->SetBranchAddress("bb.tr.n", &bb_tr_n_in);

	// computed variables MC. branches
	double MC_mc_sigmaold_in, MC_mc_sigma_in, MC_mc_sigmaPol_in, MC_mc_fnucl_in, MC_mc_THETA_in, MC_mc_BETA_in, MC_mc_omega_in;
	inputTree->SetBranchAddress("MC.mc_sigmaold", &MC_mc_sigmaold_in);
	inputTree->SetBranchAddress("MC.mc_sigma", &MC_mc_sigma_in);
	inputTree->SetBranchAddress("MC.mc_sigmaPol", &MC_mc_sigmaPol_in);
	inputTree->SetBranchAddress("MC.mc_omega", &MC_mc_omega_in);
	inputTree->SetBranchAddress("MC.mc_fnucl", &MC_mc_fnucl_in);
	inputTree->SetBranchAddress("MC.mc_THETA", &MC_mc_THETA_in);
	inputTree->SetBranchAddress("MC.mc_BETA", &MC_mc_BETA_in);

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
	    
            int int_sbs_hcal_nclus_in = (int)sbs_hcal_nclus_in;
	    sbs_hcal_nclus = int_sbs_hcal_nclus_in;
	    for (int k = 0; k < sbs_hcal_nclus; ++k){
	         sbs_hcal_clus_blk_id[k] = sbs_hcal_clus_blk_id_in[k];
	    }
	    
	    sbs_hcal_e = sbs_hcal_e_in;
	    sbs_hcal_x = sbs_hcal_x_in;
	    sbs_hcal_y = sbs_hcal_y_in;
	    sbs_hcal_rowblk = sbs_hcal_rowblk_in;
	    sbs_hcal_colblk = sbs_hcal_colblk_in;
	    sbs_hcal_idblk = sbs_hcal_idblk_in;

	    bb_sh_e = bb_sh_e_in;
	    bb_ps_e = bb_ps_e_in;
	    bb_tr_n = bb_tr_n_in;

	    MC_mc_sigmaold = MC_mc_sigmaold_in;
	    MC_mc_sigma = MC_mc_sigma_in;
	    MC_mc_sigmaPol = MC_mc_sigmaPol_in;
	    MC_mc_omega = MC_mc_omega_in;
	    MC_mc_weight = MC_mc_omega * MC_mc_sigma;
            MC_mc_fnucl = MC_mc_fnucl_in;
	    MC_mc_THETA = MC_mc_THETA_in;
	    MC_mc_BETA = MC_mc_BETA_in;

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

	    e_kine_W2 = Pprime.Mag2();
	    e_kine_Q2 = -q.Mag2();

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
