//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
//   Created by Jacob Koenemann & ChatGPT
//   contact: bxy3zr@virginia.edu
//                                                                     
//   Last Modified April 26, 2025   
//
//      
//   The purpose of this script is to compute dx     
//   and dy and add these variables to the trimmed
//   ROOT file from either data or simu.
//                                                    
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////  
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"

// Need for reading .cfg file
#include <iostream>
#include <cstdlib>  // for std::exit and EXIT_FAILURE
#include <fstream>

void Computedxdy(std::string config_file){

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
    
    std::string dir_path = "/lustre24/expphy/volatile/halla/sbs/koeneman/data/raw/" + pass +"/test/try7" + "/"  + config + "_" + target + "/";
    std::string file_name = "QE_test_data_" + exp_name + "_sbs100p_nucleon_np" + ".root";
    std::string file_path = dir_path + file_name;

    // Open ROOT file
    TFile* file = new TFile(file_path.c_str(),"UPDATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
	return;
    }

    // Get tree "Tout"
    TTree* tree = (TTree*)file->Get("Tout");
    if (!tree) {
        std::cerr << "Tree 'Tout' not found!" << std::endl;
	return;
    }

    // Declare variables for branches

    // Epics data
    double ebeam_epics;

    tree->SetBranchAddress("HALLA_p", &ebeam_epics);
    
    // BigBite variables
    std::vector<double>* px = nullptr;
    std::vector<double>* py = nullptr;
    std::vector<double>* pz = nullptr;
    double eprime_she, eprime_pse, e_over_p;
    std::vector<double>* tx = nullptr;
    std::vector<double>* ty = nullptr;
    std::vector<double>* tz = nullptr;

    tree->SetBranchAddress("bb.tr.px", &px);
    tree->SetBranchAddress("bb.tr.py", &py);
    tree->SetBranchAddress("bb.tr.pz", &pz);
    //tree->SetBranchAddress("bb.ps.e", &eprime_she);
    //tree->SetBranchAddress("bb.sh.e", &eprime_pse);
    tree->SetBranchAddress("bb.tr.vx", &tx);
    tree->SetBranchAddress("bb.tr.vy", &ty);
    tree->SetBranchAddress("bb.tr.vz", &tz);

    // Hcal variables
    double Hcalx, Hcaly, Hcale;

    tree->SetBranchAddress("sbs.hcal.x", &Hcalx);
    tree->SetBranchAddress("sbs.hcal.y", &Hcaly);
    //tree->SetBranchAddress("sbs.hcal.e", &Hcale);

    // Setting sbs.hcal.x_exp, sbs.hcal.y_exp, dx, and dy Branch
    double Hcalx_exp, Hcaly_exp, dx, dy;
    TBranch* sbs_hcal_x_exp = tree->Branch("sbs.hcal.x_exp", &Hcalx_exp, "sbs.hcal.x_exp/D");
    TBranch* sbs_hcal_y_exp = tree->Branch("sbs.hcal.y_exp", &Hcaly_exp, "sbs.hcal.y_exp/D");
    TBranch* sbs_dx = tree->Branch("dx", &dx, "dx/D");
    TBranch* sbs_dy = tree->Branch("dy", &dy, "dy/D");

    double keprime_x, keprime_y, keprime_z, target_x, target_y, target_z;
			     
    //Loop over entries of the tree
    Long64_t Nentries = tree->GetEntries();
    int last_percent = -1;
    for (Long64_t i = 0; i < Nentries; ++i) {
      tree->GetEntry(i);

      if (px->size() > 0 && py->size() > 0 && pz->size() > 0) {
	  keprime_x = px->at(0);
	  keprime_y = py->at(0);
	  keprime_z = pz->at(0);
	  target_x = tx->at(0);
	  target_y = ty->at(0);
	  target_z = tz->at(0);
      }

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

      Hcaly_exp = D_vec.Dot(HCAL_unitvector_y);
      Hcalx_exp = D_vec.Dot(HCAL_unitvector_x);

      dx = Hcalx - Hcalx_exp;
      dy = Hcaly - Hcaly_exp;

      sbs_hcal_x_exp->Fill();
      sbs_hcal_y_exp->Fill();
      sbs_dx->Fill();
      sbs_dy->Fill();

      int percent = static_cast<int>(100.0 * i / Nentries);
      if (percent != last_percent) {
	 std::cout << "\rProcessing: " << percent << "% completed" << std::flush;
	 last_percent = percent;
      }

      
    }
    std::cout << "\rProcessing: 100% completed" << std::endl;
    tree->Write("",TObject::kOverwrite);
    file->Close();

    
}
