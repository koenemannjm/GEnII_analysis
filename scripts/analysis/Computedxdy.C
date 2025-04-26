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

#include <iostream>

void Computedxdy(std::string config, std::string target, std::string pass){
    // Validate config
    std::vector<std::string> valid_configs = {"kin2", "kin3", "kin4a", "kin4b"};
    if (std::find(valid_configs.begin(), valid_configs.end(), config) == valid_configs.end()) {
        std::cerr << "Error: Invalid configuration '" << config << "'.\n";
        std::cerr << "Allowed configurations are: kin2, kin3, kin4a, kin4b.\n";
        return;
    }

    // Validate target
    std::vector<std::string> valid_targets = {"He3", "H"};
    if (std::find(valid_targets.begin(), valid_targets.end(), target) == valid_targets.end()) {
        std::cerr << "Error: Invalid target '" << target << "'.\n";
        std::cerr << "Allowed configurations are: He3, H.\n";
        return;
    }

    // Validate target
    std::vector<std::string> valid_passes = {"pass1", "pass2"};
    if (std::find(valid_passes.begin(), valid_passes.end(), pass) == valid_passes.end()) {
        std::cerr << "Error: Invalid target '" << pass << "'.\n";
        std::cerr << "Allowed configurations are: pass1, pass2.\n";
        return;
    }
  
    double ebeam;
    double MN = 0.9385;
    double HCAL_angle;
    double HCAL_distance = 17.0;

    // Set experiment name
    std::string exp_name;
    if (config == "kin2") {
        exp_name = "GEN2";
	ebeam = 4.291;
	HCAL_angle = 34.7 * TMath::DegToRad();
    }
    else if (config == "kin3") {
         exp_name = "GEN3";
	 ebeam = 6.373;
	 HCAL_angle = 21.6 * TMath::DegToRad();
    }
    else if (config == "kin4a") {
         exp_name = "GEN4";
	 ebeam = 8.448;
	 HCAL_angle = 18.0 * TMath::DegToRad();
    }
    else if (config == "kin4b") {
         exp_name = "GEN4b";
	 ebeam = 8.448;
	 HCAL_angle = 18.0 * TMath::DegToRad();
    }

    TVector3 HCAL_vector(-HCAL_distance*TMath::Sin(HCAL_angle), 0.0, HCAL_distance*TMath::Cos(HCAL_angle));
    TVector3 HCAL_unitvector_z(TMath::Sin(HCAL_angle), 0.0, TMath::Cos(HCAL_angle));
    TVector3 HCAL_unitvector_y(TMath::Sin(TMath::Pi() - HCAL_angle), 0.0, TMath::Cos(TMath::Pi() - HCAL_angle));
    TVector3 HCAL_unitvector_x = HCAL_unitvector_y.Cross(HCAL_unitvector_z);
    
    std::string dir_path = "../../data/raw/" + pass + "/"  + config + "_" + target + "/";
    std::string file_name = "QE_data_" + exp_name + "_sbs100p_nucleon_np" + ".root";
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
    double keprime_x, keprime_y, keprime_z, eprime_she, eprime_pse, e_over_p, target_x, target_y, target_z;

    tree->SetBranchAddress("bb.tr.px", &keprime_x);
    tree->SetBranchAddress("bb.tr.py", &keprime_y);
    tree->SetBranchAddress("bb.tr.pz", &keprime_z);
    //tree->SetBranchAddress("bb.ps.e", &eprime_she);
    //tree->SetBranchAddress("bb.sh.e", &eprime_pse);
    tree->SetBranchAddress("bb.tr.vx", &target_x);
    tree->SetBranchAddress("bb.tr.vy", &target_y);
    tree->SetBranchAddress("bb.tr.vz", &target_z);

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
			     
    //Loop over entries of the tree
    Long64_t Nentries = tree->GetEntries();
    int last_percent = -1;
    for (Long64_t i = 0; i < Nentries; ++i) {
      tree->GetEntry(i);

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
