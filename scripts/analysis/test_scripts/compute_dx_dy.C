//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
//   Created by Jacob Koenemann & ChatGPT
//   contact: bxy3zr@virginia.edu
//                                                                     
//   Last Modified April 6, 2025   
//
//      
//   The purpose of this script is to test out how         
//   to compute dx and dy and how to structure the 
//   a header file that I will use in the future.
//                                                    
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////  
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <iostream>

void compute_dx_dy(std::string config, std::string target, std::string pass){
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
    
    std::string dir_path = "../../../data/raw/" + pass + "/"  + config + "_" + target + "/";
    std::string file_name = "QE_data_" + exp_name + "_sbs100p_nucleon_np" + ".root";
    std::string file_path = dir_path + file_name;

    // Open ROOT file
    TFile* file = TFile::Open(file_path.c_str());
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
    tree->SetBranchAddress("bb.ps.e", &eprime_she);
    tree->SetBranchAddress("bb.sh.e", &eprime_pse);
    tree->SetBranchAddress("bb.etot_over_p", &e_over_p);
    tree->SetBranchAddress("bb.tr.vx", &target_x);
    tree->SetBranchAddress("bb.tr.vy", &target_y);
    tree->SetBranchAddress("bb.tr.vz", &target_z);

    // Hcal variables
    double Hcalx, Hcaly, Hcale;

    tree->SetBranchAddress("sbs.hcal.x", &Hcalx);
    tree->SetBranchAddress("sbs.hcal.y", &Hcaly);
    tree->SetBranchAddress("sbs.hcal.e", &Hcale);

    

    // Pre-Defining Histograms for Calculations results
    TH2D* h_dx_dy = new TH2D("h_dx_dy", "dx vs dy;dy;dx", 100, -11.0, 2.0, 100, -4.5, 2.0);
    TH2D* h_W2_dy = new TH2D("h_W2_dy", "W2 vs dy;dy;W2", 100, -11.0, 2.0, 100, -2.0, 9.0);
    TH2D* h_dx_dy_cut = new TH2D("h_dx_dy_cut", "dx vs dy", 100, -1.0, 1.0, 100, -4.5, 2.0);
    TH1D* h_dx = new TH1D("h_dx", "dx", 100, -4.5, 2.0);
    TH1D* h_dx_cut = new TH1D("h_dx_cut", "dx", 100, -4.5, 2.0);
			     
    //Loop over entries of the tree
    Long64_t Nentries = tree->GetEntries();
    int last_percent = -1;
    for (Long64_t i = 0; i < Nentries; ++i) {
      tree->GetEntry(i);
      if (TMath::Abs(target_z) <= 0.27 && TMath::Abs(target_x) <= 0.0115 && TMath::Abs(target_y) <= 0.0115 && Hcale > 0.025 && eprime_pse > 0.2 && ebeam_epics/1000 > 0) {
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

	  double nu;
	  nu = q.E();


	  TVector3 Pprime_vec = Pprime.Vect();
	  TVector3 Pprime_unitvec = Pprime_vec.Unit();

          double w = (HCAL_vector - target_vec).Dot(HCAL_unitvector_z) / (Pprime_unitvec.Dot(HCAL_unitvector_z));
	  TVector3 w_vec = target_vec + w*Pprime_unitvec;
	  TVector3 D_vec = w_vec - HCAL_vector;

	  double Hcaly_exp = D_vec.Dot(HCAL_unitvector_y);
	  double Hcalx_exp = D_vec.Dot(HCAL_unitvector_x);

	  double dx = Hcalx - Hcalx_exp;
	  double dy = Hcaly - Hcaly_exp;

	  

	  double W2 = Pprime.Mag2();
	  double Q2 = -q.Mag2();

	  if (e_over_p >= 0.85 && e_over_p <= 1.15 && Q2 >= 2.0 && W2 < 2.0 && TMath::Abs(dy) < 0.5) {
	    h_dx_dy_cut->Fill(dy,dx);
	    h_dx_cut->Fill(dx);
	  }

	  h_dx_dy->Fill(dy,dx);
	  h_dx->Fill(dx);
	  h_W2_dy->Fill(dy,W2);
      

	  int percent = static_cast<int>(100.0 * i / Nentries);
	  if (percent != last_percent) {
	    std::cout << "\rProcessing: " << percent << "% completed" << std::flush;
	    last_percent = percent;
	  }
      }

      
    }
    std::cout << "\rProcessing: 100% completed" << std::endl;

    // Draw Histograms

    TCanvas* c1 = new TCanvas("c1", "dx vs dy", 900, 500);
    c1->Divide(2,1);
    c1->cd(1); h_dx_dy->Draw("COLZ");
    c1->cd(2); h_dx_dy_cut->Draw("COLZ");

    TCanvas* c2 = new TCanvas("c2", "dx vs dy cut", 600, 500);
    h_W2_dy->Draw("COLZ");

    TCanvas* c3 = new TCanvas("c3", "dx", 900, 500);
    c3->Divide(2,1);
    c3->cd(1); h_dx->Draw();
    c3->cd(2); h_dx_cut->Draw();
    
}
