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

void compute_variables(std::string config, std::string target, std::string pass){
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

    // Set experiment name
    std::string exp_name;
    if (config == "kin2") {
        exp_name = "GEN2";
	ebeam = 4.30;
    }
    else if (config == "kin3") {
         exp_name = "GEN3";
	 ebeam = 6.40;
    }
    else if (config == "kin4a") {
         exp_name = "GEN4";
	 ebeam = 8.50;
    }
    else if (config == "kin4b") {
         exp_name = "GEN4b";
	 ebeam = 8.50;
    }
    
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
    double keprime_x, keprime_y, keprime_z, eprime_she, eprime_pse, e_over_p, tr_vz;

    tree->SetBranchAddress("bb.tr.px", &keprime_x);
    tree->SetBranchAddress("bb.tr.py", &keprime_y);
    tree->SetBranchAddress("bb.tr.pz", &keprime_z);
    tree->SetBranchAddress("bb.ps.e", &eprime_she);
    tree->SetBranchAddress("bb.sh.e", &eprime_pse);
    tree->SetBranchAddress("bb.tr.vz", &tr_vz);
    tree->SetBranchAddress("bb.etot_over_p", &e_over_p);

    // Hcal variables
    double Hcalx, Hcaly, Hcale;

    tree->SetBranchAddress("sbs.hcal.x", &Hcalx);
    tree->SetBranchAddress("sbs.hcal.y", &Hcaly);
    tree->SetBranchAddress("sbs.hcal.e", &Hcale);

    // Pre-Defining Histograms for Calculations results
    TH1D* h_keprime_mag = new TH1D("h_keprime_mag", "Track Momentum p;p;Counts", 100, 0.5, 3.5);
    TH1D* h_etheta = new TH1D("h_etheta", "etheta ;etheta;Counts", 100, 0.4, 0.7);
    TH1D* h_ephi = new TH1D("h_ephi", "ephi;ephi;Counts", 100, -0.6, 0.6);
    TH1D* h_W2 = new TH1D("h_W2", "W2;W2;Counts", 100, -6.0, 13.0);
    TH1D* h_Q2 = new TH1D("h_Q2", "Q2;Q2;Counts", 100, -6.0, 9.0);
    TH1D* h_W2_cut = new TH1D("h_W2_cut", "W2_cut;W2;Counts", 100, -6.0, 13.0);
    TH1D* h_Q2_cut = new TH1D("h_Q2_cut", "Q2_cut;Q2;Counts", 100, -6.0, 9.0);
    TH1D* h_W2_tr = new TH1D("h_W2_tr", "W2_tr;W2;Counts", 100, -6.0, 13.0);
    TH1D* h_Q2_tr = new TH1D("h_Q2_tr", "Q2_tr;Q2;Counts", 100, -6.0, 9.0);
    TH1D* h_W2_tr_cut = new TH1D("h_W2_tr_cut", "W2_tr;W2;Counts", 100, -6.0, 13.0);
    TH1D* h_Q2_tr_cut = new TH1D("h_Q2_tr_cut", "Q2_tr;Q2;Counts", 100, -6.0, 9.0);
    TH2D* h_mag_etheta = new TH2D("h_mag_etheta", "p vs etheta;etheta,p", 100, 0.4, 0.7, 100, 0.5, 3.5);
    

    //Loop over entries of the tree
    Long64_t Nentries = tree->GetEntries();
    int last_percent = -1;
    for (Long64_t i = 0; i < Nentries; ++i) {
      tree->GetEntry(i);
      if (abs(tr_vz) <= 0.27 && Hcale > 0.025 && eprime_pse > 0.2 && ebeam_epics/1000 > 0) {
	  // Defining momentum 3-vector
	  TVector3 keprime_vec(keprime_x,keprime_y,keprime_z);

	  double keprime_mag = keprime_vec.Mag();
	  double etheta = keprime_vec.Theta();
	  double ephi = keprime_vec.Phi();
	  double eprime = eprime_she + eprime_pse;

	  TLorentzVector keprime(keprime_x,keprime_y,keprime_z,eprime);
	  TLorentzVector keprime_tr(keprime_x,keprime_y,keprime_z,keprime_vec.Mag());
	  TLorentzVector ke(0.0,0.0,ebeam_epics/1000,ebeam_epics/1000);
	  TLorentzVector P(0.0,0.0,0.0,MN);
	  TLorentzVector Pprime;
	  TLorentzVector q;
	  TLorentzVector Pprime_tr;
	  TLorentzVector q_tr;
      
	  q = ke - keprime;
	  Pprime = q + P;

	  q_tr = ke - keprime_tr;
	  Pprime_tr = q_tr + P;

	  double W2 = Pprime.Mag2();
	  double Q2 = -q.Mag2();

	  double W2_tr = Pprime_tr.Mag2();
	  double Q2_tr = -q_tr.Mag2();

	  if (e_over_p >= 0.85 && e_over_p <= 1.15 && Q2_tr >= 8.0) {
	    h_W2_cut->Fill(W2);
	    h_Q2_cut->Fill(Q2);
	    h_W2_tr_cut->Fill(W2_tr);
	    h_Q2_tr_cut->Fill(Q2_tr);
	  }
      

	  h_keprime_mag->Fill(keprime_mag);
	  h_etheta->Fill(etheta);
	  h_ephi->Fill(ephi);
	  h_mag_etheta->Fill(etheta,eprime);
	  h_W2->Fill(W2);
	  h_Q2->Fill(Q2);
	  h_W2_tr->Fill(W2_tr);
	  h_Q2_tr->Fill(Q2_tr);

	  int percent = static_cast<int>(100.0 * i / Nentries);
	  if (percent != last_percent) {
	    std::cout << "\rProcessing: " << percent << "% completed" << std::flush;
	    last_percent = percent;
	  }
      }

      
    }
    std::cout << "\rProcessing: 100% completed" << std::endl;

    // Draw Histograms
    TCanvas* c1 = new TCanvas("c1", "Track Information", 1200, 400);
    c1->Divide(3,1);
    c1->cd(1); h_keprime_mag->Draw();
    c1->cd(2); h_etheta->Draw();
    c1->cd(3); h_ephi->Draw();

    TCanvas* c2 = new TCanvas("c2", "track p vs etheta", 600, 500);
    h_mag_etheta->Draw("COLZ");

    TCanvas* c3 = new TCanvas("c3", "W2 and Q2", 1200, 500);
    c3->Divide(2,1);
    c3->cd(1);
    h_W2_tr_cut->SetLineColor(kRed);
    h_W2_tr_cut->SetLineStyle(2);
    h_W2_tr_cut->SetLineWidth(2);
    h_W2_tr_cut->SetStats(0);
    h_W2_tr_cut->Draw();
    h_W2_cut->SetLineColor(kGreen);
    h_W2_cut->SetLineWidth(2);
    h_W2_cut->SetStats(0);
    h_W2_cut->Draw("SAME");
    h_W2_tr->SetLineColor(kOrange);
    h_W2_tr->SetLineStyle(2);
    h_W2_tr->SetLineWidth(2);
    h_W2_tr->SetStats(0);
    h_W2_tr->Draw("SAME");
    h_W2->SetLineColor(kBlue);
    h_W2->SetLineWidth(2);
    h_W2->SetStats(0);
    h_W2->Draw("SAME");
    TLegend *leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg1->AddEntry(h_W2, "h_W2", "l");
    leg1->AddEntry(h_W2_cut, "h_W2_cut", "l");
    leg1->AddEntry(h_W2_tr, "h_W2_tr", "l");
    leg1->AddEntry(h_W2_tr_cut, "h_W2_tr_cut", "l");
    leg1->Draw();
    
    c3->cd(2);
    h_Q2_tr->SetLineColor(kRed);
    h_Q2_tr->SetLineStyle(2);
    h_Q2_tr->SetLineWidth(2);
    h_Q2_tr->SetStats(0);
    h_Q2_tr->Draw();
    h_Q2_tr_cut->SetLineColor(kOrange);
    h_Q2_tr_cut->SetLineStyle(2);
    h_Q2_tr_cut->SetLineWidth(2);
    h_Q2_tr_cut->SetStats(0);
    h_Q2_tr_cut->Draw("SAME");
    h_Q2->SetLineColor(kBlue);
    h_Q2->SetLineWidth(2);
    h_Q2->SetStats(0);
    h_Q2->Draw("SAME");
    h_Q2_cut->SetLineColor(kGreen);
    h_Q2_cut->SetLineWidth(2);
    h_Q2_cut->SetStats(0);
    h_Q2_cut->Draw("SAME");
    TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg2->AddEntry(h_Q2, "h_Q2", "l");
    leg2->AddEntry(h_Q2_cut, "h_Q2_cut", "l");
    leg2->AddEntry(h_Q2_tr, "h_Q2_tr", "l");
    leg2->AddEntry(h_Q2_tr_cut, "h_Q2_tr_cut", "l");
    leg2->Draw();
}
