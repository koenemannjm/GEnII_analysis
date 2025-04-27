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

void compute_cointime(std::string config, std::string target, std::string pass){
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
    double Q2_cut, W2_cut;

    // Set experiment name
    std::string exp_name;
    if (config == "kin2") {
        exp_name = "GEN2";
	ebeam = 4.291;
	Q2_cut = 2.0;
	W2_cut = 2.0;
    }
    else if (config == "kin3") {
         exp_name = "GEN3";
	 ebeam = 6.373;
	 Q2_cut = 5.2;
	 W2_cut = 2.0;
    }
    else if (config == "kin4a") {
         exp_name = "GEN4";
	 ebeam = 8.448;
	 Q2_cut = 8.0;
	 W2_cut = 2.0;
    }
    else if (config == "kin4b") {
         exp_name = "GEN4b";
	 ebeam = 8.448;
	 Q2_cut = 8.0;
	 W2_cut = 2.0;
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
    double keprime_x, keprime_y, keprime_z, eprime_she, eprime_pse, e_over_p, target_x, target_y, target_z, bbcal_time;

    tree->SetBranchAddress("bb.tr.px", &keprime_x);
    tree->SetBranchAddress("bb.tr.py", &keprime_y);
    tree->SetBranchAddress("bb.tr.pz", &keprime_z);
    tree->SetBranchAddress("bb.ps.e", &eprime_she);
    tree->SetBranchAddress("bb.sh.e", &eprime_pse);
    tree->SetBranchAddress("bb.tr.vx", &target_x);
    tree->SetBranchAddress("bb.tr.vy", &target_y);
    tree->SetBranchAddress("bb.tr.vz", &target_z);
    tree->SetBranchAddress("bb.etot_over_p", &e_over_p);
    tree->SetBranchAddress("bb.sh.atimeblk", &bbcal_time);
   

    // Hcal variables
    double Hcalx, Hcaly, Hcale, dx, dy, Hcal_time;

    tree->SetBranchAddress("sbs.hcal.x", &Hcalx);
    tree->SetBranchAddress("sbs.hcal.y", &Hcaly);
    tree->SetBranchAddress("sbs.hcal.e", &Hcale);
    tree->SetBranchAddress("dx",&dx);
    tree->SetBranchAddress("dy",&dy);
    tree->SetBranchAddress("sbs.hcal.atimeblk", &Hcal_time);

    // Pre-Defining Histograms for Calculations results
    TH1D* h_cointime = new TH1D("h_cointime", "coin.time = sbs.hcal.atimeblk - bb.sh.atimeblk ;coin.time (ns);Counts", 500, -10.0, 220.0);
    TH1D* h_cointime_cut = new TH1D("h_cointime_cut", "coin.time = sbs.hcal.atimeblk - bb.sh.atimeblk ;coin.time (ns);Counts", 300, 50.0, 150.0);

    //Loop over entries of the tree
    Long64_t Nentries = tree->GetEntries();
    int last_percent = -1;
    for (Long64_t i = 0; i < Nentries; ++i) {
      tree->GetEntry(i);
      if (TMath::Abs(target_z) <= 0.27 && TMath::Abs(target_x) <= 0.0115 && TMath::Abs(target_y) <= 0.0115 && Hcale > 0.025 && eprime_pse > 0.2) {
	  // Defining momentum 3-vector
	  TVector3 keprime_vec(keprime_x,keprime_y,keprime_z);

	  double keprime_mag = keprime_vec.Mag();
	  double etheta = keprime_vec.Theta();
	  double ephi = keprime_vec.Phi();
	  double eprime = keprime_mag;

	  TLorentzVector keprime(keprime_x,keprime_y,keprime_z,eprime);
	  TLorentzVector ke(0.0,0.0,ebeam,ebeam);
	  TLorentzVector P(0.0,0.0,0.0,MN);
	  TLorentzVector Pprime;
	  TLorentzVector q;
      
	  q = ke - keprime;
	  Pprime = q + P;

	  double W2 = Pprime.Mag2();
	  double Q2 = -q.Mag2();

	  double cointime = Hcal_time - bbcal_time;

	  if (e_over_p >= 0.85 && e_over_p <= 1.15 && Q2 >= Q2_cut && W2 <= W2_cut && TMath::Abs(dy)<0.5) {
	    h_cointime_cut->Fill(cointime);
	  }
      

	  h_cointime->Fill(cointime);

	  int percent = static_cast<int>(100.0 * i / Nentries);
	  if (percent != last_percent) {
	    std::cout << "\rProcessing: " << percent << "% completed" << std::flush;
	    last_percent = percent;
	  }
      }

      
    }
    std::cout << "\rProcessing: 100% completed" << std::endl;

    // Fitting cointime peak to get mean and sigma

    int peakcointimeBin = h_cointime_cut->GetMaximumBin();
    double peakcointime = h_cointime_cut->GetXaxis()->GetBinCenter(peakcointimeBin);

    double fitcointime_left = peakcointime - 2.0;
    double fitcointime_right = peakcointime + 2.0;

    TF1 *fitFunc = new TF1("fitFunc", "gaus",  fitcointime_left,  fitcointime_right);

    h_cointime_cut->Fit(fitFunc, "R");

    // Draw Histograms
    TCanvas* c1 = new TCanvas("c1", "Track Information", 1200, 600);
    c1->Divide(2,1);
    c1->cd(1);
    h_cointime_cut->SetLineColor(kBlue);
    h_cointime_cut->SetLineStyle(2);
    h_cointime_cut->SetLineWidth(2);
    h_cointime_cut->Draw("E");
    fitFunc->Draw("SAME");
    c1->cd(2);
    h_cointime->SetLineColor(kBlack);
    h_cointime->SetLineWidth(2);
    h_cointime->Draw();
      
}
