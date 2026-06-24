//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
//   Created by Jacob Koenemann & ChatGPT
//   contact: bxy3zr@virginia.edu
//                                                                     
//   Last Modified April 26, 2025   
//
//      
//   The purpose of this script is to compute the
//   cointime between the bbcal and hcal arms. This
//   then adds
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

void Computecointime(std::string config, std::string target, std::string pass){
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

    // Set experiment name
    std::string exp_name;
    if (config == "kin2") {
        exp_name = "GEN2";
    }
    else if (config == "kin3") {
         exp_name = "GEN3";
    }
    else if (config == "kin4a") {
         exp_name = "GEN4";
    }
    else if (config == "kin4b") {
         exp_name = "GEN4b";
    }
    
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
    
    // BigBite variables
    double bbcal_time;
    tree->SetBranchAddress("bb.sh.atimeblk", &bbcal_time);
   
    // Hcal variables
    double Hcal_time;
    tree->SetBranchAddress("sbs.hcal.atimeblk", &Hcal_time);

    double cointime;
    TBranch* atimeblk_cointime = tree->Branch("atimeblk.cointime", &cointime, "atimeblk.cointime/D");

    //Loop over entries of the tree
    Long64_t Nentries = tree->GetEntries();
    int last_percent = -1;
    for (Long64_t i = 0; i < Nentries; ++i) {
      tree->GetEntry(i);
      
      cointime = Hcal_time - bbcal_time;

      atimeblk_cointime->Fill();

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
