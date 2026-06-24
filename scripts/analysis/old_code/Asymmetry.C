//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//
//   Created by Jacob Koenemann & ChatGPT
//   contact: bxy3zr@virginia.edu
//                                                                     
//   Last Modified April 26, 2025   
//
//      
//   The purpose of this script to compute the full
//   physics Asymmetry for the quasi-elastic neutron
//   and extract GEn/GMn and possibly GEn.
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

void Asymmetry(std::string config_file){

    // Reading in config file
    std::string config_path = "../../config/" + config_file;
    std::ifstream file(config_path);
    if (!file) {
       std::cerr << "Error: Config file does not exist at " << config_path << std::endl;
       std::exit(EXIT_FAILURE);
    }
    
    std::ifstream config(config_path);
    std::string line;
    std::map<std::string, std::string> settings;

    while (std::getline(config, line)) {
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
    


}
