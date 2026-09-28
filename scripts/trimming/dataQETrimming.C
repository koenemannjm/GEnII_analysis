#include "../../include/configParser.C"

#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TTreeFormula.h"

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

void dataQETrimming(const std::string& config_filename) {

  readConfig(config_filename);

  TString expConfig = getConfigString("config");
  TString rootDir = getConfigString("output_dir");
  TString QErootFile = getConfigString("QE_output_filename");
  TString QEecut = getConfigString("goode_cut");
  
  TString inputFiles = "*"+expConfig+"*.root";
  TString inputPath = rootDir+inputFiles;
  TString outputPath = rootDir+QErootFile;

  TChain C("T");
  C.Add(inputPath.Data());

  TFile f(outputPath, "recreate");
  TTree *T = C.CloneTree(0);

  TTreeFormula globalCut_expression("cut", QEecut, &C);

  std::cout << "Starting Trimming Script..." << std::endl;
  std::cout << "output path: " << rootDir << std::endl;
  std::cout << "filename: " << QErootFile << std::endl;

  std::cout << "Getting number of events to analyze..." << std::endl;
  Long64_t totEntries = C.GetEntries();
  std::cout << "Total Events: " << totEntries << std::endl;

  std::cout << std::endl;
  int currentTree = -1;
  Long64_t finalEntries = 0;
  for (Long64_t event = 0; event < totEntries; event++) {

    Long64_t local_entry = C.LoadTree(event);
    if (local_entry < 0) break;

    if (C.GetTreeNumber() != currentTree) {
      currentTree = C.GetTreeNumber();
      globalCut_expression.UpdateFormulaLeaves();
    }

    Long64_t entryLoading = C.GetEntry(event);
    if (entryLoading <= 0) break;

    if (event % 50000 == 0) {
      double percent = event * 100.0 / totEntries;
      std::cout << "\rProgress: " << std::fixed << std::setprecision(3)
		<< percent << "%"<< std::flush;
    }

    if (globalCut_expression.EvalInstance() == 0) continue;
    T->Fill();
    finalEntries++;
    
  }

  std::cout << std::endl;

  T->Write();
  f.Close();
  
  std::cout << "Trimmed QE rootfile created!" << std::endl;
  std::cout << "Events Passed: " << finalEntries << "/" << totEntries << " events" << std::endl;

  clearConfig();
  delete T;
}
