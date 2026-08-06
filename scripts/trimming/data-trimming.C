#include "../../include/configParser.C"

#include <ROOT/RDataFrame.hxx>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <TObjArray.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TString.h>
#include <TRegexp.h>
#include <TMath.h>
#include "TChainElement.h"

#include <string>
#include <vector>
#include <iostream>

void raw_trimming_pass3_3(std::string config_filename){

  std::string config_filepath = "/work/halla/sbs/koeneman/GEnII/GEnII_analysis/config/";
  std::string config_file = config_filepath + config_filename;
  
  
  readConfig(config_file);

  TString output_rootFileName = getConfigString("output_filename");
  TString output_rootDir = getConfigString("output_dir");
  TString input_rootDir = getConfigString("input_dir");
  TString globalCut = getConfigString("global_cut");

  TChain C("T");

  C.Add(input_rootDir+"*.root");

  C.SetBranchStatus("*", 0);

  C.SetBranchStatus("bb.sh.*",1);
  C.SetBranchStatus("bb.ps.*",1);
  C.SetBranchStatus("sbs.hcal.*",1);
  C.SetBranchStatus("bb.tr.*",1);
  C.SetBranchStatus("bb.hodotdc.*",1);
  C.SetBranchStatus("e.kine.W2.*",1);
  
  TFile output_rootFile(output_rootDir+output_rootFileName, "recreate");
  TTree *output_rootTree = C.CloneTree(0);

  TTreeFormula globalCut_expression("cut", globalCut, &C);

  Long64_t event = 0;
  int currentTree = -1;
  while (C.GetEntry(event)) {

    Long64_t entryInCurrentTree = C.LoadTree(event);
    if (entryInCurrentTree < 0) break;

    if (C.GetTreeNumber() != currentTree) {
      currentTree = C.GetTreeNumber();
      globalCut_expression.UpdateFormulaLeaves();
    }

    event++;

    if (globalCut_expression.EvalInstance() == 0) continue;

    output_rootTree->Fill();
    
  }

  output_rootTree->Write();
  output_rootFile.Close();
  
}
