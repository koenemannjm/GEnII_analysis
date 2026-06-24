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


void data_trimming(const std::string& config_filename) {

  readConfig(config_filename);

  TString output_rootFile = getConfigString("output_filename");
  TString output_rootDir = getConfigString("output_dir");
  TString globalCut = getConfigString("global_cut");

  std::vector<std::string> listoffiles;
  void* dir = gSystem->OpenDirectory(file_path.c_str());
  const char* entry;
  while ((entry = gSystem->GetDirEntry(dir))) {
    std::string name = entry;
    if (name.find(".root") != std::string::npos) {
      listoffiles.push_back(file_path + "/" + name);
    }
  }

  gSystem->FreeDirectory(dir);

  std::cout << "Found " << listoffiles.size() << " ROOT files" << std::endl;

  ROOT::RDataFrame df("T", listoffiles);

  auto filter = df.Filter(globalCut);

  //filter.Snapshot("T", root_file_name.c_str());

  std::cout << "\nDone! Output file: " << root_file_name << std::endl;
  std::cout << "Raw number of events: " << df.Count().GetValue() << std::endl;
  std::cout << "cut: " << cut << std::endl;
  std::cout << "Output number of events: " << filter.Count().GetValue() << std::endl;
}
