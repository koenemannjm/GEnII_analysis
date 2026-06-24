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

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

std::vector<std::string> GetBranchNames(TChain &chain, const std::vector<std::string> &variables){
  std::vector<std::string> result_branches;

  TTree *tree = chain.GetTree();
  if(!tree){
    int dummy = 0;
    chain.GetEntry(dummy);
    tree = chain.GetTree();
  }
  if(!tree){
    std::cerr << "ERROR IN FINDING BRANCHES FROM TCHAIN!" << "\n";
    std::string entry = "nope";
    result_branches.push_back(entry);
    return result_branches;
  }

  TObjArray *branches = tree->GetListOfBranches();
  Long64_t number_brances = branches->GetEntriesFast();
  int number_variables = variables.size();

  for(int i = 0; i < number_variables; ++i){

    const std::string pat = variables[i];
    const std::string prefix = pat.substr(0, pat.find('*'));
    

    for(int j = 0; j < number_brances; ++j){
      auto *br = static_cast<TBranch*>(branches->At(j));

      if(!br) continue;
      const char *cname = br->GetName();
      if (!cname) continue;

      const std::string name = cname;
      
      if(name.compare(0, prefix.size(), prefix) == 0){
	result_branches.push_back(name);
      }
    }
  }

  return result_branches;
}

void raw_trimming_pass3(const std::string &filepath,const std::string &outfile_name){
  TChain C("T");
  C.Add(filepath.c_str());

  std::cout << "Starting Files: " << C.GetListOfFiles()->GetEntries() << std::endl;

  TChain C_good("T");

  TObjArray* list_of_files = C.GetListOfFiles();
  int n_files = list_of_files->GetEntries();

  for (int i = 0; i < n_files; i++) {

    TChainElement* element = (TChainElement*)list_of_files->At(i);

    TString fpath = element->GetTitle();

    TFile* f = TFile::Open(fpath);
    if (!f || f->IsZombie()) {
      continue;
    }

    TTree* t = (TTree*)f->Get("T");
    if (t && t->GetBranch("T")) {
      C_good.Add(fpath);
      std::cout << "Added: " << fpath << std::endl; 
    }

    f->Close();
    
  }

  std::cout << "Good Files: " << C_good.GetListOfFiles()->GetEntries() << std::endl;
  
  std::string outdir_path = "/work/halla/sbs/koeneman/GEnII/outdir/outfiles/trimmed/";

  std::string outfile_path = outdir_path + outfile_name;

  const char* cut = "bb.tr.n>0&&bb.sh.e>0&&bb.ps.e>0&&sbs.hcal.e>0";

  std::vector<std::string> variables = {
    "bb.ps.*",
    "Ndata.bb.ps*",
    "bb.sh.*",
    "Ndata.bb.sh*",
    "bb.tr.*",
    "Ndata.bb.tr.*",
    "bb.etot*",
    "Ndata.bb.etot*",
    "sbs.hcal.*",
    "Ndata.sbs.hcal.*",
    "g.*",
    "bb.hodo*",
    "Ndata.bb.hodo*",
    "bb.grinch_tdc*",
    "Ndata.bb.grinch_tdc*",
    "e.kine.*"
  };

  auto branch_variables = GetBranchNames(C_good, variables);
  
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df(C_good);
  auto df_f = df.Filter(cut, "global cuts");

  ROOT::RDF::RSnapshotOptions opts;
  opts.fCompressionAlgorithm = ROOT::kZLIB;
  opts.fCompressionLevel = 6;
  
  df_f.Snapshot("TrimmedFile",outfile_name.c_str(), branch_variables, opts);
  
}
