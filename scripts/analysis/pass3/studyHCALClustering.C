#include "../../../include/configParser.C"
//#include "gen_SBSOFF.C"

void studyHCALClustering(const std::string& config_filename){
  
  readConfig(config_filename);

  double beam_energy = getConfigDouble("ebeam");
  TString target = getConfigString("target");
  TString rootFile = getConfigString("output_filename");
  TString rootDir = getConfigString("output_dir");
  TString globalCut = getConfigString("global_cut");

  TString rootPath = rootDir + rootFile;

  std::cout << "ROOT file: " << rootFile << std::endl;
  std::cout << "gloablCut: " << globalCut << std::endl;

  TChain C("T");

  C.Add(rootPath);

  TTreeFormula cutFormula("cutFormula", globalCut.Data(), &C);

  Long64_t event = 0;
  while (C.GetEntry(event) > 0) {

    C.LoadTree(event);
    event++;

    if (cutFormula.EvalInstance() == 0) continue;




  }

  
}
