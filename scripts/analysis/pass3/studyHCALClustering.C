#include "../../../include/configParser.C"
//#include "gen_SBSOFF.C"
//#include "gen_SBSON.C"

void studyHCALClustering(const std::string& config_filename) {
  
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

    double bbsh_atimeblk = T->bb_sh_atimeblk;

    // Filling the bbsh_atimeblk - hcal_adctimei vs hcal cluster index
    int number_hcal_clusters = T->Ndata_sbs_hcal_clus_adctime;
    double hcal_adctimei, hcal_ei;
    for (int i=0; i<number_hcal_clusters; i++) {
      hcal_adctimei = T->sbs_hcal_clus_adctime[i];
      hcal_ei = T->sbs_hcal_clus_e[i];
    }

    int number_hcal_goodblocks = T->Ndata_sbs_hcal_goodblock_atime;
    double hcal_gb_adctimei, hcal_gb_ei, hcal_gb_coli, hcal_gb_rowi, hcal_gb_cidi;
    for (int i=0; i<number_hcal_goodblocks; i++) {
      hcal_gb_adctimei = T->sbs_hcal_goodblock_atime[i];
      hcal_gb_ei = T->sbs_hcal_goodblock_e[i];
      hcal_gb_coli = T->sbs_hcal_goodblock_col[i];
      hcal_gb_rowi = T->sbs_hcal_goodblock_row[i];
      hcal_gb_cidi = T->sbs_hcal_goodblock_cid[i];
    }

  }

  
}
