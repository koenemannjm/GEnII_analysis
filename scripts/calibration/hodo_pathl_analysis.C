
void hodo_analysis() {

  // Defining path to ROOT file
  std::string path, filename, filepath;
  path = "../../outfiles/rootfiles/";
  filename = "QE_pass2_try8_v4_data_GEN3_sbs100p_nucleon_np.root";
  filepath = path + filename;
  
  // Opening and verifying ROOT file
  TFile *file = TFile::Open(filepath.c_str());
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file " << filepath << std::endl;
    return;
  }
  
  std::string treename = "Tout";
  TTree *tree = (TTree*)file->Get(treename.c_str());
  if (!tree) {
    std::cerr << "Tree " << treename << "not found in file." << std::endl;
    return;
  }

}
