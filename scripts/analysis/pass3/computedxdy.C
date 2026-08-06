#include "../../../include/configParser.C"
//#include "gen_SBSOFF.C"
//#include "gen_SBSON.C"

void computerdxdy(const std::string& config_filename) {
  
  readConfig(config_filename);

  double beam_energy = getConfigDouble("ebeam");
  TString target = getConfigString("target");
  TString rootFile = getConfigString("output_filename");
  TString rootDir = getConfigString("output_dir");
  TString globalCut = getConfigString("global_cut");

}