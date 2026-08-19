#include "../../include/configParser.C"
#include "../../include/computeKineVariables.C"

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
#include "TTreeFormula.h"

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

void data_trimming(const std::string& config_filename){
  
  readConfig(config_filename);

  TString output_rootFileName = getConfigString("output_filename");
  TString output_rootDir = getConfigString("output_dir");
  TString input_rootDir = getConfigString("input_dir");
  TString globalCut = getConfigString("global_cut");

  double beam_energy = getConfigDouble("ebeam");
  TString target = getConfigString("target");
  double hcal_angle = getConfigDouble("hcal_angle");
  double hcal_distance = getConfigDouble("hcal_distance");

  TChain C("T");

  C.Add(input_rootDir+"*.root");

  C.SetBranchStatus("*", 0);

  C.SetBranchStatus("bb.sh.e",1);
  C.SetBranchStatus("bb.sh.atimeblk",1);
  C.SetBranchStatus("bb.sh.nclus",1);
  C.SetBranchStatus("bb.sh.x",1);
  C.SetBranchStatus("bb.sh.y",1);
  C.SetBranchStatus("bb.sh.colblk",1);
  C.SetBranchStatus("bb.sh.rowblk",1);
  C.SetBranchStatus("bb.sh.idblk",1);
  C.SetBranchStatus("bb.sh.nblk",1);
  C.SetBranchStatus("bb.sh.clus_blk.*");
  
  C.SetBranchStatus("bb.ps.e",1);
  C.SetBranchStatus("bb.ps.atimeblk",1);
  C.SetBranchStatus("bb.ps.nclus",1);
  C.SetBranchStatus("bb.ps.x",1);
  C.SetBranchStatus("bb.ps.y",1);
  C.SetBranchStatus("bb.ps.colblk",1);
  C.SetBranchStatus("bb.ps.rowblk",1);
  C.SetBranchStatus("bb.ps.idblk",1);
  C.SetBranchStatus("bb.ps.nblk",1);
  C.SetBranchStatus("bb.ps.clus_blk.*");
  
  C.SetBranchStatus("sbs.hcal.e",1);
  C.SetBranchStatus("sbs.hcal.atimeblk",1);
  C.SetBranchStatus("sbs.hcal.nclus",1);
  C.SetBranchStatus("sbs.hcal.x",1);
  C.SetBranchStatus("sbs.hcal.y",1);
  C.SetBranchStatus("sbs.hcal.colblk",1);
  C.SetBranchStatus("sbs.hcal.rowblk",1);
  C.SetBranchStatus("sbs.hcal.idblk",1);
  C.SetBranchStatus("sbs.hcal.nblk",1);
  C.SetBranchStatus("sbs.hcal.clus_blk.*");
  C.SetBranchStatus("sbs.hcal.clus.*");
  C.SetBranchStatus("sbs.hcal.goodblock.*",1);
  
  C.SetBranchStatus("bb.tr.v*",1);
  C.SetBranchStatus("bb.tr.p*",1);
  C.SetBranchStatus("bb.etot_over_p",1);
  
  C.SetBranchStatus("bb.hodotdc.clus.bar.tdc.*",1);
  C.SetBranchStatus("bb.hodotdc.clus.tfinal",1);
  C.SetBranchStatus("bb.hodotdc.clus.tmean",1);
  C.SetBranchStatus("bb.hodotdc.clus.tmeanRFcorr",1);
  C.SetBranchStatus("bb.hodotdc.clus.xmean",1);
  C.SetBranchStatus("bb.hodotdc.clus.ymean",1);
  C.SetBranchStatus("bb.hodotdc.clus.id",1);

  C.SetBranchStatus("e.kine.*",1);
  
  TFile output_rootFile(output_rootDir+output_rootFileName, "recreate");
  TTree *output_rootTree = C.CloneTree(0);
  output_rootTree->SetAutoFlush(0);

  double sbs_hcal_dx, sbs_hcal_dy, sbs_hcal_x_exp, sbs_hcal_y_exp;

  output_rootTree->Branch("sbs.hcal.dx", &sbs_hcal_dx, "sbs.hcal.dx/D");
  output_rootTree->Branch("sbs.hcal.dy", &sbs_hcal_dy, "sbs.hcal.dy/D");
  output_rootTree->Branch("sbs.hcal.x_exp", &sbs_hcal_x_exp, "sbs.hcal.x_exp/D");
  output_rootTree->Branch("sbs.hcal.y_exp", &sbs_hcal_y_exp, "sbs.hcal.y_exp/D");

  double sbs_hcal_x, sbs_hcal_y;
  C.SetBranchAddress("sbs.hcal.x", &sbs_hcal_x);
  C.SetBranchAddress("sbs.hcal.y", &sbs_hcal_y);

  double bb_tr_px[100], bb_tr_py[100], bb_tr_pz[100], bb_tr_p[100], bb_tr_vx[100], bb_tr_vy[100], bb_tr_vz[100];
  C.SetBranchAddress("bb.tr.px", bb_tr_px);
  C.SetBranchAddress("bb.tr.py", bb_tr_py);
  C.SetBranchAddress("bb.tr.pz", bb_tr_pz);
  C.SetBranchAddress("bb.tr.vx", bb_tr_vx);
  C.SetBranchAddress("bb.tr.vy", bb_tr_vy);
  C.SetBranchAddress("bb.tr.vz", bb_tr_vz);

  TTreeFormula globalCut_expression("cut", globalCut, &C);

  std::cout << "Starting Trimming Script..." << std::endl;
  std::cout << "output path: " << output_rootDir << std::endl;
  std::cout << "filename: " << output_rootFileName << std::endl;

  std::cout << "Getting number of events to analyze..." << std::endl;
  Long64_t totEntries = C.GetEntries();
  std::cout << "Total Events: " << totEntries << std::endl;

  std::cout << std::endl;
  int currentTree = -1;
  Long64_t finalEntries = 0;
  for (Long64_t event = 0; event < totEntries; event++) {

    Long64_t entryLoading = C.GetEntry(event);
    if (entryLoading <= 0) break;

    if (C.GetTreeNumber() != currentTree) {
      currentTree = C.GetTreeNumber();
      globalCut_expression.UpdateFormulaLeaves();
    }

    if (event % 50000 == 0) {
      double percent = event * 100.0 / totEntries;
      std::cout << "\rProgress: " << std::fixed << std::setprecision(3)
		<< percent << "%"<< std::flush;
    }

    if (globalCut_expression.EvalInstance() == 0) continue;

    TVector3 kf(bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0]);
    TVector3 v(bb_tr_vx[0], bb_tr_vy[0], bb_tr_vz[0]);

    std::vector<double> dxdy = computeDxDy(target,
					   beam_energy,
					   hcal_angle,
					   hcal_distance,
					   kf,
					   v,
					   sbs_hcal_x,
					   sbs_hcal_y);

    sbs_hcal_dx = dxdy[0];
    sbs_hcal_dy = dxdy[1];

    sbs_hcal_x_exp = sbs_hcal_x - sbs_hcal_dx;
    sbs_hcal_y_exp = sbs_hcal_y - sbs_hcal_dy;
    
    output_rootTree->Fill();
    finalEntries++;
    
  }

  std::cout << std::endl;

  output_rootTree->Write();
  output_rootFile.Close();

  std::cout << "Trimmed rootfile created!" << std::endl;
  std::cout << "Events Passed: " << finalEntries << "/" << totEntries << " events" << std::endl;
  
}
