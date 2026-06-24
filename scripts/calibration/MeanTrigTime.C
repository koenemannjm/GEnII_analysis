#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
//#include "gen_tree.C"
#include "gen_tree_old.C"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"

#include <iostream>
#include <cstdlib>
#include <fstream>

struct KineInfo_t {
  double min_runnum;
  double max_runnum;
  double bin_runnum;
};

KineInfo_t KineInfo(const std::string& kine){
  KineInfo_t info;
  if(kine=="GEN2"){
    info.min_runnum = 1997.5;
    info.max_runnum = 2323.5;
    info.bin_runnum = 326;
  }
  if(kine=="GEN3"){
    info.min_runnum = 2463.5;
    info.max_runnum = 3265.5;
    info.bin_runnum = 802;
  }
  if(kine=="GEN4"){
    info.min_runnum = 3448.5;
    info.max_runnum = 4587.5;
    info.bin_runnum = 1139;
  }
  if(kine=="GEN4b"){
    info.min_runnum = 4982.5;
    info.max_runnum = 6086.5;
    info.bin_runnum = 1104;
  }
  return info;
}

void MeanTrigTime(std::string root_file_path, std::string fig_title, std::string kine_name){

  auto kine = KineInfo(kine_name);
  double min_runnum = kine.min_runnum;
  double max_runnum = kine.max_runnum;
  double bin_runnum = kine.bin_runnum;

  TChain *C = new TChain("T");
  C->Add(root_file_path.c_str());
  int numtrees = C->GetNtrees();
  std::cout << "Number of Trees Added: " << numtrees << std::endl;

  TCut cut = "";

  cut += "bb.ps.e>0.2&&sbs.hcal.e>0.02&&fabs(bb.tr.vz[0])<0.27&&fabs(bb.etot_over_p[0]-1.0)<0.1&&g.trigbits==4";

  TH2D *hHODOtrigtime_runnum = new TH2D("hHODOtrigtime_runnum",(fig_title + ";run number;t^{trigtime}_{HODO} (ns)").c_str(),bin_runnum,min_runnum,max_runnum, 300, 320, 380);
  TH2D *hHODOrftime_runnum = new TH2D("hHODOrftime_runnum",(fig_title + "; run number;t^{rftime}_{HODO} (ns)").c_str(),bin_runnum,min_runnum,max_runnum,300,-105,105);
  
  TGraphErrors *gmeanHODOtrigtime_runnum = new TGraphErrors(bin_runnum);
  TGraphErrors *gmeanHODOrftime_runnum = new TGraphErrors(bin_runnum);

  //gen_tree *T = new gen_tree(C);
  gen_tree_old *T = new gen_tree_old(C);

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("e.kine.W2",1);
  C->SetBranchStatus("bb.*",1);
  C->SetBranchStatus("sbs.hcal.*",1);
  C->SetBranchStatus("bb.hodotdc.*",1);
  C->SetBranchStatus("g.*",1);

   TTreeFormula *cutFormula = new TTreeFormula("cut",cut,C);

  int nevent = 0;
  int treenum = -1;
  int oldtreenum = -1;
  while(T->GetEntry(nevent)){
    if(nevent % 100000 == 0){
      std::cout << "Event number: " << nevent << std::endl;
    }

    treenum = C->GetTreeNumber();
    if( nevent == 0 || treenum != oldtreenum ){
      oldtreenum = treenum;
      cutFormula->UpdateFormulaLeaves();
    }
    nevent++;

    if(cutFormula->EvalInstance(0)==0) continue;

    double hodotrigtime = T->bb_hodotdc_trigtime;
    double hodorftime = T->bb_hodotdc_rftime;
    double runnum = T->g_runnum;

    hHODOtrigtime_runnum->Fill(runnum,hodotrigtime);
    hHODOrftime_runnum->Fill(runnum,hodorftime);

  }
  
  for (int xi = 1; xi<=bin_runnum; ++xi){

    TH1D *projTT = hHODOtrigtime_runnum->ProjectionY(Form("Runnum=%d",xi), xi, xi);
    TH1D *projRF = hHODOrftime_runnum->ProjectionY(Form("Runnum=%d",xi), xi, xi);

    double meanTT = 0.0;
    double meanRF = 0.0;
    
    double stdTT = 1.0;
    double stdRF = 1.0;
    if(projTT->GetEntries() > 0){
      meanTT = projTT->GetMean();
      //stdY = projY->GetStdDev()/sqrt(nentries);
      //stdY = projY->GetStdDev();
    }
    if(projRF->GetEntries() > 0){
      meanRF = projRF->GetMean();
    }

    double x = hHODOtrigtime_runnum->GetXaxis()->GetBinCenter(xi);

    gmeanHODOtrigtime_runnum->SetPoint(xi-1, x, meanTT);
    gmeanHODOtrigtime_runnum->SetPointError(xi-1, 0.0, stdTT);

    gmeanHODOrftime_runnum->SetPoint(xi-1, x, meanRF);
    gmeanHODOrftime_runnum->SetPointError(xi-1, 0.0, stdRF);

    projTT->SetDirectory(0);
    projRF->SetDirectory(0);
    
    delete projTT;
    delete projRF;
  }

  gmeanHODOtrigtime_runnum->SetMarkerStyle(20);
  gmeanHODOtrigtime_runnum->SetMarkerColor(kRed);
  gmeanHODOtrigtime_runnum->SetLineColor(kRed);
  gmeanHODOtrigtime_runnum->SetTitle("Mean HODO trigger ref. time vs Run number; run number; t^{trigtime}_{HODO} (ns)");

  gmeanHODOrftime_runnum->SetMarkerStyle(20);
  gmeanHODOrftime_runnum->SetMarkerColor(kRed);
  gmeanHODOrftime_runnum->SetLineColor(kRed);
  gmeanHODOrftime_runnum->SetTitle("Mean HODO RF time vs Run number; run number; t^{rftime}_{HODO} (ns)");

  std::string base_dir = "/work/halla/sbs/koeneman/GEnII/";
  std::string out_dir = base_dir + "outdir/figures/MeanTrigTime/";
  std::string out_hist_name = "MeanTrigTime_" + fig_title + ".pdf";
  std::string out_hist_path = out_dir + out_hist_name;
  
  TCanvas *c = new TCanvas("c","c",800,600);

  hHODOtrigtime_runnum->Draw("colz");
  c->Print((out_hist_path + "(").c_str());

  c->Clear();
  hHODOrftime_runnum->Draw("colz");
  c->Print(out_hist_path.c_str());
  
  c->Clear();
  gmeanHODOtrigtime_runnum->Draw("AP");
  c->Print(out_hist_path.c_str());

  c->Clear();
  gmeanHODOrftime_runnum->Draw("AP");
  c->Print(out_hist_path.c_str());

  c->Clear();
  hHODOtrigtime_runnum->Draw("colz");
  gmeanHODOtrigtime_runnum->Draw("P SAME");
  c->Print(out_hist_path.c_str());

  c->Clear();
  hHODOrftime_runnum->Draw("colz");
  gmeanHODOrftime_runnum->Draw("P SAME");
  c->Print((out_hist_path+")").c_str());

  delete c;
  delete hHODOtrigtime_runnum;
  delete gmeanHODOtrigtime_runnum;
}
