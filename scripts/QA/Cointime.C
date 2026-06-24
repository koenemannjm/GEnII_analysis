#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TH2.h"
#include "TH1.h"
#include "TMath.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "gen_tree.C"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <vector>

// Physical constants:
const double C_M_PER_NS = 0.299792458;
const double MP = 0.938272;
const double MN = 0.939565;


////////////////////////////////////////////////////////////
//// NEED TO UPDATE THESE PARAMETERS FOR EACH KINEMATIC ////

struct KineInfo_t {
  double E_BEAM;
  double HCAL_DIST;
  double HCAL_ANGLE;
  double min_runnum;
  double max_runnum;
  double bin_runnum;
  double dx0p;
  double dy0p;
  double dx0n;
  double dy0n;
  double sigma_dx;
  double sigma_dy;
  double nsigma_dxdy;
};

KineInfo_t KineInfo(const std::string& kine){

  KineInfo_t info;
  if(kine=="GEN2"){
    // ++++ GEN2 +++++
    info.E_BEAM = 4.291;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 34.7;
    info.min_runnum = 1997.5;
    info.max_runnum = 2323.5;
    info.bin_runnum = 326;
    info.dx0p = 0.0;
    info.dy0p = 0.0;
    info.dx0n = 0.0;
    info.dy0n = 0.0;
    info.sigma_dx = 0.3;
    info.sigma_dy = 0.3;
    info.nsigma_dxdy = 1.0;
  }
  if(kine=="GEN3"){
    // ++++ GEN3 +++++
    info.E_BEAM = 6.373;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 21.6;
    info.min_runnum = 2463.5;
    info.max_runnum = 3265.5;
    info.bin_runnum = 802;
    info.dx0p = -1.6;
    info.dy0p = -0.1;
    info.dx0n = 0.0;
    info.dy0n = -0.1;
    info.sigma_dx = 0.3;
    info.sigma_dy = 0.3;
    info.nsigma_dxdy = 1.0;
  }
  if(kine=="GEN4"){
    // ++++ GEN4 +++++
    info.E_BEAM = 8.448;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 18.0;
    info.min_runnum = 3448.5;
    info.max_runnum = 4587.5;
    info.bin_runnum = 1139;
    info.dx0p = -1.1;
    info.dy0p = -0.05;
    info.dx0n = 0.0;
    info.dy0n = -0.05;
    info.sigma_dx = 0.3;
    info.sigma_dy = 0.3;
    info.nsigma_dxdy = 1.0;
  }
  if(kine=="GEN4b"){
    // ++++ GEN4b +++++
    info.E_BEAM = 8.448;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 18.0;
    info.min_runnum = 4982.5;
    info.max_runnum = 6086.5;
    info.bin_runnum = 1104;
    info.dx0p = -1.1;
    info.dy0p = -0.1;
    info.dx0n = 0.0;
    info.dy0n = -0.1;
    info.sigma_dx = 0.3;
    info.sigma_dy = 0.3;
    info.nsigma_dxdy = 1.0;
  }
  return info;
}

TGraphErrors* ComputeMeanAndStdDev(TH2D *h2){
  int nbinsx = h2->GetNbinsX();
  TGraphErrors* g = new TGraphErrors(nbinsx);
  for (int i = 1; i<=nbinsx; ++i){
    TH1D *proj = h2->ProjectionY(Form("proj_%d",i), i, i);
    proj->SetDirectory(nullptr);
    double nentries = proj->GetEntries();
    double mean = proj->GetMean();
    double std = proj->GetStdDev();
    std /= sqrt(nentries);
    double x = h2->GetXaxis()->GetBinCenter(i);
    g->SetPoint(i-1, x, mean);
    g->SetPointError(i-1, 0, std);
    
    delete proj;
  }
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kRed);
  g->SetLineColor(kBlack);
  return g;
}

TFitResultPtr FitPeak(TH1D *h, double threshold = 0.7){
  int nbins = h->GetXaxis()->GetNbins();

  double highest_sum = -1;
  int binmax = 2;
  for(int i = 2; i < nbins; i++){
    double sum = h->GetBinContent(i-1) + h->GetBinContent(i) + h->GetBinContent(i);
    if( sum > highest_sum ){
      highest_sum = sum;
      binmax = i;
    }
  }

  int binlow = binmax-1;
  int binhigh = binmax+1;
  double peakheight = h->GetBinContent(binmax);

  while(binlow > 1 && h->GetBinContent(binlow) >= threshold * peakheight){
    binlow--;
  }
  while(binhigh < nbins && h->GetBinContent(binhigh) >= threshold * peakheight){
    binhigh++;
  }

  double xlow = h->GetBinLowEdge(binlow);
  double xhigh = h->GetBinLowEdge(binhigh);

  TFitResultPtr fr = h->Fit("gaus","SQ0","",xlow,xhigh);
  
  return fr;
}

TGraphErrors* FitMeanAndStdDev(TH2D *h2, double threshold = 0.7){
  int nbinsx = h2->GetXaxis()->GetNbins();
  TGraphErrors* g = new TGraphErrors(nbinsx);

  g->SetName(Form("g_%s_mean", h2->GetName()) );
  g->SetTitle(Form("Mean vs %s; %s; %s",
		   h2->GetXaxis()->GetTitle(),
		   h2->GetXaxis()->GetTitle(),
		   h2->GetYaxis()->GetTitle()) );
  
  for(int i = 1; i<=nbinsx; ++i){
    TH1D* proj = h2->ProjectionY("_projy", i, i);
    proj->SetDirectory(nullptr);

    int nbinsy = proj->GetNbinsX();
    int nentries = proj->GetEntries();

    if(nentries == 0){
      double x = h2->GetXaxis()->GetBinCenter(i);
      g->SetPoint(i-1, x, 0);
      g->SetPointError(i-1, 0, 0);
      delete proj;
      continue;
    }
    
    double highest_sum =-1;
    int binmax = 2;
    for(int j = 2; j < nbinsy; j++){
      double sum = proj->GetBinContent(j-1) + proj->GetBinContent(j) + proj->GetBinContent(j+1);
      if( sum > highest_sum ){
	highest_sum = sum;
	binmax = j;
      }
    }

    int binlow = binmax-1;
    int binhigh = binmax+1;
    double peakheight = proj->GetBinContent(binmax);

    while(binlow > 1 && proj->GetBinContent(binlow) >= threshold * peakheight){
      binlow--;
    }
    while(binhigh < nbinsy && proj->GetBinContent(binhigh) >= threshold * peakheight){
      binhigh++;
    }

    double xlow = proj->GetBinLowEdge(binlow);
    double xhigh = proj->GetBinLowEdge(binhigh);

    double mean=0.0, sigma=0.0;
    if(nentries>150 && xhigh>xlow){
      TFitResultPtr fr = proj->Fit("gaus","SQ0","",xlow,xhigh);
      if( fr && fr->IsValid() ){
	mean = fr->Parameter(1);
	sigma = fr->Parameter(2);
      }
      else{
	mean = proj->GetMean();
	sigma = proj->GetRMS();
      }
    }
    else{
      mean = proj->GetMean();
      sigma = proj->GetRMS();
    }

    double x = h2->GetXaxis()->GetBinCenter(i);
    g->SetPoint(i-1, x, mean);
    g->SetPointError(i-1, 0, sigma);

    delete proj;
  }
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kRed);
  g->SetLineColor(kBlack);
  return g;
}

TH2D* RemoveRunnumGap(TH2D* h2){
  double ymin = h2->GetYaxis()->GetXmin();
  double ymax = h2->GetYaxis()->GetXmax();
  int nbinsx = h2->GetNbinsX();
  int nbinsy = h2->GetNbinsY();
  
  vector<TH1D*> goodProjections;
  for (int i = 1; i<=nbinsx; ++i){
    TH1D *proj = h2->ProjectionY(Form("proj_%d",i), i, i);
    proj->SetDirectory(nullptr);
    if(proj->GetEntries() > 0){
      goodProjections.push_back(proj);
    }
    else{
      delete proj;
    }
  }

  int new_binsx = goodProjections.size();
  double min_x_bin_center = 0.0;
  double max_x_bin_center = new_binsx - 1.0;
  double new_bin_widthx = (max_x_bin_center - min_x_bin_center)/(new_binsx-1);
  double new_xmin = min_x_bin_center - new_bin_widthx/2.0;
  double new_xmax = max_x_bin_center + new_bin_widthx/2.0;

  TH2D* h2_new = new TH2D(Form("%s_nogap",h2->GetName()),h2->GetTitle(),new_binsx,new_xmin,new_xmax,nbinsy,ymin,ymax);

  for (int i = 0; i<new_binsx; ++i){
    TH1D* proj = goodProjections[i];
    for (int j = 1; j<=nbinsy; ++j){
      double content = proj->GetBinContent(j);
      h2_new->SetBinContent(i+1, j, content);
    }
    delete proj;
  }

  return h2_new;
}

void Cointime(std::vector<std::string> root_file_path, std::string fig_title, std::string kine_name){

  gStyle->SetPalette(kRainbow);
  gStyle->SetOptFit(1);
  gStyle->SetGridStyle(1);
  gStyle->SetGridColor(kBlack);
  gStyle->SetGridWidth(1);

  // Constants //

  auto kine = KineInfo(kine_name);
  double E_BEAM = kine.E_BEAM;
  double HCAL_DIST = kine.HCAL_DIST;
  double HCAL_ANGLE = kine.HCAL_ANGLE;
  double min_runnum = kine.min_runnum;
  double max_runnum = kine.max_runnum;
  double bin_runnum = kine.bin_runnum;
  double dx0p = kine.dx0p;
  double dy0p = kine.dy0p;
  double dx0n = kine.dx0n;
  double dy0n = kine.dy0n;
  double sigma_dx = kine.sigma_dx;
  double sigma_dy = kine.sigma_dy;
  double nsigma_dxdy = kine.nsigma_dxdy;
  double HCAL_THETA = HCAL_ANGLE*TMath::Pi()/180.0;

  TVector3 z_HCAL(-sin(HCAL_THETA),0.,cos(HCAL_THETA));
  TVector3 x_HCAL(0.,-1.,0.) ;
  TVector3 y_HCAL = (z_HCAL.Cross(x_HCAL)).Unit();

  TVector3 HCAL_origin = HCAL_DIST*z_HCAL;

  double Pprime_central = 2.0*MN*E_BEAM*(MN + E_BEAM)*cos(HCAL_THETA)/(pow(MN,2) + 2.0*MN*E_BEAM + pow(E_BEAM*sin(HCAL_THETA),2));
  double Eprime_central = sqrt(pow(Pprime_central,2) + pow(MN,2));
  double beta_central = Pprime_central/Eprime_central;
  double TOF_HCAL_central = HCAL_DIST/(beta_central*C_M_PER_NS);


  TChain *C = new TChain("T");
  for(const std::string& file : root_file_path){
    C->Add(file.c_str());
  }
  
  int numtrees = C->GetNtrees();
  std::cout << "Number of Trees Added: " << numtrees << std::endl;

  TH1D *hdt_BBSH_HCAL = new TH1D("hdt_BBSH_HCAL","BBSH - HCAL;t_{BBSH}^{FADC} - t_{HCAL}^{FADC} (ns);Counts",200,-20,20);
  TH1D *hdt_HODO_HCAL = new TH1D("hdt_HODO_HCAL","HODO - HCAL;t_{HODO}^{tfinal} - t_{HCAL}^{FADC} (ns);Counts",200,-20,20);
  TH2D *hdt_BBSH_BBPS_HODO_HCAL = new TH2D("hdt_BBSH_BBPS_HODO_HCAL","AVG of HCAL Coincidences;HCAL ID;(#Delta t^{HODO}_{HCAL} + #Delta t^{BBSH}_{HCAL} + #Delta t^{BBPS}_{HCAL})/3 (ns)",288,0.5,288.5,100,-20,20);

  TH2D *hdt_HODO_tfinal_IDHODO = new TH2D("hdt_HODO_tfinal_IDHODO","HODO tfinal vs ID;HODO ID;t_{HODO}^{tfinal} (ns)",90,-0.5,89.5,300,-20,20);
  TH2D *hdt_HODO_RFcorr_IDHODO = new TH2D("hdt_HODO_RFCorr_IDHODO","HODO tmeanRFcorr vs ID;HODO ID;t_{HODO} - t_{RF} (ns)",90,-0.5,89.5,500,-30,30);

  TH2D *hdt_HODO_HCAL_IDBLK = new TH2D("hdt_HODO_HCAL_IDBLK","HODO - HCAL vs IDBLK;HCAL ID;t_{HODO}^{tfinal} - t_{HCAL}^{FADC} (ns)",288,0.5,288.5,200,-20,20);
  TH2D *hdt_HODO_BBSH_IDBLK = new TH2D("hdt_HODO_BBSH_IDBLK","HODO - BBSH vs IDBLK;BBSH ID;t_{HODO}^{tfinal} - t_{BBSH}^{FADC} (ns)",189,-0.5,188.5,200,-20,20);
  TH2D *hdt_HODO_BBPS_IDBLK = new TH2D("hdt_HODO_BBPS_IDBLK","HODO - BBPS vs IDBLK;BBPS ID;t_{HODO}^{tfinal} - t_{BBPS}^{FADC} (ns)",52,-0.5,51.5,200,-20,20);
  TH2D *hdt_HODO_GRINCH_PMTNUM = new TH2D("hdt_HODO_GRINCH_PMTNUM","HODO - GRINCH vs PMTNUM;PMTNUM;t_{HODO}^{tfinal} - t_{GRINCH}^{hit-time} (ns)",512,-0.5,511.5,200,-30,30);

  TH2D *hHCAL_IDBLK = new TH2D("hHCAL_IDBLK","HCAL vs IDBLK;HCAL ID; t_{HCAL}^{FADC} (ns)",288,0.5,288.5,200,-20,20);
  TH2D *hBBSH_IDBLK = new TH2D("hBBSH_IDBLK","BBSH;BBSH ID; t_{BBSH}^{FADC} (ns)",189,-0.5,188.5,200,-20,20);
  TH2D *hBBPS_IDBLK = new TH2D("hBBPS_IDBLK",";BBPS ID; t_{BBPS}^{FADC} (ns)",52,-0.5,51.5,200,-20,20);
  TH2D *hGRINCH_PMTNUM = new TH2D("hGRINCH_PMTNUM","GRINCH vs PMTNUM;PMTNUM;t_{HODO}^{tfinal} - t_{GRINCH}^{hit-time} (ns)",512,-0.5,511.5,200,-30,30);

  TH1D *hdt_cluster_BBSH = new TH1D("hdt_cluster_BBSH","BBSH resolution;(bb.sh.atimeblk - bb.sh.clus_blk.atime[j]) (ns);Counts",200,-20,20);
  TH1D *hdt_cluster_BBPS = new TH1D("hdt_cluster_BBPS","BBPS resolution;(bb.ps.atimeblk - bb.ps.clus_blk.atime[j]) (ns);Counts",200,-20,20);
  TH1D *hdt_cluster_HCAL = new TH1D("hdt_cluster_HCAL","HCAL resolution;(sbs.hcal.atimeblk - sbs.hcal.clus_blk.atime[j]) (ns);Counts",200,-20,20);

  TH2D *hdt_cluster_HODO_BAR = new TH2D("hdt_cluster_HODO_BAR","HODO resolution vs Bar;HODO BAR ID; bb.hodotdc.clus.bar.tdc.tfinal[0] - bb.hodotdc.clus.bar.tdc.tfinal[j] (ns)",90,-0.5,89.5,300,-20,20);

  TH2D *hdt_cluster_BBSH_COL = new TH2D("hdt_cluster_BBSH_COL","BBSH resolution vs Column;BBSH COL ID;(bb.sh.atimeblk - bb.sh.clus_blk.atime[j]) (ns)",7,-0.5,6.5,200,-20,20);
  TH2D *hdt_cluster_BBPS_COL = new TH2D("hdt_cluster_BBPS_COL","BBPS resolution vs Column;BBPS COL ID;(bb.ps.atimeblk - bb.ps.clus_blk.atime[j]) (ns)",2,-0.5,1.5,200,-20,20);
  TH2D *hdt_cluster_HCAL_COL = new TH2D("hdt_cluster_HCAL_COL","HCAL resolution vs Column;HCAL COL ID;(sbs.hcal.atimeblk - sbs.hcal.clus_blk.atime[j]) (ns)",12,-0.5,11.5,200,-20,20);

  TH2D *hdt_cluster_BBSH_ROW = new TH2D("hdt_cluster_BBSH_ROW","BBSH resolution vs Row;BBSH ROW ID;(bb.sh.atimeblk - bb.sh.clus_blk.atime[j]) (ns)",26,-0.5,25.5,200,-20,20);
  TH2D *hdt_cluster_BBPS_ROW = new TH2D("hdt_cluster_BBPS_ROW","BBPS resolution vs Row;BBPS ROW ID;(bb.ps.atimeblk - bb.ps.clus_blk.atime[j]) (ns)",24,-0.5,23.5,200,-20,20);
  TH2D *hdt_cluster_HCAL_ROW = new TH2D("hdt_cluster_HCAL_ROW","HCAL resolution vs Row;HCAL ROW ID;(sbs.hcal.atimeblk - sbs.hcal.clus_blk.atime[j]) (ns)",24,-0.5,23.5,200,-20,20);

  TH2D *hdt_cluster_GRINCH_X = new TH2D("hdt_cluster_GRINCH_X","GRINCH resolution vs GRINCH X; GRINCH X (m); t_{GRINCH}^{tdcmean} - t_{GRINCH}^{hit-time[i]} (ns)",200,-1.,1.,200,-20,20);
  TH2D *hdt_cluster_GRINCH_Y = new TH2D("hdt_cluster_GRINCH_Y","GRINCH resolution vs GRINCH Y; GRINCH Y (m); t_{GRINCH}^{tdcmean} - t_{GRINCH}^{hit-time[i]} (ns)",200,-.15,.15,200,-20,20);

  TH2D *hdt_avg_BBSH_HCAL = new TH2D("hdt_avg_BBSH_HCAL","BBSH - HCAL vs IDBLK;HCAL ID;<bb.sh.clus_blk> - <sbs.hcal.clus_blk> (ns)",288,0.5,288.5,200,-20,20);
  TH2D *hdt_avg_BBPS_HCAL = new TH2D("hdt_avg_BBPS_HCAL","BBPS - HCAL vs IDBLK;HCAL ID;<bb.ps.clus_blk> - <sbs.hcal.clus_blk> (ns)",288,0.5,288.5,200,-20,20);
  TH2D *hdt_avg_BBSH_BBPS = new TH2D("hdt_avg_BBSH_BBPS","BBSH - BBPS vs IDBLK;BBSH ID;<bb.sh.clus_blk> - <bb.ps.clus_blk> (ns)",189,-0.5,188.5,200,-20,20);
  TH2D *hdt_avg_HODO_HCAL = new TH2D("hdt_avg_HODO_HCAL","HODO - HCAL vs IDBLK;HCAL ID;<bb.hodotdc.clus> - <sbs.hcal.clus_blk> (ns)",288,0.5,288.5,200,-20,20);

  TH2D *hdxdy = new TH2D("hdxdy","dx vs dy;dy (m);dx (m)",300,-4,4,300,-4,4);
  TH2D *hdtBBSH_HCAL_dx = new TH2D("hdtBBSH_HCAL_dx","dx vs BBSH - HCAL;dx (m);t_{BBSH}^{FADC} - t_{HCAL}^{FADC} (ns)",300,-4,4,300,-20,20);
  TH2D *hdtBBSH_HCAL_dy = new TH2D("hdtBBSH_HCAL_dy","dy vs BBSH - HCAL;dy (m);t_{BBSH}^{FADC} - t_{HCAL}^{FADC} (ns)",300,-4,4,300,-20,20);
  TH2D *hdtHODO_HCAL_dx = new TH2D("hdtHODO_HCAL_dx","dx vs HODO - HCAL;dx (m);t_{HODO} - t_{HCAL}^{FADC} (ns)",300,-4,4,300,-20,20);
  TH2D *hdtHODO_HCAL_dy = new TH2D("hdtHODO_HCAL_dy","dy vs HODO - HCAL;dy (m);t_{HODO} - t_{HCAL}^{FADC} (ns)",300,-4,4,300,-20,20);
  TH2D *hdtBBSH_HCAL_W2 = new TH2D("hdtBBSH_HCAL_W2","W2 vs BBSH - HCAL;W^{2} (GeV^{2});t_{BBSH}^{FADC} - t_{HCAL}^{FADC} (ns)",300,-1.0,6.0,300,-20,20);
  TH2D *hdtHODO_HCAL_W2 = new TH2D("hdtHODO_HCAL_W2","W2 vs HODO - HCAL;W^{2} (GeV^{2});t_{HODO} - t_{HCAL}^{FADC} (ns)",300,-1.0,6.0,300,-20,20);

  TH2D *hdtHODO_HCAL_runnum = new TH2D("hdtHODO_HCAL_runnum","HODO - HCAL vs run number;run number; t_{HODO} - t_{HCAL}^{FADC}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hdtHODO_BBSH_runnum = new TH2D("hdtHODO_BBSH_runnum","HODO - BBSH vs run number;run number; t_{HODO} - t_{BBSH}^{FADC}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hdtHODO_BBPS_runnum = new TH2D("hdtHODO_BBPS_runnum","HODO - BBPS vs run number;run number; t_{HODO} - t_{BBPS}^{FADC}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hdtHODO_GRINCH_runnum = new TH2D("hdtHODO_GRINCH_runnum","HODO - GRINCH vs run number;run number; t_{HODO} - t_{GRINCH}^{TDCmean}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  
  TH2D *hdtBBSH_HCAL_runnum = new TH2D("hdtBBSH_HCAL_runnum","BBSH - HCAL vs run number;run number; t_{BBSH}^{FADC} - t_{HCAL}^{FADC}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hdtBBSH_BBPS_runnum = new TH2D("hdtBBSH_BBPS_runnum","BBSH - BBPS vs run number;run number; t_{BBSH}^{FADC} - t_{BBPS}^{FADC}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hdtBBPS_HCAL_runnum = new TH2D("hdtBBPS_HCAL_runnum","BBPS - HCAL vs run number;run number; t_{BBPS}^{FADC} - t_{HCAL}^{FADC}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hdtBBSH_GRINCH_runnum = new TH2D("hdtBBSH_GRINCH_runnum","BBSH - GRINCH vs run number;run number; t_{BBSH}^{FADC} - t_{GRINCH}^{TDCmean}",bin_runnum,min_runnum,max_runnum,300,-20,20);

  TH2D *hHCAL_runnum = new TH2D("hHCAL_runnum","HCAL vs run number;run number; t_{HCAL}^{FADC}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hHODO_runnum = new TH2D("hHODO_runnum","HODO vs run number;run number; t_{HODO}^{tfinal}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hBBPS_runnum = new TH2D("hBBPS_runnum","BBPS vs run number;run number; t_{BBPS}^{FADC}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hBBSH_runnum = new TH2D("hBBSH_runnum","BBSH vs run number;run number; t_{BBSH}^{FADC}",bin_runnum,min_runnum,max_runnum,300,-20,20);
  TH2D *hGRINCH_runnum = new TH2D("hGRINCH_runnum","GRINCH vs run number;run number; t_{GRINCH}^{TDCmean}",bin_runnum,min_runnum,max_runnum,300,-20,20);

  TH2D *hdt_HODO_BBPS_trX = new TH2D("hdt_HODO_BBPS_trX","trX vs HODO - BBPS;track X (m); t_{HODO}^{tfinal} - t_{BBPS}^{FADC}",300,-0.6,0.6,300,-20,20);
  TH2D *hdt_HODO_BBPS_trY = new TH2D("hdt_HODO_BBPS_trY","trY vs HODO - BBPS;track Y (m); t_{HODO}^{tfinal} - t_{BBPS}^{FADC}",300,-0.2,0.2,300,-20,20);
  TH2D *hdt_HODO_BBPS_trPh = new TH2D("hdt_HODO_BBPS_trPh","trPh vs HODO - BBPS;track #phi; t_{HODO}^{tfinal} - t_{BBPS}^{FADC}",300,-0.1,0.1,300,-20,20);
  TH2D *hdt_HODO_BBPS_trTh = new TH2D("hdt_HODO_BBPS_trTh","trTh vs HODO - BBPS;track #theta; t_{HODO}^{tfinal} - t_{BBPS}^{FADC}",300,-0.2,0.2,300,-20,20);

  TH2D *hdt_HODO_BBSH_trX = new TH2D("hdt_HODO_BBSH_trX","trX vs HODO - BBSH;track X (m); t_{HODO}^{tfinal} - t_{BBSH}^{FADC}",300,-0.6,0.6,300,-20,20);
  TH2D *hdt_HODO_BBSH_trY = new TH2D("hdt_HODO_BBSH_trY","trY vs HODO - BBSH;track Y (m); t_{HODO}^{tfinal} - t_{BBSH}^{FADC}",300,-0.2,0.2,300,-20,20);
  TH2D *hdt_HODO_BBSH_trPh = new TH2D("hdt_HODO_BBSH_trPh","trPh vs HODO - BBSH;track #phi; t_{HODO}^{tfinal} - t_{BBSH}^{FADC}",300,-0.1,0.1,300,-20,20);
  TH2D *hdt_HODO_BBSH_trTh = new TH2D("hdt_HODO_BBSH_trTh","trTh vs HODO - BBSH;track #theta; t_{HODO}^{tfinal} - t_{BBSH}^{FADC}",300,-0.2,0.2,300,-20,20);

  TH2D *hdt_HODO_HCAL_trX = new TH2D("hdt_HODO_HCAL_trX","trX vs HODO - HCAL;track X (m); t_{HODO}^{tfinal} - t_{HCAL}^{FADC}",300,-0.6,0.6,300,-20,20);
  TH2D *hdt_HODO_HCAL_trY = new TH2D("hdt_HODO_HCAL_trY","trY vs HODO - HCAL;track Y (m); t_{HODO}^{tfinal} - t_{HCAL}^{FADC}",300,-0.2,0.2,300,-20,20);
  TH2D *hdt_HODO_HCAL_trPh = new TH2D("hdt_HODO_HCAL_trPh","trPh vs HODO - HCAL;track #phi; t_{HODO}^{tfinal} - t_{HCAL}^{FADC}",300,-0.1,0.1,300,-20,20);
  TH2D *hdt_HODO_HCAL_trTh = new TH2D("hdt_HODO_HCAL_trTh","trTh vs HODO - HCAL;track #theta; t_{HODO}^{tfinal} - t_{HCAL}^{FADC}",300,-0.2,0.2,300,-20,20);

  gen_tree *T = new gen_tree(C);

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("e.kine.W2",1);
  C->SetBranchStatus("bb.tr*",1);
  C->SetBranchStatus("bb.sh*",1);
  C->SetBranchStatus("bb.ps*",1);
  C->SetBranchStatus("bb.hodo*",1);
  C->SetBranchStatus("bb.etot*",1);
  C->SetBranchStatus("sbs.hcal.*",1);
  C->SetBranchStatus("bb.hodotdc.*",1);
  C->SetBranchStatus("bb.grinch_tdc.*",1);
  C->SetBranchStatus("g.*",1);

  int nevent = 0;
  Long64_t nentries = C->GetEntries();
  for(Long64_t i = 0; i < nentries; i++){
    T->GetEntry(i);
    if(nevent % 50000 == 0){
      std::cout << "Event number: " << nevent << '\n';
    }

    bool first_cut = (T->bb_tr_n > 0);
    if(!first_cut) continue;
    /*
    bool cutFormula = (T->bb_ps_e>0.2) &&
      (std::fabs(T->bb_tr_vz[0])<0.27) &&
      (std::fabs(T->bb_etot_over_p[0]-1.0)<0.3) &&
      (T->sbs_hcal_e>0.02) &&
      (T->g_trigbits>3) &&
      (T->g_trigbits<5) &&
      (T->bb_tr_n>0);
    */
    bool cutFormula = (T->bb_ps_e>0.2) &&
      (std::fabs(T->bb_tr_vz[0])<0.27) &&
      (std::fabs(T->bb_etot_over_p[0]-1.0)<0.3) &&
      (T->sbs_hcal_e>0.02) &&
      (T->bb_tr_n>0);

    nevent++;

    if(!cutFormula) continue;

    TLorentzVector kprime(T->bb_tr_px[0],T->bb_tr_py[0],T->bb_tr_pz[0],T->bb_tr_p[0]);
    TLorentzVector k(0.,0.,E_BEAM,E_BEAM);
    TLorentzVector P(0.,0.,0.,MN);
    TLorentzVector q = k - kprime;
    TLorentzVector Pprime = q + P;

    TVector3 vertex(T->bb_tr_vx[0],T->bb_tr_vy[0],T->bb_tr_vz[0]);
    TVector3 hcal_vect = T->sbs_hcal_x*x_HCAL + T->sbs_hcal_y*y_HCAL;
    TVector3 Phat = Pprime.Vect().Unit();

    double s_intersect = (HCAL_origin - vertex).Dot(z_HCAL)/(Phat.Dot(z_HCAL));
    TVector3 HCAL_intersect = vertex + s_intersect*Phat;

    double xHCAL_exp = (HCAL_intersect - HCAL_origin).Dot(x_HCAL);
    double yHCAL_exp = (HCAL_intersect - HCAL_origin).Dot(y_HCAL);

    double dx = T->sbs_hcal_x - xHCAL_exp;
    double dy = T->sbs_hcal_y - yHCAL_exp;

    double hodo_tfinal = T->bb_hodotdc_clus_tfinal[0];
    double hodo_tmeanRFcorr = T->bb_hodotdc_clus_tmeanRFcorr[0];
    double hodo_tmean = T->bb_hodotdc_clus_tmean[0];
    double hodo_id = T->bb_hodotdc_clus_id[0];

    double sh_adctime = T->bb_sh_clus_adctime[0];
    double ps_adctime = T->bb_ps_clus_adctime[0];
    double hcal_adctime = T->sbs_hcal_clus_adctime[0];
    double grinch_tdcmean = T->bb_grinch_tdc_clus_t_mean_corr;

    int hcal_idblk = int(T->sbs_hcal_idblk);
    int sh_idblk = int(T->bb_sh_idblk);
    int ps_idblk = int(T->bb_ps_idblk);

    double bb_tr_x = T->bb_tr_x[0];
    double bb_tr_y = T->bb_tr_y[0];
    double bb_tr_th = T->bb_tr_th[0];
    double bb_tr_ph = T->bb_tr_ph[0];

    double W2 = T->e_kine_W2;

    hdtBBSH_HCAL_W2->Fill(W2,sh_adctime - hcal_adctime);
    hdtHODO_HCAL_W2->Fill(W2,hodo_tfinal - hcal_adctime);

    hdt_HODO_BBPS_trX->Fill(bb_tr_x, hodo_tfinal - ps_adctime);
    hdt_HODO_BBPS_trY->Fill(bb_tr_y, hodo_tfinal - ps_adctime);
    hdt_HODO_BBPS_trPh->Fill(bb_tr_ph, hodo_tfinal - ps_adctime);
    hdt_HODO_BBPS_trTh->Fill(bb_tr_th, hodo_tfinal - ps_adctime);

    hdt_HODO_BBSH_trX->Fill(bb_tr_x, hodo_tfinal - sh_adctime);
    hdt_HODO_BBSH_trY->Fill(bb_tr_y, hodo_tfinal - sh_adctime);
    hdt_HODO_BBSH_trPh->Fill(bb_tr_ph, hodo_tfinal - sh_adctime);
    hdt_HODO_BBSH_trTh->Fill(bb_tr_th, hodo_tfinal - sh_adctime);

    hdt_HODO_HCAL_trX->Fill(bb_tr_x, hodo_tfinal - hcal_adctime);
    hdt_HODO_HCAL_trY->Fill(bb_tr_y, hodo_tfinal - hcal_adctime);
    hdt_HODO_HCAL_trPh->Fill(bb_tr_ph, hodo_tfinal - hcal_adctime);
    hdt_HODO_HCAL_trTh->Fill(bb_tr_th, hodo_tfinal - hcal_adctime);

    if(T->g_trigbits<5&&T->g_trigbits>3){
      hdt_HODO_RFcorr_IDHODO->Fill(hodo_id, hodo_tmeanRFcorr);
    }

    bool good_W2 = (W2>0.0)&&(W2<1.6);
    if(!good_W2) continue;

    double dt, sum_te, sum_e;
    int nclus;

    double bb_sh_e = T->bb_sh_e;
    double bb_sh_col = T->bb_sh_colblk;
    double bb_sh_row = T->bb_sh_rowblk;
    double bb_sh_t = T->bb_sh_atimeblk;
    sum_te = 0.0;
    sum_e = 0.0;
    nclus = int(T->bb_sh_clus_nblk[0]);
    for(int i=0; i<nclus; i++){
      double sh_ei = T->bb_sh_clus_blk_e[i];
      double sh_ti = T->bb_sh_clus_blk_atime[i];
      double sh_idi = T->bb_sh_clus_blk_id[i];
      if( (sh_ei<0.1*bb_sh_e) ) continue;
      sum_te += sh_ti*sh_ei;
      sum_e += T->bb_sh_clus_blk_e[i];
      hdt_HODO_BBSH_IDBLK->Fill(sh_idi,hodo_tfinal-sh_ti);
      dt = bb_sh_t - sh_ti;
      hBBSH_IDBLK->Fill(sh_idi,sh_ti);
      if( i>0 ){
	      hdt_cluster_BBSH->Fill(dt);
	      hdt_cluster_BBSH_COL->Fill(bb_sh_col, dt);
	      hdt_cluster_BBSH_ROW->Fill(bb_sh_row, dt);
      }
    }

    double avg_bb_sh_time = sum_te / sum_e;

    double bb_ps_e = T->bb_ps_e;
    double bb_ps_col = T->bb_ps_colblk;
    double bb_ps_row = T->bb_ps_rowblk;
    double bb_ps_t = T->bb_ps_atimeblk;
    sum_te = 0.0;
    sum_e = 0.0;
    nclus = int(T->bb_ps_clus_nblk[0]);
    for(int i=0; i<nclus; i++){
      double ps_ei = T->bb_ps_clus_blk_e[i];
      double ps_ti = T->bb_ps_clus_blk_atime[i];
      double ps_idi = T->bb_ps_clus_blk_id[i];
      if( (ps_ei<0.1*bb_ps_e) ) continue;
      sum_te += ps_ti*ps_ei;
      sum_e += ps_ei;
      hdt_HODO_BBPS_IDBLK->Fill(ps_idi,hodo_tfinal-ps_ti);
      dt = bb_ps_t - ps_ti;
      hBBPS_IDBLK->Fill(ps_idi,ps_ti);
      if( i>0 ){
	      hdt_cluster_BBPS->Fill(dt);
	      hdt_cluster_BBPS_COL->Fill(bb_ps_col, dt);
	      hdt_cluster_BBPS_ROW->Fill(bb_ps_row, dt);
      }
    }

    double avg_bb_ps_time = sum_te / sum_e;

    double sbs_hcal_e = T->sbs_hcal_e;
    double sbs_hcal_col = T->sbs_hcal_colblk;
    double sbs_hcal_row = T->sbs_hcal_rowblk;
    double sbs_hcal_t = T->sbs_hcal_atimeblk;
    sum_te = 0.0;
    sum_e = 0.0;
    nclus = int(T->sbs_hcal_clus_nblk[0]);
    for(int i=0; i<nclus; i++){
      double hcal_ei = T->sbs_hcal_clus_blk_e[i];
      double hcal_ti = T->sbs_hcal_clus_blk_atime[i];
      double hcal_idi = T->sbs_hcal_clus_blk_id[i];
      if( (hcal_ei<0.1*sbs_hcal_e) ) continue;
      sum_te += hcal_ti*hcal_ei;
      sum_e += hcal_ei;
      hdt_HODO_HCAL_IDBLK->Fill(hcal_idi,hodo_tfinal-hcal_ti);
      dt = sbs_hcal_t - hcal_ti;
      hHCAL_IDBLK->Fill(hcal_idi,hcal_ti);
      if( i>0 ){
	      hdt_cluster_HCAL->Fill(dt);
	      hdt_cluster_HCAL_COL->Fill(sbs_hcal_col, dt);
	      hdt_cluster_HCAL_ROW->Fill(sbs_hcal_row, dt);
      }
    }

    double avg_hcal_time = sum_te / sum_e;

    int nhits = int(T->bb_grinch_tdc_ngoodhits);
    double grinchx = T->bb_grinch_tdc_clus_x_mean;
    double grinchy = T->bb_grinch_tdc_clus_y_mean;
    for(int i=0; i<nhits; i++){
      double grinch_ti = T->bb_grinch_tdc_hit_time_corr[i];
      double grinch_pmti = T->bb_grinch_tdc_hit_pmtnum[i];
      if((T->bb_grinch_tdc_hit_trackindex[i]==0) && (T->bb_grinch_tdc_hit_clustindex[i]==T->bb_grinch_tdc_bestcluster)){
	hdt_HODO_GRINCH_PMTNUM->Fill(grinch_pmti, hodo_tfinal - grinch_ti);
	hGRINCH_PMTNUM->Fill(grinch_pmti, grinch_ti);
	if( (i>0) && (T->bb_grinch_tdc_clus_size>1) ){
	  double dt = grinch_tdcmean - grinch_ti;
	  hdt_cluster_GRINCH_X->Fill(grinchx, dt);
	  hdt_cluster_GRINCH_Y->Fill(grinchy, dt);
	}
      }
    }

    int nbars = int(T->Ndata_bb_hodotdc_clus_bar_tdc_tfinal);
    if(nbars>1){
      double hodo_t = T->bb_hodotdc_clus_bar_tdc_tfinal[0];
      double barid = T->bb_hodotdc_clus_id[0];
      for(int i=0; i<nbars; i++){
	if(i>0){
	  double hodo_ti = T->bb_hodotdc_clus_bar_tdc_tfinal[i];
	  double dt = hodo_t - hodo_ti;
	  hdt_cluster_HODO_BAR->Fill(barid, dt);
	}
      }
    }

    hdxdy->Fill(dy,dx);
    hdtBBSH_HCAL_dx->Fill(dx,sh_adctime - hcal_adctime);
    hdtBBSH_HCAL_dy->Fill(dy,sh_adctime - hcal_adctime);
    hdtHODO_HCAL_dx->Fill(dx,hodo_tfinal - hcal_adctime);
    hdtHODO_HCAL_dy->Fill(dy,hodo_tfinal - hcal_adctime);

    int runnum = int(T->g_runnum);
    hdtHODO_HCAL_runnum->Fill(runnum,hodo_tfinal-hcal_adctime);
    hdtHODO_BBSH_runnum->Fill(runnum,hodo_tfinal-sh_adctime);
    hdtHODO_BBPS_runnum->Fill(runnum,hodo_tfinal-ps_adctime);
    hdtHODO_GRINCH_runnum->Fill(runnum,hodo_tfinal - grinch_tdcmean);
    hdtBBSH_HCAL_runnum->Fill(runnum,sh_adctime-hcal_adctime);
    hdtBBSH_BBPS_runnum->Fill(runnum,sh_adctime-ps_adctime);
    hdtBBPS_HCAL_runnum->Fill(runnum,ps_adctime-hcal_adctime);
    hdtBBSH_GRINCH_runnum->Fill(runnum,sh_adctime - grinch_tdcmean);

    hHCAL_runnum->Fill(runnum,hcal_adctime); 
    hHODO_runnum->Fill(runnum,hodo_tfinal);
    hBBPS_runnum->Fill(runnum,ps_adctime);
    hBBSH_runnum->Fill(runnum,sh_adctime);
    hGRINCH_runnum->Fill(runnum,grinch_tdcmean);

    if(T->g_trigbits<5&&T->g_trigbits>3){
      hdt_HODO_tfinal_IDHODO->Fill(hodo_id, hodo_tfinal);
    }
    
    bool good_dxdy_n = (pow((dx-dx0n)/sigma_dx,2) + pow((dy-dy0n)/sigma_dy,2))<nsigma_dxdy;
    bool good_dxdy_p = (pow((dx-dx0p)/sigma_dx,2) + pow((dy-dy0p)/sigma_dy,2))<nsigma_dxdy;
    //bool good_dxdy = (fabs(dy)<0.5) && (fabs(dx)<0.5);
    if(!good_dxdy_n) continue;

    hdt_BBSH_HCAL->Fill(sh_adctime-hcal_adctime);
    hdt_HODO_HCAL->Fill(hodo_tfinal-hcal_adctime);

    double hodo_dt = hodo_tfinal - avg_hcal_time;
    double bbsh_dt = avg_bb_sh_time - avg_hcal_time;
    double bbps_dt = avg_bb_ps_time - avg_hcal_time;

    double avgdt = (hodo_dt + bbsh_dt + bbps_dt) / 3.0;

    hdt_avg_BBSH_HCAL->Fill( hcal_idblk, bbsh_dt );
    hdt_avg_BBPS_HCAL->Fill( hcal_idblk, bbps_dt );
    hdt_avg_BBSH_BBPS->Fill( sh_idblk, avg_bb_sh_time - avg_bb_ps_time );
    hdt_avg_HODO_HCAL->Fill( hcal_idblk, hodo_dt );
    hdt_BBSH_BBPS_HODO_HCAL->Fill( hcal_idblk, avgdt );
    
  }

  std::string base_dir = "/work/halla/sbs/koeneman/GEnII/";
  std::string out_dir = base_dir + "outdir/outfiles/Cointime/";
  std::string out_hist_name_root = "Cointime_" + fig_title + ".root";
  std::string out_hist_name_pdf = "Cointime_" + fig_title + ".pdf";
  std::string out_hist_path_root = out_dir + out_hist_name_root;
  std::string out_hist_path_pdf = out_dir + out_hist_name_pdf;
  TFile *out_hist_file_root = TFile::Open(out_hist_path_root.c_str(),"RECREATE");
  TCanvas *c = new TCanvas("c","c",800,600);
  c->cd();
  TString tempname = out_hist_path_pdf + "(";
  

  std::string fitfunction = "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2) + [6]";

  std::string fitfunction1 = "[0]*exp(-0.5*((x-[1])/[2])^2)";

  double Amp, mean, sigma;
  TFitResultPtr fr;

  gPad->SetGridx();
  TF1 *func0 = new TF1("fit0",fitfunction.c_str(),hdt_HODO_HCAL->GetXaxis()->GetXmin(),hdt_HODO_HCAL->GetXaxis()->GetXmax());
  func0->SetParameters(1000,0,2.0,500,-2.0,5.0,10);
  func0->SetParNames("Constant0","Mean0","Sigma0","Constant1","Mean1","Sigma1","Accidentals");
  hdt_HODO_HCAL->SetLineColor(kBlack);
  hdt_HODO_HCAL->SetMarkerStyle(20);
  hdt_HODO_HCAL->SetMarkerSize(0.8);
  hdt_HODO_HCAL->Draw("E");
  hdt_HODO_HCAL->Fit(func0,"RQ");
  func0->SetLineColor(kRed);
  func0->SetLineWidth(2);
  func0->Draw("SAME");
  gPad->Modified();
  gPad->Update();
  hdt_HODO_HCAL->Write();
  c->Print(tempname.Data());
  delete func0;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  hdt_HODO_HCAL->SetLineColor(kBlack);
  hdt_HODO_HCAL->SetMarkerStyle(20);
  hdt_HODO_HCAL->SetMarkerSize(0.8);
  hdt_HODO_HCAL->Draw("E");
  fr = FitPeak(hdt_HODO_HCAL);
  Amp = fr->Parameter(0);
  mean = fr->Parameter(1);
  sigma = fr->Parameter(2);
  TF1 *func = new TF1("fit","gaus",mean-0.7*sigma,mean+0.7*sigma);
  func->SetParameters(Amp,mean,sigma);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->Draw("SAME");
  c->Print(out_hist_path_pdf.c_str());
  delete func;
  delete hdt_HODO_HCAL;
    
  c->cd();
  c->Clear();
  gPad->SetGridx();
  TF1 *func1 = new TF1("fit1",fitfunction.c_str(),hdt_BBSH_HCAL->GetXaxis()->GetXmin(),hdt_BBSH_HCAL->GetXaxis()->GetXmax());
  func1->SetParameters(1000,0,2.0,500,-2.0,5.0,10);
  func1->SetParNames("Constant0","Mean0","Sigma0","Constant1","Mean1","Sigma1","Accidentals");
  hdt_BBSH_HCAL->SetLineColor(kBlack);
  hdt_BBSH_HCAL->SetMarkerStyle(20);
  hdt_BBSH_HCAL->SetMarkerSize(0.8);
  hdt_BBSH_HCAL->Draw("E");
  hdt_BBSH_HCAL->Fit(func1,"RQ");
  func1->SetLineColor(kRed);
  func1->SetLineWidth(2);
  func1->Draw("SAME");
  gPad->Modified();
  gPad->Update();
  hdt_BBSH_HCAL->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete func1;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  hdt_BBSH_HCAL->SetLineColor(kBlack);
  hdt_BBSH_HCAL->SetMarkerStyle(20);
  hdt_BBSH_HCAL->SetMarkerSize(0.8);
  hdt_BBSH_HCAL->Draw("E");
  fr = FitPeak(hdt_BBSH_HCAL);
  Amp = fr->Parameter(0);
  mean = fr->Parameter(1);
  sigma = fr->Parameter(2);
  func = new TF1("fit","gaus",mean-0.7*sigma,mean+0.7*sigma);
  func->SetParameters(Amp,mean,sigma);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->Draw("SAME");
  c->Print(out_hist_path_pdf.c_str());
  delete func;  
  delete hdt_BBSH_HCAL;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  TH1D *hdt_avg_BBSH_HCAL_proj = hdt_avg_BBSH_HCAL->ProjectionY("hdt_avg_BBSH_HCAL_proj");
  TF1 *func2 = new TF1("fit2",fitfunction.c_str(),hdt_avg_BBSH_HCAL_proj->GetXaxis()->GetXmin(),hdt_avg_BBSH_HCAL_proj->GetXaxis()->GetXmax());
  func2->SetParameters(1000,0,2.0,500,-2.0,5.0,10);
  func2->SetParNames("Constant0","Mean0","Sigma0","Constant1","Mean1","Sigma1","Accidentals");
  hdt_avg_BBSH_HCAL_proj->SetLineColor(kBlack);
  hdt_avg_BBSH_HCAL_proj->SetMarkerStyle(20);
  hdt_avg_BBSH_HCAL_proj->SetMarkerSize(0.8);
  hdt_avg_BBSH_HCAL_proj->Draw("E");
  hdt_avg_BBSH_HCAL_proj->Fit(func2,"RQ");
  func2->SetLineColor(kRed);
  func2->SetLineWidth(2);
  func2->Draw("SAME");
  gPad->Modified();
  gPad->Update();
  hdt_avg_BBSH_HCAL_proj->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete func2;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  hdt_avg_BBSH_HCAL_proj->SetLineColor(kBlack);
  hdt_avg_BBSH_HCAL_proj->SetMarkerStyle(20);
  hdt_avg_BBSH_HCAL_proj->SetMarkerSize(0.8);
  hdt_avg_BBSH_HCAL_proj->Draw("E");
  fr = FitPeak(hdt_avg_BBSH_HCAL_proj);
  Amp = fr->Parameter(0);
  mean = fr->Parameter(1);
  sigma = fr->Parameter(2);
  func = new TF1("fit","gaus",mean-0.7*sigma,mean+0.7*sigma);
  func->SetParameters(Amp,mean,sigma);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->Draw("SAME");
  c->Print(out_hist_path_pdf.c_str());
  delete func;
  delete hdt_avg_BBSH_HCAL_proj;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  TH1D *hdt_BBSH_BBPS_HODO_HCAL_proj = hdt_BBSH_BBPS_HODO_HCAL->ProjectionY("hdt_BBSH_BBPS_HODO_HCAL_proj");
  TF1 *func3 = new TF1("fit3",fitfunction.c_str(),hdt_BBSH_BBPS_HODO_HCAL_proj->GetXaxis()->GetXmin(),hdt_BBSH_BBPS_HODO_HCAL_proj->GetXaxis()->GetXmax());
  func3->SetParameters(1000,0,2.0,500,-2.0,5.0,10);
  func3->SetParNames("Constant0","Mean0","Sigma0","Constant1","Mean1","Sigma1","Accidentals");
  hdt_BBSH_BBPS_HODO_HCAL_proj->SetLineColor(kBlack);
  hdt_BBSH_BBPS_HODO_HCAL_proj->SetMarkerStyle(20);
  hdt_BBSH_BBPS_HODO_HCAL_proj->SetMarkerSize(0.8);
  hdt_BBSH_BBPS_HODO_HCAL_proj->Draw("E");
  hdt_BBSH_BBPS_HODO_HCAL_proj->Fit(func3,"RQ");
  func3->SetLineColor(kRed);
  func3->SetLineWidth(2);
  func3->Draw("SAME");
  gPad->Modified();
  gPad->Update();
  hdt_BBSH_BBPS_HODO_HCAL_proj->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete func3;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  hdt_BBSH_BBPS_HODO_HCAL_proj->SetLineColor(kBlack);
  hdt_BBSH_BBPS_HODO_HCAL_proj->SetMarkerStyle(20);
  hdt_BBSH_BBPS_HODO_HCAL_proj->SetMarkerSize(0.8);
  hdt_BBSH_BBPS_HODO_HCAL_proj->Draw("E");
  fr = FitPeak(hdt_BBSH_BBPS_HODO_HCAL_proj);
  Amp = fr->Parameter(0);
  mean = fr->Parameter(1);
  sigma = fr->Parameter(2);
  func = new TF1("fit","gaus",mean-0.7*sigma,mean+0.7*sigma);
  func->SetParameters(Amp,mean,sigma);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->Draw("SAME");
  c->Print(out_hist_path_pdf.c_str());
  delete func;
  delete hdt_BBSH_BBPS_HODO_HCAL_proj;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_BBSH_BBPS_HODO_HCAL->Draw("colz");
  gPad->Update();
  TGraphErrors *g_BBSH_BBPS_HODO_HCAL_mean = FitMeanAndStdDev(hdt_BBSH_BBPS_HODO_HCAL);
  g_BBSH_BBPS_HODO_HCAL_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_BBSH_BBPS_HODO_HCAL->Write();
  g_BBSH_BBPS_HODO_HCAL_mean->Write("g_BBSH_BBPS_HODO_HCAL_mean");
  c->Print(out_hist_path_pdf.c_str());
  delete hdt_BBSH_BBPS_HODO_HCAL;
  delete g_BBSH_BBPS_HODO_HCAL_mean;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_RFcorr_IDHODO->Draw("colz");
  hdt_HODO_RFcorr_IDHODO->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdt_HODO_RFcorr_IDHODO;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_tfinal_IDHODO->Draw("colz");
  hdt_HODO_tfinal_IDHODO->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdt_HODO_tfinal_IDHODO;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_HCAL_IDBLK->Draw("colz");
  gPad->Update();
  TGraphErrors *g_HODO_HCAL_IDBLK_mean = FitMeanAndStdDev(hdt_HODO_HCAL_IDBLK);
  g_HODO_HCAL_IDBLK_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_HODO_HCAL_IDBLK->Write();
  g_HODO_HCAL_IDBLK_mean->Write("g_HODO_HCAL_IDBLK_mean");
  c->Print(out_hist_path_pdf.c_str());

  delete g_HODO_HCAL_IDBLK_mean;
  delete hdt_HODO_HCAL_IDBLK;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBSH_IDBLK->Draw("colz");
  gPad->Update();
  TGraphErrors *g_HODO_BBSH_IDBLK_mean = FitMeanAndStdDev(hdt_HODO_BBSH_IDBLK);
  g_HODO_BBSH_IDBLK_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_HODO_BBSH_IDBLK->Write();
  g_HODO_BBSH_IDBLK_mean->Write("g_HODO_BBSH_IDBLK_mean");
  c->Print(out_hist_path_pdf.c_str());

  delete g_HODO_BBSH_IDBLK_mean;
  delete hdt_HODO_BBSH_IDBLK;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBPS_IDBLK->Draw("colz");
  gPad->Update();
  TGraphErrors *g_HODO_BBPS_IDBLK_mean = FitMeanAndStdDev(hdt_HODO_BBPS_IDBLK);
  g_HODO_BBPS_IDBLK_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_HODO_BBPS_IDBLK->Write();
  g_HODO_BBPS_IDBLK_mean->Write("g_HODO_BBPS_IDBLK_mean");
  c->Print(out_hist_path_pdf.c_str());

  delete g_HODO_BBPS_IDBLK_mean;
  delete hdt_HODO_BBPS_IDBLK;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_GRINCH_PMTNUM->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_HODO_GRINCH_PMTNUM_mean = FitMeanAndStdDev(hdt_HODO_GRINCH_PMTNUM);
  g_hdt_HODO_GRINCH_PMTNUM_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_HODO_GRINCH_PMTNUM->Write();
  g_hdt_HODO_GRINCH_PMTNUM_mean->Write("g_hdt_HODO_GRINCH_PMTNUM_mean");
  c->Print(out_hist_path_pdf.c_str());

  delete g_hdt_HODO_GRINCH_PMTNUM_mean;
  delete hdt_HODO_GRINCH_PMTNUM;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  TH1D *hdt_cluster_HODO = hdt_cluster_HODO_BAR->ProjectionY("hdt_cluster_HODO");
  hdt_cluster_HODO->SetLineColor(kBlack);
  hdt_cluster_HODO->SetMarkerStyle(20);
  hdt_cluster_HODO->SetMarkerSize(0.8);
  hdt_cluster_HODO->Draw("E");
  fr = FitPeak(hdt_cluster_HODO);
  Amp = fr->Parameter(0);
  mean = fr->Parameter(1);
  sigma = fr->Parameter(2);
  func = new TF1("fit","gaus",mean-0.7*sigma,mean+0.7*sigma);
  func->SetParameters(Amp,mean,sigma);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->Draw("SAME");
  hdt_cluster_HODO->Write();
  gPad->Modified();
  gPad->Update();
  c->Print(out_hist_path_pdf.c_str());
  delete func;
  delete hdt_cluster_HODO;

  c->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_cluster_HODO_BAR->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_cluster_HODO_BAR_mean = FitMeanAndStdDev(hdt_cluster_HODO_BAR);
  g_hdt_cluster_HODO_BAR_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_cluster_HODO_BAR->Write();
  g_hdt_cluster_HODO_BAR_mean->Write("g_hdt_cluster_HODO_BAR_mean");
  c->Print(out_hist_path_pdf.c_str());
  delete g_hdt_cluster_HODO_BAR_mean;
  delete hdt_cluster_HODO_BAR;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  hdt_cluster_BBSH->SetLineColor(kBlack);
  hdt_cluster_BBSH->SetMarkerStyle(20);
  hdt_cluster_BBSH->SetMarkerSize(0.8);
  hdt_cluster_BBSH->Draw("E");
  fr = FitPeak(hdt_cluster_BBSH);
  Amp = fr->Parameter(0);
  mean = fr->Parameter(1);
  sigma = fr->Parameter(2);
  func = new TF1("fit","gaus",mean-0.7*sigma,mean+0.7*sigma);
  func->SetParameters(Amp,mean,sigma);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->Draw("SAME");
  hdt_cluster_BBSH->Write();
  gPad->Modified();
  gPad->Update();
  c->Print(out_hist_path_pdf.c_str());
  delete func;
  delete hdt_cluster_BBSH;

  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_cluster_BBSH_COL->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_cluster_BBSH_COL_mean = FitMeanAndStdDev(hdt_cluster_BBSH_COL);
  g_hdt_cluster_BBSH_COL_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_cluster_BBSH_COL->Write();
  g_hdt_cluster_BBSH_COL_mean->Write("g_hdt_cluster_BBSH_COL_mean");
  
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_cluster_BBSH_ROW->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_cluster_BBSH_ROW_mean = FitMeanAndStdDev(hdt_cluster_BBSH_ROW);
  g_hdt_cluster_BBSH_ROW_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_cluster_BBSH_ROW->Write();
  g_hdt_cluster_BBSH_ROW_mean->Write("g_hdt_cluster_BBSH_ROW_mean");
  c->Print(out_hist_path_pdf.c_str());
  delete g_hdt_cluster_BBSH_COL_mean;
  delete g_hdt_cluster_BBSH_ROW_mean;
  delete hdt_cluster_BBSH_COL;
  delete hdt_cluster_BBSH_ROW;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  hdt_cluster_BBPS->SetLineColor(kBlack);
  hdt_cluster_BBPS->SetMarkerStyle(20);
  hdt_cluster_BBPS->SetMarkerSize(0.8);
  hdt_cluster_BBPS->Draw("E");
  fr = FitPeak(hdt_cluster_BBPS);
  Amp = fr->Parameter(0);
  mean = fr->Parameter(1);
  sigma = fr->Parameter(2);
  func = new TF1("fit","gaus",mean-0.7*sigma,mean+0.7*sigma);
  func->SetParameters(Amp,mean,sigma);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->Draw("SAME");
  hdt_cluster_BBPS->Write();
  gPad->Modified();
  gPad->Update();
  c->Print(out_hist_path_pdf.c_str());
  delete func;
  delete hdt_cluster_BBPS;

  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_cluster_BBPS_COL->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_cluster_BBPS_COL_mean = FitMeanAndStdDev(hdt_cluster_BBPS_COL);
  g_hdt_cluster_BBPS_COL_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_cluster_BBPS_COL->Write();
  g_hdt_cluster_BBPS_COL_mean->Write("g_hdt_cluster_BBPS_COL_mean");
  
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_cluster_BBPS_ROW->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_cluster_BBPS_ROW_mean = FitMeanAndStdDev(hdt_cluster_BBPS_ROW);
  g_hdt_cluster_BBPS_ROW_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_cluster_BBPS_ROW->Write();
  g_hdt_cluster_BBPS_ROW_mean->Write("g_hdt_cluster_BBPS_ROW_mean");
  c->Print(out_hist_path_pdf.c_str());
  delete g_hdt_cluster_BBPS_COL_mean;
  delete g_hdt_cluster_BBPS_ROW_mean;
  delete hdt_cluster_BBPS_COL;
  delete hdt_cluster_BBPS_ROW;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  hdt_cluster_HCAL->SetLineColor(kBlack);
  hdt_cluster_HCAL->SetMarkerStyle(20);
  hdt_cluster_HCAL->SetMarkerSize(0.8);
  hdt_cluster_HCAL->Draw("E");
  fr = FitPeak(hdt_cluster_HCAL);
  Amp = fr->Parameter(0);
  mean = fr->Parameter(1);
  sigma = fr->Parameter(2);
  func = new TF1("fit","gaus",mean-0.7*sigma,mean+0.7*sigma);
  func->SetParameters(Amp,mean,sigma);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->Draw("SAME");
  hdt_cluster_HCAL->Write();
  gPad->Modified();
  gPad->Update();
  c->Print(out_hist_path_pdf.c_str());
  delete func;
  delete hdt_cluster_HCAL;

  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_cluster_HCAL_COL->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_cluster_HCAL_COL_mean = FitMeanAndStdDev(hdt_cluster_HCAL_COL);
  g_hdt_cluster_HCAL_COL_mean->SetMarkerStyle(20);
  g_hdt_cluster_HCAL_COL_mean->SetMarkerColor(kRed);
  g_hdt_cluster_HCAL_COL_mean->SetLineColor(kBlack);
  g_hdt_cluster_HCAL_COL_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_cluster_HCAL_COL->Write();
  g_hdt_cluster_HCAL_COL_mean->Write("g_hdt_cluster_HCAL_COL_mean");
  
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_cluster_HCAL_ROW->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_cluster_HCAL_ROW_mean = FitMeanAndStdDev(hdt_cluster_HCAL_ROW);
  g_hdt_cluster_HCAL_ROW_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_cluster_HCAL_ROW->Write();
  g_hdt_cluster_HCAL_ROW_mean->Write("g_hdt_cluster_HCAL_ROW_mean");
  c->Print(out_hist_path_pdf.c_str());
  delete g_hdt_cluster_HCAL_COL_mean;
  delete g_hdt_cluster_HCAL_ROW_mean;
  delete hdt_cluster_HCAL_COL;
  delete hdt_cluster_HCAL_ROW;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  TH1D *hdt_cluster_GRINCH = hdt_cluster_GRINCH_X->ProjectionY("hdt_cluster_GRINCH");
  hdt_cluster_GRINCH->SetLineColor(kBlack);
  hdt_cluster_GRINCH->SetMarkerStyle(20);
  hdt_cluster_GRINCH->SetMarkerSize(0.8);
  hdt_cluster_GRINCH->Draw("E");
  fr = FitPeak(hdt_cluster_GRINCH);
  Amp = fr->Parameter(0);
  mean = fr->Parameter(1);
  sigma = fr->Parameter(2);
  func = new TF1("fit","gaus",mean-0.7*sigma,mean+0.7*sigma);
  func->SetParameters(Amp,mean,sigma);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->Draw("SAME");
  hdt_cluster_GRINCH->Write();
  gPad->Modified();
  gPad->Update();
  c->Print(out_hist_path_pdf.c_str());
  delete func;
  delete hdt_cluster_GRINCH;
  
  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_cluster_GRINCH_X->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_cluster_GRINCH_X_mean = FitMeanAndStdDev(hdt_cluster_GRINCH_X);
  g_hdt_cluster_GRINCH_X_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_cluster_GRINCH_X->Write();
  g_hdt_cluster_GRINCH_X_mean->Write("g_hdt_cluster_GRINCH_X_mean");

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_cluster_GRINCH_Y->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdt_cluster_GRINCH_Y_mean = FitMeanAndStdDev(hdt_cluster_GRINCH_Y);
  g_hdt_cluster_GRINCH_Y_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_cluster_GRINCH_Y->Write();
  g_hdt_cluster_GRINCH_Y_mean->Write("g_hdt_cluster_GRINCH_Y_mean");
  c->Print(out_hist_path_pdf.c_str());
  delete g_hdt_cluster_GRINCH_X_mean;
  delete g_hdt_cluster_GRINCH_Y_mean;
  delete hdt_cluster_GRINCH_X;
  delete hdt_cluster_GRINCH_Y;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_avg_BBPS_HCAL->Draw("colz");
  gPad->Update();
  TGraphErrors *g_avg_BBPS_HCAL_mean = FitMeanAndStdDev(hdt_avg_BBPS_HCAL);
  g_avg_BBPS_HCAL_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_avg_BBPS_HCAL->Write();
  g_avg_BBPS_HCAL_mean->Write("g_avg_BBPS_HCAL_mean");
  c->Print(out_hist_path_pdf.c_str());

  delete g_avg_BBPS_HCAL_mean;
  delete hdt_avg_BBPS_HCAL;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_avg_BBSH_BBPS->Draw("colz");
  gPad->Update();
  TGraphErrors *g_avg_BBSH_BBPS_mean = FitMeanAndStdDev(hdt_avg_BBSH_BBPS);
  g_avg_BBSH_BBPS_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_avg_BBSH_BBPS->Write();
  g_avg_BBSH_BBPS_mean->Write("g_avg_BBSH_BBPS_mean");
  c->Print(out_hist_path_pdf.c_str());

  delete g_avg_BBSH_BBPS_mean;
  delete hdt_avg_BBSH_BBPS;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_avg_HODO_HCAL->Draw("colz");
  gPad->Update();
  TGraphErrors *g_avg_HODO_HCAL_mean = FitMeanAndStdDev(hdt_avg_HODO_HCAL);
  g_avg_HODO_HCAL_mean->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdt_avg_HODO_HCAL->Write();
  g_avg_HODO_HCAL_mean->Write("g_avg_HODO_HCAL_mean");
  c->Print(out_hist_path_pdf.c_str());

  delete g_avg_HODO_HCAL_mean;
  delete hdt_avg_HODO_HCAL;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hHCAL_IDBLK->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hHCAL_IDBLK = FitMeanAndStdDev(hHCAL_IDBLK);
  g_hHCAL_IDBLK->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hHCAL_IDBLK->Write();
  g_hHCAL_IDBLK->Write("g_hHCAL_IDBLK");
  c->Print(out_hist_path_pdf.c_str());
  delete hHCAL_IDBLK;
  delete g_hHCAL_IDBLK;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hBBSH_IDBLK->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hBBSH_IDBLK = FitMeanAndStdDev(hBBSH_IDBLK);
  g_hBBSH_IDBLK->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hBBSH_IDBLK->Write();
  g_hBBSH_IDBLK->Write("g_hBBSH_IDBLK");
  c->Print(out_hist_path_pdf.c_str());
  delete hBBSH_IDBLK;
  delete g_hBBSH_IDBLK;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hBBPS_IDBLK->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hBBPS_IDBLK = FitMeanAndStdDev(hBBPS_IDBLK);
  g_hBBPS_IDBLK->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hBBPS_IDBLK->Write();
  g_hBBPS_IDBLK->Write("g_hBBPS_IDBLK");
  c->Print(out_hist_path_pdf.c_str());
  delete hBBPS_IDBLK;
  delete g_hBBPS_IDBLK;

  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hGRINCH_PMTNUM->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hGRINCH_PMTNUM = FitMeanAndStdDev(hGRINCH_PMTNUM);
  g_hGRINCH_PMTNUM->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hGRINCH_PMTNUM->Write();
  g_hGRINCH_PMTNUM->Write("g_hGRINCH_PMTNUM");
  c->Print(out_hist_path_pdf.c_str());
  delete hGRINCH_PMTNUM;
  delete g_hGRINCH_PMTNUM;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdxdy->Draw("colz");
  hdxdy->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdxdy;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdtBBSH_HCAL_dx->Draw("colz");
  hdtBBSH_HCAL_dx->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdtBBSH_HCAL_dx;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdtBBSH_HCAL_dy->Draw("colz");
  hdtBBSH_HCAL_dy->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdtBBSH_HCAL_dy;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdtHODO_HCAL_dx->Draw("colz");
  hdtHODO_HCAL_dx->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdtHODO_HCAL_dx;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdtHODO_HCAL_dy->Draw("colz");
  hdtHODO_HCAL_dy->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdtHODO_HCAL_dy;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdtBBSH_HCAL_W2->Draw("colz");
  hdtBBSH_HCAL_W2->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdtBBSH_HCAL_W2;
  
  c->cd();
  c->Clear();
  gPad->SetGridx();
  gPad->SetGridy();
  hdtHODO_HCAL_W2->Draw("colz");
  hdtHODO_HCAL_W2->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdtHODO_HCAL_W2;

  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBPS_trX->Draw("colz");
  hdt_HODO_BBPS_trX->Write();
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBPS_trY->Draw("colz");
  hdt_HODO_BBPS_trY->Write();
  c->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBPS_trPh->Draw("colz");
  hdt_HODO_BBPS_trPh->Write();
  c->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBPS_trTh->Draw("colz");
  hdt_HODO_BBPS_trTh->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdt_HODO_BBPS_trX;
  delete hdt_HODO_BBPS_trY;
  delete hdt_HODO_BBPS_trPh;
  delete hdt_HODO_BBPS_trTh;

  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBSH_trX->Draw("colz");
  hdt_HODO_BBSH_trX->Write();
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBSH_trY->Draw("colz");
  hdt_HODO_BBSH_trY->Write();
  c->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBSH_trPh->Draw("colz");
  hdt_HODO_BBSH_trPh->Write();
  c->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_BBSH_trTh->Draw("colz");
  hdt_HODO_BBSH_trTh->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdt_HODO_BBSH_trX;
  delete hdt_HODO_BBSH_trY;
  delete hdt_HODO_BBSH_trPh;
  delete hdt_HODO_BBSH_trTh;

  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_HCAL_trX->Draw("colz");
  hdt_HODO_HCAL_trX->Write();
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_HCAL_trY->Draw("colz");
  hdt_HODO_HCAL_trY->Write();
  c->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_HCAL_trPh->Draw("colz");
  hdt_HODO_HCAL_trPh->Write();
  c->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  hdt_HODO_HCAL_trTh->Draw("colz");
  hdt_HODO_HCAL_trTh->Write();
  c->Print(out_hist_path_pdf.c_str());
  delete hdt_HODO_HCAL_trX;
  delete hdt_HODO_HCAL_trY;
  delete hdt_HODO_HCAL_trPh;
  delete hdt_HODO_HCAL_trTh;

  TH2D* hHCAL_runnum_nogap = RemoveRunnumGap(hHCAL_runnum);
  TH2D* hHODO_runnum_nogap = RemoveRunnumGap(hHODO_runnum);
  TH2D* hBBPS_runnum_nogap = RemoveRunnumGap(hBBPS_runnum);
  TH2D* hBBSH_runnum_nogap = RemoveRunnumGap(hBBSH_runnum);
  TH2D* hGRINCH_runnum_nogap = RemoveRunnumGap(hGRINCH_runnum);

  hHCAL_runnum->Write();
  hHODO_runnum->Write();
  hBBPS_runnum->Write();
  hBBSH_runnum->Write();
  delete hHCAL_runnum;
  delete hHODO_runnum;
  delete hBBPS_runnum;
  delete hBBSH_runnum;
  
  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hHCAL_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hHCAL_runnum_nogap = FitMeanAndStdDev(hHCAL_runnum_nogap);
  g_hHCAL_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hHCAL_runnum_nogap->Write();
  g_hHCAL_runnum_nogap->Write("g_hHCAL_runnum_nogap");
  
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hHODO_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hHODO_runnum_nogap = FitMeanAndStdDev(hHODO_runnum_nogap);
  g_hHODO_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hHODO_runnum_nogap->Write();
  g_hHODO_runnum_nogap->Write("g_hHODO_runnum_nogap");
  c->Print(out_hist_path_pdf.c_str());
  
  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hBBPS_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hBBPS_runnum_nogap = FitMeanAndStdDev(hBBPS_runnum_nogap);
  g_hBBPS_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hBBPS_runnum_nogap->Write();
  g_hBBPS_runnum_nogap->Write("g_hBBPS_runnum_nogap");
  
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hBBSH_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hBBSH_runnum_nogap = FitMeanAndStdDev(hBBSH_runnum_nogap);
  g_hBBSH_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hBBSH_runnum_nogap->Write();
  g_hBBSH_runnum_nogap->Write("g_hBBSH_runnum_nogap");
  c->Print(out_hist_path_pdf.c_str());

  c->Clear();
  c->cd();
  gPad->SetGridx();
  gPad->SetGridy();
  hGRINCH_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hGRINCH_runnum_nogap = FitMeanAndStdDev(hGRINCH_runnum_nogap);
  g_hGRINCH_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hGRINCH_runnum_nogap->Write();
  g_hGRINCH_runnum_nogap->Write("g_hGRINCH_runnum_nogap");
  c->Print(out_hist_path_pdf.c_str());
  
  delete g_hHCAL_runnum_nogap;
  delete g_hHODO_runnum_nogap;
  delete g_hBBPS_runnum_nogap;
  delete g_hBBSH_runnum_nogap;
  delete g_hGRINCH_runnum_nogap;
  delete hHCAL_runnum_nogap;
  delete hHODO_runnum_nogap;
  delete hBBPS_runnum_nogap;
  delete hBBSH_runnum_nogap;
  delete hGRINCH_runnum_nogap;
  

  TH2D* hdtHODO_HCAL_runnum_nogap = RemoveRunnumGap(hdtHODO_HCAL_runnum);
  TH2D* hdtHODO_BBSH_runnum_nogap = RemoveRunnumGap(hdtHODO_BBSH_runnum);
  TH2D* hdtHODO_BBPS_runnum_nogap = RemoveRunnumGap(hdtHODO_BBPS_runnum);
  TH2D* hdtHODO_GRINCH_runnum_nogap = RemoveRunnumGap(hdtHODO_GRINCH_runnum);

  hdtHODO_HCAL_runnum->Write();
  hdtHODO_BBSH_runnum->Write();
  hdtHODO_BBPS_runnum->Write();
  hdtHODO_GRINCH_runnum->Write();
  delete hdtHODO_HCAL_runnum;
  delete hdtHODO_BBSH_runnum;
  delete hdtHODO_BBPS_runnum;
  delete hdtHODO_GRINCH_runnum;
  
  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdtHODO_HCAL_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdtHODO_HCAL_runnum_nogap = FitMeanAndStdDev(hdtHODO_HCAL_runnum_nogap);
  g_hdtHODO_HCAL_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdtHODO_HCAL_runnum_nogap->Write();
  g_hdtHODO_HCAL_runnum_nogap->Write("g_hdtHODO_HCAL_runnum_nogap");
  
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdtHODO_BBSH_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdtHODO_BBSH_runnum_nogap = FitMeanAndStdDev(hdtHODO_BBSH_runnum_nogap);
  g_hdtHODO_BBSH_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdtHODO_BBSH_runnum_nogap->Write();
  g_hdtHODO_BBSH_runnum_nogap->Write("g_hdtHODO_BBSH_runnum_nogap");
  c->Print(out_hist_path_pdf.c_str());

  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdtHODO_BBPS_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdtHODO_BBPS_runnum_nogap = FitMeanAndStdDev(hdtHODO_BBPS_runnum_nogap);
  g_hdtHODO_BBPS_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdtHODO_BBPS_runnum_nogap->Write();
  g_hdtHODO_BBPS_runnum_nogap->Write("g_hdtHODO_BBPS_runnum_nogap");
  
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdtHODO_GRINCH_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdtHODO_GRINCH_runnum_nogap = FitMeanAndStdDev(hdtHODO_GRINCH_runnum_nogap);
  g_hdtHODO_GRINCH_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdtHODO_GRINCH_runnum_nogap->Write();
  g_hdtHODO_GRINCH_runnum_nogap->Write("g_hdtHODO_GRINCH_runnum_nogap");
  c->Print(out_hist_path_pdf.c_str());
  
  delete g_hdtHODO_HCAL_runnum_nogap;
  delete g_hdtHODO_BBSH_runnum_nogap;
  delete g_hdtHODO_BBPS_runnum_nogap;
  delete g_hdtHODO_GRINCH_runnum_nogap;
  delete hdtHODO_HCAL_runnum_nogap;
  delete hdtHODO_BBSH_runnum_nogap;
  delete hdtHODO_BBPS_runnum_nogap;
  delete hdtHODO_GRINCH_runnum_nogap;

  TH2D* hdtBBSH_HCAL_runnum_nogap = RemoveRunnumGap(hdtBBSH_HCAL_runnum);
  TH2D* hdtBBSH_BBPS_runnum_nogap = RemoveRunnumGap(hdtBBSH_BBPS_runnum);
  TH2D* hdtBBPS_HCAL_runnum_nogap = RemoveRunnumGap(hdtBBPS_HCAL_runnum);
  TH2D* hdtBBSH_GRINCH_runnum_nogap = RemoveRunnumGap(hdtBBSH_GRINCH_runnum);

  hdtBBSH_HCAL_runnum->Write();
  hdtBBSH_BBPS_runnum->Write();
  hdtBBPS_HCAL_runnum->Write();
  hdtBBSH_GRINCH_runnum->Write();
  delete hdtBBSH_HCAL_runnum;
  delete hdtBBSH_BBPS_runnum;
  delete hdtBBPS_HCAL_runnum;
  delete hdtBBSH_GRINCH_runnum;

  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdtBBSH_HCAL_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdtBBSH_HCAL_runnum_nogap = FitMeanAndStdDev(hdtBBSH_HCAL_runnum_nogap);
  g_hdtBBSH_HCAL_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdtBBSH_HCAL_runnum_nogap->Write();
  g_hdtBBSH_HCAL_runnum_nogap->Write("g_hdtBBSH_HCAL_runnum_nogap");
  
  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdtBBSH_BBPS_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdtBBSH_BBPS_runnum_nogap = FitMeanAndStdDev(hdtBBSH_BBPS_runnum_nogap);
  g_hdtBBSH_BBPS_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdtBBSH_BBPS_runnum_nogap->Write();
  g_hdtBBSH_BBPS_runnum_nogap->Write("g_hdtBBSH_BBPS_runnum_nogap");
  c->Print(out_hist_path_pdf.c_str());

  c->Clear();
  c->Divide(1,2);
  c->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hdtBBPS_HCAL_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdtBBPS_HCAL_runnum_nogap = FitMeanAndStdDev(hdtBBPS_HCAL_runnum_nogap);
  g_hdtBBPS_HCAL_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdtBBPS_HCAL_runnum_nogap->Write();
  g_hdtBBPS_HCAL_runnum_nogap->Write("g_hdtBBPS_HCAL_runnum_nogap");

  c->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  hdtBBSH_GRINCH_runnum_nogap->Draw("colz");
  gPad->Update();
  TGraphErrors *g_hdtBBSH_GRINCH_runnum_nogap = FitMeanAndStdDev(hdtBBSH_GRINCH_runnum_nogap);
  g_hdtBBSH_GRINCH_runnum_nogap->Draw("P SAME");
  gPad->Modified();
  gPad->Update();
  hdtBBSH_GRINCH_runnum_nogap->Write();
  g_hdtBBSH_GRINCH_runnum_nogap->Write("g_hdtBBSH_GRINCH_runnum_nogap");
  tempname = out_hist_path_pdf + ")";
  c->Print(tempname.Data());
  delete g_hdtBBSH_HCAL_runnum_nogap;
  delete g_hdtBBSH_BBPS_runnum_nogap;
  delete g_hdtBBPS_HCAL_runnum_nogap;
  delete g_hdtBBSH_GRINCH_runnum_nogap;
  delete hdtBBSH_HCAL_runnum_nogap;
  delete hdtBBSH_BBPS_runnum_nogap;
  delete hdtBBPS_HCAL_runnum_nogap;
  delete hdtBBSH_GRINCH_runnum_nogap;
  out_hist_file_root->Close();
  
}
