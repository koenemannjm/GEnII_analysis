#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "gen_tree.C"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>
#include <cstdlib>
#include <fstream>

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
  double t0;
};

KineInfo_t KineInfo(const std::string& kine){

  KineInfo_t info;
  if(kine=="GEN2"){
    // ++++ GEN2 +++++
    info.E_BEAM = 4.291;
    info.HCAL_DIST = 17.0;
    //info.HCAL_ANGLE = 34.7;
    info.HCAL_ANGLE = 34.0;
    info.min_runnum = 1997.5;
    info.max_runnum = 2323.5;
    info.bin_runnum = 326;
    info.t0 = 130.0;
  }
  if(kine=="GEN3"){
    // ++++ GEN3 +++++
    info.E_BEAM = 6.373;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 21.6;
    info.min_runnum = 2463.5;
    info.max_runnum = 3265.5;
    info.bin_runnum = 802;
    info.t0 = 119.0;
  }
  if(kine=="GEN4" || kine=="GEN4b"){
    // ++++ GEN4 +++++
    info.E_BEAM = 8.448;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 18.0;
    info.min_runnum = 3448.5;
    info.max_runnum = 4587.5;
    info.bin_runnum = 1139;
    info.t0 = 122.0;
  }
  if(kine=="GEN4b"){
    // ++++ GEN4b +++++
    info.E_BEAM = 8.448;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 18.0;
    info.min_runnum = 4982.5;
    info.max_runnum = 6086.5;
    info.bin_runnum = 1104;
    info.t0 = 185.0;
  }
  return info;
}


void Cointime_pass2(std::vector<std::string> root_file_path, std::string fig_title, std::string kine_name){

  // Constants //

  auto kine = KineInfo(kine_name);
  double E_BEAM = kine.E_BEAM;
  double HCAL_DIST = kine.HCAL_DIST;
  double HCAL_ANGLE = kine.HCAL_ANGLE;
  double min_runnum = kine.min_runnum;
  double max_runnum = kine.max_runnum;
  double bin_runnum = kine.bin_runnum;
  double t0HCAL = kine.t0;
  
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

  TCut cut = "";

  cut += "bb.ps.e>0.2&&sbs.hcal.e>0.02&&fabs(bb.tr.vz[0])<0.27&&fabs(bb.etot_over_p[0]-1.0)<0.3";


  TH1D *hdt_BBSH_HCAL = new TH1D("hdt_BBSH_HCAL",(fig_title + ";t_{BBSH}^{FADC} - t_{HCAL}^{FADC} (ns);Counts").c_str(),500,-50,50);
  TH2D *hdt_BBSH_BBPS_HODO_HCAL = new TH2D("hdt_BBSH_BBPS_HODO_HCAL",(fig_title + ";HCAL ID;(#Delta t^{HODO}_{HCAL} + #Delta t^{BBSH}_{HCAL} + #Delta t^{BBPS}_{HCAL})/3 (ns)").c_str(),288,0.5,288.5,500,-50,50);

  TH2D *hdt_HODO_tfinal_IDHODO = new TH2D("hdt_HODO_tfinal_IDHODO",(fig_title + ";HODO ID;t_{HODO} (ns)").c_str(),90,-0.5,89.5,500,-20,20);
  TH2D *hdt_HODO_RFcorr_IDHODO = new TH2D("hdt_HODO_RFCorr_IDHODO",(fig_title + ";HODO ID;t_{HODO} - t_{RF} (ns)").c_str(),90,-0.5,89.5,500,-30,30);

  TH2D *hdt_HODO_HCAL = new TH2D("hdt_HODO_HCAL",(fig_title + ";HCAL ID;t_{HODO} - t_{HCAL}^{FADC} (ns)").c_str(),288,0.5,288.5,500,-50,50);
  TH2D *hdt_HODO_BBSH = new TH2D("hdt_HODO_BBSH",(fig_title + ";BBSH ID;t_{HODO} - t_{BBSH}^{FADC} (ns)").c_str(),189,-0.5,188.5,500,-50,50);
  TH2D *hdt_HODO_BBPS = new TH2D("hdt_HODO_BBPS",(fig_title + ";BBPS ID;t_{HODO} - t_{BBPS}^{FADC} (ns)").c_str(),52,-0.5,51.5,500,-50,50);

  TH1D *hdt_cluster_BBSH = new TH1D("hdt_cluster_BBSH",(fig_title + ";(bb.sh.clus_blk.atimet[i] - bb.sh.clus_blk.atime[j]) (ns);Counts").c_str(),500,-50,50);
  TH1D *hdt_cluster_BBPS = new TH1D("hdt_cluster_BBPS",(fig_title + ";(bb.ps.clus_blk.atime[i] - bb.ps.clus_blk.atime[j]) (ns);Counts").c_str(),500,-50,50);
  TH1D *hdt_cluster_HCAL = new TH1D("hdt_cluster_HCAL",(fig_title + ";(sbs.hcal.clus_blk.atime[i] - sbs.hcal.clus_blk.atime[j]) (ns);Counts").c_str(),500,-50,50);

  TH2D *hdt_avg_BBSH_HCAL = new TH2D("hdt_avg_BBSH_HCAL",(fig_title + ";HCAL ID;<bb.sh.clus_blk> - <sbs.hcal.clus_blk> (ns)").c_str(),288,0.5,288.5,500,-50,50);
  TH2D *hdt_avg_BBPS_HCAL = new TH2D("hdt_avg_BBPS_HCAL",(fig_title + ";HCAL ID;<bb.ps.clus_blk> - <sbs.hcal.clus_blk> (ns)").c_str(),288,0.5,288.5,500,-50,50);
  TH2D *hdt_avg_BBSH_BBPS = new TH2D("hdt_avg_BBSH_BBPS",(fig_title + ";BBSH ID;<bb.sh.clus_blk> - <bb.ps.clus_blk> (ns)").c_str(),189,-0.5,188.5,500,-50,50);
  TH2D *hdt_avg_HODO_HCAL = new TH2D("hdt_avg_HODO_HCAL",(fig_title + ";HCAL ID;<bb.hodotdc.clus> - <sbs.hcal.clus_blk> (ns)").c_str(),288,0.5,288.5,500,-50,50);

  TH2D *hdxdy = new TH2D("hdxdy",(fig_title + ";dy (m);dx (m)").c_str(),300,-4,4,300,-4,4);
  TH2D *hdtBBSH_HCAL_dx = new TH2D("hdtBBSH_HCAL_dx",(fig_title + ";dx (m);t_{BBSH}^{FADC} - t_{HCAL}^{FADC} (ns)").c_str(),300,-4,4,300,-50,50);
  TH2D *hdtBBSH_HCAL_dy = new TH2D("hdtBBSH_HCAL_dy",(fig_title + ";dy (m);t_{BBSH}^{FADC} - t_{HCAL}^{FADC} (ns)").c_str(),300,-4,4,300,-50,50);
  TH2D *hdtHODO_HCAL_dx = new TH2D("hdtHODO_HCAL_dx",(fig_title + ";dx (m);t_{HODO} - t_{HCAL}^{FADC} (ns)").c_str(),300,-4,4,300,-50,50);
  TH2D *hdtHODO_HCAL_dy = new TH2D("hdtHODO_HCAL_dy",(fig_title + ";dy (m);t_{HODO} - t_{HCAL}^{FADC} (ns)").c_str(),300,-4,4,300,-50,50);
  TH2D *hdtBBSH_HCAL_W2 = new TH2D("hdtBBSH_HCAL_W2",(fig_title + ";W^{2} (GeV^{2});t_{BBSH}^{FADC} - t_{HCAL}^{FADC} (ns)").c_str(),300,-1.0,6.0,300,-50,50);
  TH2D *hdtHODO_HCAL_W2 = new TH2D("hdtHODO_HCAL_W2",(fig_title + ";W^{2} (GeV^{2});t_{HODO} - t_{HCAL}^{FADC} (ns)").c_str(),300,-1.0,6.0,300,-50,50);

  TH2D *hdtHODO_HCAL_runnum = new TH2D("hdtHODO_HCAL_runnum",(fig_title + ";run number; t_{HODO} - t_{HCAL}^{FADC}").c_str(),bin_runnum,min_runnum,max_runnum,300,-50,50);
  TH2D *hdtHODO_BBSH_runnum = new TH2D("hdtHODO_BBSH_runnum",(fig_title + ";run number; t_{HODO} - t_{BBSH}^{FADC}").c_str(),bin_runnum,min_runnum,max_runnum,300,-50,50);
  TH2D *hdtHODO_BBPS_runnum = new TH2D("hdtHODO_BBPS_runnum",(fig_title + ";run number; t_{HODO} - t_{BBPS}^{FADC}").c_str(),bin_runnum,min_runnum,max_runnum,300,-50,50);
  TH2D *hdtBBSH_HCAL_runnum = new TH2D("hdtBBSH_HCAL_runnum",(fig_title + ";run number; t_{BBSH}^{FADC} - t_{HCAL}^{FADC}").c_str(),bin_runnum,min_runnum,max_runnum,300,-50,50);
  TH2D *hdtBBSH_BBPS_runnum = new TH2D("hdtBBSH_BBPS_runnum",(fig_title + ";run number; t_{BBSH}^{FADC} - t_{BBPS}^{FADC}").c_str(),bin_runnum,min_runnum,max_runnum,300,-50,50);
  TH2D *hdtBBPS_HCAL_runnum = new TH2D("hdtBBPS_HCAL_runnum",(fig_title + ";run number; t_{BBPS}^{FADC} - t_{HCAL}^{FADC}").c_str(),bin_runnum,min_runnum,max_runnum,300,-50,50);


  gen_tree *T = new gen_tree(C);

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
    if(nevent % 10000 == 0){
      std::cout << "Event number: " << nevent << std::endl;
    }

    treenum = C->GetTreeNumber();
    if( nevent == 0 || treenum != oldtreenum ){
      oldtreenum = treenum;
      cutFormula->UpdateFormulaLeaves();
    }
    nevent++;

    if(cutFormula->EvalInstance(0)==0) continue;

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
    TVector3 HCAL_intersect_actual = HCAL_origin + hcal_vect - vertex;

    double L_path_HCAL_expect = (HCAL_intersect - vertex).Mag();
    double L_path_HCAL_actual = HCAL_intersect_actual.Mag();
    double beta = Pprime.Vect().Mag() / Pprime.E();

    double TOF_HCAL_expect = L_path_HCAL_expect/(beta*C_M_PER_NS);
    double TOF_HCAL_actual = L_path_HCAL_actual/(beta*C_M_PER_NS);
    double HCAL_TOF_corr = TOF_HCAL_expect - TOF_HCAL_central;

    double xHCAL_exp = (HCAL_intersect - HCAL_origin).Dot(x_HCAL);
    double yHCAL_exp = (HCAL_intersect - HCAL_origin).Dot(y_HCAL);

    double dx = T->sbs_hcal_x - xHCAL_exp;
    double dy = T->sbs_hcal_y - yHCAL_exp;

    //double hodo_tfinal = T->bb_hodotdc_clus_tfinal[0];
    double hodo_tfinal = T->bb_hodotdc_clus_tmean[0];
    double hodo_tmeanRFcorr = T->bb_hodotdc_clus_tmeanRFcorr[0];
    double hodo_id = T->bb_hodotdc_clus_id[0];

    double sh_adctime = T->bb_sh_clus_adctime[0];
    double ps_adctime = T->bb_ps_clus_adctime[0];
    double hcal_adctime = T->sbs_hcal_clus_adctime[0] - t0HCAL;
    //hcal_adctime -= HCAL_TOF_corr;

    int hcal_idblk = int(T->sbs_hcal_idblk);
    int sh_idblk = int(T->bb_sh_idblk);
    int ps_idblk = int(T->bb_ps_idblk);

    double W2 = T->e_kine_W2;

    hdtBBSH_HCAL_W2->Fill(W2,sh_adctime - hcal_adctime);
    hdtHODO_HCAL_W2->Fill(W2,hodo_tfinal - hcal_adctime);

    bool good_W2 = (W2>0.0)&&(W2<1.6);
    if(!good_W2) continue;

    hdxdy->Fill(dy,dx);
    hdtBBSH_HCAL_dx->Fill(dx,sh_adctime - hcal_adctime);
    hdtBBSH_HCAL_dy->Fill(dy,sh_adctime - hcal_adctime);
    hdtHODO_HCAL_dx->Fill(dx,hodo_tfinal - hcal_adctime);
    hdtHODO_HCAL_dy->Fill(dy,hodo_tfinal - hcal_adctime);

    int runnum = int(T->g_runnum);
    hdtHODO_HCAL_runnum->Fill(runnum,hodo_tfinal-hcal_adctime);
    hdtHODO_BBSH_runnum->Fill(runnum,hodo_tfinal-sh_adctime);
    hdtHODO_BBPS_runnum->Fill(runnum,hodo_tfinal-ps_adctime);
    hdtBBSH_HCAL_runnum->Fill(runnum,sh_adctime-hcal_adctime);
    hdtBBSH_BBPS_runnum->Fill(runnum,sh_adctime-ps_adctime);
    hdtBBPS_HCAL_runnum->Fill(runnum,ps_adctime-hcal_adctime);

    hdt_HODO_HCAL->Fill(hcal_idblk,hodo_tfinal-hcal_adctime);
    hdt_HODO_BBSH->Fill(sh_idblk,hodo_tfinal-sh_adctime);
    hdt_HODO_BBPS->Fill(ps_idblk,hodo_tfinal-ps_adctime);

    hdt_HODO_tfinal_IDHODO->Fill(hodo_id, hodo_tfinal);
    hdt_HODO_RFcorr_IDHODO->Fill(hodo_id, hodo_tmeanRFcorr);

    // bool good_dxdy = (pow((dx+2.814)/0.3,2) + pow(dy/0.4,2))<1.;
    bool good_dxdy = (fabs(dy)<0.5) && (fabs(dx)<0.5);
    if(!good_dxdy) continue;

    hdt_BBSH_HCAL->Fill(sh_adctime-hcal_adctime);
    
    double dt, sum_te, sum_e;
    int nclus;

    sum_te = 0.0;
    sum_e = 0.0;
    nclus = int(T->bb_sh_clus_nblk[0]);
    for(int i=0; i<nclus; i++){
      if(fabs(T->bb_sh_clus_blk_atime[i]-T->bb_sh_atimeblk)<5.0){
        sum_te += T->bb_sh_clus_blk_atime[i]*T->bb_sh_clus_blk_e[i];
        sum_e += T->bb_sh_clus_blk_e[i];
      }
      for(int j=i+1; j<nclus; j++){
        dt = T->bb_sh_clus_blk_atime[i] - T->bb_sh_clus_blk_atime[j];
        hdt_cluster_BBSH->Fill(dt);
      }
    }

    double avg_bb_sh_time = sum_te / sum_e;

    sum_te = 0.0;
    sum_e = 0.0;
    nclus = int(T->bb_ps_clus_nblk[0]);
    for(int i=0; i<nclus; i++){
      if(fabs(T->bb_ps_clus_blk_atime[i]-T->bb_ps_atimeblk)<5.0){
        sum_te += T->bb_ps_clus_blk_atime[i]*T->bb_ps_clus_blk_e[i];
        sum_e += T->bb_ps_clus_blk_e[i];
      }
      for(int j=i+1; j<nclus; j++){
        dt = T->bb_ps_clus_blk_atime[i] - T->bb_ps_clus_blk_atime[j];
        hdt_cluster_BBPS->Fill(dt);
      }
    }

    double avg_bb_ps_time = sum_te / sum_e;

    sum_te = 0.0;
    sum_e = 0.0;
    nclus = int(T->sbs_hcal_clus_nblk[0]);
    for(int i=0; i<nclus; i++){
      if(fabs(T->sbs_hcal_clus_blk_atime[i]-T->sbs_hcal_atimeblk)<5.0){
        sum_te += (T->sbs_hcal_clus_blk_atime[i]-t0HCAL)*T->sbs_hcal_clus_blk_e[i];
        sum_e += T->sbs_hcal_clus_blk_e[i];
      }
      for(int j=i+1; j<nclus; j++){
        dt = T->sbs_hcal_clus_blk_atime[i] - T->sbs_hcal_clus_blk_atime[j];
        hdt_cluster_HCAL->Fill(dt);
      }
    }

    double avg_hcal_time = sum_te / sum_e;

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
  std::string out_hist_name = "Cointime_" + fig_title + ".root";
  std::string out_hist_path = out_dir + out_hist_name;
  TFile *out_hist_file = TFile::Open(out_hist_path.c_str(),"RECREATE");
  gStyle->SetPalette(kRainbow);
  gStyle->SetOptFit(1111);
  gStyle->SetGridStyle(1);
  gStyle->SetGridColor(kBlack);
  gStyle->SetGridWidth(1);

  std::string fitfunction = "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]";

  TF1 *func1 = new TF1("fit",fitfunction.c_str(),hdt_BBSH_HCAL->GetXaxis()->GetXmin(),hdt_BBSH_HCAL->GetXaxis()->GetXmax());
  func1->SetParameters(1000,0,2.0,500,-2.0,5.0,10);
  func1->SetParNames("Constant","Mean","Sigma","Accidentals");
  hdt_BBSH_HCAL->Fit(func1,"RQ");
  hdt_BBSH_HCAL->Write();
  hdt_HODO_RFcorr_IDHODO->Write();
  hdt_HODO_tfinal_IDHODO->Write();
  hdt_HODO_HCAL->Write();
  hdt_HODO_BBSH->Write();
  hdt_HODO_BBPS->Write();
  hdt_cluster_BBSH->Write();
  hdt_cluster_BBPS->Write();
  hdt_cluster_HCAL->Write();
  hdt_avg_BBSH_HCAL->Write();
  TH1D *hdt_avg_BBSH_HCAL_proj = hdt_avg_BBSH_HCAL->ProjectionY("hdt_avg_BBSH_HCAL_proj");
  TF1 *func2 = new TF1("fit",fitfunction.c_str(),hdt_avg_BBSH_HCAL_proj->GetXaxis()->GetXmin(),hdt_avg_BBSH_HCAL_proj->GetXaxis()->GetXmax());
  func2->SetParameters(1000,0,2.0,500,-2.0,5.0,10);
  func2->SetParNames("Constant","Mean","Sigma","Accidentals");
  hdt_avg_BBSH_HCAL_proj->Fit(func2,"RQ");
  hdt_avg_BBSH_HCAL_proj->Write();
  hdt_avg_BBPS_HCAL->Write();
  hdt_avg_BBSH_BBPS->Write();
  hdxdy->Write();
  hdtBBSH_HCAL_dx->Write();
  hdtBBSH_HCAL_dy->Write();
  hdtHODO_HCAL_dx->Write();
  hdtHODO_HCAL_dy->Write();
  hdtBBSH_HCAL_W2->Write();
  hdtHODO_HCAL_W2->Write();
  hdt_avg_HODO_HCAL->Write();
  hdt_BBSH_BBPS_HODO_HCAL->Write();
  TH1D *hdt_BBSH_BBPS_HODO_HCAL_proj = hdt_BBSH_BBPS_HODO_HCAL->ProjectionY("hdt_BBSH_BBPS_HODO_HCAL_proj");
  TF1 *func3 = new TF1("fit",fitfunction.c_str(),hdt_BBSH_BBPS_HODO_HCAL_proj->GetXaxis()->GetXmin(),hdt_BBSH_BBPS_HODO_HCAL_proj->GetXaxis()->GetXmax());
  func3->SetParameters(1000,0,2.0,500,-2.0,5.0,10);
  func3->SetParNames("Constant","Mean","Sigma","Accidentals");
  hdt_BBSH_BBPS_HODO_HCAL_proj->Fit(func3,"RQ");
  func1->Write();
  func2->Write();
  func3->Write();
  hdt_BBSH_BBPS_HODO_HCAL_proj->Write();
  hdtHODO_HCAL_runnum->Write();
  hdtHODO_BBSH_runnum->Write();
  hdtHODO_BBPS_runnum->Write();
  hdtBBSH_HCAL_runnum->Write();
  hdtBBSH_BBPS_runnum->Write();
  hdtBBPS_HCAL_runnum->Write();
  out_hist_file->Close();

}
