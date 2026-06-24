#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "gen_tree.C"

#include <iostream>
#include <cstdlib>
#include <fstream>

// Physical constants:
const double C_M_PER_NS = 0.299792458;
const double MP = 0.938272;

////////////////////////////////////////////////////////////
//// NEED TO UPDATE THESE PARAMETERS FOR EACH KINEMATIC ////

struct KineInfo_t {
  double E_BEAM;
  double HCAL_DIST;
  double HCAL_ANGLE;
};

KineInfo_t KineInfo(const std::string& kine){

  KineInfo_t info;
  if(kine=="GEN2"){
    // ++++ GEN2 +++++
    info.E_BEAM = 4.291;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 34.7;
  }
  if(kine=="GEN3"){
    // ++++ GEN3 +++++
    info.E_BEAM = 6.373;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 21.6;
  }
  if(kine=="GEN4" || kine=="GEN4b"){
    // ++++ GEN4 +++++
    info.E_BEAM = 8.448;
    info.HCAL_DIST = 17.0;
    info.HCAL_ANGLE = 18.0;
  }
  return info;
}


void SBSbbcal(std::string root_file_path, std::string fig_title, std::string kine_name){

  // Constants //

  auto kine = KineInfo(kine_name);
  double E_BEAM = kine.E_BEAM;
  double HCAL_DIST = kine.HCAL_DIST;
  double HCAL_ANGLE = kine.HCAL_ANGLE;
  double HCAL_THETA = HCAL_ANGLE*TMath::Pi()/180.0;

  TVector3 z_HCAL(-sin(HCAL_THETA),0.,cos(HCAL_THETA));
  TVector3 x_HCAL(0.,-1.,0.) ;
  TVector3 y_HCAL = (z_HCAL.Cross(x_HCAL)).Unit();

  TVector3 HCAL_origin = HCAL_DIST*z_HCAL;

  TChain *C = new TChain("T");
  C->Add(root_file_path.c_str());
  int numtrees = C->GetNtrees();
  std::cout << "Number of Trees Added: " << numtrees << std::endl;

  TCut cut = "";

  cut += "bb.sh.nblk>0&&sbs.hcal.nblk>0&&e.kine.W2<2.0&&e.kine.W2>0.0&&bb.ps.e>0.2&&sbs.hcal.e>0.02&&fabs(bb.tr.vz[0])<0.27&&fabs(bb.etot_over_p[0]-1.0)<0.3&&fabs(bb.sh.atimeblk - sbs.hcal.atimeblk)<3.0";

  TH2D *hdxdy = new TH2D("hdxdy",(fig_title + ";dy (m);dx (m)").c_str(),300,-4,4,300,-4,4);
  TH2D *hPexp_over_Pmeas_HCAL = new TH2D("hPexp_over_Pmeas_HCAL",(fig_title + ";HCAL BLOCK ID;2E^{clus}_{HCAL}M_{p}/Q^{2}").c_str(),288,0.5,288.5,100,-0.1,0.4);
  TH2D *hPexp_over_Pmeas_colHCAL = new TH2D("hPexp_over_Pmeas_colHCAL",(fig_title + ";HCAL col (m);2E^{clus}_{HCAL}M_{p}/Q^{2}").c_str(),24,0.5,24.5,100,-0.1,0.4);
  TH2D *hPexp_over_Pmeas_rowHCAL = new TH2D("hPexp_over_Pmeas_rowHCAL",(fig_title + ";HCAL row (m);2E^{clus}_{HCAL}M_{p}/Q^{2}").c_str(),12,0.5,12.5,100,-0.1,0.4);

  TH2D *hHCALe_vs_clusindex = new TH2D("hHCALe_vs_clusindex",(fig_title + ";HCAL Clus Index;HCAL E^{clus} (GeV)").c_str(),50,-0.5,48.5,200,0,2.0);
  TH1D *hHCALnclus = new TH1D("hHCALnclus",(fig_title + ";sbs.hcal.nclus;Counts").c_str(),50, -0.5, 48.5);
  gen_tree *T = new gen_tree(C);

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("e.kine.W2",1);
  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("bb.ps.*",1);
  C->SetBranchStatus("bb.sh.*",1);
  C->SetBranchStatus("sbs.hcal.*",1);
  C->SetBranchStatus("bb.tr.v*",1);
  C->SetBranchStatus("bb.tr.p*",1);

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
    TLorentzVector P(0.,0.,0.,MP);
    TLorentzVector q = k - kprime;
    TLorentzVector Pprime = q + P;

    TVector3 vertex(T->bb_tr_vx[0],T->bb_tr_vy[0],T->bb_tr_vz[0]);
    TVector3 hcal_vect = T->sbs_hcal_x*x_HCAL + T->sbs_hcal_y*y_HCAL;
    TVector3 Phat = Pprime.Vect().Unit();
    double Pprime_mag2 = Pprime.Vect().Mag2();

    double s_intersect = (HCAL_origin - vertex).Dot(z_HCAL)/(Phat.Dot(z_HCAL));
    TVector3 HCAL_intersect = vertex + s_intersect*Phat;
    TVector3 HCAL_intersect_actual = HCAL_origin + hcal_vect - vertex;

    double xHCAL_exp = (HCAL_intersect - HCAL_origin).Dot(x_HCAL);
    double yHCAL_exp = (HCAL_intersect - HCAL_origin).Dot(y_HCAL);

    double dx = T->sbs_hcal_x - xHCAL_exp;
    double dy = T->sbs_hcal_y - yHCAL_exp;

    hdxdy->Fill(dy,dx);

    int HCAL_nclus = int(T->sbs_hcal_nclus);
    hHCALnclus->Fill(HCAL_nclus);
    // int HCAL_nblks = int(T->sbs_hcal_nblk);
    for(int i=0; i<HCAL_nclus; i++){

        double Tmeasi = T->sbs_hcal_clus_e[i];
        double cointimei = T->bb_sh_atimeblk - T->sbs_hcal_clus_atimeblk[i];
        hHCALe_vs_clusindex->Fill(i, Tmeasi);
        if(Tmeasi>0.2*T->sbs_hcal_e && fabs(cointimei)<2.0){

            int HCAL_blocki = int(T->sbs_hcal_clus_id[i]);
            double colHCAL_measi = T->sbs_hcal_clus_col[i];
            double rowHCAL_measi = T->sbs_hcal_clus_row[i];

            double ratioi = 2*Tmeasi*MP / Pprime_mag2;

            hPexp_over_Pmeas_HCAL->Fill(HCAL_blocki, ratioi);
            hPexp_over_Pmeas_colHCAL->Fill(colHCAL_measi, ratioi);
            hPexp_over_Pmeas_rowHCAL->Fill(rowHCAL_measi, ratioi);
        }
    }

  }

  std::string base_dir = "/work/halla/sbs/koeneman/GEnII/";
  std::string out_dir = base_dir + "./../../../outdir/outfiles/Hcal/";
  std::string out_hist_name = "Hcal_" + fig_title + ".root";
  std::string out_hist_path = out_dir + out_hist_name;
  TFile *out_hist_file = TFile::Open(out_hist_path.c_str(),"RECREATE");

  hPexp_over_Pmeas_HCAL->Write();
  hPexp_over_Pmeas_XHCAL->Write();
  hPexp_over_Pmeas_YHCAL->Write();
  hdxdy->Write();
  out_hist_file->Close();
}