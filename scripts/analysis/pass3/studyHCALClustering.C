#include "../../../include/configParser.C"

#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTreeFormula.h"
#include "TChain.h"

void studyHCALClustering(const std::string& config_filename) {
  
  readConfig(config_filename);

  double beam_energy = getConfigDouble("ebeam");
  TString target = getConfigString("target");
  TString configKine = getConfigString("config");
  TString passKine = getConfigString("pass");
  TString sbsConfigKine = getConfigString("sbs_config");
  TString rootFile = getConfigString("output_filename");
  TString rootDir = getConfigString("output_dir");
  TString goodeCut = getConfigString("goode_cut");
  double hcal_angle = getConfigDouble("hcal_angle");
  double hcal_distance = getConfigDouble("hcal_distance");
  double de_offset = getConfigDouble("de_offset");


  TString rootFileAll = rootFile(0, rootFile.Length() - 5) + "*";
  TString rootPath = rootDir + rootFileAll;

  std::cout << "ROOT file: " << rootFile << std::endl;
  std::cout << "goodeCut: " << goodeCut << std::endl;

  TChain C("T");

  C.Add(rootPath);

  TTreeFormula cutFormula("cutFormula", goodeCut.Data(), &C);

  TH2D *h2d_gb_cluster_0 = new TH2D("h2d_gb_cluster_0", "Primary Cluster GoodBlock; HCAL Y; HCAL X",14,-1.00711,1.00711,26,-2.734375,1.234375);
  TH2D *h2d_gb_cluster_1 = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_1");
  TH2D *h2d_gb_cluster_2 = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_2");
  TH2D *h2d_gb_cluster_3 = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_3");

  TH2D *h2d_gb_cluster_0_tcut = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_0_tcut");
  TH2D *h2d_gb_cluster_1_tcut = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_1_tcut");
  TH2D *h2d_gb_cluster_2_tcut = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_2_tcut");
  TH2D *h2d_gb_cluster_3_tcut = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_3_tcut");

  TH2D *h2d_gb_cluster_0_tweight = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_0_tweight");
  TH2D *h2d_gb_cluster_1_tweight = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_1_tweight");
  TH2D *h2d_gb_cluster_2_tweight = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_2_tweight");
  TH2D *h2d_gb_cluster_3_tweight = (TH2D*) h2d_gb_cluster_0->Clone("h2d_gb_cluster_3_tweight");

  TH2D *h2d_gb_dxblock_dtblocki = new TH2D("h2d_gb_dxblock_dtblocki",";t_{HCAL}^{main}-t_{HCAL}^{i};X_{HCAL}^{main}-X_{HCAL}^{i}",100,-20,20,39,-3.05288462,3.05288462);
  TH2D *h2d_gb_dyblock_dtblocki = new TH2D("h2d_gb_dyblock_dtblocki",";t_{HCAL}^{main}-t_{HCAL}^{i};Y_{HCAL}^{main}-Y_{HCAL}^{i}",100,-20,20,19,-1.4387286,1.4387286);
  TH2D *h2d_gb_dxblock_dtblock = new TH2D("h2d_gb_dxblock_dtblock",";t_{HCAL}^{main}-t_{HCAL}^{i};X_{HCAL}^{main}-X_{HCAL}^{i}",100,-20,20,39,-3.05288462,3.05288462);
  TH2D *h2d_gb_dyblock_dtblock = new TH2D("h2d_gb_dyblock_dtblock",";t_{HCAL}^{main}-t_{HCAL}^{i};Y_{HCAL}^{main}-Y_{HCAL}^{i}",100,-20,20,19,-1.4387286,1.4387286);

  TH2D *h2d_clust_survive_si = new TH2D("h2d_clust_survive_si","; Cluster id; S",10,-0.5,9.5,300,-0.1,1.1);
  TH2D *h2d_gb_dt_e = new TH2D("h2d_gb_dt_e", "Good Block Cointime vs Block Energy; E_{HCAL,gb}; t_{HCAL,gb} - t_{BBSH}",300,0,2.0,300,-15,15);
  TH2D *h2d_gb_dxdy = new TH2D("h2d_gb_dxdy", "Good Block position weighted=e^2 (weight>1); dy; dx",300,-2.,2.,300,-3.,1.5);
  TH1D *h1d_gb_dx = h2d_gb_dxdy->ProjectionY("h1d_gb_dx");
  TH1D *h1d_gb_dy = h2d_gb_dxdy->ProjectionX("h1d_gb_dy");
  TH2D *h2d_dxdy = new TH2D("h2d_dxdy", "Good Block position weighted=e^2 (weight>1); dy; dx",300,-2.,2.,300,-3.,1.5);
  TH1D *h1d_dx = h2d_dxdy->ProjectionY("h1d_dx");
  TH1D *h1d_dy = h2d_dxdy->ProjectionX("h1d_dy");
  TH1D *h1d_W2 = new TH1D("h1d_W2", "W2 with abs(dx_gb)<0.3 abs(dy_gb)<0.2; W^{2}; Counts",300,-1.,2.0);
  TH1D *h1d_W2_1 = (TH1D*) h1d_W2->Clone("h1d_W2_1");

  TGraph *gr = new TGraph();
  TGraph *gr2 = new TGraph();
  
  TString outpdf_file = configKine + "_" + target + "_" + passKine + "_" + sbsConfigKine + "_" + "studyHCALClustering" + ".pdf";
  TCanvas *c3 = new TCanvas("c3", "c3", 1200, 800);
  c3->Divide(3,1);
  TCanvas *c2 = new TCanvas("c2", "c2", 1200, 800);
  c2->Divide(2,1);

  double sbs_hcal_goodblock_atime[300], sbs_hcal_goodblock_e[300], sbs_hcal_goodblock_col[300], sbs_hcal_goodblock_row[300], sbs_hcal_goodblock_cid[300];
  double sbs_hcal_goodblock_x[300], sbs_hcal_goodblock_y[300], sbs_hcal_goodblock_id[300];
  double bb_sh_atimeblk, sbs_hcal_nclus;
  double sbs_hcal_x, sbs_hcal_y, sbs_hcal_e;
  double sbs_hcal_x_exp, sbs_hcal_y_exp;
  double sbs_hcal_dx, sbs_hcal_dy;
  double e_kine_W2, e_kine_Q2;
  int Ndata_sbs_hcal_goodblock_atime;
  
  C.SetBranchAddress("sbs.hcal.goodblock.atime", sbs_hcal_goodblock_atime);
  C.SetBranchAddress("sbs.hcal.goodblock.e", sbs_hcal_goodblock_e);
  C.SetBranchAddress("Ndata.sbs.hcal.goodblock.atime", &Ndata_sbs_hcal_goodblock_atime);
  C.SetBranchAddress("sbs.hcal.goodblock.col", sbs_hcal_goodblock_col);
  C.SetBranchAddress("sbs.hcal.goodblock.row", sbs_hcal_goodblock_row);
  C.SetBranchAddress("sbs.hcal.goodblock.cid", sbs_hcal_goodblock_cid);
  C.SetBranchAddress("sbs.hcal.goodblock.x", sbs_hcal_goodblock_x);
  C.SetBranchAddress("sbs.hcal.goodblock.y", sbs_hcal_goodblock_y);
  C.SetBranchAddress("sbs.hcal.goodblock.id", sbs_hcal_goodblock_id);
  C.SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk);
  C.SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  C.SetBranchAddress("sbs.hcal.x", &sbs_hcal_x);
  C.SetBranchAddress("sbs.hcal.y", &sbs_hcal_y);
  C.SetBranchAddress("sbs.hcal.e", &sbs_hcal_e);
  C.SetBranchAddress("sbs.hcal.x_exp", &sbs_hcal_x_exp);
  C.SetBranchAddress("sbs.hcal.y_exp", &sbs_hcal_y_exp);
  C.SetBranchAddress("sbs.hcal.dx", &sbs_hcal_dx);
  C.SetBranchAddress("sbs.hcal.dy", &sbs_hcal_dy);
  C.SetBranchAddress("e.kine.W2", &e_kine_W2);
  C.SetBranchAddress("e.kine.Q2", &e_kine_Q2);

  double coin_time_resolution = 0.77;
  double confidence_coin_time = 3.0;
  double dt_conf = coin_time_resolution*confidence_coin_time;
  double hcal_pos_resolution = 0.161;
  double confidence_hcal_pos = 3.0;
  double dr_conf = hcal_pos_resolution*confidence_hcal_pos;
  double hcal_e_resolution = 0.6;
  double confidence_hcal_e = 3.0;
  double de_conf = hcal_e_resolution*confidence_hcal_e;
  double si_threshold = 0.55;

  double confidence_weight = (dt_conf)*(dr_conf);

  double weight_dt_min = 0.5*dt_conf/(dt_conf*dt_conf + 0.25*dt_conf*dt_conf)/3.14159265359;
  double weight_dr_min = 0.5*dr_conf/(dr_conf*dr_conf + 0.25*dr_conf*dr_conf)/3.14159265359;

  //double weight_min = weight_dt_min*weight_dr_min;
  double weight_min = 1.0;

  Long64_t max_event = C.GetEntries();
  int goodblock_tracker = 0;
  int max_goodblock_tracker = 100;
  int currentTree = -1;
  
  for (Long64_t event=0; event<max_event; event++) {

    Long64_t local_entry = C.LoadTree(event);
    if (local_entry < 0) break;
    
    if (C.GetTreeNumber() != currentTree) {
      currentTree = C.GetTreeNumber();
      cutFormula.UpdateFormulaLeaves();
    }

    Long64_t entryLoading = C.GetEntry(event);
    if (entryLoading <=0) break;

    if (event % 50000 == 0) {
      std::cout << "\rProgress: " << event << std::flush;
    }
    
    if (cutFormula.EvalInstance() == 0) continue;

    int number_hcal_goodblocks = Ndata_sbs_hcal_goodblock_atime;
    if (number_hcal_goodblocks < 3) continue;

    bool bool_event_display = (goodblock_tracker < max_goodblock_tracker)&&(sbs_hcal_nclus>3);

    if (bool_event_display) {
      h2d_gb_cluster_0->Reset();
      h2d_gb_cluster_0->SetTitle(Form("Event %lld", event));
      h2d_gb_cluster_0->SetStats(0);
      h2d_gb_cluster_0->SetFillColor(kBlue);
      
      h2d_gb_cluster_1->Reset();
      h2d_gb_cluster_1->SetTitle(Form("Event %lld", event));
      h2d_gb_cluster_1->SetStats(0);
      h2d_gb_cluster_1->SetFillColor(kMagenta);
      
      h2d_gb_cluster_2->Reset();
      h2d_gb_cluster_2->SetTitle(Form("Event %lld", event));
      h2d_gb_cluster_2->SetStats(0);
      h2d_gb_cluster_2->SetFillColor(kOrange);
      
      h2d_gb_cluster_3->Reset();
      h2d_gb_cluster_3->SetTitle(Form("Event %lld", event));
      h2d_gb_cluster_3->SetStats(0);
      h2d_gb_cluster_3->SetFillColor(kBlack);
      
      h2d_gb_cluster_0_tweight->Reset();
      h2d_gb_cluster_0_tweight->SetTitle(Form("Coin Weighted Event %lld", event));
      h2d_gb_cluster_0_tweight->SetStats(0);
      h2d_gb_cluster_0_tweight->SetFillColor(kBlue);
      
      h2d_gb_cluster_1_tweight->Reset();
      h2d_gb_cluster_1_tweight->SetTitle(Form("Coin Weighted Event %lld", event));
      h2d_gb_cluster_1_tweight->SetStats(0);
      h2d_gb_cluster_1_tweight->SetFillColor(kMagenta);
      
      h2d_gb_cluster_2_tweight->Reset();
      h2d_gb_cluster_2_tweight->SetTitle(Form("Coin Weighted Event %lld", event));
      h2d_gb_cluster_2_tweight->SetStats(0);
      h2d_gb_cluster_2_tweight->SetFillColor(kOrange);
      
      h2d_gb_cluster_3_tweight->Reset();
      h2d_gb_cluster_3_tweight->SetTitle(Form("Coin Weighted Event %lld", event));
      h2d_gb_cluster_3_tweight->SetStats(0);
      h2d_gb_cluster_3_tweight->SetFillColor(kBlack);
      
      h2d_gb_cluster_0_tcut->Reset();
      h2d_gb_cluster_0_tcut->SetTitle(Form("Coin Cut (5 ns) Event %lld", event));
      h2d_gb_cluster_0_tcut->SetStats(0);
      h2d_gb_cluster_0_tcut->SetFillColor(kBlue);
      
      h2d_gb_cluster_1_tcut->Reset();
      h2d_gb_cluster_1_tcut->SetTitle(Form("Coin Cut (5 ns) Event %lld", event));
      h2d_gb_cluster_1_tcut->SetStats(0);
      h2d_gb_cluster_1_tcut->SetFillColor(kMagenta);
      
      h2d_gb_cluster_2_tcut->Reset();
      h2d_gb_cluster_2_tcut->SetTitle(Form("Coin Cut (5 ns) Event %lld", event));
      h2d_gb_cluster_2_tcut->SetStats(0);
      h2d_gb_cluster_2_tcut->SetFillColor(kOrange);
      
      h2d_gb_cluster_3_tcut->Reset();
      h2d_gb_cluster_3_tcut->SetTitle(Form("Coin Cut (5 ns) Event %lld", event));
      h2d_gb_cluster_3_tcut->SetStats(0);
      h2d_gb_cluster_3_tcut->SetFillColor(kBlack);

      h2d_gb_dxblock_dtblocki->Reset();
      h2d_gb_dxblock_dtblocki->SetTitle(Form("X dist. from best block vs relative time btw. best blocks | Event %lld", event));
      h2d_gb_dxblock_dtblocki->SetStats(0);

      h2d_gb_dyblock_dtblocki->Reset();
      h2d_gb_dyblock_dtblocki->SetTitle(Form("Y dist. from best block vs relative time btw. best blocks | Event %lld", event));
      h2d_gb_dyblock_dtblocki->SetStats(0);

      gr->Set(0);
      gr2->Set(0);
      
    }
    
    double hcal_gb_adctimei, hcal_gb_ei, hcal_gb_coli, hcal_gb_rowi, hcal_gb_cidi, hcal_gb_xi, hcal_gb_yi;
    double dti, tdiff, weight;
    double weight_dt, weight_dr, weight_de;
    double dxi, dyi, dri;
    double dei;
    double si;
    double si_total = 0.0;
    double si_max = -1.0;
    int id_max;
    double xgb_max = 0.0;
    double ygb_max = 0.0;
    double tgb_max = 0.0;
    // computing the total si = Ei(conf_dt/dt)^2(conf_dr/dr)^2
    for (int i=0; i<number_hcal_goodblocks; i++) {
      hcal_gb_adctimei = sbs_hcal_goodblock_atime[i];
      hcal_gb_ei = sbs_hcal_goodblock_e[i];
      hcal_gb_coli = sbs_hcal_goodblock_col[i];
      hcal_gb_rowi = sbs_hcal_goodblock_row[i];
      hcal_gb_cidi = sbs_hcal_goodblock_cid[i];
      hcal_gb_xi = sbs_hcal_goodblock_x[i];
      hcal_gb_yi = sbs_hcal_goodblock_y[i];

      dxi = hcal_gb_xi - sbs_hcal_x_exp;
      dyi = hcal_gb_yi - sbs_hcal_y_exp;
      
      dri = sqrt(dxi*dxi + dyi*dyi);
      
      dti = hcal_gb_adctimei-bb_sh_atimeblk;

      dei = 2.0*(hcal_gb_ei*hcal_gb_ei + 0.939*hcal_gb_ei) - e_kine_Q2 - de_offset;
      
      weight_dt = pow(dt_conf/dti,2);
      weight_dr = pow(dr_conf/dri,2);
      //weight_dr = 1.0;
      //weight_de = 1.0;
      weight_de = pow(de_conf/dei,2);

      //weight_dt = exp(-0.5*pow(dti/dt_conf,2))/sqrt(2.0*3.14159265359)/dt_conf;
      //weight_de = exp(-0.5*pow(dei/de_conf,2))/sqrt(2.0*3.14159265359)/de_conf;
      //weight_dr = exp(-0.5*pow(dri/dr_conf,2))/sqrt(2.0*3.14159265359)/dr_conf;
      
      //weight_dt = 0.5*dt_conf/(dti*dti + 0.25*dt_conf*dt_conf)/3.14159265359;
      //weight_dr = 0.5*dr_conf/(dri*dri + 0.25*dr_conf*dr_conf)/3.14159265359;
      
      weight = weight_dt*weight_dr*weight_de;

      si = hcal_gb_ei*weight;
      si_total += si;

      if (si>si_max) {
	si_max = si;
	id_max = int(sbs_hcal_goodblock_id[i]);
	xgb_max = hcal_gb_xi;
	ygb_max = hcal_gb_yi;
	tgb_max = hcal_gb_adctimei;
      }
    }
    
    double sum_e = 0.0;
    double sum_t = 0.0;
    int sum_blocks = 0;
    double sum_weight_posx = 0.0;
    double sum_weight_posy = 0.0;
    double sum_weight = 0.0;
    double ratioi, gbdti, gbdxi, gbdyi;
    for (int i=0; i<number_hcal_goodblocks; i++) {
      hcal_gb_adctimei = sbs_hcal_goodblock_atime[i];
      hcal_gb_ei = sbs_hcal_goodblock_e[i];
      hcal_gb_coli = sbs_hcal_goodblock_col[i];
      hcal_gb_rowi = sbs_hcal_goodblock_row[i];
      hcal_gb_cidi = sbs_hcal_goodblock_cid[i];
      hcal_gb_xi = sbs_hcal_goodblock_x[i];
      hcal_gb_yi = sbs_hcal_goodblock_y[i];
      
      dxi = hcal_gb_xi - sbs_hcal_x_exp;
      dyi = hcal_gb_yi - sbs_hcal_y_exp;
      
      dri = sqrt(dxi*dxi + dyi*dyi);
      
      dti = hcal_gb_adctimei-bb_sh_atimeblk;
      tdiff = abs(dti);

      dei = 2.0*(hcal_gb_ei*hcal_gb_ei + 0.939*hcal_gb_ei) - e_kine_Q2 - de_offset;
      
      weight_dt = pow(dt_conf/dti,2);
      weight_dr = pow(dr_conf/dri,2);
      //weight_dr = 1.0;
      //weight_de = 1.0;
      weight_de = pow(de_conf/dei,2);

      //weight_dt = exp(-0.5*pow(dti/dt_conf,2))/sqrt(2.0*3.14159265359)/dt_conf;
      //weight_de = exp(-0.5*pow(dei/de_conf,2))/sqrt(2.0*3.14159265359)/de_conf;
      //weight_dr = exp(-0.5*pow(dri/dr_conf,2))/sqrt(2.0*3.14159265359)/dr_conf;
      
      //weight_dt = 0.5*dt_conf/(dti*dti + 0.25*dt_conf*dt_conf)/3.14159265359;
      //weight_dr = 0.5*dr_conf/(dri*dri + 0.25*dr_conf*dr_conf)/3.14159265359;
      
      weight = weight_dt*weight_dr*weight_de;

      si = hcal_gb_ei*weight;
      si /= si_total;
      si_max /= si_total;

      ratioi = si/si_max;
      gbdti = tgb_max-hcal_gb_adctimei;
      gbdxi = xgb_max-hcal_gb_xi;
      gbdyi = ygb_max-hcal_gb_yi;

      if ((abs(gbdti)<5.5)&&(sqrt(gbdxi*gbdxi + gbdyi*gbdyi)<0.6)) {
	sum_blocks += 1;
	sum_e += hcal_gb_ei;
	sum_t += hcal_gb_adctimei*pow(hcal_gb_ei,2);;
	sum_weight_posx += hcal_gb_xi*pow(hcal_gb_ei,2);
	sum_weight_posy += hcal_gb_yi*pow(hcal_gb_ei,2);
	sum_weight += pow(hcal_gb_ei,2);
	if (int(sbs_hcal_goodblock_id[i])==id_max) {
	  h2d_clust_survive_si->Fill(hcal_gb_cidi,si);
	}
      }

      if (int(sbs_hcal_goodblock_id[i])!=id_max) {
	h2d_gb_dxblock_dtblock->Fill(gbdti,gbdxi);
	h2d_gb_dyblock_dtblock->Fill(gbdti,gbdyi);
      }
      
      if (bool_event_display) {
	
	if (int(hcal_gb_cidi) == 0) {
	  h2d_gb_cluster_0->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  h2d_gb_cluster_0_tweight->Fill(hcal_gb_yi,hcal_gb_xi,si);
	  if (tdiff<5.0) {
	    h2d_gb_cluster_0_tcut->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  }
	}
	else if (int(hcal_gb_cidi) == 1) {
	  h2d_gb_cluster_1->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  h2d_gb_cluster_1_tweight->Fill(hcal_gb_yi,hcal_gb_xi,si);
	  if (tdiff<5.0) {
	    h2d_gb_cluster_1_tcut->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  }
	}
	else if (int(hcal_gb_cidi) == 2) {
	  h2d_gb_cluster_2->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  h2d_gb_cluster_2_tweight->Fill(hcal_gb_yi,hcal_gb_xi,si);
	  if (tdiff<5.0) {
	    h2d_gb_cluster_2_tcut->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  }
	}
	else {
	  h2d_gb_cluster_3->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  h2d_gb_cluster_3_tweight->Fill(hcal_gb_yi,hcal_gb_xi,si);
	  if (tdiff<5.0) {
	    h2d_gb_cluster_3_tcut->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  }
	}

	if (int(sbs_hcal_goodblock_id[i])!=id_max) {
	  h2d_gb_dxblock_dtblocki->Fill(gbdti,gbdxi,ratioi);
	  h2d_gb_dyblock_dtblocki->Fill(gbdti,gbdyi,ratioi);
	}
      }
    }

    if (sum_blocks>1) {
      sum_weight_posx /= sum_weight;
      sum_weight_posy /= sum_weight;
      sum_t /= sum_weight;
      h2d_gb_dt_e->Fill(sum_e,sum_t - bb_sh_atimeblk);
      if ((abs(sum_t - bb_sh_atimeblk)<2.0)&&(sum_e>0.02)) {
	h2d_gb_dxdy->Fill(sum_weight_posy-sbs_hcal_y_exp,sum_weight_posx-sbs_hcal_x_exp);
	h2d_dxdy->Fill(sbs_hcal_dy,sbs_hcal_dx);
	h1d_gb_dy->Fill(sum_weight_posy-sbs_hcal_y_exp);
	h1d_dy->Fill(sbs_hcal_dy);
	if (abs(sum_weight_posy-sbs_hcal_y_exp)<0.2) {
	  h1d_gb_dx->Fill(sum_weight_posx-sbs_hcal_x_exp);
	  if (abs(sum_weight_posx-sbs_hcal_x_exp)<0.3) {
	    h1d_W2->Fill(e_kine_W2);
	  }
	}
	if (abs(sbs_hcal_dy)<0.2) {
	  h1d_dx->Fill(sbs_hcal_dx);
	  if (abs(sbs_hcal_dx)<0.3) {
	    h1d_W2_1->Fill(e_kine_W2);
	  }
	}
      }
    }
    else {
      sum_weight_posx = 0.0;
      sum_weight_posy = 0.0;
    }

    if (bool_event_display) {
      
      gr->SetPoint(0, sbs_hcal_y_exp, sbs_hcal_x_exp);
      gr->SetMarkerStyle(30);
      gr->SetMarkerSize(2);
      gr->SetMarkerColor(kRed);

      gr2->SetPoint(0, sum_weight_posy, sum_weight_posx);
      gr2->SetMarkerStyle(30);
      gr2->SetMarkerSize(2);
      gr2->SetMarkerColor(kBlack);

      c3->cd(1);
      gPad->Clear();
      h2d_gb_cluster_0->Draw("box");
      h2d_gb_cluster_1->Draw("box sames");
      h2d_gb_cluster_2->Draw("box sames");
      h2d_gb_cluster_3->Draw("box sames");
      gr->Draw("P sames");

      c3->Update();

      c3->cd(2);
      gPad->Clear();
      h2d_gb_cluster_0_tweight->Draw("box");
      h2d_gb_cluster_1_tweight->Draw("box sames");
      h2d_gb_cluster_2_tweight->Draw("box sames");
      h2d_gb_cluster_3_tweight->Draw("box sames");
      gr->Draw("P sames");

      c3->Update();

      c3->cd(3);
      gPad->Clear();
      h2d_gb_cluster_0_tcut->Draw("box");
      h2d_gb_cluster_1_tcut->Draw("box sames");
      h2d_gb_cluster_2_tcut->Draw("box sames");
      h2d_gb_cluster_3_tcut->Draw("box sames");
      gr->Draw("P sames");
      gr2->Draw("P sames");
      
      c3->Update();

      if (goodblock_tracker == 0) {
	c3->Print((outpdf_file + "(").Data());
      }
      else {
	c3->Print(outpdf_file.Data());
      }
      
      c2->cd(1);
      gPad->Clear();
      h2d_gb_dyblock_dtblocki->Draw("colz");
      c2->Update();

      c2->cd(2);
      gPad->Clear();
      h2d_gb_dxblock_dtblocki->Draw("colz");
      c2->Update();

      c2->Print(outpdf_file.Data());
      goodblock_tracker++;
    }

  }

  TCanvas *c = new TCanvas("c", "c", 1200, 800);

  c->cd();
  h2d_gb_dt_e->Draw("colz");
  c->Print(outpdf_file.Data());

  TH1D *h1d_gb_e = h2d_gb_dt_e->ProjectionX("h1d_gb_e");
  TH1D *h1d_gb_t = h2d_gb_dt_e->ProjectionY("h1d_gb_t");

  c->Clear();

  c->Divide(2,1);
  c->cd(1);
  h1d_gb_e->Draw();
  c->Update();
  
  c->cd(2);
  h1d_gb_t->Draw();
  c->Update();
  c->Print(outpdf_file.Data());

  c->cd();
  h2d_clust_survive_si->Draw("colz");
  c->Print(outpdf_file.Data());

  TH1D *h1d_clust_survive = h2d_clust_survive_si->ProjectionX("h1d_clust_survive");

  c->Clear();

  c->cd();
  h1d_clust_survive->Draw();
  c->Print(outpdf_file.Data());

  c->Clear();
  c->Divide(2,1);
  c->cd(1);
  h2d_gb_dxblock_dtblock->Draw("colz");
  c->Update();
  
  c->cd(2);
  h2d_gb_dyblock_dtblock->Draw("colz");
  c->Update();
  c->Print(outpdf_file.Data());

  c->Clear();

  c->Divide(2,1);
  c->cd(1);
  h2d_gb_dxdy->Draw("colz");
  c->Update();

  c->cd(2);
  h2d_dxdy->Draw("colz");
  c->Update();
  c->Print(outpdf_file.Data());

  c->Clear();
  
  c->Divide(2,1);
  c->cd(1);
  h1d_gb_dx->SetLineColor(kRed);
  h1d_gb_dx->SetLineWidth(2);
  h1d_gb_dx->Draw();
  h1d_dx->SetLineColor(kBlue);
  h1d_dx->SetLineWidth(2);
  h1d_dx->Draw("same");
  c->Update();

  c->cd(2);
  h1d_gb_dy->SetLineColor(kRed);
  h1d_gb_dy->SetLineWidth(2);
  h1d_gb_dy->Draw();
  h1d_dy->SetLineColor(kBlue);
  h1d_dy->SetLineWidth(2);
  h1d_dy->Draw("same");
  c->Update();
  c->Print(outpdf_file.Data());

  c->cd();
  h1d_W2->SetLineColor(kRed);
  h1d_W2->SetLineWidth(2);
  h1d_W2->Draw();
  h1d_W2_1->SetLineColor(kBlue);
  h1d_W2_1->SetLineWidth(2);
  h1d_W2_1->Draw("same");

  c->Print((outpdf_file + ")").Data());

  delete c;
  delete c2;
  delete c3;
  delete gr;
  delete gr2;
  delete h2d_gb_cluster_0;
  delete h2d_gb_cluster_1;
  delete h2d_gb_cluster_2;
  delete h2d_gb_cluster_3;
  delete h2d_gb_cluster_0_tweight;
  delete h2d_gb_cluster_1_tweight;
  delete h2d_gb_cluster_2_tweight;
  delete h2d_gb_cluster_3_tweight;
  delete h2d_gb_cluster_0_tcut;
  delete h2d_gb_cluster_1_tcut;
  delete h2d_gb_cluster_2_tcut;
  delete h2d_gb_cluster_3_tcut;
  delete h2d_gb_dt_e;
  delete h1d_gb_e;
  delete h1d_gb_t;
  delete h2d_gb_dxdy;
  delete h2d_dxdy;
  delete h1d_gb_dx;
  delete h1d_gb_dy;
  delete h1d_dx;
  delete h1d_dy;
  delete h1d_W2;
  delete h1d_W2_1;
  delete h2d_gb_dxblock_dtblock;
  delete h2d_gb_dyblock_dtblock;
  delete h2d_gb_dyblock_dtblocki;
  delete h2d_gb_dxblock_dtblocki;
  delete h2d_clust_survive_si;
  delete h1d_clust_survive;
}
