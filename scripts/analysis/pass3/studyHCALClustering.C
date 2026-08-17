#include "../../../include/configParser.C"
#include "../../../include/computeKineVariables.C"
//#include "gen_SBSOFF.C"
//#include "gen_SBSON.C"

#include "TH2D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TGraph.h"

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

  TString outpdf_file = configKine + "_" + target + "_" + passKine + "_" + sbsConfigKine + "_" + "studyHCALClustering" + ".pdf";
  TCanvas *c = new TCanvas("c", "c", 1200, 800);
  

  double sbs_hcal_goodblock_atime[256], sbs_hcal_goodblock_e[256], sbs_hcal_goodblock_col[256], sbs_hcal_goodblock_row[256], sbs_hcal_goodblock_cid[256];
  double sbs_hcal_goodblock_x[256], sbs_hcal_goodblock_y[256];
  double bb_tr_px[100], bb_tr_py[100], bb_tr_pz[100], bb_tr_p[100], bb_tr_vx[100], bb_tr_vy[100], bb_tr_vz[100];
  double bb_sh_atimeblk, sbs_hcal_nclus;
  double sbs_hcal_x, sbs_hcal_y, sbs_hcal_e;
  int Ndata_sbs_hcal_goodblock_atime;
  
  C.SetBranchAddress("sbs.hcal.goodblock.atime", sbs_hcal_goodblock_atime);
  C.SetBranchAddress("sbs.hcal.goodblock.e", sbs_hcal_goodblock_e);
  C.SetBranchAddress("Ndata.sbs.hcal.goodblock.atime", &Ndata_sbs_hcal_goodblock_atime);
  C.SetBranchAddress("sbs.hcal.goodblock.col", sbs_hcal_goodblock_col);
  C.SetBranchAddress("sbs.hcal.goodblock.row", sbs_hcal_goodblock_row);
  C.SetBranchAddress("sbs.hcal.goodblock.cid", sbs_hcal_goodblock_cid);
  C.SetBranchAddress("sbs.hcal.goodblock.x", sbs_hcal_goodblock_x);
  C.SetBranchAddress("sbs.hcal.goodblock.y", sbs_hcal_goodblock_y);
  C.SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk);
  C.SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  C.SetBranchAddress("bb.tr.px", bb_tr_px);
  C.SetBranchAddress("bb.tr.py", bb_tr_py);
  C.SetBranchAddress("bb.tr.pz", bb_tr_pz);
  C.SetBranchAddress("bb.tr.p", bb_tr_p);
  C.SetBranchAddress("bb.tr.vx", bb_tr_vx);
  C.SetBranchAddress("bb.tr.vy", bb_tr_vy);
  C.SetBranchAddress("bb.tr.vz", bb_tr_vz);
  C.SetBranchAddress("sbs.hcal.x", &sbs_hcal_x);
  C.SetBranchAddress("sbs.hcal.y", &sbs_hcal_y);
  C.SetBranchAddress("sbs.hcal.e", &sbs_hcal_e);

  double coin_time_resolution = 1.11;
  double confidence_coin_time = 5.0;
  double hcal_pos_resolution = 0.161;
  double confidence_hcal_pos = 5.0;

  double confidence_weight = (confidence_coin_time*coin_time_resolution)*(confidence_hcal_pos*hcal_pos_resolution);
  
  Long64_t event = 0;
  int goodblock_tracker = 0;
  int max_goodblock_tracker = 100;
  while (C.GetEntry(event) > 0) {

    C.LoadTree(event);
    event++;

    if (cutFormula.EvalInstance() == 0) continue;

    TVector3 kf(bb_tr_px[0], bb_tr_py[0], bb_tr_pz[0]);
    TVector3 v(bb_tr_vx[0], bb_tr_vy[0], bb_tr_vz[0]);

    std::vector<double> dxdy_temp;
    double hcal_x_exp, hcal_y_exp;
    dxdy_temp = computeDxDy(target,
			    beam_energy,
			    hcal_angle,
			    hcal_distance,
			    kf,
			    v,
			    sbs_hcal_x,
			    sbs_hcal_y);

    hcal_x_exp = sbs_hcal_x - dxdy_temp[0];
    hcal_y_exp = sbs_hcal_y - dxdy_temp[1];

    Double_t gr_x[1], gr_y[1];
    gr_x[0] = hcal_y_exp;
    gr_y[0] = hcal_x_exp;

    std::vector<double> dxdyi;
    double dxi, dyi, dri;

    if ((event % 10 == 0) && (goodblock_tracker < max_goodblock_tracker) && (sbs_hcal_nclus>4)) {

      h2d_gb_cluster_0->Reset();
      h2d_gb_cluster_0->SetTitle(Form("Event %lld", event));
      h2d_gb_cluster_0->SetFillColor(kBlue);

      h2d_gb_cluster_1->Reset();
      h2d_gb_cluster_1->SetTitle(Form("Event %lld", event));
      h2d_gb_cluster_1->SetFillColor(kMagenta);

      h2d_gb_cluster_2->Reset();
      h2d_gb_cluster_2->SetTitle(Form("Event %lld", event));
      h2d_gb_cluster_2->SetFillColor(kOrange);

      h2d_gb_cluster_3->Reset();
      h2d_gb_cluster_3->SetTitle(Form("Event %lld", event));
      h2d_gb_cluster_3->SetFillColor(kBlack);


      h2d_gb_cluster_0_tweight->Reset();
      h2d_gb_cluster_0_tweight->SetTitle(Form("Coin Weighted Event %lld", event));
      h2d_gb_cluster_0_tweight->SetFillColor(kBlue);

      h2d_gb_cluster_1_tweight->Reset();
      h2d_gb_cluster_1_tweight->SetTitle(Form("Coin Weighted Event %lld", event));
      h2d_gb_cluster_1_tweight->SetFillColor(kMagenta);

      h2d_gb_cluster_2_tweight->Reset();
      h2d_gb_cluster_2_tweight->SetTitle(Form("Coin Weighted Event %lld", event));
      h2d_gb_cluster_2_tweight->SetFillColor(kOrange);

      h2d_gb_cluster_3_tweight->Reset();
      h2d_gb_cluster_3_tweight->SetTitle(Form("Coin Weighted Event %lld", event));
      h2d_gb_cluster_3_tweight->SetFillColor(kBlack);

      h2d_gb_cluster_0_tcut->Reset();
      h2d_gb_cluster_0_tcut->SetTitle(Form("Coin Cut (5 ns) Event %lld", event));
      h2d_gb_cluster_0_tcut->SetFillColor(kBlue);

      h2d_gb_cluster_1_tcut->Reset();
      h2d_gb_cluster_1_tcut->SetTitle(Form("Coin Cut (5 ns) Event %lld", event));
      h2d_gb_cluster_1_tcut->SetFillColor(kMagenta);

      h2d_gb_cluster_2_tcut->Reset();
      h2d_gb_cluster_2_tcut->SetTitle(Form("Coin Cut (5 ns) Event %lld", event));
      h2d_gb_cluster_2_tcut->SetFillColor(kOrange);

      h2d_gb_cluster_3_tcut->Reset();
      h2d_gb_cluster_3_tcut->SetTitle(Form("Coin Cut (5 ns) Event %lld", event));
      h2d_gb_cluster_3_tcut->SetFillColor(kBlack);

      TGraph *gr = new TGraph(1, gr_x, gr_y);
      
      int number_hcal_goodblocks = Ndata_sbs_hcal_goodblock_atime;
      double hcal_gb_adctimei, hcal_gb_ei, hcal_gb_coli, hcal_gb_rowi, hcal_gb_cidi, hcal_gb_xi, hcal_gb_yi;
      for (int i=0; i<number_hcal_goodblocks; i++) {
	hcal_gb_adctimei = sbs_hcal_goodblock_atime[i];
	hcal_gb_ei = sbs_hcal_goodblock_e[i];
	hcal_gb_coli = sbs_hcal_goodblock_col[i];
	hcal_gb_rowi = sbs_hcal_goodblock_row[i];
	hcal_gb_cidi = sbs_hcal_goodblock_cid[i];
	hcal_gb_xi = sbs_hcal_goodblock_x[i];
	hcal_gb_yi = sbs_hcal_goodblock_y[i];
	

	dxdyi = computeDxDy(target,
			    beam_energy,
			    hcal_angle,
			    hcal_distance,
			    kf,
			    v,
			    hcal_gb_xi,
			    hcal_gb_yi);

	dxi = dxdyi[0];
	dyi = dxdyi[1];

	dri = sqrt(dxi*dxi + dyi*dyi);

	double tdiff = abs(hcal_gb_adctimei-bb_sh_atimeblk);
	double weight = confidence_weight/(tdiff*dri);

	if (int(hcal_gb_cidi) == 0) {
	  h2d_gb_cluster_0->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  h2d_gb_cluster_0_tweight->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei*weight);
	  if (tdiff<1) {
	    h2d_gb_cluster_0_tcut->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  }
	}
	else if (int(hcal_gb_cidi) == 1) {
	  h2d_gb_cluster_1->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  h2d_gb_cluster_1_tweight->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei*weight);
	  if (tdiff<1) {
	    h2d_gb_cluster_1_tcut->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  }
	}
	else if (int(hcal_gb_cidi) == 2) {
	  h2d_gb_cluster_2->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  h2d_gb_cluster_2_tweight->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei*weight);
	  if (tdiff<1) {
	    h2d_gb_cluster_2_tcut->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  }
	}
	else {
	  h2d_gb_cluster_3->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  h2d_gb_cluster_3_tweight->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei*weight);
	  if (tdiff<1) {
	    h2d_gb_cluster_3_tcut->Fill(hcal_gb_yi,hcal_gb_xi,hcal_gb_ei);
	  }
	}
	
      }

      gr->SetMarkerStyle(30);
      gr->SetMarkerSize(2);
      gr->SetMarkerColor(kRed);

      c->Divide(3,1);
      c->cd(1);
      
      h2d_gb_cluster_0->Draw("box");
      h2d_gb_cluster_1->Draw("box sames");
      h2d_gb_cluster_2->Draw("box sames");
      h2d_gb_cluster_3->Draw("box sames");
      gr->Draw("P sames");

      c->Update();

      c->cd(2);
      h2d_gb_cluster_0_tweight->Draw("box");
      h2d_gb_cluster_1_tweight->Draw("box sames");
      h2d_gb_cluster_2_tweight->Draw("box sames");
      h2d_gb_cluster_3_tweight->Draw("box sames");
      gr->Draw("P sames");

      c->Update();

      c->cd(3);
      h2d_gb_cluster_0_tcut->Draw("box");
      h2d_gb_cluster_1_tcut->Draw("box sames");
      h2d_gb_cluster_2_tcut->Draw("box sames");
      h2d_gb_cluster_3_tcut->Draw("box sames");
      gr->Draw("P sames");
      
      c->Update();

      if (goodblock_tracker == 0) {
	c->Print((outpdf_file + "(").Data());
      }
      else {
	c->Print(outpdf_file.Data());
      }
      
      c->Clear();
      delete gr;

      goodblock_tracker++;
    }

    if (goodblock_tracker == max_goodblock_tracker) break;

  }

  c->Print((outpdf_file + ")").Data());

  delete c;
  delete h2d_gb_cluster_0;
  delete h2d_gb_cluster_1;
  delete h2d_gb_cluster_2;
  delete h2d_gb_cluster_3;
  
}
