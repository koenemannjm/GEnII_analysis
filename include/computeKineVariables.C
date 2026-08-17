#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <vector>


std::vector<double> computeDxDy(TString target,
				double beam_energy,
				double hcal_angle,
				double hcal_distance,
				TVector3 kf,
				TVector3 v,
				double hcalx,
				double hcaly) {

  // Some constant(s)
  double Mp = 0.938272088; // PDG 2025
  double Mn = 0.939565420; // PDG 2025
  double MN = 0.5*(Mp + Mn);
  double pi = 3.1415926535;

  TVector3 hcal_unit_z(-TMath::Sin(hcal_angle*pi/180.0),
		       0.0,
		       TMath::Cos(hcal_angle*pi/180.0));
  TVector3 hcal_unit_x(0,-1.0,0);
  TVector3 hcal_unit_y = hcal_unit_z.Cross(hcal_unit_x).Unit();
  TVector3 hcal_vector = hcal_distance*hcal_unit_z;

  TVector3 ki(0.0,0.0,beam_energy);
  TVector3 q = ki - kf;
  TVector3 q_unit = q.Unit();
  
  // computing dx and dy

  double w = (hcal_vector - v).Dot(hcal_unit_z) / (q_unit.Dot(hcal_unit_z));
  TVector3 W = v + w*q_unit;
  TVector3 D = W - hcal_vector;

  double hcalx_exp = D.Dot(hcal_unit_x);
  double hcaly_exp = D.Dot(hcal_unit_y);

  double dx = hcalx - hcalx_exp;
  double dy = hcaly - hcaly_exp;

  return {dx, dy};

}


double computePseudoMissingMass(TString target,
				double beam_energy,
				double hcal_angle,
				double hcal_distance,
				TVector3 kf,
				TVector3 v,
				double hcalx,
				double hcaly) {

  // Some constant(s)
  double Mp = 0.938272088; // PDG 2025
  double Mn = 0.939565420; // PDG 2025
  double MN = 0.5*(Mp + Mn);
  double pi = 3.1415926535;

  TVector3 hcal_unit_z(-TMath::Sin(hcal_angle*pi/180.0),
		       0.0,
		       TMath::Cos(hcal_angle*pi/180.0));
  TVector3 hcal_unit_x(0,-1.0,0);
  TVector3 hcal_unit_y = hcal_unit_z.Cross(hcal_unit_x).Unit();
  TVector3 hcal_vector = hcal_distance*hcal_unit_z;

  TVector3 ki(0.0,0.0,beam_energy);
  TLorentzVector pi_3he4(0.0, 0.0, 0.0, 2*Mp + Mn);
  TVector3 q = ki - kf;
  TLorentzVector q4(q, beam_energy + MN);
  TVector3 q_unit = q.Unit();

  // computing Missing Mass sq

  TVector3 hcal_hit = hcalx*hcal_unit_x + hcaly*hcal_unit_y;
  TVector3 pf_hit_unit = (hcal_hit + hcal_vector - v).Unit();
  TVector3 pf_hit = pf_hit_unit*(q.Mag());
  double pf_e = sqrt(pow(pf_hit.Mag(),2) + MN*MN);
  TLorentzVector pf_hit4(pf_hit, pf_e);

  double pf_para = q_unit.Dot(q - pf_hit);
  double pf_perp = (q - pf_hit - q_unit*pf_para).Mag();
  double pf_missing_mass_sq = (pi_3he4 + q4 - pf_hit4).M2();
  
  return pf_missing_mass_sq;
}

std::vector<double> computePseudoMissingMommentum(TString target,
						  double beam_energy,
						  double hcal_angle,
						  double hcal_distance,
						  TVector3 kf,
						  TVector3 v,
						  double hcalx,
						  double hcaly) {

  // Some constant(s)
  double Mp = 0.938272088; // PDG 2025
  double Mn = 0.939565420; // PDG 2025
  double MN = 0.5*(Mp + Mn);
  double pi = 3.1415926535;

  TVector3 hcal_unit_z(-TMath::Sin(hcal_angle*pi/180.0),
		       0.0,
		       TMath::Cos(hcal_angle*pi/180.0));
  TVector3 hcal_unit_x(0,-1.0,0);
  TVector3 hcal_unit_y = hcal_unit_z.Cross(hcal_unit_x).Unit();
  TVector3 hcal_vector = hcal_distance*hcal_unit_z;

  TVector3 ki(0.0,0.0,beam_energy);
  TVector3 q = ki - kf;
  TVector3 q_unit = q.Unit();

  // computing Missing Mass sq

  TVector3 hcal_hit = hcalx*hcal_unit_x + hcaly*hcal_unit_y;
  TVector3 pf_hit_unit = (hcal_hit + hcal_vector - v).Unit();
  TVector3 pf_hit = pf_hit_unit*(q.Mag());
  double pf_e = sqrt(pow(pf_hit.Mag(),2) + MN*MN);

  double pf_para = q_unit.Dot(q - pf_hit);
  double pf_perp = (q - pf_hit - q_unit*pf_para).Mag();
  
  return {pf_para, pf_perp};
}
