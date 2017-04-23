#ifndef _LCORE_H_
#define _LCORE_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "Math/Functor.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/Interpolator.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"

#include "Lparticle.h"

using namespace std;

const double Mp = PARTICLE::proton.M();
const double MJpsi = PARTICLE::Jpsi.M();

namespace MODEL{
  
  const double Eb = 0.004;//binding energy
  const double Mass = Mp + MJpsi - Eb;//bound state mass
  const double Width = 0.003627;//bound state width
  const double FractionNJpsi = 0.58;//NJpsi component fraction
  const double BrNJpsi = 0.0149;//branch ratio of decay from NJpsi component

  ROOT::Math::Interpolator Rr_INTER(300, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator Ur_INTER(300, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator Veff_INTER(100, ROOT::Math::Interpolation::kCSPLINE);

  int SetVeff(){//
    double x[100], y[100];
    ifstream infile("wave/V_r.dat");
    for (int i = 0; i < 100; i++)
      infile >> x[i] >> y[i];
    Veff_INTER.SetData(100, x, y);
    infile.close();
    return 0;
  }

  int SetUr(){//
    double x[300], y[300], z[300];
    ifstream infile("wave/wf_5cc_1.dat");
    for (int i = 0; i < 300; i++)
      infile >> x[i] >> y[i] >> z[i];
    Rr_INTER.SetData(300, x, y);
    Ur_INTER.SetData(300, x, z);
    infile.close();
    return 0;
  }

  double Veff(const double r){//effective potential
    double rfm = r * Phys::hbar;//convert GeV^-1 to fm
    if (rfm < 0 || rfm > 5.0)
      return 0;
    return Veff_INTER.Eval(rfm) / 1000.0;//in GeV
  }

  double Rr(const double r){//radial wave function
    double N = 12.1729018;//Normalization factor
    double rfm = r * Phys::hbar;//convert GeV^-1 to fm
    if (rfm < 0 || rfm > 15.0)
      return 0;
    return N * Rr_INTER.Eval(rfm);//in GeV^3/2
  }

  double Ur(const double r){//radial wave function
    double N = 61.68892;//Normalization factor
    double rfm = r * Phys::hbar;//convert GeV^-1 to fm
    if (rfm < 0 || rfm > 15.0)
      return 0;
    return N * Ur_INTER.Eval(rfm);//in GeV^1/2
  }

  double FQ_integrand(const double r, void * par){
    double * k = (double *) par;
    if (k[0] == 0)
      return -r * Veff(r) * Ur(r);
    return -sin(k[0] * r) / k[0] * Veff(r) * Ur(r);
  }

  double FQk(double k){
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-4);
    ig.SetFunction(&FQ_integrand, &k);
    double result = ig.Integral(0.0, 30.0);
    return result / (4.0 * M_PI);
  }

  int CalculateFQ(){
    FILE * fp = fopen("wave/FQ0.dat", "w");
    double k, fk;
    for (int i = 0; i < 500; i++){
      k = i * 0.002;
      fk = FQk(k);
      cout << k << "   " << fk << endl;
      fprintf(fp, "%.6E\t%.6E\n", k, fk);
    }
    fclose(fp);
    return 0;
  }
      
  ROOT::Math::Interpolator FQ_INTER(500, ROOT::Math::Interpolation::kCSPLINE);
  int SetFQ(){
    ifstream infile("wave/FQ0.dat");
    double x[500], y[500];
    for (int i = 0; i < 500; i++)
      infile >> x[i] >> y[i];
    FQ_INTER.SetData(500, x, y);
    infile.close();
    return 0;
  }

  double FQ(const double k){
    if (k < 0.0 || k > 0.99)
      return 0;
    return FQ_INTER.Eval(k);
  }
  
  double BreitWigner(const double * E, const double * par){
    return 1.0 / (pow(E[0] * E[0] - Mass * Mass, 2) + E[0] * E[0] * Width * Width) / 25.8032;
  }

  TF1 TF_fMass("fM", BreitWigner, Mass - 5.0 * Width, Mass + 5.0 * Width, 0);

  int SetMODEL(){
    SetVeff();
    SetUr();
    SetFQ();
    TF_fMass.SetNpx(1000);
    return 0;
  }

}

namespace GOLD{

  const double NA = 197.0;
  const double ProtonDensity = 79.0 / (4.0 * M_PI * pow(7.3, 3) / 3.0) * pow(Phys::hbar, 3);//GeV^3

  double fMomentum(const double * p0, const double * par = 0){//non-normalized
    double p = p0[0];//nucleon momentum in Au197 in unit of GeV
    if (p < 0.0){
      std::cerr << "Unphysical momentum value in GoldMomentum!" << std::endl;
      return -1.0;
    }
    double A0 = 58.3382;
    double A2 = 69.2938;
    double B2 = 7.82756;
    double result = (A0 + pow(A2 * p, 2)) * exp(-pow(B2 * p, 2));
    return p * p * result / 0.162508;//Normalized momentum distribution
  }
  
  double fEnergy(const double * E0, const double * par = 0){
    double E = E0[0];//nucleon missing energy in Au197 in unit of GeV
    if (E <= 0.0){
      std::cerr << "Unphysical energy value in GoldEnergy!" << std::endl;
      return -1.0;
    }
    double A1 = 1.73622;
    double a1 = 3.07375;
    double b1 = 0.645561;
    double A2 = 14.1433;
    double a2 = 0.795058;
    double result = A1 * atan(A2 * pow(E/0.01, a1)) * exp(-b1 * pow(E/0.01, a2));
    return result / 0.0433967;//Normalized missing energy distribution
  }

  TF1 TF_fMomentum("fp", fMomentum, 0.0, 1.0, 0);
  TF1 TF_fEnergy("fE", fEnergy, 0.0, 0.3, 0);
 
  int SetGOLD(){
    TF_fMomentum.SetNpx(1000);
    TF_fEnergy.SetNpx(1000);
    return 0;
  }

}

namespace GENERATE{

  TRandom3 random(0);
  double Weight = 0.0;
  
  int NucleonGold(TLorentzVector * P){
    double p = GOLD::TF_fMomentum.GetRandom();
    double cth = random.Uniform(-1.0, 1.0);
    double phi = random.Uniform(-M_PI, M_PI);
    double dE = GOLD::TF_fEnergy.GetRandom();
    P->SetXYZT(p * sqrt(1.0 - cth * cth) * cos(phi), p * sqrt(1.0 - cth * cth) * sin(phi), p * cth, sqrt(p * p + Mp * Mp) - dE);
    return 0;
  }

  double dSigmaJpsi2g(const double * t, const double * s){
    if (s[0] <= pow(Mp + MJpsi, 2))
	return 0;
    const double N2g = 7.5671e-4 / pow(Phys::hbar, 4);//GeV^-4
    const double RM = MJpsi / Phys::hbar;//unit 1
    const double x = (2.0 * Mp * MJpsi + MJpsi * MJpsi) / (s[0] - Mp * Mp);
    const double FF = exp(1.13 * t[0]);
    return N2g * pow(1.0 - x, 2) * FF / (16.0 * M_PI * RM * RM);//GeV^-4
  }

  double JpsiPhotoproduction(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma, N; kf: Jpsi, N'
    TLorentzVector Pout = ki[0] + ki[1];//Total
    const double s = Pout.M2();//c.m. energy square
    if (s <= pow(MJpsi + Mp, 2)){
      weight[0] = 0;
      return 0;
    }
    const double q = sqrt( (s - pow(ki[0].M() + ki[1].M(), 2)) * (s - pow(ki[0].M() - ki[1].M(), 2))) / (2.0 * Pout.M());
    const double Q = sqrt( (s - pow(MJpsi + Mp, 2)) * (s - pow(MJpsi - Mp, 2))) / (2.0 * Pout.M());
    const double t0 = pow(ki[0].M2() - ki[1].M2() - MJpsi * MJpsi + Mp * Mp, 2) / (4.0 * s) - pow(q - Q, 2);
    const double t1 = pow(ki[0].M2() - ki[1].M2() - MJpsi * MJpsi + Mp * Mp, 2) / (4.0 * s) - pow(q + Q, 2);
    const double t = random.Uniform(t1, t0);
    weight[0] = dSigmaJpsi2g(&t, &s) * (t0 - t1);
    double theta = asin(sqrt((t0 - t) / (4.0 * q * Q))) * 2.0;
    double phi = random.Uniform(-M_PI, M_PI);
    kf[0].SetXYZM(Q * sin(theta) * cos(phi), Q * sin(theta) * sin(phi), Q * cos(theta), MJpsi);
    kf[0].Boost(Pout.BoostVector());//Jpsi
    kf[1] = Pout - kf[0];//N'
    return weight[0];//GeV^-2
  }

  double JpsiPhotoproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma; kf: Jpsi, N'
    TLorentzVector ki1[2];
    ki1[0] = ki[0];//photon
    NucleonGold(&ki1[1]);//off-shell nucleon
    JpsiPhotoproduction(ki1, kf, weight);
    weight[0] *= GOLD::NA;
    return weight[0];//GeV^-2 
  }
  
  double BoundStateFormationGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){//
    //ki: Jpsi; kf: d
    const double Md = MODEL::TF_fMass.GetRandom();//bound state mass
    const double dE = GOLD::TF_fEnergy.GetRandom();//missing energy
    const double Et = ki[0].E() - dE;
    const double k = ki[0].P();//Jpsi momentum
    const double MM = Md * Md + k * k - Et * Et - Mp * Mp;
    const double cth = random.Uniform(-1.0, 1.0);
    const double phi = random.Uniform(-M_PI, M_PI);
    const double a = Et * Et - k * k * cth * cth;
    const double b = -MM * k * cth;
    const double c = Et * Et * Mp * Mp - MM * MM / 4.0;
    const double DD = b * b - 4.0 * a * c;
    if (DD < 0){
      weight[0] = 0;
      return 0;
    }
    double p2;
    if (a * c < 0)
      p2 = (-b + sqrt(DD)) / (2.0 * a);
    else
      p2 = (-b - sqrt(DD)) / (2.0 * a);
    if (p2 < 0){
      weight[0] = 0;
      return 0;
    }
    weight[0] = GOLD::fMomentum(&p2);
    TLorentzVector P2;
    P2.SetXYZT(p2 * sqrt(1.0 - cth * cth) * cos(phi), p2 * sqrt(1.0 - cth * cth) * sin(phi), p2 * cth, sqrt(Mp * Mp + p2 * p2) - dE);
    P2.RotateY(ki[0].Theta());
    P2.RotateZ(ki[0].Phi());
    const double Q = sqrt( (Md * Md - pow(ki[0].M() + P2.M(), 2)) * (Md * Md - pow(ki[0].M() - P2.M(), 2))) / (2.0 * Md);
    weight[0] *= pow(MODEL::FQ(Q), 2) * MODEL::FractionNJpsi;
    kf[0] = ki[0] + P2;
    weight[0] *= GOLD::ProtonDensity * ki[0].Gamma() / PARTICLE::Jpsi.Gamma();
    return weight[0];
  }

  double BoundStatePhotoproductionGold(const TLorentzVector * ki, TLorentzVector * kf, double * weight = &Weight){
    //ki: gamma; kf: N', d
    TLorentzVector kf1[2];
    double weight1;
    JpsiPhotoproductionGold(ki, kf1, &weight1);//produce Jpsi
    if (weight1 == 0){
      weight[0] = 0;
      return 0;
    }
    kf[0] = kf1[1];//N'
    double weight2;
    BoundStateFormationGold(&kf1[0], &kf[1], &weight2);//form bound state
    weight[0] = weight1 * weight2;
    return weight[0];
  }

    
}





#endif
