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
  const double m0 = Mp * MJpsi / (Mp + MJpsi);//reduced mass

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


  int test(){
    return 0;
  }


}









#endif
