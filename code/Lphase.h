#ifndef _LPHASE_H_
#define _LPHASE_H_

#include <cmath>
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/Interpolator.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"

namespace LPHASE{
  double VolumePHS2(const double E, const double mass[2]);
  double VolumePHS3(const double E, const double mass[3]);
  double E1E2Range3(const double * x, const double * par);
  double VolumePHS(const double E, const double * mass, const int Nf);

}

double LPHASE::VolumePHS2(const double E, const double mass[2]){
  if (E <= mass[0] + mass[1])
    return 0;
  double p = sqrt( (E * E - pow(mass[0] + mass[1], 2)) * (E * E - pow(mass[0] - mass[1], 2))) / (2.0 * E);
  return M_PI * p / E;
}

double LPHASE::E1E2Range3(const double * x, const double * par){
  double E1 = x[0];
  double E = par[0];
  double m1 = par[1];
  double m2 = par[2];
  double m3 = par[3];
  if (E - E1 < m2 + m3 || E1 < m1)
    return 0;
  double M23 = sqrt(E * E - 2.0 * E * E1 + m1 * m1);
  double k = sqrt( (M23 * M23 - pow(m2 + m3, 2)) * (M23 * M23 - pow(m2 - m3, 2))) / (2.0 * M23);
  double gamma = (E - E1) / M23;
  double beta = sqrt(E1 * E1 - m1 * m1) / (E - E1);
  double E2max = gamma * (sqrt(m2 * m2 + k * k) + beta * k);
  double E2min = m2;
  if (beta > k / sqrt(m2 * m2 + k * k))
    E2min = gamma * (sqrt(m2 * m2 + k * k) - beta * k);
  return E2max - E2min;
}

double LPHASE::VolumePHS3(const double E, const double mass[3]){
  if (E <= mass[0] + mass[1] + mass[2])
    return 0;
  double par[4] = {E, mass[0], mass[1], mass[2]};
  double E1min = mass[0];
  double E1max = (E * E + mass[0] * mass[0] - pow(mass[1] + mass[2], 2)) / (2.0 * E);
  TF1 f0("E1 integrand", LPHASE::E1E2Range3, E1min, E1max, 4);
  f0.SetParameters(par);
  ROOT::Math::WrappedTF1 wf0(f0);
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 0.0, 1.0e-6, 1000);
  ig.SetFunction(wf0);
  double result = ig.Integral(E1min, E1max);
  return M_PI * M_PI * result;
}

double LPHASE::VolumePHS(const double E, const double * mass, const int Nf){
  if (Nf == 2)
    return LPHASE::VolumePHS2(E, mass);
  if (Nf == 3)
    return LPHASE::VolumePHS3(E, mass);
  return 0;
}


#endif
