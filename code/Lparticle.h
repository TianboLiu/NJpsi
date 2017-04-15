/* A class defining particles
   Declare some particles
*/

#ifndef _LPARTICLE_H_
#define _LPARTICLE_H_

#include <cmath>

#include "TRandom.h"
#include "TRandom3.h"

/* physics constant */
namespace Phys{
  const double hbar = 0.1973269718;//GeV fm
  const double c = 299792458.0;// m/s
}

/* class */
class Lparticle{
 private:
  double kmass;//Breit-Wigner mass
  double kwidth;//Breit-Wigner width
  double klife;//mean life
  int kpid;//particle id
  double kj;//spin
  TRandom3 * kran;//a pointer to random number
 public:
  Lparticle();
  Lparticle(const double mass, const double width);
  Lparticle(const double mass, const double life, const int ctr);
  double Gamma();
  double GetLife();
  double GetMass();
  int GetPID();
  double GetSpin();
  double GetWidth();
  double J();
  double M();
  double RandomM();
  double RandomM(const double Mmin, const double Mmax);
  int SetMass(const double mass);
  int SetPID(const int pid);
  int SetSpin(const double j);
  int SetWidth(const double width);
  double Tau();
};

/* member functions */
Lparticle::Lparticle(){
  kmass = 0;
  kwidth = 0;
  klife = 1.0;
  kran = 0;
  kj = 0;
}

Lparticle::Lparticle(const double mass, const double width){
  kmass = mass;
  kwidth = width;
  if (width > 0)
    klife = Phys::hbar * 1.0e-15 / width / Phys::c;
  else
    klife = 1.0;
  kran = 0;
  kj = 0;
}

Lparticle::Lparticle(const double mass, const double life, const int ctr){
  kmass = mass;
  klife = life;
  kwidth = Phys::hbar * 1.0e-15 / life / Phys::c;
  kran = 0;
  kj = 0;
}

double Lparticle::Gamma(){
  return GetWidth();
}

double Lparticle::GetLife(){
  return klife;
}

double Lparticle::GetMass(){
  return kmass;
}

int Lparticle::GetPID(){
  return kpid;
}

double Lparticle::GetSpin(){
  return kj;
}

double Lparticle::GetWidth(){
  return kwidth;
}

double Lparticle::J(){
  return GetSpin();
}

double Lparticle::M(){
  return GetMass();
}

double Lparticle::RandomM(){
  double Mmin = kmass - 5.0 * kwidth;
  double Mmax = kmass + 5.0 * kwidth;
  return RandomM(Mmin, Mmax);
}

double Lparticle::RandomM(const double Mmin, const double Mmax){
  if (kwidth == 0)
    return kmass;
  if (kran == 0)
    kran = new TRandom3(0);
  double mass = kmass;
  do {
    mass = kran->Uniform(Mmin, Mmax);
  } while (kran->Uniform(0.0, 1.0) > pow(kmass * kwidth, 2) / (pow(mass * mass - kmass * kmass, 2) + pow(kmass * kwidth, 2)));
  return mass;
}

int Lparticle::SetMass(const double mass){
  kmass = mass;
  return 1;
}

int Lparticle::SetPID(const int pid){
  kpid = pid;
  return 1;
}

int Lparticle::SetSpin(const double j){
  kj = j;
  return 1;
}

int Lparticle::SetWidth(const double width){
  kwidth = width;
  return 1;
}

double Lparticle::Tau(){
  return GetLife();
}




/* particles */


#endif
