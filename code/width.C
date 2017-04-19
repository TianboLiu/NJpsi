#include <iostream>
#include <fstream>

#include "TString.h"

#include "Lparticle.h"
#include "Lphase.h"

using namespace std;
using namespace PARTICLE;
using namespace LPHASE;

Lparticle SetParticle(TString str);

int main(const int argc, const char * argv[]){

  if (argc < 2){

    return 0;
  }

  int opt = atoi(argv[1]);

  if (opt == 1){
    ifstream input(argv[2]);
    TString str0, str1, str2, str3, str4;
    double E, mass[3];
    double br, w;
    double free, bind;
    double dE = 0.004;
    Lparticle parent, child1, child2, child3;
    int Nc;
    FILE * fp = fopen("bindingwidth.dat", "w");
    while (input >> str0 >> str1 >> str2 >> str3 >> str4){
      parent = SetParticle(str0);
      child1 = SetParticle(str1);
      child2 = SetParticle(str2);
      child3 = SetParticle(str3);
      br = str4.Atof() / 100.0;
      w = parent.Gamma() * br;
      E = parent.M();
      mass[0] = child1.M();
      mass[1] = child2.M();
      mass[2] = child3.M();
      if (str3.EqualTo("--")) Nc = 2;
      else Nc = 3;
      free = VolumePHS(E, mass, Nc);
      bind = VolumePHS(E - dE, mass, Nc);
      fprintf(fp, "%s\t\t%s\t%s\t%s\t%.3f\t%.6f\t%.6f\n", str0.Data(), str1.Data(), str2.Data(), str3.Data(), br * 100.0, w * 1000.0, w * 1000.0 * bind / free);
    }
    fclose(fp);
    input.close();
  }




  return 0;
}

Lparticle SetParticle(TString str){
  if (str.EqualTo("e+")) return e$p;
  if (str.EqualTo("e-")) return e$m;
  if (str.EqualTo("mu+")) return mu$p;
  if (str.EqualTo("mu-")) return mu$m;
  if (str.EqualTo("pi")) return pi;
  if (str.EqualTo("pi+")) return pi$p;
  if (str.EqualTo("pi-")) return pi$m;
  if (str.EqualTo("pi0")) return pi$0;
  if (str.EqualTo("rho")) return rho;
  if (str.EqualTo("K+")) return K$p;
  if (str.EqualTo("K-")) return K$m;
  if (str.EqualTo("K0")) return K$S;
  if (str.EqualTo("KS")) return K$S;
  if (str.EqualTo("KL")) return K$L;
  if (str.EqualTo("K+*")) return K892$p;
  if (str.EqualTo("K-*")) return K892$m;
  if (str.EqualTo("p")) return proton;
  if (str.EqualTo("n")) return neutron;
  if (str.EqualTo("Jpsi")) return Jpsi;
  if (str.EqualTo("Lambda")) return Lambda;
  if (str.EqualTo("Lambdac")) return Lambdac;
  if (str.EqualTo("D+")) return D$p;
  if (str.EqualTo("D-")) return D$m;
  if (str.EqualTo("D0")) return D$0;
  if (str.EqualTo("D+*")) return D2010$p;
  if (str.EqualTo("D-*")) return D2010$m;
  if (str.EqualTo("D0*")) return D2007;
  if (str.EqualTo("Sigma+")) return Sigma$p;
  if (str.EqualTo("Sigma0")) return Sigma$0;
  if (str.EqualTo("Sigma-")) return Sigma$m;
  if (str.EqualTo("Sigmac++")) return Sigmac$pp;
  if (str.EqualTo("Sigmac+")) return Sigmac$p;
  if (str.EqualTo("Sigmac0")) return Sigmac$0;
  if (str.EqualTo("Sigmac++*")) return Sigmac2520$pp;
  if (str.EqualTo("Sigmac+*")) return Sigmac2520$p;
  if (str.EqualTo("Sigmac0*")) return Sigmac2520$0;
  if (str.EqualTo("g")) return pi$0;
  if (str.EqualTo("gamma")) return photon;


  return Lparticle(0.0, 0.0);
}
