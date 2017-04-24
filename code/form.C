#include "Lcore.h"

int main(const int argc, const char * argv[]){
  
  if (argc < 2) return 0;

  double pz = atof(argv[1]);

  MODEL::SetMODEL();
  GOLD::SetGOLD();

  TLorentzVector ki[2], kf[2];
  ki[0].SetXYZM(0, 0, pz, MJpsi);
  
  Long64_t Nsim = 10000000;

  TH1D * h0 = new TH1D("h0", "", 1, 0.0, 2.0);

  double weight;

  for (Long64_t i = 0; i < Nsim; i++){

    weight = GENERATE::BoundStateFormationGold(ki, kf);
    if (weight > 0){
      h0->Fill(1.0, weight);
    }
  }

  cout << pz << "  " << h0->Integral(1, -1) / Nsim << endl;
  
  return 0;
}
