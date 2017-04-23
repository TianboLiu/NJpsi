#include <iostream>
#include <fstream>

#include "Lcore.h" 

int main(){

  //MODEL::SetVeff();
  //MODEL::SetUr();

  //MODEL::CalculateFQ();

  //MODEL::SetFQ();
  //cout << MODEL::FQ(0.0) << endl;
  //cout << MODEL::FQ(0.96) << endl;

  //cout << GOLD::TF_fEnergy.Integral(0.0, 0.3) << endl;

  //cout << MODEL::TF_fMass.Integral(MODEL::TF_fMass.GetXmin(), MODEL::TF_fMass.GetXmax()) << endl;

  MODEL::SetMODEL();
  GOLD::SetGOLD();

  TLorentzVector ki[2], kf[2];
  ki[0].SetXYZT(0.0, 0.0, 8.8, 8.8);
  ki[1].SetXYZT(0.0, 0.0, 0.0, Mp);

  TLorentzVector kii[2];
  kii[0].SetXYZM(0, 0, 9.0, MJpsi);


  //GENERATE::JpsiPhotoproduction(ki, kf);

  double wt;
  for (int i = 0; i < 1000; i++){
    //cout << GENERATE::BoundStateFormationGold(kii, kf) << "\t"  << kf[0].M() << endl;
    //cout << GENERATE::JpsiPhotoproductionGold(ki, kf) << "\t"  << kf[0].Z() << endl;
    wt = GENERATE::BoundStatePhotoproductionGold(ki, kf);
    if (wt > 0) cout << wt << " " << kf[0].M() << " " << kf[1].M() << endl;
  }

  return 0;
}
