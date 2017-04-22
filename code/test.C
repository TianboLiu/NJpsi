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
  ki[0].SetXYZM(0.0, 0.0, 1.0, PARTICLE::Jpsi.M());
  ki[1].SetXYZM(0.0, 0.0, 0.0, 196.0 * Mp);

  for (int i = 0; i < 10; i++){
    cout << GENERATE::BoundStateFormation(ki, kf) << endl;
  }

  return 0;
}
