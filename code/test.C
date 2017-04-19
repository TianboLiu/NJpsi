#include <iostream>
#include <fstream>

#include "Lparticle.h"
#include "Lphase.h"

using namespace std;
using namespace PARTICLE;

int main(){

  double mass[3] = {0.14, 1.0, 0.14};
  cout << LPHASE::VolumePHS3(2.6, mass) << endl;

  cout << pi.M() << endl;

  return 0;
}
