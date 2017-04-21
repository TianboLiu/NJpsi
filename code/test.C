#include <iostream>
#include <fstream>

#include "Lparticle.h"
#include "Lphase.h"

using namespace std;
using namespace PARTICLE;

int main(){

  double mass[7] = {0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14};
  cout << LPHASE::VolumePHS3(4.6, mass) << endl;
  cout << LPHASE::VPHS(3.097, mass, 7) << endl;

  cout << pi.M() << endl;

  return 0;
}
