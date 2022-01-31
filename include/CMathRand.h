#ifndef CMATH_RAND_H
#define CMATH_RAND_H

#include <string>

namespace CMathRand {
  // seed
  void seedRand(int s);
  void timeSeedRand();

  // int
  int randInRange(int rmin, int rmax);

  // real
  double unitRand();
  double symmetricUnitRand();
  double randInRange(double rmin, double rmax);

  // string
  std::string randString();
}

#endif
