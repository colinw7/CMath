#ifndef CMATH_RAND_H
#define CMATH_RAND_H

namespace CMathRand {
  void   seedRand(int s);
  int    randInRange(int rmin, int rmax);
  double unitRand();
  double symmetricUnitRand();
  double randInRange(double rmin, double rmax);
}

#endif
