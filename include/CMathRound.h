#ifndef CMATH_ROUND_H
#define CMATH_ROUND_H

namespace CMathRound {
  enum Rounding {
    ROUND_DOWN,
    ROUND_UP,
    ROUND_NEAREST
  };

  int Round(double x, Rounding rounding=ROUND_NEAREST);

  int RoundNearest(double x);
  int RoundUp     (double x);
  int RoundDown   (double x);

  double RoundF(double x, Rounding rounding=ROUND_NEAREST);

  double RoundNearestF(double x);
  double RoundUpF     (double x);
  double RoundDownF   (double x);
}

#endif
