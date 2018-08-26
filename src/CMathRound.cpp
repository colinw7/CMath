#include <CMathRound.h>
#include <cmath>
#include <climits>
#include <cerrno>

namespace CMathRound {

double EPSILON_E6 = 1E-6;

int Round(double x, Rounding rounding) {
  switch (rounding) {
    case ROUND_UP  : return RoundUp     (x);
    case ROUND_DOWN: return RoundDown   (x);
    default        : return RoundNearest(x);
  }
}

int RoundNearest(double x) {
  double x1;

  if (x <= 0.0)
    x1 = (x - 0.499999);
  else
    x1 = (x + 0.500001);

  if (x1 < INT_MIN || x1 > INT_MAX)
    errno = ERANGE;

  return int(x1);
}

int RoundUp(double x) {
  double x1;

  if (x <= 0.0)
    x1 = (x       - EPSILON_E6);
  else
    x1 = (x + 1.0 - EPSILON_E6);

  if (x1 < INT_MIN || x1 > INT_MAX)
    errno = ERANGE;

  return int(x1);
}

int RoundDown(double x) {
  double x1;

  if (x >= 0.0)
    x1 = (x       + EPSILON_E6);
  else
    x1 = (x - 1.0 + EPSILON_E6);

  if (x1 < INT_MIN || x1 > INT_MAX)
    errno = ERANGE;

  return int(x1);
}

double RoundF(double x, Rounding rounding) {
  switch (rounding) {
    case ROUND_UP  : return RoundUpF     (x);
    case ROUND_DOWN: return RoundDownF   (x);
    default        : return RoundNearestF(x);
  }
}

double RoundNearestF(double x) {
  double x1;

  if (x <= 0.0)
    x1 = (x - 0.499999);
  else
    x1 = (x + 0.500001);

  return std::trunc(x1);
}

double RoundUpF(double x) {
  double x1;

  if (x <= 0.0)
    x1 = (x       - EPSILON_E6);
  else
    x1 = (x + 1.0 - EPSILON_E6);

  return std::trunc(x1);
}

double RoundDownF(double x) {
  double x1;

  if (x >= 0.0)
    x1 = (x       + EPSILON_E6);
  else
    x1 = (x - 1.0 + EPSILON_E6);

  return std::trunc(x1);
}

}
