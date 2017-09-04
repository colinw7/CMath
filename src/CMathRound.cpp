#include <CMathRound.h>
#include <climits>
#include <cerrno>

namespace {
  double EPSILON_E6 = 1E-6;
}

int
CMathRound::
Round(double x, Rounding rounding)
{
  switch (rounding) {
    case ROUND_UP  : return RoundUp(x);
    case ROUND_DOWN: return RoundDown(x);
    default        : return RoundNearest(x);
  }
}

int
CMathRound::
RoundNearest(double x)
{
  double x1;

  if (x <= 0.0)
    x1 = (x - 0.499999);
  else
    x1 = (x + 0.500001);

  if (x1 < INT_MIN || x1 > INT_MAX)
    errno = ERANGE;

  return int(x1);
}

int
CMathRound::
RoundUp(double x)
{
  double x1;

  if (x <= 0.0)
    x1 = (x       - EPSILON_E6);
  else
    x1 = (x + 1.0 - EPSILON_E6);

  if (x1 < INT_MIN || x1 > INT_MAX)
    errno = ERANGE;

  return int(x1);
}

int
CMathRound::
RoundDown(double x)
{
  double x1;

  if (x >= 0.0)
    x1 = (x       + EPSILON_E6);
  else
    x1 = (x - 1.0 + EPSILON_E6);

  if (x1 < INT_MIN || x1 > INT_MAX)
    errno = ERANGE;

  return int(x1);
}
