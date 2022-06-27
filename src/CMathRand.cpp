#include <CMathRand.h>
#include <COSRand.h>

#include <set>
#include <time.h>

static std::set<std::string> s_rand_strings;

void
CMathRand::
timeSeedRand()
{
  COSRand::srand(int(time(nullptr)));
}

void
CMathRand::
seedRand(int s)
{
  COSRand::srand(s);
}

int
CMathRand::
randInRange(int imin, int imax)
{
  return COSRand::randIn(imin, imax);
}

long
CMathRand::
randInRange(long lmin, long lmax)
{
  return COSRand::randIn(lmin, lmax);
}

double
CMathRand::
unitRand()
{
  return COSRand::randIn(0.0, 1.0);
}

double
CMathRand::
symmetricUnitRand()
{
  return COSRand::randIn(-1.0, 1.0);
}

double
CMathRand::
randInRange(double rmin, double rmax)
{
  return COSRand::randIn(rmin, rmax);
}

std::string
CMathRand::
randString()
{
  int n = randInRange(3, 6);

  std::string s;

  while (true) {
    s = "";

    // A-Z (65 to 90)
    // a-z (97 to 122)
    for (int i = 0; i < n; ++i) {
      auto n = randInRange(0, 51);

      auto c = char(n >= 26 ? n + 71 : n + 65);

      s += c;
    }

    if (s_rand_strings.find(s) == s_rand_strings.end()) {
      s_rand_strings.insert(s);
      break;
    }
  }

  return s;
}
