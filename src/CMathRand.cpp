#include <CMathRand.h>
#include <COSRand.h>

void
CMathRand::
seedRand(int s)
{
  COSRand::srand(s);
}

int
CMathRand::
randInRange(int rmin, int rmax)
{
  return COSRand::randIn(rmin, rmax);
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
