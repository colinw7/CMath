#include <CMathGen.h>
#include <CMathMacros.h>
#include <CThrow.h>

#ifdef CMATH_NAN
#include <COSNaN.h>
#endif

#include <climits>
#include <cassert>

double CMathGen::EPSILON_E4 = 1E-4;
double CMathGen::EPSILON_E5 = 1E-5;
double CMathGen::EPSILON_E6 = 1E-6;

float  CMathGen::EPSILON_E4_F = 1E-4f;
float  CMathGen::EPSILON_E5_F = 1E-5f;
float  CMathGen::EPSILON_E6_F = 1E-6f;

double
CMathGen::
DSin(double angle)
{
  return sin(DEG_TO_RAD(angle));
}

double
CMathGen::
DCos(double angle)
{
  return cos(DEG_TO_RAD(angle));
}

double
CMathGen::
DTan(double angle)
{
  return tan(DEG_TO_RAD(angle));
}

//------

int
CMathGen::
fastDistance(int x, int y)
{
  x = abs(x);
  y = abs(y);

  int mn = std::min(x, y);

  return (x + y - (mn >> 1) - (mn >> 2) + (mn >> 4));
}

double
CMathGen::
fastDistance(double fx, double fy, double fz)
{
  int x = int(fabs(fx) * 1024);
  int y = int(fabs(fy) * 1024);
  int z = int(fabs(fz) * 1024);

  if (y < x) std::swap(x, y);
  if (z < y) std::swap(y, z);
  if (y < x) std::swap(x, y);

  int dist = (z + 11*(y >> 5) + (x >> 2));

  return ((double)(dist >> 10));
}

//------

double
CMathGen::
pow(double real1, double real2)
{
  double real;

#ifdef CMATH_NAN
  if (COSNaN::is_nan(real1) || COSNaN::is_nan(real2)) {
    COSNaN::set_nan(&real);
    return real;
  }
#endif

  if (real1 < 0.0 && ! isInteger(real2))
    CTHROW("Non Integer Power of Negative");

  if (real1 == 0.0 && real2 < 0.0)
    CTHROW("Zero to Negative Power is Undefined");

  errno = 0;

  if (real2 < 0.0)
    real = 1.0/::pow(real1, -real2);
  else
    real = ::pow(real1, real2);

  if (errno != 0)
    CTHROW("Power Failed");

  return real;
}

double
CMathGen::
pow(long integer1, long integer2)
{
  if (integer1 == 0 && integer2 < 0)
    CTHROW("Zero to Negative Power is Undefined");

  errno = 0;

  double real;

  if (integer2 < 0.0)
    real = 1.0/::pow((double) integer1, (double) -integer2);
  else
    real = ::pow((double) integer1, (double) integer2);

  if (errno != 0)
    CTHROW("Power Failed");

  return real;
}

double
CMathGen::
ipow(double value, uint power)
{
  if (power == 0)
    return 1.0;

  double pvalue = value;

  for (uint i = 1; i < power; ++i)
    pvalue *= value;

  return pvalue;
}

//------

double
CMathGen::
log(double real)
{
  double result;

#ifdef CMATH_NAN
  if (COSNaN::is_nan(real)) {
    COSNaN::set_nan(&result);
    return real;
  }
#endif

  if (real <= 0.0)
    CTHROW("Negative or Zero Value for Log");

  result = ::log(real);

  return result;
}

double
CMathGen::
logN(int n, double real)
{
  double result;

#ifdef CMATH_NAN
  if (COSNaN::is_nan(real)) {
    COSNaN::set_nan(&result);
    return real;
  }
#endif

  if (real <= 0.0)
    CTHROW("Negative or Zero Value for Log");

  result = ::log(real)/log(n);

  return result;
}

double
CMathGen::
log10(double real)
{
  double result;

#ifdef CMATH_NAN
  if (COSNaN::is_nan(real)) {
    COSNaN::set_nan(&result);
    return real;
  }
#endif

  if (real <= 0.0)
    CTHROW("Negative or Zero Value for Log");

  result = ::log10(real);

  return result;
}

//------

bool
CMathGen::
isPowerOf(uint base, uint value)
{
  assert(base > 1);

  if (base == 2) return isPowerOf2(value);

  if (value == 0) return false;

  double n = logN(base, value);

  double r = n - int(n);

  return (r < 1E-10);
}

bool
CMathGen::
isPowerOf2(uint x)
{
  // extra methods at:
  //  http://www.exploringbinary.com/ten-ways-to-check-if-an-integer-is-a-power-of-two-in-c

  // Decrement and Compare
  return ((x != 0) && !(x & (x - 1)));

  // Complement and Compare
  //return ((x != 0) && ((x & (~x + 1)) == x));
}

//------

int
CMathGen::
numSetBits(uint x)
{
  uint count  = 0;

  while (x) {
    x = clearLastBit(x);

    ++count;
  }

  return count;
}

//------

uint
CMathGen::
clearLastBit(uint x)
{
  return (x & (x - 1));
}

//------

double
CMathGen::
sqrt(double real)
{
  double result;

#ifdef CMATH_NAN
  if (COSNaN::is_nan(real)) {
    COSNaN::set_nan(&result);
    return result;
  }
#endif

  if (real < 0.0)
    CTHROW("Negative Value for Sqrt");

  result = ::sqrt(real);

  return result;
}

//------

double
CMathGen::
acos(double x)
{
  if      (x >=  1.0)
    return 0.0;
  else if (x <= -1.0)
    return M_PI;
  else
    return ::acos(x);
}

double
CMathGen::
asin(double x)
{
  if      (x >=  1.0)
    return  PI_DIV_2;
  else if (x <= -1.0)
    return -PI_DIV_2;
  else
    return ::asin(x);
}

// note: usual argument order of atan2 flipped - ugh
float
CMathGen::
atan2(float x, float y)
{
  if      (fabsf(x) < EPSILON_E6_F) {
    if      (fabsf(y) < EPSILON_E6_F)
      return 0.0f;
    else if (y > 0.0f)
      return PI_DIV_2_F;
    else
      return PI3_DIV_2_F;
  }
  else if (fabsf(y) < EPSILON_E6_F) {
    if      (fabsf(x) < EPSILON_E6_F)
      return 0.0f;
    else if (x > 0.0f)
      return 0.0f;
    else
      return PI_F;
  }
  else {
    float angle = ::atan2(y, x);

    if (angle >= 0.0f)
      return angle;
    else
      return (PI2_F + angle);
  }
}

// note: usual argument order of atan2 flipped - ugh
double
CMathGen::
atan2(double x, double y)
{
  if      (fabs(x) < EPSILON_E6) {
    if      (fabs(y) < EPSILON_E6)
      return 0.0;
    else if (y > 0.0)
      return PI_DIV_2;
    else
      return PI3_DIV_2;
  }
  else if (fabs(y) < EPSILON_E6) {
    if      (fabs(x) < EPSILON_E6)
      return 0.0;
    else if (x > 0.0)
      return 0.0;
    else
      return M_PI;
  }
  else {
    double angle = ::atan2(y, x);

    if (angle >= 0.0)
      return angle;
    else
      return (PI2 + angle);
  }
}

//------

int
CMathGen::
sign(long x)
{
  if      (x > 0)
    return 1;
  else if (x < 0)
    return -1;
  else
    return 0;
}

int
CMathGen::
sign(double x)
{
  if      (x >  EPSILON_E6)
    return 1;
  else if (x < -EPSILON_E6)
    return -1;
  else
    return 0;
}

//------

long
CMathGen::
modulus(long integer1, long integer2)
{
  long result = 0;

  if (integer2 == 0.0)
    CTHROW("Divide by Zero");

  result = integer1 % integer2;

  return result;
}

double
CMathGen::
modulus(double real1, double real2)
{
  double result = 0;

#ifdef CMATH_NAN
  if (COSNaN::is_nan(real1) || COSNaN::is_nan(real2)) {
    COSNaN::set_nan(&result);
    return result;
  }
#endif

  if (real2 == 0.0)
    CTHROW("Divide by Zero");

  int factor = (int) (real1/real2);

  result = real1 - (real2*factor);

  return result;
}

//------

uint
CMathGen::
GetHCF(uint a, uint b)
{
  if (b == 0)
    return 0;

  if (a < b)
    std::swap(a, b);

  uint remainder = a % b;

  while (remainder != 0) {
    a = b;
    b = remainder;

    remainder = a % b;
  }

  return b;
}

ulong
CMathGen::
GetHCF(ulong a, ulong b)
{
  if (b == 0)
    return 0;

  if (a < b)
    std::swap(a, b);

  ulong remainder = a % b;

  while (remainder != 0) {
    a = b;
    b = remainder;

    remainder = a % b;
  }

  return b;
}

uint
CMathGen::
GetLCM(uint a, uint b)
{
  uint hcf = GetHCF(a, b);

  if (hcf == 0)
    return 0;

  return (a/hcf)*b;
}

ulong
CMathGen::
GetLCM(ulong a, ulong b)
{
  ulong hcf = GetHCF(a, b);

  if (hcf == 0)
    return 0;

  return (a/hcf)*b;
}

//------

double
CMathGen::
DegToRad(double deg)
{
  return deg*M_PI/180.0;
}

double
CMathGen::
RadToDeg(double rad)
{
  return rad*180.0/M_PI;
}

//------

bool
CMathGen::
isInteger(double r)
{
  return (long(r) == r);
}

//------

bool
CMathGen::
cmp(double a, double b, double prec)
{
  return (std::fabs(a - b) < prec);
}

//------

uint
CMathGen::
twosCompliment(int value)
{
  return ~((uint) -value) + 1;
}

ushort
CMathGen::
twosCompliment(short value)
{
  return ~((ushort) -value) + 1;
}

//------

bool
CMathGen::
range(double *x, double *y, uint num_xy,
      double *xmin, double *ymin, double *xmax, double *ymax)
{
  if (num_xy < 1)
    return false;

  *xmin = x[0];
  *ymin = y[0];
  *xmax = x[0];
  *ymax = y[0];

  for (uint i = 1; i < num_xy; ++i) {
    *xmin = std::min(*xmin, x[i]);
    *ymin = std::min(*ymin, y[i]);
    *xmax = std::max(*xmax, x[i]);
    *ymax = std::max(*ymax, y[i]);
  }

  return true;
}

//------

// Calculate p = sqrt(a*a + b*b)
//
// Uses most accurate computer calculation method
double
CMathGen::
pythag(double a, double b)
{
  double p = 0.0;

  double at = fabs(a);
  double bt = fabs(b);

  if      (at > bt) {
    double ct = bt/at;

    p = at*sqrt(1.0 + ct*ct);
  }
  else if (bt > 0.0) {
    double ct = at/bt;

    p = bt*sqrt(1.0 + ct*ct);
  }

  return p;
}

//------

float
CMathGen::
FastInvSqrtF(float fValue)
{
  float fHalf = 0.5f*fValue;

  int i = *(int *) &fValue;

  i = 0x5f3759df - (i >> 1);

  fValue = *(float*) &i;

  fValue = fValue*(1.5f - fHalf*fValue*fValue);

  return fValue;
}

#if 0
double
CMathGen::
FastInvSqrt(double dValue)
{
  double dHalf = 0.5*dValue;
  long long i = *(long long*) &dValue;
  i = 0x5fe6ec85e7de30da - (i >> 1);
  dValue = *(double*) &i;
  dValue = dValue*(1.5 - dHalf*dValue*dValue);
  return dValue;
}
#endif

double
CMathGen::
FastSin0(double fAngle)
{
  double fASqr = fAngle*fAngle;

  double fResult = (double)7.61e-03;
  fResult *= fASqr;

  fResult -= (double)1.6605e-01;
  fResult *= fASqr;

  fResult += (double)1.0;
  fResult *= fAngle;

  return fResult;
}

double
CMathGen::
FastSin1(double fAngle)
{
  double fASqr = fAngle*fAngle;

  double fResult = -(double)2.39e-08;
  fResult *= fASqr;

  fResult += (double)2.7526e-06;
  fResult *= fASqr;

  fResult -= (double)1.98409e-04;
  fResult *= fASqr;

  fResult += (double)8.3333315e-03;
  fResult *= fASqr;

  fResult -= (double)1.666666664e-01;
  fResult *= fASqr;

  fResult += (double)1.0;
  fResult *= fAngle;

  return fResult;
}

double
CMathGen::
FastCos0(double fAngle)
{
  double fASqr = fAngle*fAngle;

  double fResult = (double)3.705e-02;
  fResult *= fASqr;

  fResult -= (double)4.967e-01;
  fResult *= fASqr;

  fResult += (double)1.0;

  return fResult;
}

double
CMathGen::
FastCos1(double fAngle)
{
  double fASqr = fAngle*fAngle;

  double fResult = -(double)2.605e-07;
  fResult *= fASqr;

  fResult += (double)2.47609e-05;
  fResult *= fASqr;

  fResult -= (double)1.3888397e-03;
  fResult *= fASqr;

  fResult += (double)4.16666418e-02;
  fResult *= fASqr;

  fResult -= (double)4.999999963e-01;
  fResult *= fASqr;

  fResult += (double)1.0;

  return fResult;
}

double
CMathGen::
FastTan0(double fAngle)
{
  double fASqr = fAngle*fAngle;

  double fResult = (double)2.033e-01;
  fResult *= fASqr;

  fResult += (double)3.1755e-01;
  fResult *= fASqr;

  fResult += (double)1.0;
  fResult *= fAngle;

  return fResult;
}

double
CMathGen::
FastTan1(double fAngle)
{
  double fASqr = fAngle*fAngle;

  double fResult = (double)9.5168091e-03;
  fResult *= fASqr;

  fResult += (double)2.900525e-03;
  fResult *= fASqr;

  fResult += (double)2.45650893e-02;
  fResult *= fASqr;

  fResult += (double)5.33740603e-02;
  fResult *= fASqr;

  fResult += (double)1.333923995e-01;
  fResult *= fASqr;

  fResult += (double)3.333314036e-01;
  fResult *= fASqr;

  fResult += (double)1.0;
  fResult *= fAngle;

  return fResult;
}

double
CMathGen::
FastInvSin0(double fValue)
{
  double fRoot = sqrt(((double)1.0)-fValue);

  double fResult = -(double)0.0187293;
  fResult *= fValue;

  fResult += (double)0.0742610;
  fResult *= fValue;

  fResult -= (double)0.2121144;
  fResult *= fValue;

  fResult += (double)1.5707288;

  fResult = PI_DIV_2 - fRoot*fResult;

  return fResult;
}

double
CMathGen::
FastInvSin1(double fValue)
{
  double fRoot = sqrt(fabs(((double)1.0)-fValue));

  double fResult = -(double)0.0012624911;
  fResult *= fValue;

  fResult += (double)0.0066700901;
  fResult *= fValue;

  fResult -= (double)0.0170881256;
  fResult *= fValue;

  fResult += (double)0.0308918810;
  fResult *= fValue;

  fResult -= (double)0.0501743046;
  fResult *= fValue;

  fResult += (double)0.0889789874;
  fResult *= fValue;

  fResult -= (double)0.2145988016;
  fResult *= fValue;

  fResult += (double)1.5707963050;

  fResult = PI_DIV_2 - fRoot*fResult;

  return fResult;
}

double
CMathGen::
FastInvCos0(double fValue)
{
  double fRoot = sqrt(((double)1.0)-fValue);

  double fResult = -(double)0.0187293;
  fResult *= fValue;

  fResult += (double)0.0742610;
  fResult *= fValue;

  fResult -= (double)0.2121144;
  fResult *= fValue;

  fResult += (double)1.5707288;
  fResult *= fRoot;

  return fResult;
}

double
CMathGen::
FastInvCos1(double fValue)
{
  double fRoot = sqrt(fabs(((double)1.0)-fValue));

  double fResult = -(double)0.0012624911;
  fResult *= fValue;

  fResult += (double)0.0066700901;
  fResult *= fValue;

  fResult -= (double)0.0170881256;
  fResult *= fValue;

  fResult += (double)0.0308918810;
  fResult *= fValue;

  fResult -= (double)0.0501743046;
  fResult *= fValue;

  fResult += (double)0.0889789874;
  fResult *= fValue;

  fResult -= (double)0.2145988016;
  fResult *= fValue;

  fResult += (double)1.5707963050;
  fResult *= fRoot;

  return fResult;
}

double
CMathGen::
FastInvTan0(double fValue)
{
  double fVSqr = fValue*fValue;

  double fResult = (double)0.0208351;
  fResult *= fVSqr;

  fResult -= (double)0.085133;
  fResult *= fVSqr;

  fResult += (double)0.180141;
  fResult *= fVSqr;

  fResult -= (double)0.3302995;
  fResult *= fVSqr;

  fResult += (double)0.999866;
  fResult *= fValue;

  return fResult;
}

double
CMathGen::
FastInvTan1(double fValue)
{
  double fVSqr = fValue*fValue;

  double fResult = (double)0.0028662257;
  fResult *= fVSqr;

  fResult -= (double)0.0161657367;
  fResult *= fVSqr;

  fResult += (double)0.0429096138;
  fResult *= fVSqr;

  fResult -= (double)0.0752896400;
  fResult *= fVSqr;

  fResult += (double)0.1065626393;
  fResult *= fVSqr;

  fResult -= (double)0.1420889944;
  fResult *= fVSqr;

  fResult += (double)0.1999355085;
  fResult *= fVSqr;

  fResult -= (double)0.3333314528;
  fResult *= fVSqr;

  fResult += (double)1.0;
  fResult *= fValue;

  return fResult;
}

//------

double
CMathGen::
mapToReal(char i)
{
  static double factor1 = -1.0/CHAR_MIN;
  static double factor2 =  1.0/CHAR_MAX;

  if (i < 0)
    return 0.5*(i*factor1 + 1.0);
  else
    return 0.5*(i*factor2 + 1.0);
}

double
CMathGen::
mapToReal(unsigned char i)
{
  static double factor = 1.0/UCHAR_MAX;

  return i*factor;
}

double
CMathGen::
mapToReal(short i)
{
  static double factor1 = -1.0/SHRT_MIN;
  static double factor2 =  1.0/SHRT_MAX;

  if (i < 0)
    return 0.5*(i*factor1 + 1.0);
  else
    return 0.5*(i*factor2 + 1.0);
}

double
CMathGen::
mapToReal(ushort i)
{
  static double factor = 1.0/USHRT_MAX;

  return i*factor;
}

double
CMathGen::
mapToReal(int i)
{
  static double factor1 = -1.0/INT_MIN;
  static double factor2 =  1.0/INT_MAX;

  if (i < 0)
    return 0.5*(i*factor1 + 1.0);
  else
    return 0.5*(i*factor2 + 1.0);
}

double
CMathGen::
mapToReal(uint i)
{
  static double factor = 1.0/UINT_MAX;

  return i*factor;
}

double
CMathGen::
mapToReal(long i)
{
  static double factor1 = -1.0/LONG_MIN;
  static double factor2 =  1.0/LONG_MAX;

  if (i < 0)
    return 0.5*(i*factor1 + 1.0);
  else
    return 0.5*(i*factor2 + 1.0);
}

double
CMathGen::
mapToReal(ulong i)
{
  static double factor = 1.0/ULONG_MAX;

  return i*factor;
}

//------

ulong
CMathGen::
factorial(ushort f)
{
  enum { MAX_FACTORIAL = 21 };

  assert(f <= MAX_FACTORIAL);

  static ushort computed_value;
  static ulong  values[MAX_FACTORIAL + 1];

  if (! f) return 0;

  if (f > computed_value) {
    if (computed_value == 0) {
      values[1] = 1;

      computed_value = 1;
    }

    ushort last_value;

    while (computed_value < f) {
      last_value = computed_value++;

      values[computed_value] = computed_value*values[last_value];
    }
  }

  return values[f];
}

double
CMathGen::
factorial(double f)
{
  enum { MAX_FACTORIAL = 99 };

  uint fi = uint(f);

  assert(fi <= MAX_FACTORIAL);

  static ushort computed_value;
  static double values[MAX_FACTORIAL + 1];

  if (! fi) return 0;

  if (fi > computed_value) {
    if (computed_value == 0) {
      values[1] = 1;

      computed_value = 1;
    }

    ushort last_value;

    while (computed_value < fi) {
      last_value = computed_value++;

      values[computed_value] = computed_value*values[last_value];
    }
  }

  return values[fi];
}

//------

ulong
CMathGen::
binomialCoeff(uint n, uint k)
{
  enum { MAX_COEFF = 28 };

  assert(n <= MAX_COEFF && k <= n);

  if (! n) return 0;

  static ushort computed_value;
  static ulong  values[MAX_COEFF + 1][MAX_COEFF + 1];

  if (n > computed_value) {
    if (computed_value == 0) {
      values[1][0] = 1;
      values[1][1] = 1;

      computed_value = 1;
    }

    uint n2 = n/2;

    for (uint k1 = 0; k1 <= n2; ++k1) {
      ulong t = 1;
      ulong b = 1;

      for (uint i = 0; i < k1; ++i) {
        t *= n - i;
        b *= i + 1;
      }

      ulong tb = t/b;

      values[n][    k1] = tb;
      values[n][n - k1] = tb;
    }

    if (n & 1) {
      ulong t = 1;
      ulong b = 1;

      uint k1 = n2 + 1;

      for (uint i = 0; i < k1; ++i) {
        t *= n - i;
        b *= i + 1;
      }

      ulong tb = t/b;

      values[n][k1] = tb;
    }
  }

  return values[n][k];
}

//------

void
CMathGen::
getXYVals(const CIPoint2D *points, uint npoints,
          double **xvals, int *num_xvals, double **yvals, int *num_yvals)
{
  double *x = new double [npoints];
  double *y = new double [npoints];

  for (uint i = 0; i < npoints; ++i) {
    x[i] = points[i].x;
    y[i] = points[i].y;
  }

  getXYVals(x, y, npoints, xvals, num_xvals, yvals, num_yvals);

  delete [] x;
  delete [] y;
}

void
CMathGen::
getXYVals(double *x, double *y, int num_xy,
          double **xvals, int *num_xvals, double **yvals, int *num_yvals)
{
  /* Allocate return value for maximum possible size */

  *xvals = new double [num_xy];
  *yvals = new double [num_xy];

  /*------*/

  *num_xvals = 0;
  *num_yvals = 0;

  for (int i = 0; i < num_xy; ++i) {
    /* Search for matching x value (within tolerance)
       (could use sorted property to speed up search) */

    int j = 0;

    for ( ; j < *num_xvals; ++j)
      if (REAL_EQ((*xvals)[j], x[i]))
        break;

    /* No value found so add to array */

    if (j >= *num_xvals) {
      /* Find sorted position for add */
      int j = 0;

      for ( ; j < *num_xvals; ++j)
        if (x[i] < (*xvals)[j])
          break;

      /* Make room */

      for (int k = *num_xvals - 1; k >= j; --k)
        (*xvals)[k + 1] = (*xvals)[k];

      /* Add to array */

      (*xvals)[j] = x[i];

      ++(*num_xvals);
    }

    /*------*/

    /* Search for matching y value (within tolerance)
       (could used sorted property to speed up search) */

    j = 0;

    for ( ; j < *num_yvals; ++j)
      if (REAL_EQ((*yvals)[j], y[i]))
        break;

    /* No value found so add to array */

    if (j >= *num_yvals) {
      /* Find sorted position for add */

      int j = 0;

      for ( ; j < *num_yvals; ++j)
        if (y[i] < (*yvals)[j])
          break;

      /* Make room */

      for (int k = *num_yvals - 1; k >= j; --k)
        (*yvals)[k + 1] = (*yvals)[k];

      /* Add to array */

      (*yvals)[j] = y[i];

      ++(*num_yvals);
    }
  }
}

double
CMathGen::
getNaN()
{
  double r;

  SetNaN(r);

  return r;
}

static double
noiseFade(double t)
{
  return t * t * t * (t * (t * 6 - 15) + 10);
}

static double
noiseLerp(double t, double a, double b)
{
  return a + t * (b - a);
}

static double
noiseGrad(int hash, double x, double y, double z)
{
  int h = hash & 15;                      // CONVERT LO 4 BITS OF HASH CODE

  double u = h<8 ? x : y,                 // INTO 12 GRADIENT DIRECTIONS.
         v = h<4 ? y : h==12||h==14 ? x : z;

  return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
}

static int permutation[] =
{ 151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,
142,8,99,37,240,21,10,23,190,6,148,247,120,234,75,0,26,197,62,94,252,219,
203,117,35,11,32,57,177,33,88,237,149,56,87,174,20,125,136,171,168,68,175,
74,165,71,134,139,48,27,166,77,146,158,231,83,111,229,122,60,211,133,230,
220,105,92,41,55,46,245,40,244,102,143,54,65,25,63,161,1,216,80,73,209,76,
132,187,208,89,18,169,200,196,135,130,116,188,159,86,164,100,109,198,173,
186,3,64,52,217,226,250,124,123,5,202,38,147,118,126,255,82,85,212,207,
206,59,227,47,16,58,17,182,189,28,42,223,183,170,213,119,248,152,2,44,
154,163,70,221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,
113,224,232,178,185,112,104,218,246,97,228,251,34,242,193,238,210,144,
12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,214,31,181,199,
106,157,184,84,204,176,115,121,50,45,127,4,150,254,138,236,205,93,222,
114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
};

static int p[512];

static void
noiseInit()
{
  static bool initialized;

  if (! initialized) {
    for (int i = 0; i < 256 ; i++)
      p[256 + i] = p[i] = permutation[i];

    initialized = true;
  }
}

double
CMathGen::
noise(double x, double y, double z)
{
  noiseInit();

  int X = int(floor(x)) & 255, // FIND UNIT CUBE THAT
      Y = int(floor(y)) & 255, // CONTAINS POINT.
      Z = int(floor(z)) & 255;

  x -= floor(x); // FIND RELATIVE X,Y,Z
  y -= floor(y); // OF POINT IN CUBE.
  z -= floor(z);

  double u = noiseFade(x), // COMPUTE FADE CURVES
         v = noiseFade(y), // FOR EACH OF X,Y,Z.
         w = noiseFade(z);

  int A = p[X  ]+Y, AA = p[A]+Z, AB = p[A+1]+Z, // HASH COORDINATES OF
      B = p[X+1]+Y, BA = p[B]+Z, BB = p[B+1]+Z; // THE 8 CUBE CORNERS,

  return noiseLerp(w,
    noiseLerp(v, noiseLerp(u, noiseGrad(p[AA  ], x  , y  , z   ),  // AND ADD
                              noiseGrad(p[BA  ], x-1, y  , z   )), // BLENDED
                 noiseLerp(u, noiseGrad(p[AB  ], x  , y-1, z   ),  // RESULTS
                              noiseGrad(p[BB  ], x-1, y-1, z   ))),// FROM  8
    noiseLerp(v, noiseLerp(u, noiseGrad(p[AA+1], x  , y  , z-1 ),  // CORNERS
                              noiseGrad(p[BA+1], x-1, y  , z-1 )), // OF CUBE
                 noiseLerp(u, noiseGrad(p[AB+1], x  , y-1, z-1 ),
                              noiseGrad(p[BB+1], x-1, y-1, z-1 ))));
}

bool
CMathGen::
realEq(double v1, double v2, double tol)
{
  return (fabs(v1 - v2) < tol);
}
