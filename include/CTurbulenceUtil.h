#ifndef CTurbulenceUtil_H
#define CTurbulenceUtil_H

class CTurbulenceUtil {
 private:
  static const int RAND_m = 2147483647; /* 2**31 - 1 */
  static const int RAND_a = 16807; /* 7**5; primitive root of m */
  static const int RAND_q = 127773; /* m / a */
  static const int RAND_r = 2836; /* m % a */

  static const int BSize   = 0x100;
  static const int BM      = 0xff;
  static const int PerlinN = 0x1000;
  static const int NP      = 12; /* 2^PerlinN */
  static const int NM      = 0xfff;

  //------

  int    latticeSelector[BSize + BSize + 2];
  double gradient[4][BSize + BSize + 2][2];

  //------

 public:
  CTurbulenceUtil(int seed) {
    init(seed);
  }

  /* Produces results in the range [1, 2**31 - 2].

     Algorithm is: r = (a * r) mod m
     where a = 16807 and m = 2**31 - 1 = 2147483647

     See [Park & Miller], CACM vol. 31 no. 10 p. 1195, Oct. 1988

     To test: the algorithm should produce the result 1043618065
     as the 10,000th generated number if the original seed is 1.
  */

 private:
  int initSeed(int lSeed) {
    if (lSeed <= 0)
      lSeed = -(lSeed % (RAND_m - 1)) + 1;

    if (lSeed > RAND_m - 1)
      lSeed = RAND_m - 1;

    return lSeed;
  }

  int random(int lSeed) {
    int result = RAND_a*(lSeed % RAND_q) - RAND_r*(lSeed / RAND_q);

    if (result <= 0)
      result += RAND_m;

    return result;
  }

  void init(int lSeed) {
    double s;

    lSeed = initSeed(lSeed);

    for (int k = 0; k < 4; k++) {
      for (int i = 0; i < BSize; i++) {
        latticeSelector[i] = i;

        for (int j = 0; j < 2; j++) {
          lSeed = random(lSeed);

          gradient[k][i][j] =
            double((lSeed) % (BSize + BSize) - BSize)/BSize;
        }

        s = sqrt(gradient[k][i][0]*gradient[k][i][0] +
                 gradient[k][i][1]*gradient[k][i][1]);

        gradient[k][i][0] /= s;
        gradient[k][i][1] /= s;
      }
    }

    int i = BSize;

    while (--i) {
      int k = latticeSelector[i];

      lSeed = random(lSeed);

      int j = lSeed % BSize;

      latticeSelector[i] = latticeSelector[j];
      latticeSelector[j] = k;
    }

    for (int i = 0; i < BSize + 2; i++) {
      latticeSelector[BSize + i] = latticeSelector[i];

      for (int k = 0; k < 4; k++)
        for (int j = 0; j < 2; j++)
          gradient[k][BSize + i][j] = gradient[k][i][j];
    }
  }

  double s_curve(double t) {
    return (t*t*(3.0 - 2.0*t));
  }

  double lerp(double t, double a, double b) {
    return (a + t*(b - a));
  }

  double noise2(int channelNum, double vec[2]) {
    double t = vec[0] + PerlinN;

    int bx0 = int(t)    & BM;
    int bx1 = (bx0 + 1) & BM;

    double rx0 = t - int(t);
    double rx1 = rx0 - 1.0;

    t = vec[1] + PerlinN;

    int by0 = int(t)    & BM;
    int by1 = (by0 + 1) & BM;

    double ry0 = t - int(t);
    double ry1 = ry0 - 1.0;

    int i = latticeSelector[bx0];
    int j = latticeSelector[bx1];

    int b00 = latticeSelector[i + by0];
    int b10 = latticeSelector[j + by0];
    int b01 = latticeSelector[i + by1];
    int b11 = latticeSelector[j + by1];

    double sx = s_curve(rx0);
    double sy = s_curve(ry0);

    double u, v, *q;

    q = gradient[channelNum][b00]; u = rx0*q[0] + ry0*q[1];
    q = gradient[channelNum][b10]; v = rx1*q[0] + ry0*q[1];

    double a = lerp(sx, u, v);

    q = gradient[channelNum][b01]; u = rx0*q[0] + ry1*q[1];
    q = gradient[channelNum][b11]; v = rx1*q[0] + ry1*q[1];

    double b = lerp(sx, u, v);

    return lerp(sy, a, b);
  }

  // Returns 'turbFunctionResult'

 public:
  double turbulence(int channelNum, double *point, double baseFreqX,
                    double baseFreqY, int numOctaves, bool fractal) {
    double vec[2];

    double sum        = 0;
    double frequencyX = baseFreqX;
    double frequencyY = baseFreqY;

    for (int nOctave = 0; nOctave < numOctaves; nOctave++) {
      vec[0] = frequencyX*point[0];
      vec[1] = frequencyY*point[1];

      double amplitudeX = baseFreqX/frequencyX;
    //double amplitudeY = baseFreqY/frequencyY;

      double sum1 = noise2(channelNum, vec)*amplitudeX;

      if (fractal)
        sum += sum1;
      else
        sum += fabs(sum1);

      frequencyX *= 2;
    //frequencyY *= 2;
    }

    return sum;
  }
};

#endif
