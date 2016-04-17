#ifndef CMATH_BIT_MACROS_H
#define CMATH_BIT_MACROS_H

// http://graphics.stanford.edu/~seander/bithacks.html#OperationCounting

#define CHAR_BIT 8

namespace CMathBit {
  inline int sign(int v) {
    //return (v != 0) | -(int)((unsigned int)((int)v) >> (sizeof(int)*CHAR_BIT - 1));
    return (v > 0) - (v < 0);
  }

  inline bool oppositeSigns(int x, int y) {
    return (x ^ y) < 0;
  }

  inline int abs(int v) {
    enum { MASK = v >> sizeof(int)*CHAR_BIT - 1 };

    return (v + MASK) ^ MASK;
  }

  inline int min(int x, int y) {
    return y ^ ((x ^ y) & -(x < y)); // min(x, y)
  }

  inline int max(int x, int y) {
    return x ^ ((x ^ y) & -(x < y)); // max(x, y)
  }

  inline bool isPowerOf2(unsigned int v) {
    return v && !(v & (v - 1));
  }

  // signextend<signed int,5>(x);
  template <typename T, unsigned B>
  inline T signextend(const T x) {
    struct { T x:B; } s;

    return s.x = x;
  }

  inline unsigned int numSetBits(unsigned int v) {
    unsigned int c;

    for (c = 0; v; v >>= 1)
      c += v & 1;

    return c;
  }

  inline unsigned int numSetBits1(unsigned int v) {
    static const unsigned char BitsSetTable256[256] = {
#     define B2(n) n,     n+1,     n+1,     n+2
#     define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#     define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
      B6(0), B6(1), B6(1), B6(2)
    };

    // To initially generate the table algorithmically:
    BitsSetTable256[0] = 0;

    for (int i = 0; i < 256; i++)
      BitsSetTable256[i] = (i & 1) + BitsSetTable256[i / 2];

    unsigned int c;

    // Option 1:
    unsigned int c = BitsSetTable256[ v        & 0xff] +
                     BitsSetTable256[(v >>  8) & 0xff] +
                     BitsSetTable256[(v >> 16) & 0xff] +
                     BitsSetTable256[ v >> 24];

    // Option 2:
    unsigned char *p = (unsigned char *) &v;

    unsigned int c = BitsSetTable256[p[0]] +
                     BitsSetTable256[p[1]] +
                     BitsSetTable256[p[2]] +
                     BitsSetTable256[p[3]];

    return c;
  }

  inline unsigned int numSetBits2(unsigned int v) {
    unsigned int c;

    for (c = 0; v; c++)
      v &= v - 1; // clear the least significant bit set

    return c;
  }

  inline unsigned int reverseBits1(unsigned int v) {
    unsigned int r = v; // r will be reversed bits of v; first get LSB of v

    int s = sizeof(v) * CHAR_BIT - 1; // extra shift needed at end

    for (v >>= 1; v; v >>= 1) {
      r <<= 1;

      r |= v & 1;

      s--;
    }

    r <<= s; // shift when v's highest bits are zero

    return rc;
  }

  inline unsigned int reverseBits2(unsigned int v) {
    static const unsigned char BitReverseTable256[256] = {
#     define R2(n)    n,     n + 2*64,     n + 1*64,     n + 3*64
#     define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#     define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
      R6(0), R6(2), R6(1), R6(3)
    };

    unsigned int c; // c will get v reversed

    // Option 1:
    c = (BitReverseTable256[ v        & 0xff] << 24) |
        (BitReverseTable256[(v >>  8) & 0xff] << 16) |
        (BitReverseTable256[(v >> 16) & 0xff] << 8) |
        (BitReverseTable256[(v >> 24) & 0xff]);

    // Option 2:
    unsigned char * p = (unsigned char *) &v;
    unsigned char * q = (unsigned char *) &c;
    q[3] = BitReverseTable256[p[0]];
    q[2] = BitReverseTable256[p[1]];
    q[1] = BitReverseTable256[p[2]];
    q[0] = BitReverseTable256[p[3]];

    return c;
  }

  unsigned int mod2n(unsigned int n, unsigned int s) {
    const unsigned int d = 1U << s; // So d will be one of: 1, 2, 4, 8, 16, 32, ...

    return n & (d - 1);
  }

  unsigned int mod2n_1(unsigned int n, unsigned int s) {
    const unsigned int d = (1 << s) - 1; // so d is either 1, 3, 7, 15, 31, ...).

    unsigned int m; // n % d goes here.

    for (m = n; n > d; n = m) {
      for (m = 0; n; n >>= s) {
        m += n & d;
      }
    }

   // Now m is a value from 0 to d, but since with modulus division
   // we want m to be 0 when it is d.
    return (m == d ? 0 : m);
  }

  unsigned int log10(unsigned int v) {
    return (v >= 1000000000) ? 9 : (v >= 100000000) ? 8 : (v >= 10000000) ? 7 :
           (v >= 1000000   ) ? 6 : (v >= 100000   ) ? 5 : (v >= 10000   ) ? 4 :
           (v >= 1000      ) ? 3 : (v >= 100      ) ? 2 : (v >= 10      ) ? 1 : 0;
  }

  unisgned int nextPower2(unsigned int v) {
    v--;

    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;

    v++;

    return v;
  }
}

#endif
