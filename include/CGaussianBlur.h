#ifndef CGAUSSIAN_BLUR_H
#define CGAUSSIAN_BLUR_H

#include <CRGBA.h>
#include <vector>

template<typename IMAGE>
class CGaussianBlur {
 public:
  CGaussianBlur() { }

 ~CGaussianBlur() { }

  bool blur(const IMAGE &src, IMAGE &dst, double bx, double by, int nx, int ny) {
    if (bx <= 0 && by <= 0)
      return false;

    // init matrix (if needed)
    if (nx == 0) {
      nx = int(6*bx + 1);

      if (nx > 4) nx = 4;
    }

    if (ny == 0) {
      ny = int(6*by + 1);

      if (ny > 4) ny = 4;
    }

    setMatrix(bx, by, nx, ny);

    //------

    int nx1 = -nx_/2;
    int nx2 =  nx_/2;
    int ny1 = -ny_/2;
    int ny2 =  ny_/2;

    int dx = nx2 - nx1 + 1;
    int dy = ny2 - ny1 + 1;

    assert(int(m_.size()) == dx && int(m_[0].size()) == dy);

    //---

    int px1, py1, px2, py2;

    src.getPixelRange(&px1, &py1, &px2, &py2);

    int wx1, wy1, wx2, wy2;

    src.getWindow(&wx1, &wy1, &wx2, &wy2);

    //---

    for (int y1 = wy1 + ny1, y2 = wy1, y3 = wy1 + ny2; y2 <= wy2; ++y1, ++y2, ++y3) {
      if (y2 < py1 || y2 > py2) continue;

      for (int x1 = wx1 + nx1, x2 = wx1, x3 = wx1 + nx2; x2 <= wx2; ++x1, ++x2, ++x3) {
        if (x2 < px1 || x2 > px2) continue;

        double r = 0.0;
        double g = 0.0;
        double b = 0.0;
        double a = 0.0;
        int    n = 0;

        for (int i = 0, x = x1; i < dx; ++i, ++x) {
          if (x < px1 || x > px2) continue;

          for (int j = 0, y = y1; j < dy; ++j, ++y) {
            if (y < py1 || y > py2) continue;

            double r1, g1, b1, a1;

            src.getRGBA(x, y, &r1, &g1, &b1, &a1);

            double f = m_[i][j]/sm_;

            r += r1*a1*f;
            g += g1*a1*f;
            b += b1*a1*f;

            a += a1;

            ++n;
          }
        }

        a /= n;

        if (a < 1E-3)
          continue;

        CRGBA rgba1(r/a, g/a, b/a, a);
        CRGBA rgba2(0, 0, 0, 0);

        rgba2.combine(rgba1);

        dst.setRGBA(x2, y2, rgba2.getRed(), rgba2.getGreen(), rgba2.getBlue(), rgba2.getAlpha());
      }
    }

    return true;
  }

  void setMatrix(double bx, double by, int nx, int ny) {
    if (realEq(bx, bx_) && realEq(by_, by) && nx == nx_ && ny == ny_)
      return;

    int nx1 = -nx/2;
    int nx2 =  nx/2;
    int ny1 = -ny/2;
    int ny2 =  ny/2;

    int dx = nx2 - nx1 + 1;
    int dy = ny2 - ny1 + 1;

    m_.resize(dx);

    for (int i = 0; i < dx; ++i)
      m_[i].resize(dy);

    double sm = 0.0;

    // x,y blur
    if      (bx > 0 && by > 0) {
      double bxy   = bx*by;
      double bxy2  = 2*bxy;
      double ibxy2 = 1.0/(M_PI*bxy2);

      for (int i = 0, i1 = nx1; i < dx; ++i, ++i1) {
        int x2 = i1*i1;

        for (int j = 0, j1 = ny1; j < dy; ++j, ++j1) {
          int y2 = j1*j1;

          m_[i][j] = ibxy2*exp(-(x2 + y2)/bxy2);

          sm += m_[i][j];
        }
      }
    }
    // x blur
    else if (bx > 0) {
      double bx2  = bx*bx;
      double bx22 = 2*bx*bx;
      double ibx2 = 1.0/sqrt(M_PI*bx2);

      for (int i = 0, i1 = nx1; i < dx; ++i, ++i1) {
        int x2 = i1*i1;

        for (int j = 0, j1 = ny1; j < dy; ++j, ++j1) {
          m_[i][j] = ibx2*exp(-x2/bx22);

          sm += m_[i][j];
        }
      }
    }
    // y blur
    else {
      double by2  = by*by;
      double by22 = 2*by*by;
      double iby2 = 1.0/sqrt(M_PI*by2);

      for (int i = 0, i1 = nx1; i < dx; ++i, ++i1) {
        for (int j = 0, j1 = ny1; j < dy; ++j, ++j1) {
          int y2 = j1*j1;

          m_[i][j] = iby2*exp(-y2/by22);

          sm += m_[i][j];
        }
      }
    }

    //---

    bx_ = bx;
    by_ = by;
    nx_ = nx;
    ny_ = ny;
    sm_ = sm;
  }

 private:
  static bool realEq(double r1, double r2) {
    return (std::fabs(r1 - r2) < 1E-6);
  }

 private:
  typedef std::vector<double> Reals;
  typedef std::vector<Reals>  RealsArray;

  RealsArray m_;
  double     bx_ { 0 }, by_ { 0 };
  int        nx_ { 0 }, ny_ { 0 };
  double     sm_ { 0 };
};

#endif
