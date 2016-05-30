#ifndef CGAUSSIAN_BLUR_H
#define CGAUSSIAN_BLUR_H

#include <vector>

template<typename IMAGE>
class CGaussianBlur {
 public:
  CGaussianBlur() { }

 ~CGaussianBlur() { }

  bool blur(const IMAGE &src, IMAGE &dst, double bx, double by, int nx, int ny) {
    if (std::min(bx, by) <= 0)
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

    assert(m_.size() == dx && m_[0].size() == dy);

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

            r += r1*f;
            g += g1*f;
            b += b1*f;
            a += a1;

            ++n;
          }
        }

        a /= n;

        r = std::min(r, 1.0);
        g = std::min(g, 1.0);
        b = std::min(b, 1.0);
        a = std::min(a, 1.0);

        dst.setRGBA(x2, y2, r, g, b, a);
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

    double bxy  = bx*by;
    double bxy1 = 2*bxy;
    double bxy2 = 1.0/sqrt(M_PI*bxy1);

    double sm = 0.0;

    for (int i = 0, i1 = nx1; i < dx; ++i, ++i1) {
      int i2 = i1*i1;

      for (int j = 0, j1 = ny1; j < dy; ++j, ++j1) {
        int j2 = j1*j1;

        m_[i][j] = bxy2*exp(-(i2 + j2)/bxy1);

        sm += m_[i][j];
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
