#ifndef C3BEZIER_2D_H
#define C3BEZIER_2D_H

#include <CPoint2D.h>
#include <C2Bezier2D.h>

class C3Bezier2D {
 public:
  C3Bezier2D() :
   p1_(), p2_(), p3_(), p4_() {
  }

  C3Bezier2D(double x1, double y1, double x2, double y2,
             double x3, double y3, double x4, double y4) :
   p1_(x1, y1), p2_(x2, y2), p3_(x3, y3), p4_(x4, y4) {
  }

  C3Bezier2D(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3, const CPoint2D &p4) :
   p1_(p1), p2_(p2), p3_(p3), p4_(p4) {
  }

  C3Bezier2D(const C3Bezier2D &bezier) :
    p1_(bezier.p1_), p2_(bezier.p2_), p3_(bezier.p3_), p4_(bezier.p4_) {
  }

  // create order 3 bezier from order 2 using 'degree elevation'
  C3Bezier2D(const C2Bezier2D &bezier2) {
    p1_ = bezier2.getFirstPoint();
    p4_ = bezier2.getLastPoint ();

    const CPoint2D &p = bezier2.getControlPoint();

    p2_ = (p1_ + 2*p)/3;
    p3_ = (2*p + p4_)/3;
  }

  C3Bezier2D &operator=(const C3Bezier2D &bezier) {
    p1_ = bezier.p1_;
    p2_ = bezier.p2_;
    p3_ = bezier.p3_;
    p4_ = bezier.p4_;

    return *this;
  }

  const CPoint2D &getFirstPoint   () const { return p1_; }
  const CPoint2D &getControlPoint1() const { return p2_; }
  const CPoint2D &getControlPoint2() const { return p3_; }
  const CPoint2D &getLastPoint    () const { return p4_; }

  void setFirstPoint   (const CPoint2D &p1) { p1_ = p1; }
  void setControlPoint1(const CPoint2D &p2) { p2_ = p2; };
  void setControlPoint2(const CPoint2D &p3) { p3_ = p3; };
  void setLastPoint    (const CPoint2D &p4) { p4_ = p4; };

  void getFirstPoint   (double *x, double *y) const { *x = p1_.x; *y = p1_.y; }
  void getControlPoint1(double *x, double *y) const { *x = p2_.x; *y = p2_.y; }
  void getControlPoint2(double *x, double *y) const { *x = p3_.x; *y = p3_.y; }
  void getLastPoint    (double *x, double *y) const { *x = p4_.x; *y = p4_.y; }

  void setFirstPoint   (double x, double y) { setFirstPoint   (CPoint2D(x, y)); }
  void setControlPoint1(double x, double y) { setControlPoint1(CPoint2D(x, y)); }
  void setControlPoint2(double x, double y) { setControlPoint2(CPoint2D(x, y)); }
  void setLastPoint    (double x, double y) { setLastPoint    (CPoint2D(x, y)); }

  void setPoints(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
    setPoints(CPoint2D(x1, y1), CPoint2D(x2, y2), CPoint2D(x3, y3), CPoint2D(x4, y4));
  }

  void setPoints(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3, const CPoint2D &p4) {
    p1_ = p1; p2_ = p2; p3_ = p3; p4_ = p4;
  }

  void getPoints(CPoint2D &p1, CPoint2D &p2, CPoint2D &p3, CPoint2D &p4) const {
    p1 = p1_; p2 = p2_; p3 = p3_; p4 = p4_;
  }

  void calc(double t, double *x, double *y) const {
    CPoint2D p;

    calc(t, p);

    *x = p.x;
    *y = p.y;
  }

  void calc(double t, CPoint2D &p) const {
    p = calc(t);
  }

  CPoint2D calc(double t) const {
    double u = (1.0 - t);

    double tt  = t*t;
    double ttt = tt*t;

    double uu  = u*u;
    double uuu = uu*u;

    return p1_*uuu + 3.0*p2_*t*uu + 3.0*p3_*tt*u + p4_*ttt;
  }

  bool interp(double x, double y, double *t) const {
    return interp(CPoint2D(x, y), t);
  }

  bool interp(const CPoint2D &p, double *t) const {
    double t1 = (::fabs(p.x   - p1_.x) + ::fabs(p.y   - p1_.y))/
           (::fabs(p4_.x - p1_.x) + ::fabs(p4_.y - p1_.y));

    CPoint2D pp;

    calc(t1, pp);

    double dx1 = ::fabs(p.x - pp.x);
    double dy1 = ::fabs(p.y - pp.y);

    while (dx1 > 1E-5 || dy1 > 1E-5) {
      if ((pp.x < p.x && pp.x < p4_.x) || (pp.x > p.x && pp.x > p4_.x)) {
        if (pp.x != p4_.x)
          t1 = (1.0 - t1)*(p.x - pp.x)/(p4_.x - pp.x ) + t1;
        else
          t1 = 1.0;
      }
      else {
        if (pp.x != p1_.x)
          t1 = t1*(p.x - p1_.x)/(pp.x - p1_.x);
        else
          t1 = 0.0;
      }

      calc(t1, pp);

      double dx2 = ::fabs(p.x - pp.x);
      double dy2 = ::fabs(p.y - pp.y);

      if (dx2 < dx1 && dy2 < dy1)
        goto next;

      if ((pp.y < p.y && pp.y < p4_.y) || (pp.y > p.y && pp.y > p4_.y)) {
        if (pp.y != p4_.y)
          t1 = (1.0 - t1)*(p.y - pp.y )/(p4_.y - pp.y) + t1;
        else
          t1 = 1.0;
      }
      else {
        if (pp.y != p1_.y)
          t1 = t1*(p.y - p1_.y)/(pp.y - p1_.y);
        else
          t1 = 0.0;
      }

      calc(t1, pp);

      dx2 = ::fabs(p.x - pp.x);
      dy2 = ::fabs(p.y - pp.y);

      if (dx2 >= dx1 || dy2 >= dy1)
        return false;

   next:
      dx1 = dx2;
      dy1 = dy2;
    }

    *t = t1;

    return true;
  }

  double gradientStart() const {
    //return CMathGen::atan2(p2_.x - p1_.x, p2_.y - p1_.y);
    return atan2(p2_.y - p1_.y, p2_.x - p1_.x);
  }

  double gradientEnd() const {
    //return CMathGen::atan2(p4_.x - p3_.x, p4_.y - p3_.y);
    return atan2(p4_.y - p3_.y, p4_.x - p3_.x);
  }

  double gradient(double t) const {
    double u = 1.0 - t;

    double tt = t*t;
    double uu = u*u;
    double tu = t*u;

    CPoint2D p = (p2_ - p1_)*uu + 2.0*(p3_ - p2_)*tu + (p4_ - p3_)*tt;

    //double g = CMathGen::atan2(p.x, p.y);
    double g = atan2(p.y, p.x);

    return g;
  }

  void getHullPolygon(std::vector<CPoint2D> &points) const {
    points.push_back(p1_);
    points.push_back(p2_);
    points.push_back(p3_);
    points.push_back(p4_);
  }

  void getHullBBox(CBBox2D &bbox) const {
    bbox.reset();

    bbox.add(p1_);
    bbox.add(p2_);
    bbox.add(p3_);
    bbox.add(p4_);
  }

  void split(C3Bezier2D &bezier1, C3Bezier2D &bezier2) const {
    split(0.5, bezier1, bezier2);
  }

  void split(double t, C3Bezier2D &bezier1, C3Bezier2D &bezier2) const {
    // split at t (0 - 1) of curve
    double u = 1.0 - t;

    CPoint2D p11 = u*p1_ + t*p2_;
    CPoint2D p12 = u*p2_ + t*p3_;
    CPoint2D p13 = u*p3_ + t*p4_;

    CPoint2D p21 = u*p11 + t*p12;
    CPoint2D p22 = u*p12 + t*p13;

    CPoint2D p31 = u*p21 + t*p22;

    bezier1 = C3Bezier2D(p1_, p11, p21, p31);
    bezier2 = C3Bezier2D(p31, p22, p13, p4_);
  }

  bool split(const CPoint2D &p, C3Bezier2D &bezier1, C3Bezier2D &bezier2) const {
    double t;

    if (! interp(p, &t)) return false;

    split(t, bezier1, bezier2);

    return true;
  }

  CPoint2D deCasteljauInterp(double t) const {
    double u = 1.0 - t;

    CPoint2D p11 = u*p1_ + t*p2_;
    CPoint2D p12 = u*p2_ + t*p3_;
    CPoint2D p13 = u*p3_ + t*p4_;

    CPoint2D p21 = u*p11 + t*p12;
    CPoint2D p22 = u*p12 + t*p13;

    CPoint2D p31 = u*p21 + t*p22;

    return p31;
  }

  template<typename FUNC>
  static C3Bezier2D bestFit(FUNC f, double x1, double y1, double g1,
                            double x2, double y2, double g2, int steps=50) {
    double c1 = y1 - g1*x1;
    double c2 = y2 - g2*x2;

    double y1e = g1*x2 + c1;
    double y2e = g2*x1 + c2;

    double s = 1.0/steps;

    double x11 = x1;
    double y11 = y1;

    double x21 = x2;
    double y21 = y2;

    int minI1 = -1;
    int minI2 = -1;

    bool changed1 = true;
    bool changed2 = true;

    int max_iter = 100;

    while (max_iter && (changed1 || changed2)) {
      int    oldI1 = minI1;
      double minD1 = 1E50;

      for (int i = 0; i <= 2*steps; ++i) {
        double x11 = (x2  - x1)*i*s + x1;
        double y11 = (y1e - y1)*i*s + y1;

        C3Bezier2D b(x1, y1, x11, y11, x21, y21, x2, y2);

        double bx, by;

        b.calc(0.33, &bx, &by);

        double d = fabs(by - f(bx));

        if (d < minD1) {
          minD1 = d;
          minI1 = i;
        }
        else
          break;
      }

      x11 = (x2  - x1)*minI1*s + x1;
      y11 = (y1e - y1)*minI1*s + y1;

      changed1 = (minI1 != oldI1);

      //----

      int    oldI2 = minI2;
      double minD2 = 1E50;

      for (int i = 0; i <= 2*steps; ++i) {
        double x21 = (x1  - x2)*i*s + x2;
        double y21 = (y2e - y2)*i*s + y2;

        C3Bezier2D b(x1, y1, x11, y11, x21, y21, x2, y2);

        double bx, by;

        b.calc(0.66, &bx, &by);

        double d = fabs(by - f(bx));

        if (d < minD2) {
          minD2 = d;
          minI2 = i;
        }
        else
          break;
      }

      x21 = (x1  - x2)*minI2*s + x2;
      y21 = (y2e - y2)*minI2*s + y2;

      changed2 = (minI2 != oldI2);

      --max_iter;
    }

    return C3Bezier2D(x1, y1, x11, y11, x21, y21, x2, y2);
  }

  template<typename FUNC>
  static C3Bezier2D bestParamFit(FUNC f, double t1=0, double t2=1, int steps=50) {
    assert(t2 > t1);

    CPoint2D p1 = f(t1);
    CPoint2D p2 = f(t1 + 0.01);
    CPoint2D p3 = f(t2 - 0.01);
    CPoint2D p4 = f(t2);

    double x1 = p1.x, y1 = p1.y;
    double x2 = p4.x, y2 = p4.y;

    double g1 = (p2.y - p1.y)/(p2.x - p1.x);
    double g2 = (p4.y - p3.y)/(p4.x - p3.x);

    double x11 = x1, y11 = y1;
    double x21 = x2, y21 = y2;

    double dx = (x2 - x1)/(steps - 1);

    int minI1 = -1;
    int minI2 = -1;

    bool changed1 = true;
    bool changed2 = true;

    int max_iter = 100;

    while (max_iter && (changed1 || changed2)) {
      int    oldI1 = minI1;
      double minD1 = 1E50;

      for (int i = 1; i < 2*steps; ++i) {
        double x11 = i*dx + x1;
        double y11 = g1*(x11 - x1) + y1;

        C3Bezier2D b(x1, y1, x11, y11, x21, y21, x2, y2);

        double bx, by;

        double it = 0.33*(t2 - t1) + t1;

        b.calc(it, &bx, &by);

        CPoint2D pi = f(it);

        double d = std::max(fabs(bx - pi.x), fabs(by - pi.y));

        if (d < minD1) {
          minD1 = d;
          minI1 = i;
        }
      }

      x11 = minI1*dx + x1;
      y11 = g1*(x11 - x1) + y1;

      changed1 = (minI1 != oldI1);

      //----

      int    oldI2 = minI2;
      double minD2 = 1E50;

      for (int i = steps - 1; i > -steps; --i) {
        double x21 = i*dx + x1;
        double y21 = g2*(x21 - x2) + y2;

        C3Bezier2D b(x1, y1, x11, y11, x21, y21, x2, y2);

        double bx, by;

        double it = 0.66*(t2 - t1) + t1;

        b.calc(it, &bx, &by);

        CPoint2D pi = f(it);

        double d = std::max(fabs(bx - pi.x), fabs(by - pi.y));

        if (d < minD2) {
          minD2 = d;
          minI2 = i;
        }
      }

      x21 = minI2*dx + x1;
      y21 = g2*(x21 - x2) + y2;

      changed2 = (minI2 != oldI2);

      --max_iter;
    }

    return C3Bezier2D(x1, y1, x11, y11, x21, y21, x2, y2);
  }

  template<typename FUNC>
  static C3Bezier2D bestPolarFit(FUNC f, double a1, double r1, double g1,
                                 double a2, double r2, double g2, int steps=50) {
    double x1 = r1*cos(a1);
    double y1 = r1*sin(a1);

    double x2 = r2*cos(a2);
    double y2 = r2*sin(a2);

    double c1 = y1 - g1*x1;
    double c2 = y2 - g2*x2;

    double dx = (x2 - x1)/steps;
    double dy = (y2 - y1)/steps;

    double x11, y11, x21, y21;

    bool use_dx = (fabs(g1) < 4 && fabs(g2) < 4);
    //bool use_dx = false;
    //bool use_dx = true;

    if (use_dx) {
      x11 = x1 + dx; y11 = g1*x11 + c1;
      x21 = x2 - dx; y21 = g2*x21 + c2;
    }
    else {
      y11 = y1 + dy; x11 = (y11 - c1)/g1;
      y21 = y2 - dy; x21 = (y21 - c2)/g2;
    }

    int minI1 = -1;
    int minI2 = -1;

    bool changed1 = true;
    bool changed2 = true;

    while (changed1 || changed2) {
      int    oldI1 = minI1;
      double minD1 = 1E50;

      for (int i = 1; i < steps; ++i) {
        double x11, y11;

        if (use_dx) {
          x11 = x1 + i*dx;
          y11 = g1*x11 + c1;
        }
        else {
          y11 = y1 + i*dy;
          x11 = (y11 - c1)/g1;
        }

        C3Bezier2D b(x1, y1, x11, y11, x21, y21, x2, y2);

        double bx, by;

        b.calc(0.33, &bx, &by);

        double ba = a1 + 0.33*(a2 - a1);
        double br = f(ba);

        double d;

        if  (use_dx)
          d = fabs(by - br*sin(ba));
        else
          d = fabs(bx - br*cos(ba));

        if (d < minD1) {
          minD1 = d;
          minI1 = i;
        }
      }

      if  (use_dx) {
        x11 = x1 + minI1*dx;
        y11 = g1*x11 + c1;
      }
      else {
        y11 = y1 + minI1*dy;
        x11 = (y11 - c1)/g1;
      }

      changed1 = (minI1 != oldI1);

      //----

      int    oldI2 = minI2;
      double minD2 = 1E50;

      for (int i = 1; i < steps; ++i) {
        double x21, y21;

        if (use_dx) {
          x21 = x2 - i*dx;
          y21 = g2*x21 + c2;
        }
        else {
          y21 = y2 - i*dy;
          x21 = (y21 - c2)/g2;
        }

        C3Bezier2D b(x1, y1, x11, y11, x21, y21, x2, y2);

        double bx, by;

        b.calc(0.66, &bx, &by);

        double ba = a1 + 0.66*(a2 - a1);
        double br = f(ba);

        double d;

        if  (use_dx)
          d = fabs(by - br*sin(ba));
        else
          d = fabs(bx - br*cos(ba));

        if (d < minD2) {
          minD2 = d;
          minI2 = i;
        }
      }

      if  (use_dx) {
        x21 = x2 - minI2*dx;
        y21 = g2*x21 + c2;
      }
      else {
        y21 = y2 - minI2*dy;
        x21 = (y21 - c2)/g2;
      }

      changed2 = (minI2 != oldI2);
    }

    return C3Bezier2D(x1, y1, x11, y11, x21, y21, x2, y2);
  }

  double arcLength(double tol=1E-3) const {
    double l1 = p1_.distanceTo(p4_);
    double l2 = p1_.distanceTo(p2_) + p2_.distanceTo(p3_) + p3_.distanceTo(p4_);

    if (fabs(l2 - l1) < tol)
      return l1;

    C3Bezier2D bezier1, bezier2;

    split(bezier1, bezier2);

    return bezier1.arcLength(tol) + bezier2.arcLength(tol);
  }

  void print(std::ostream &os) const {
    os << "[[" << p1_.x << ", " << p1_.y << "] [" <<
                  p2_.x << ", " << p2_.y << "] [" <<
                  p3_.x << ", " << p3_.y << "] [" <<
                  p4_.x << ", " << p4_.y << "]]";
  }

  friend std::ostream &operator<<(std::ostream &os, const C3Bezier2D &bezier) {
    bezier.print(os);

    return os;
  }

 private:
  CPoint2D p1_, p2_, p3_, p4_;
};

#endif
