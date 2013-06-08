#ifndef C3BEZIER_2D_H
#define C3BEZIER_2D_H

#include <CPoint2D.h>
#include <CMathGen.h>
#include <C2Bezier2D.h>

template<typename T>
class C3Bezier2DT {
 private:
  typedef CPoint2DT<T>   Point;
  typedef C3Bezier2DT<T> Bezier;

 public:
  C3Bezier2DT() :
   p1_(), p2_(), p3_(), p4_() {
  }

  C3Bezier2DT(T x1, T y1, T x2, T y2, T x3, T y3, T x4, T y4) :
   p1_(x1, y1), p2_(x2, y2), p3_(x3, y3), p4_(x4, y4) {
  }

  C3Bezier2DT(const Point &p1, const Point &p2, const Point &p3, const Point &p4) :
   p1_(p1), p2_(p2), p3_(p3), p4_(p4) {
  }

  C3Bezier2DT(const C3Bezier2DT<T> &bezier) :
    p1_(bezier.p1_), p2_(bezier.p2_), p3_(bezier.p3_), p4_(bezier.p4_) {
  }

  // create order 3 bezier from order 2 using 'degree elevation'
  C3Bezier2DT(const C2Bezier2DT<T> &bezier2) {
    p1_ = bezier2.getFirstPoint();
    p4_ = bezier2.getLastPoint ();

    const Point &p = bezier2.getControlPoint();

    p2_ = (p1_ + 2*p)/3;
    p3_ = (2*p + p4_)/3;
  }

  C3Bezier2DT &operator=(const C3Bezier2DT<T> &bezier) {
    p1_ = bezier.p1_;
    p2_ = bezier.p2_;
    p3_ = bezier.p3_;
    p4_ = bezier.p4_;

    return *this;
  }

  const Point &getFirstPoint   () const { return p1_; }
  const Point &getControlPoint1() const { return p2_; }
  const Point &getControlPoint2() const { return p3_; }
  const Point &getLastPoint    () const { return p4_; }

  void setFirstPoint   (const Point &p1) { p1_ = p1; }
  void setControlPoint1(const Point &p2) { p2_ = p2; };
  void setControlPoint2(const Point &p3) { p3_ = p3; };
  void setLastPoint    (const Point &p4) { p4_ = p4; };

  void getFirstPoint   (T *x, T *y) const { *x = p1_.x; *y = p1_.y; }
  void getControlPoint1(T *x, T *y) const { *x = p2_.x; *y = p2_.y; }
  void getControlPoint2(T *x, T *y) const { *x = p3_.x; *y = p3_.y; }
  void getLastPoint    (T *x, T *y) const { *x = p4_.x; *y = p4_.y; }

  void setFirstPoint   (T x, T y) { setFirstPoint   (Point(x, y)); }
  void setControlPoint1(T x, T y) { setControlPoint1(Point(x, y)); }
  void setControlPoint2(T x, T y) { setControlPoint2(Point(x, y)); }
  void setLastPoint    (T x, T y) { setLastPoint    (Point(x, y)); }

  void setPoints(T x1, T y1, T x2, T y2, T x3, T y3, T x4, T y4) {
    setPoints(Point(x1, y1), Point(x2, y2), Point(x3, y3), Point(x4, y4));
  }

  void setPoints(const Point &p1, const Point &p2, const Point &p3, const Point &p4) {
    p1_ = p1; p2_ = p2; p3_ = p3; p4_ = p4;
  }

  void getPoints(Point &p1, Point &p2, Point &p3, Point &p4) const {
    p1 = p1_; p2 = p2_; p3 = p3_; p4 = p4_;
  }

  void calc(T t, T *x, T *y) const {
    Point p;

    calc(t, p);

    *x = p.x;
    *y = p.y;
  }

  void calc(T t, Point &p) const {
    p = calc(t);
  }

  Point calc(T t) const {
    T u = (1.0 - t);

    T tt  = t*t;
    T ttt = tt*t;

    T uu  = u*u;
    T uuu = uu*u;

    return p1_*uuu + 3.0*p2_*t*uu + 3.0*p3_*tt*u + p4_*ttt;
  }

  bool interp(T x, T y, T *t) const {
    return interp(Point(x, y), t);
  }

  bool interp(const Point &p, T *t) const {
    T t1 = (::fabs(p.x   - p1_.x) + ::fabs(p.y   - p1_.y))/
           (::fabs(p4_.x - p1_.x) + ::fabs(p4_.y - p1_.y));

    Point pp;

    calc(t1, pp);

    T dx1 = ::fabs(p.x - pp.x);
    T dy1 = ::fabs(p.y - pp.y);

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

      T dx2 = ::fabs(p.x - pp.x);
      T dy2 = ::fabs(p.y - pp.y);

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

  T gradientStart() const {
    return atan2(p2_.y - p1_.y, p2_.x - p1_.x);
  }

  T gradientEnd() const {
    return atan2(p4_.y - p3_.y, p4_.x - p3_.x);
  }

  T gradient(T t) const {
    T u = 1.0 - t;

    T tt = t*t;
    T uu = u*u;
    T tu = t*u;

    Point p = (p2_ - p1_)*uu + 2.0*(p3_ - p2_)*tu + (p4_ - p3_)*tt;

    T g = atan2(p.y, p.x);

    return g;
  }

  void getHullPolygon(std::vector<Point> &points) const {
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

  void split(Bezier &bezier1, Bezier &bezier2) const {
    split(0.5, bezier1, bezier2);
  }

  void split(double t, Bezier &bezier1, Bezier &bezier2) const {
    // split at t (0 - 1) of curve
    T u = 1.0 - t;

    Point p11 = u*p1_ + t*p2_;
    Point p12 = u*p2_ + t*p3_;
    Point p13 = u*p3_ + t*p4_;

    Point p21 = u*p11 + t*p12;
    Point p22 = u*p12 + t*p13;

    Point p31 = u*p21 + t*p22;

    bezier1 = Bezier(p1_, p11, p21, p31);
    bezier2 = Bezier(p31, p22, p13, p4_);
  }

  bool split(const CPoint2D &p, Bezier &bezier1, Bezier &bezier2) const {
    double t;

    if (! interp(p, &t)) return false;

    split(t, bezier1, bezier2);

    return true;
  }

  Point deCasteljauInterp(double t) const {
    T u = 1.0 - t;

    Point p11 = u*p1_ + t*p2_;
    Point p12 = u*p2_ + t*p3_;
    Point p13 = u*p3_ + t*p4_;

    Point p21 = u*p11 + t*p12;
    Point p22 = u*p12 + t*p13;

    Point p31 = u*p21 + t*p22;

    return p31;
  }

  template<typename FUNC>
  static Bezier bestFit(FUNC f, double x1, double y1, double g1,
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

        Bezier b(x1, y1, x11, y11, x21, y21, x2, y2);

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

        Bezier b(x1, y1, x11, y11, x21, y21, x2, y2);

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

    return Bezier(x1, y1, x11, y11, x21, y21, x2, y2);
  }

  template<typename FUNC>
  static Bezier bestParamFit(FUNC f, double t1=0, double t2=1, int steps=50) {
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

        Bezier b(x1, y1, x11, y11, x21, y21, x2, y2);

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

        Bezier b(x1, y1, x11, y11, x21, y21, x2, y2);

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

    return Bezier(x1, y1, x11, y11, x21, y21, x2, y2);
  }

  template<typename FUNC>
  static Bezier bestPolarFit(FUNC f, double a1, double r1, double g1,
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

      for (uint i = 1; i < steps; ++i) {
        double x11, y11;

        if (use_dx) {
          x11 = x1 + i*dx;
          y11 = g1*x11 + c1;
        }
        else {
          y11 = y1 + i*dy;
          x11 = (y11 - c1)/g1;
        }

        Bezier b(x1, y1, x11, y11, x21, y21, x2, y2);

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

      for (uint i = 1; i < steps; ++i) {
        double x21, y21;

        if (use_dx) {
          x21 = x2 - i*dx;
          y21 = g2*x21 + c2;
        }
        else {
          y21 = y2 - i*dy;
          x21 = (y21 - c2)/g2;
        }

        Bezier b(x1, y1, x11, y11, x21, y21, x2, y2);

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

    return Bezier(x1, y1, x11, y11, x21, y21, x2, y2);
  }

  void print(std::ostream &os) const {
    os << "[[" << p1_.x << ", " << p1_.y << "] [" <<
                  p2_.x << ", " << p2_.y << "] [" <<
                  p3_.x << ", " << p3_.y << "] [" <<
                  p4_.x << ", " << p4_.y << "]]";
  }

  friend std::ostream &operator<<(std::ostream &os, const Bezier &bezier) {
    bezier.print(os);

    return os;
  }

 private:
  Point p1_, p2_, p3_, p4_;
};

typedef C3Bezier2DT<double> C3Bezier2D;

#endif
