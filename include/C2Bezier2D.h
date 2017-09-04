#ifndef C2BEZIER_2D_H
#define C2BEZIER_2D_H

#include <CPoint2D.h>
#include <CBBox2D.h>
#include <vector>

class C2Bezier2D {
 public:
  C2Bezier2D() :
   p1_(), p2_(), p3_() {
  }

  C2Bezier2D(double x1, double y1, double x2, double y2, double x3, double y3) :
   p1_(x1, y1), p2_(x2, y2), p3_(x3, y3) {
  }

  C2Bezier2D(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3) :
   p1_(p1), p2_(p2), p3_(p3) {
  }

  C2Bezier2D(const C2Bezier2D &bezier) :
    p1_(bezier.p1_), p2_(bezier.p2_), p3_(bezier.p3_) {
  }

  C2Bezier2D &operator=(const C2Bezier2D &bezier) {
    p1_ = bezier.p1_;
    p2_ = bezier.p2_;
    p3_ = bezier.p3_;

    return *this;
  }

  const CPoint2D &getFirstPoint  () const { return p1_; }
  const CPoint2D &getControlPoint() const { return p2_; }
  const CPoint2D &getLastPoint   () const { return p3_; }

  void setFirstPoint  (const CPoint2D &p1) { p1_ = p1; }
  void setControlPoint(const CPoint2D &p2) { p2_ = p2; };
  void setLastPoint   (const CPoint2D &p3) { p3_ = p3; };

  void getFirstPoint  (double *x, double *y) const { *x = p1_.x; *y = p1_.y; }
  void getControlPoint(double *x, double *y) const { *x = p2_.x; *y = p2_.y; }
  void getLastPoint   (double *x, double *y) const { *x = p3_.x; *y = p3_.y; }

  void setFirstPoint  (double x, double y) { setFirstPoint  (CPoint2D(x, y)); }
  void setControlPoint(double x, double y) { setControlPoint(CPoint2D(x, y)); }
  void setLastPoint   (double x, double y) { setLastPoint   (CPoint2D(x, y)); }

  void setPoints(double x1, double y1, double x2, double y2, double x3, double y3) {
    setPoints(CPoint2D(x1, y1), CPoint2D(x2, y2), CPoint2D(x3, y3));
  }

  void setPoints(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3) {
    p1_ = p1; p2_ = p2; p3_ = p3;
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

    double tt = t*t;
    double uu = u*u;

    return p1_*uu + 2.0*p2_*t*u + p3_*tt;
  }

  bool interp(double x, double y, double *t) const {
    return interp(CPoint2D(x, y), t);
  }

  bool interp(const CPoint2D &p, double *t) const {
    double t1 = (::fabs(p.x   - p1_.x) + ::fabs(p.y   - p1_.y))/
           (::fabs(p3_.x - p1_.x) + ::fabs(p3_.y - p1_.y));

    CPoint2D pp;

    calc(t1, pp);

    double dx1 = ::fabs(p.x - pp.x);
    double dy1 = ::fabs(p.y - pp.y);

    while (dx1 > 1E-5 || dy1 > 1E-5) {
      if ((pp.x < p.x && pp.x < p3_.x) || (pp.x > p.x && pp.x > p3_.x)) {
        if (pp.x != p3_.x)
          t1 = (1.0 - t1)*(p.x - pp.x)/(p3_.x - pp.x) + t1;
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

      if ((pp.y < p.y && pp.y < p3_.y) || (pp.y > p.y && pp.y > p3_.y)) {
        if (pp.y != p3_.y)
          t1 = (1.0 - t1)*(p.y - pp.y)/(p3_.y - pp.y) + t1;
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
    //return CMathGen::atan2(p3_.x - p2_.x, p3_.y - p2_.y);
    return atan2(p3_.y - p2_.y, p3_.x - p2_.x);
  }

  double gradient(double t) const {
    double u = 1.0 - t;

    CPoint2D p = (p2_ - p1_)*u + (p3_ - p2_)*t;

    //double g = CMathGen::atan2(p.x, p.y);
    double g = atan2(p.y, p.x);

    return g;
  }

  void getHullPolygon(std::vector<CPoint2D> &points) const {
    points.push_back(p1_);
    points.push_back(p2_);
    points.push_back(p3_);
  }

  void getHullBBox(CBBox2D &bbox) const {
    bbox.reset();

    bbox.add(p1_);
    bbox.add(p2_);
    bbox.add(p3_);
  }

  void split(C2Bezier2D &bezier1, C2Bezier2D &bezier2) const {
    // split at control point
    CPoint2D p12 = (p1_ + p2_)/2.0;
    CPoint2D p23 = (p2_ + p3_)/2.0;

    CPoint2D pm = (p12 + p23)/2.0;

    bezier1 = C2Bezier2D(p1_, p12, pm );
    bezier2 = C2Bezier2D(pm , p23, p3_);
  }

  double arcLength(double tol=1E-3) const {
    double l1 = p1_.distanceTo(p3_);
    double l2 = p1_.distanceTo(p2_) + p2_.distanceTo(p3_);

    if (fabs(l2 - l1) < tol)
      return l1;

    C2Bezier2D bezier1, bezier2;

    split(bezier1, bezier2);

    return bezier1.arcLength(tol) + bezier2.arcLength(tol);
  }

  void print(std::ostream &os) const {
    os << "[[" << p1_.x << ", " << p1_.y << "] [" <<
                  p2_.x << ", " << p2_.y << "] [" <<
                  p3_.x << ", " << p3_.y << "]]";
  }

  friend std::ostream &operator<<(std::ostream &os, const C2Bezier2D &bezier) {
    bezier.print(os);

    return os;
  }

 private:
  CPoint2D p1_, p2_, p3_;
};

#endif
