#ifndef C2BEZIER_3D_H
#define C2BEZIER_3D_H

#include <CPoint3D.h>
#include <CBBox3D.h>
#include <vector>

class C2Bezier3D {
 public:
  C2Bezier3D() { }

  C2Bezier3D(const CPoint3D &p1, const CPoint3D &p2, const CPoint3D &p3) :
   p1_(p1), p2_(p2), p3_(p3) {
  }

  //---

  const CPoint3D &getFirstPoint  () const { return p1_; }
  const CPoint3D &getControlPoint() const { return p2_; }
  const CPoint3D &getLastPoint   () const { return p3_; }

  void setFirstPoint  (const CPoint3D &p1) { p1_ = p1; lengthValid_ = false; }
  void setControlPoint(const CPoint3D &p2) { p2_ = p2; lengthValid_ = false; }
  void setLastPoint   (const CPoint3D &p3) { p3_ = p3; lengthValid_ = false; }

  //---

  void setPoints(const CPoint3D &p1, const CPoint3D &p2, const CPoint3D &p3) {
    p1_ = p1; p2_ = p2; p3_ = p3;

    lengthValid_ = false;
  }

  //---

  void calc(double t, CPoint3D &p) const {
    p = calc(t);
  }

  CPoint3D calc(double t) const {
    double u = (1.0 - t);

    double tt = t*t;
    double uu = u*u;

    return p1_*uu + 2.0*p2_*t*u + p3_*tt;
  }

  //---

#if 0
  bool interp(const CPoint3D &p, double *t) const {
    double t1 = (::fabs(p.x   - p1_.x) + ::fabs(p.y   - p1_.y))/
                (::fabs(p3_.x - p1_.x) + ::fabs(p3_.y - p1_.y));

    CPoint3D pp;

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
#endif

  //---

#if 0
  double gradientStart() const {
    return atan2(p2_.y - p1_.y, p2_.x - p1_.x);
  }

  double gradientEnd() const {
    return atan2(p3_.y - p2_.y, p3_.x - p2_.x);
  }

  double gradient(double t) const {
    double u = 1.0 - t;

    CPoint3D p = (p2_ - p1_)*u + (p3_ - p2_)*t;

    double g = atan2(p.y, p.x);

    return g;
  }
#endif

  //---

  void getHullPolygon(std::vector<CPoint3D> &points) const {
    points.push_back(p1_);
    points.push_back(p2_);
    points.push_back(p3_);
  }

  void getHullBBox(CBBox3D &bbox) const {
    bbox.reset();

    bbox.add(p1_);
    bbox.add(p2_);
    bbox.add(p3_);
  }

  //---

  void split(C2Bezier3D &bezier1, C2Bezier3D &bezier2) const {
    // split at control point
    CPoint3D p12 = (p1_ + p2_)/2.0;
    CPoint3D p23 = (p2_ + p3_)/2.0;

    CPoint3D pm = (p12 + p23)/2.0;

    bezier1 = C2Bezier3D(p1_, p12, pm );
    bezier2 = C2Bezier3D(pm , p23, p3_);
  }

  //---

  double arcLength(double tol=1E-3) const {
    if (! lengthValid_) {
      double l1 = p1_.distanceTo(p3_);
      double l2 = p1_.distanceTo(p2_) + p2_.distanceTo(p3_);

      if (fabs(l2 - l1) < tol)
        return l1;

      C2Bezier3D bezier1, bezier2;
      split(bezier1, bezier2);

      length_ = bezier1.arcLength(tol) + bezier2.arcLength(tol);

      lengthValid_ = true;
    }

    return length_;
  }

  //---

  void print(std::ostream &os) const {
    os << "[[" << p1_.x << ", " << p1_.y << "] [" <<
                  p2_.x << ", " << p2_.y << "] [" <<
                  p3_.x << ", " << p3_.y << "]]";
  }

  friend std::ostream &operator<<(std::ostream &os, const C2Bezier3D &bezier) {
    bezier.print(os);

    return os;
  }

 private:
  CPoint3D p1_, p2_, p3_;

  mutable bool   lengthValid_ { false };
  mutable double length_      { 0.0 };
};

#endif
