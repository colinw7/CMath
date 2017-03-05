#ifndef CTRIANGLE_2D_H
#define CTRIANGLE_2D_H

#include <CPoint2D.h>
#include <CMathGeom2D.h>

class CTriangle2D {
 public:
  CTriangle2D() { }

  CTriangle2D(const CPoint2D &point1, const CPoint2D &point2, const CPoint2D &point3) :
   point1_(point1), point2_(point2), point3_(point3) {
  }

  const CPoint2D &getPoint1() const { return point1_; }
  const CPoint2D &getPoint2() const { return point2_; }
  const CPoint2D &getPoint3() const { return point3_; }

  CPoint2D centroid() { return 0.333333333*(point1_ + point2_ + point3_); }

  double area() { return ::fabs(0.5*area2()); }

  double area2() {
    return CMathGeom2D::TriangleArea2(point1_, point2_, point3_);
  }

  CPolygonOrientation orientation() {
    return CMathGeom2D::PolygonOrientation(point1_.x, point1_.y,
                                           point2_.x, point2_.y,
                                           point3_.x, point3_.y);
  }

  void getBarycentrics(const CPoint2D &point, double *u, double *v) const {
    CVector2D v1 = CVector2D(point  , point3_);
    CVector2D v2 = CVector2D(point3_, point1_);
    CVector2D v3 = CVector2D(point3_, point2_);

    double a, b, c, d, e, f;

    v2.getXY(&a, &d);
    v3.getXY(&b, &e);
    v1.getXY(&c, &f);

    if (a == 0 && b == 0) {
      std::swap(a, d);
      std::swap(b, e);
      std::swap(c, f);
   }

    double a1 = b*f - c*e;
    double a2 = a*e - b*d;
    double b1 = a*f - c*d;
    double b2 = b*d - a*e;

    *u = (a2 != 0 ? a1/a2 : 0.0);
    *v = (b2 != 0 ? b1/b2 : 0.0);
  }

 private:
  CPoint2D point1_;
  CPoint2D point2_;
  CPoint2D point3_;
};

#endif
