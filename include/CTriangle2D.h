#ifndef CTRIANGLE_2D_H
#define CTRIANGLE_2D_H

#include <CPoint2D.h>
#include <CMathGeom2D.h>

template<typename T>
class CTriangle2DT {
 private:
  typedef CPoint2DT<T>  Point;
  typedef CVector2DT<T> Vector;

 private:
  Point point1_;
  Point point2_;
  Point point3_;

 public:
  CTriangle2DT() { }

  CTriangle2DT(const Point &point1, const Point &point2, const Point &point3) :
   point1_(point1), point2_(point2), point3_(point3) {
  }

  const Point &getPoint1() const { return point1_; }
  const Point &getPoint2() const { return point2_; }
  const Point &getPoint3() const { return point3_; }

  Point centroid() { return 0.333333333*(point1_ + point2_ + point3_); }

  double area() { return ::fabs(0.5*area2()); }

  double area2() {
    return CMathGeom2D::TriangleArea2(point1_, point2_, point3_);
  }

  CPolygonOrientation orientation() {
    return CMathGeom2D::PolygonOrientation(point1_.x, point1_.y,
                                           point2_.x, point2_.y,
                                           point3_.x, point3_.y);
  }

  void getBarycentrics(const Point &point, T *u, T *v) const {
    Vector v1 = Vector(point  , point3_);
    Vector v2 = Vector(point3_, point1_);
    Vector v3 = Vector(point3_, point2_);

    T a, b, c, d, e, f, g, h, i;

    v2.getXY(&a, &d);
    v3.getXY(&b, &e);
    v1.getXY(&c, &f);

    if (a == 0 && b == 0) {
      swap(a, d);
      swap(b, e);
      swap(c, f);
   }

    T a1 = b*f - c*e;
    T a2 = a*e - b*d;
    T b1 = a*f - c*d;
    T b2 = b*d - a*e;

    *u = (a2 != 0 ? a1/a2 : 0.0);
    *v = (b2 != 0 ? b1/b2 : 0.0);
  }

};

typedef CTriangle2DT<double> CTriangle2D;

#endif
