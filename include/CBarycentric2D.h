#ifndef CBARYCENTRIC_2D_H
#define CBARYCENTRIC_2D_H

#include <CTriplet.h>
#include <CPoint2D.h>

// Use barycentric equation of triangle  p(b1, b2)
//
//  b0 = 1 - b1 - b2
//
//  p(b1, b2) = b0*p0 + b1*p1 + b2*p2
//
// where b1 >= 0 and b2 >= 0 and b1 + b2 <= 1
//
class CBarycentric2D {
 private:
  typedef CPoint2D  Point;
  typedef CVector2D Vector;
  typedef CTriplet  Triplet;

 public:
  CBarycentric2D(const Point &p1, const Point &p2, const Point &p3) :
   p1_(p1), p2_(p2), p3_(p3) {
  }

  Triplet fromPoint(const Point &p) {
    double area1 = triangleArea(p1_, p2_, p);
    double area2 = triangleArea(p1_, p3_, p);
    double area3 = triangleArea(p2_, p3_, p);

    double is = 1.0/(area1 + area2 + area3);

    return Triplet(area1*is, area2*is, area3*is);
  }

  Point toPoint(const Triplet &t) {
    return (p1_*t.first + p2_*t.second + p3_*t.third);
  }

 private:
  double triangleArea(const Point &p1, const Point &p2, const Point &p3) {
    Vector p21 = p2 - p1;
    Vector p31 = p3 - p1;

    return 0.5*(p21.crossProduct(p31)).length();
  }

 private:
  Point p1_, p2_, p3_;
};

#endif
