#ifndef CARC_2D_H
#define CARC_2D_H

template<typename T>
class CArc2DT {
 private:
  typedef CPoint2D<T>  Point;
  typedef CVector2D<T> Vector;

  Point center_;
  T     radius_;
  Point start_, end_;

 public:
  CArc2DT();
  CArc2DT(const Point &center, T radius, const Point &start, const Point &end) :
   center_(center), radius_(radius), start_(start), end_(end) {
  }

  bool contains(const Point &point) const {
    // Assert: |P-C| = R where P is the input point, C is the circle center,
    // and R is the circle radius.  For P to be on the arc from A to B, it
    // must be on the side of the plane containing A with normal N = Perp(B-A)
    // where Perp(u,v) = (v,-u).

    Vector pe = point - end_;
    Vector ee = end_ - start_;

    T d = pe.dotPerpendicular(ee);

    return d >= 0.0;
  }
};

typedef CArc2DT<double> CArc2D;

#endif
