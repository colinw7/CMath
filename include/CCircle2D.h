#ifndef CCIRCLE_2D_H
#define CCIRCLE_2D_H

#include <CPoint2D.h>
#include <CLine2D.h>
#include <CMathGeom2D.h>
#include <CShape2D.h>

template<typename T>
class CCircle2DT : public CShape2DT<T> {
 private:
  typedef CPoint2DT<T>  Point;
  typedef CLine2DT<T>   Line;
  typedef CCircle2DT<T> Circle;
  typedef CBBox2DT<T>   BBox;

 private:
  Point center_;
  T     radius_;

 public:
  CCircle2DT() { }

  CCircle2DT(T x, T y, T radius) :
   center_(x, y), radius_(radius) {
  }

  CCircle2DT(const Point &center, T radius) :
   center_(center), radius_(radius) {
  }

  void setCenter(const Point &center) { center_ = center; }
  void setRadius(double radius      ) { radius_ = radius; }

  const Point &getCenter() const { return center_; }
  const T      getRadius() const { return radius_; }

  double area() { return M_PI*radius_*radius_; }

  BBox getBBox() const {
    return BBox(center_.x - radius_, center_.y - radius_,
                center_.x + radius_, center_.y + radius_);
  }

  bool inside(const Point &p) const {
    T d = p.distanceTo(center_);

    return (d < radius_);
  }

  void moveBy(const Point &p) {
    center_ += p;
  }

  void resizeBy(const Point &ll, const Point &ur) {
    radius_ += max(ll.x, max(ll.y, max(ur.x, ur.y)));
  }

  void rotateBy(double, const Point &) { }

  bool lineIntersect(const Line &line, Point &point1, Point &point2) const {
    uint ni;
    T    xi1, yi1, xi2, yi2;

    if (! lineIntersect(line.start().x, line.start().y, line.end().x, line.end().y,
                        &xi1, &yi1, &xi2, &yi2, &ni))
      return false;

    point1 = Point(xi1, yi1);
    point2 = Point(xi2, yi2);

    return true;
  }

  bool lineIntersect(T lx1, T ly1, T lx2, T ly2, T *xi1, T *yi1, T *xi2, T *yi2, uint *ni) const {
    return CMathGeom2D::CircleLineIntersect(center_.x, center_.y, radius_,
                                            lx1,  ly1, lx2, ly2, xi1, yi1, xi2, yi2, ni);
  }
};

typedef CCircle2DT<double> CCircle2D;

#endif
