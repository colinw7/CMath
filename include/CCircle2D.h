#ifndef CCIRCLE_2D_H
#define CCIRCLE_2D_H

#include <CPoint2D.h>
#include <CLine2D.h>
#include <CMathGeom2D.h>
#include <CShape2D.h>

class CCircle2D : public CShape2D {
 private:
  typedef CPoint2D  Point;
  typedef CLine2D   Line;
  typedef CCircle2D Circle;
  typedef CBBox2D   BBox;

 public:
  CCircle2D() { }

  CCircle2D(double x, double y, double radius) :
   center_(x, y), radius_(radius) {
  }

  CCircle2D(const Point &center, double radius) :
   center_(center), radius_(radius) {
  }

  void setCenter(const Point &center) { center_ = center; }
  void setRadius(double radius      ) { radius_ = radius; }

  const Point& getCenter() const { return center_; }
  double       getRadius() const { return radius_; }

  double area() { return M_PI*radius_*radius_; }

  BBox getBBox() const {
    return BBox(center_.x - radius_, center_.y - radius_,
                center_.x + radius_, center_.y + radius_);
  }

  bool inside(const Point &p) const {
    double d = p.distanceTo(center_);

    return (d < radius_);
  }

  void moveBy(const Point &p) {
    center_ += p;
  }

  void resizeBy(const Point &ll, const Point &ur) {
    radius_ += std::max(ll.x, std::max(ll.y, std::max(ur.x, ur.y)));
  }

  void rotateBy(double, const Point &) { }

  bool lineIntersect(const Line &line, Point &point1, Point &point2) const {
    uint   ni;
    double xi1, yi1, xi2, yi2;

    if (! lineIntersect(line.start().x, line.start().y, line.end().x, line.end().y,
                        &xi1, &yi1, &xi2, &yi2, &ni))
      return false;

    point1 = Point(xi1, yi1);
    point2 = Point(xi2, yi2);

    return true;
  }

  bool lineIntersect(double lx1, double ly1, double lx2, double ly2,
                     double *xi1, double *yi1, double *xi2, double *yi2, uint *ni) const {
    return CMathGeom2D::CircleLineIntersect(center_.x, center_.y, radius_,
                                            lx1,  ly1, lx2, ly2, xi1, yi1, xi2, yi2, ni);
  }

 private:
  Point  center_;
  double radius_ { 1 };
};

#endif
