#ifndef CELLIPSE_2D_H
#define CELLIPSE_2D_H

#include <CPoint2D.h>
#include <CLine2D.h>
#include <CMathGeom2D.h>
#include <CShape2D.h>

class CEllipse2D : public CShape2D {
 private:
  typedef CPoint2D   Point;
  typedef CLine2D    Line;
  typedef CEllipse2D Ellipse;
  typedef CBBox2D    BBox;

 public:
  CEllipse2D() { }

  CEllipse2D(double x, double y, double xr, double yr) :
   center_(x, y), xr_(xr), yr_(yr) {
  }

  CEllipse2D(const Point &center, double xr, double yr) :
   center_(center), xr_(xr), yr_(yr) {
  }

  void setCenter(const Point &center) { center_ = center; }

  void setRadii(double xr, double yr) { xr_ = xr; yr_ = yr; }

  const Point &getCenter() const { return center_; }

  const double getXRadius() const { return xr_; }
  const double getYRadius() const { return yr_; }

  double area() { return M_PI*xr_*yr_; }

  BBox getBBox() const {
    return BBox(center_.x - xr_, center_.y - yr_, center_.x + xr_, center_.y + yr_);
  }

  bool inside(const Point &p) const {
    double x2 = p.x*p.x;
    double y2 = p.y*p.y;

    double xr2 = xr_*xr_;
    double yr2 = yr_*yr_;

    double f = x2/xr2 + y2/yr2 - 1;

    return (f <= 0);
  }

  void moveBy(const Point &p) {
    center_ += p;
  }

  void resizeBy(const Point &ll, const Point &ur) {
    center_.x += (ll.x + ur.x)/2;
    center_.y += (ll.y + ur.y)/2;

    xr_ = xr_ + (ur.x - ll.x)/2;
    yr_ = yr_ + (ur.x - ll.x)/2;
  }

  void rotateBy(double, const Point &) { }

 private:
  Point  center_;
  double xr_ { 1 }, yr_ { 1 };
};

#endif
