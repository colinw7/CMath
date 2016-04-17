#ifndef CELLIPSE_2D_H
#define CELLIPSE_2D_H

#include <CPoint2D.h>
#include <CLine2D.h>
#include <CMathGeom2D.h>
#include <CShape2D.h>

template<typename T>
class CEllipse2DT : public CShape2DT<T> {
 private:
  typedef CPoint2DT<T>   Point;
  typedef CLine2DT<T>    Line;
  typedef CEllipse2DT<T> Ellipse;
  typedef CBBox2DT<T>    BBox;

 private:
  Point center_;
  T     xr_, yr_;

 public:
  CEllipse2DT() { }

  CEllipse2DT(T x, T y, T xr, T yr) :
   center_(x, y), xr_(xr), yr_(yr) {
  }

  CEllipse2DT(const Point &center, T xr, T yr) :
   center_(center), xr_(xr), yr_(yr) {
  }

  void setCenter(const Point &center) { center_ = center; }

  void setRadii(T xr, T yr) { xr_ = xr; yr_ = yr; }

  const Point &getCenter() const { return center_; }

  const T getXRadius() const { return xr_; }
  const T getYRadius() const { return yr_; }

  T area() { return M_PI*xr_*yr_; }

  BBox getBBox() const {
    return BBox(center_.x - xr_, center_.y - yr_, center_.x + xr_, center_.y + yr_);
  }

  bool inside(const Point &p) const {
    T x2 = p.x*p.x;
    T y2 = p.y*p.y;

    T xr2 = xr_*xr_;
    T yr2 = yr_*yr_;

    T f = x2/xr2 + y2/yr2 - 1;

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

  void rotateBy(T, const Point &) { }
};

typedef CEllipse2DT<double> CEllipse2D;

#endif
