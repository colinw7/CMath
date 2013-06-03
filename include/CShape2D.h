#ifndef CSHAPE_2D_H
#define CSHAPE_2D_H

#include <CPoint2D.h>
#include <CBBox2D.h>

template<typename T>
class CShape2DT {
 private:
  typedef CPoint2DT<T> Point;
  typedef CBBox2DT<T>  BBox;

 public:
  CShape2DT() { }

  virtual ~CShape2DT() { }

  virtual BBox getBBox() const = 0;

  virtual bool inside(const Point &p) const = 0;

  virtual void moveBy(const Point &p) = 0;
  virtual void resizeBy(const Point &ll, const Point &ur) = 0;
  virtual void rotateBy(double angle, const Point &o) = 0;

 protected:
  CPoint2D rotatePoint(const CPoint2D &point, double da, const CPoint2D &o) {
    double s = sin(da), c = cos(da);

    double x1 = point.x - o.x, y1 = point.y - o.y;
    double x2 = x1*c - y1*s  , y2 = x1*s + y1*c  ;

    return CPoint2D(x2 + o.x, y2 + o.y);
  }
};

typedef CShape2DT<double> CShape2D;

#endif
