#ifndef CSHAPE_2D_H
#define CSHAPE_2D_H

#include <CPoint2D.h>
#include <CBBox2D.h>

class CShape2D {
 public:
  CShape2D() { }

  virtual ~CShape2D() { }

  virtual CBBox2D getBBox() const = 0;

  virtual bool inside(const CPoint2D &p) const = 0;

  virtual void moveBy(const CPoint2D &p) = 0;
  virtual void resizeBy(const CPoint2D &ll, const CPoint2D &ur) = 0;
  virtual void rotateBy(double angle, const CPoint2D &o) = 0;

 protected:
  CPoint2D rotatePoint(const CPoint2D &point, double da, const CPoint2D &o) {
    double s = sin(da), c = cos(da);

    double x1 = point.x - o.x, y1 = point.y - o.y;
    double x2 = x1*c - y1*s  , y2 = x1*s + y1*c  ;

    return CPoint2D(x2 + o.x, y2 + o.y);
  }
};

#endif
