#ifndef CROT_BBOX_2D_H
#define CROT_BBOX_2D_H

#include <CBBox2D.h>
#include <COptVal.h>

template<typename T>
class CRotBBox2DT {
 private:
  typedef CRotBBox2DT<T>  RotBBox;
  typedef CBBox2DT<T>     BBox;
  typedef CPoint2DT<T>    Point;
  typedef COptValT<Point> OptPoint;

 public:
  CRotBBox2DT() :
   angle_(0) {
  }

  CRotBBox2DT(const BBox &bbox, double angle=0) :
   bbox_(bbox), angle_(angle) {
  }

  CRotBBox2DT(const BBox &bbox, double angle, const Point &origin) :
   bbox_(bbox), angle_(angle), origin_(origin) {
  }

  BBox boundingBox() const {
    CPoint2D p1 = bbox_.getLL();
    CPoint2D p2 = bbox_.getLR();
    CPoint2D p3 = bbox_.getUL();
    CPoint2D p4 = bbox_.getUR();

    CPoint2D o = origin_.getValue(bbox_.getCenter());

    CPoint2D pr1 = CMathGeom2D::RotatePoint(p1, angle_, o);
    CPoint2D pr2 = CMathGeom2D::RotatePoint(p2, angle_, o);
    CPoint2D pr3 = CMathGeom2D::RotatePoint(p3, angle_, o);
    CPoint2D pr4 = CMathGeom2D::RotatePoint(p4, angle_, o);

    BBox bbox(pr1, pr2);

    bbox.add(pr3);
    bbox.add(pr4);

    return bbox;
  }

 private:
  BBox     bbox_;
  double   angle_;
  OptPoint origin_;
};

typedef CRotBBox2DT<double> CRotBBox2D;

#endif
