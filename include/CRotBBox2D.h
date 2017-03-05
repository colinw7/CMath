#ifndef CROT_BBOX_2D_H
#define CROT_BBOX_2D_H

#include <CBBox2D.h>
#include <COptVal.h>

class CRotBBox2D {
 private:
  typedef COptValT<CPoint2D> OptPoint;

 public:
  CRotBBox2D() :
   angle_(0) {
  }

  CRotBBox2D(const CBBox2D &bbox, double angle=0) :
   bbox_(bbox), angle_(angle) {
  }

  CRotBBox2D(const CBBox2D &bbox, double angle, const CPoint2D &origin) :
   bbox_(bbox), angle_(angle), origin_(origin) {
  }

  CBBox2D boundingBox() const {
    CPoint2D p1 = bbox_.getLL();
    CPoint2D p2 = bbox_.getLR();
    CPoint2D p3 = bbox_.getUL();
    CPoint2D p4 = bbox_.getUR();

    CPoint2D o = origin_.getValue(bbox_.getCenter());

    CPoint2D pr1 = CMathGeom2D::RotatePoint(p1, angle_, o);
    CPoint2D pr2 = CMathGeom2D::RotatePoint(p2, angle_, o);
    CPoint2D pr3 = CMathGeom2D::RotatePoint(p3, angle_, o);
    CPoint2D pr4 = CMathGeom2D::RotatePoint(p4, angle_, o);

    CBBox2D bbox(pr1, pr2);

    bbox.add(pr3);
    bbox.add(pr4);

    return bbox;
  }

 private:
  CBBox2D  bbox_;
  double   angle_;
  OptPoint origin_;
};

#endif
