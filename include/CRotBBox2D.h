#ifndef CROT_BBOX_2D_H
#define CROT_BBOX_2D_H

#include <CBBox2D.h>

#include <optional>

class CRotBBox2D {
 public:
  CRotBBox2D() { }

  CRotBBox2D(const CBBox2D &bbox, double angle=0.0) :
   bbox_(bbox), angle_(angle) {
  }

  CRotBBox2D(const CBBox2D &bbox, double angle, const CPoint2D &origin) :
   bbox_(bbox), angle_(angle), origin_(origin) {
  }

  CBBox2D boundingBox() const {
    auto p1 = bbox_.getLL();
    auto p2 = bbox_.getLR();
    auto p3 = bbox_.getUL();
    auto p4 = bbox_.getUR();

    auto o = origin_.value_or(bbox_.getCenter());

    auto pr1 = CMathGeom2D::RotatePoint(p1, angle_, o);
    auto pr2 = CMathGeom2D::RotatePoint(p2, angle_, o);
    auto pr3 = CMathGeom2D::RotatePoint(p3, angle_, o);
    auto pr4 = CMathGeom2D::RotatePoint(p4, angle_, o);

    CBBox2D bbox(pr1, pr2);

    bbox.add(pr3);
    bbox.add(pr4);

    return bbox;
  }

 private:
  using OptPoint = std::optional<CPoint2D>;

 private:
  CBBox2D  bbox_;
  double   angle_ { 0.0 };
  OptPoint origin_;
};

#endif
