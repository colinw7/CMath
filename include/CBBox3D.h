#ifndef CBBOX_3D_H
#define CBBOX_3D_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CVector3D.h>

class CBBox3D {
 public:
  CBBox3D() :
   pmin_(), pmax_(), set_(false) {
  }

  CBBox3D(const CBBox3D &bbox) :
   pmin_(bbox.pmin_), pmax_(bbox.pmax_), set_(bbox.set_) {
  }

  CBBox3D(const CPoint3D &point) {
    pmin_ = point;
    pmax_ = point;

    set_ = true;
  }

  CBBox3D(const CPoint3D &point1, const CPoint3D &point2) {
    pmin_ = CPoint3D::min(point1, point2);
    pmax_ = CPoint3D::max(point1, point2);

    set_ = true;
  }

  CBBox3D(double x1, double y1, double z1, double x2, double y2, double z2) {
    pmin_ = CPoint3D(std::min(x1,x2),std::min(y1,y2),std::min(z1,z2));
    pmax_ = CPoint3D(std::max(x1,x2),std::max(y1,y2),std::max(z1,z2));

    set_ = true;
  }

  void reset() { set_ = false; }

  bool isSet() const { return set_; }

  const CPoint3D &getMin() const { return pmin_; }
  const CPoint3D &getMax() const { return pmax_; }

  void setMin(const CPoint3D &pmin) { pmin_ = pmin; }
  void setMax(const CPoint3D &pmax) { pmax_ = pmax; }

  void setXMin(double x) { pmin_.x = x; }
  void setYMin(double y) { pmin_.y = y; }
  void setZMin(double z) { pmin_.z = z; }

  void setXMax(double x) { pmax_.x = x; }
  void setYMax(double y) { pmax_.y = y; }
  void setZMax(double z) { pmax_.z = z; }

  double getXMin() const { return pmin_.x; }
  double getYMin() const { return pmin_.y; }
  double getZMin() const { return pmin_.z; }

  double getXMax() const { return pmax_.x; }
  double getYMax() const { return pmax_.y; }
  double getZMax() const { return pmax_.z; }

  CPoint3D getBottomMax() const { return CPoint3D(pmax_.x, pmax_.y, pmin_.z); }
  CPoint3D getTopMin   () const { return CPoint3D(pmin_.x, pmin_.y, pmax_.z); }

  CBBox3D operator+(const CPoint3D &rhs) const {
    CBBox3D t(*this);

    t += rhs;

    return t;
  }

  CBBox3D &operator+=(const CPoint3D &rhs) {
    add(rhs);

    return *this;
  }

  CBBox3D operator+(const CVector3D &rhs) const {
    CBBox3D t(*this);

    t += rhs;

    return t;
  }

  CBBox3D &operator+=(const CVector3D &rhs) {
    add(rhs);

    return *this;
  }

  CBBox3D operator+(const CBBox3D &rhs) const {
    CBBox3D t(*this);

    t += rhs;

    return t;
  }

  CBBox3D &operator+=(const CBBox3D &rhs) {
    add(rhs);

    return *this;
  }

  CBBox3D &operator*=(double rhs) {
    scale(rhs);

    return *this;
  }

  void add(double x, double y, double z) {
    add(CPoint3D(x, y, z));
  }

  void add(const CPoint3D &point) {
    if (! set_) {
      pmin_ = point;
      pmax_ = point;

      set_ = true;
    }
    else {
      pmin_ = CPoint3D::min(pmin_, point);
      pmax_ = CPoint3D::max(pmax_, point);
    }
  }

  void add(const CVector3D &vector)  {
    add(vector.point());
  }

  void add(const CBBox3D &bbox) {
    if (! bbox.set_) return;

    if (! set_) {
      pmin_ = bbox.pmin_;
      pmax_ = bbox.pmax_;

      set_ = true;
    }
    else {
      pmin_ = CPoint3D::min(pmin_, bbox.pmin_);
      pmax_ = CPoint3D::max(pmax_, bbox.pmax_);
    }
  }

  bool overlaps(const CBBox3D &bbox) {
    return ((pmax_.x >= bbox.pmin_.x && pmin_.x <= bbox.pmax_.x) &&
            (pmax_.y >= bbox.pmin_.y && pmin_.y <= bbox.pmax_.y) &&
            (pmax_.z >= bbox.pmin_.z && pmin_.z <= bbox.pmax_.z));
  }

  bool inside(const CBBox3D &bbox) {
    if ((pmin_.x > bbox.pmax_.x || pmax_.x < bbox.pmax_.x) ||
        (pmin_.y > bbox.pmax_.y || pmax_.y < bbox.pmax_.y) ||
        (pmin_.z > bbox.pmax_.z || pmax_.z < bbox.pmax_.z))
      return false;

    return true;
  }

  bool inside(const CPoint3D &point) {
    return ((point.x >= pmin_.x && point.x <= pmax_.x) &&
            (point.y >= pmin_.y && point.y <= pmax_.y) &&
            (point.z >= pmin_.z && point.z <= pmax_.z));
  }

  void expand(double delta) {
    pmin_ += delta;
    pmax_ += delta;
  }

  double volume() const {
    CVector3D diag = getSize();

    return diag.getX()*diag.getY()*diag.getZ();
  }

  CMathGen::AxisType3D maxAxis() const {
    CVector3D diag = getSize();

    if      (diag.getX() > diag.getY() && diag.getX() > diag.getZ())
      return CMathGen::X_AXIS_3D;
    else if (diag.getY() > diag.getZ())
      return CMathGen::Y_AXIS_3D;
    else
      return CMathGen::Z_AXIS_3D;
  }

  CPoint3D getCenter() const {
    return 0.5*(pmin_ + pmax_);
  }

  CVector3D getSize() const {
    return CVector3D(pmin_, pmax_);
  }

  double getXSize() const {
    return (pmax_.x - pmin_.x);
  }

  double getYSize() const {
    return (pmax_.y - pmin_.y);
  }

  double getZSize() const {
    return (pmax_.z - pmin_.z);
  }

  double getRadius() const {
    CVector3D radius = 0.5*getSize();

    return radius.length();
  }

  void scale(double factor) {
    CPoint3D center = getCenter();

    pmin_ = center + CVector3D(center, pmin_)*factor;
    pmax_ = center + CVector3D(center, pmax_)*factor;
  }

 public:
  bool lineInteracts(const CPoint3D &l1, const CPoint3D &l2) const {
    // while line outside (no overlap)
    if (l1.x < pmin_.x && l2.x < pmin_.x) return false;
    if (l1.x > pmax_.x && l2.x > pmax_.x) return false;
    if (l1.y < pmin_.y && l2.y < pmin_.y) return false;
    if (l1.y > pmax_.y && l2.y > pmax_.y) return false;
    if (l1.z < pmin_.z && l2.z < pmin_.z) return false;
    if (l1.z > pmax_.z && l2.z > pmax_.z) return false;

    return true;
  }

  // returns true if line (l1, l2) intersects with the box (pmin_, pmax_)
  // returns intersection point in pi
  bool lineIntersect(const CPoint3D &l1, const CPoint3D &l2, CPoint3D &pi) const {
    // while line outside (no overlap)
    if (l1.x < pmin_.x && l2.x < pmin_.x) return false;
    if (l1.x > pmax_.x && l2.x > pmax_.x) return false;
    if (l1.y < pmin_.y && l2.y < pmin_.y) return false;
    if (l1.y > pmax_.y && l2.y > pmax_.y) return false;
    if (l1.z < pmin_.z && l2.z < pmin_.z) return false;
    if (l1.z > pmax_.z && l2.z > pmax_.z) return false;

    // if l1 is inside we are done ?
    if (l1.x > pmin_.x && l1.x < pmax_.x &&
        l1.y > pmin_.y && l1.y < pmax_.y &&
        l1.z > pmin_.z && l1.z < pmax_.z) {
      pi = l1;
      return true;
    }

#if 0
    // if get here either l1 is inside or l2 is inside
    if (l1.x >= pmin_.x && l1.x <= pmax_.x &&
        l1.y >= pmin_.y && l1.y <= pmax_.y &&
        l1.z >= pmin_.z && l1.z <= pmax_.z &&
        l2.x >= pmin_.x && l2.x <= pmax_.x &&
        l2.y >= pmin_.y && l2.y <= pmax_.y &&
        l2.z >= pmin_.z && l2.z <= pmax_.z) {
      pi = l1;
      return true;
    }
#endif

    // check each side for intersect
    if ((getIntersect(l1.x - pmin_.x, l2.x - pmin_.x, l1, l2, pi) && inBox(pi, 1)) ||
        (getIntersect(l1.x - pmax_.x, l2.x - pmax_.x, l1, l2, pi) && inBox(pi, 1)) ||
        (getIntersect(l1.y - pmin_.y, l2.y - pmin_.y, l1, l2, pi) && inBox(pi, 2)) ||
        (getIntersect(l1.y - pmax_.y, l2.y - pmax_.y, l1, l2, pi) && inBox(pi, 2)) ||
        (getIntersect(l1.z - pmin_.z, l2.z - pmin_.z, l1, l2, pi) && inBox(pi, 3)) ||
        (getIntersect(l1.z - pmax_.z, l2.z - pmax_.z, l1, l2, pi) && inBox(pi, 3)))
      return true;

    return false;
  }

 private:
  int getIntersect(double d1, double d2, const CPoint3D &p1, const CPoint3D &p2,
                   CPoint3D &pi) const {
    if ((d1*d2) >= 0.0) return 0;

    if (d1 == d2) return 0;

    pi = p1 + (p2 - p1)*(-d1/(d2 - d1));

    return 1;
  }

  bool inBox(const CPoint3D &pi, int axis) const {
    if      (axis == 1) {
      if (pi.z > pmin_.z && pi.z < pmax_.z && pi.y > pmin_.y && pi.y < pmax_.y) return true;
    }
    else if (axis == 2) {
      if (pi.z > pmin_.z && pi.z < pmax_.z && pi.x > pmin_.x && pi.x < pmax_.x) return true;
    }
    else if (axis == 3) {
      if (pi.x > pmin_.x && pi.x < pmax_.x && pi.y > pmin_.y && pi.y < pmax_.y) return true;
    }

    return false;
  }

 public:
  void print(std::ostream &os) const {
    if (! set_)
      os << "( not set )";
    else
      os << "(" << pmin_ << ") (" << pmax_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CBBox3D &bbox) {
    bbox.print(os);

    return os;
  }

 protected:
  CPoint3D pmin_;
  CPoint3D pmax_;
  bool     set_ { false };
};

#endif
