#ifndef CBBOX_3D_H
#define CBBOX_3D_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CVector3D.h>

template<typename T>
class CBBox3DT {
 private:
  typedef CPoint3DT<T>  Point;
  typedef CVector3DT<T> Vector;

 protected:
  Point pmin_;
  Point pmax_;
  bool  set_;

 public:
  CBBox3DT() :
   pmin_(), pmax_(), set_(false) {
  }

  CBBox3DT(const CBBox3DT &bbox) :
   pmin_(bbox.pmin_), pmax_(bbox.pmax_), set_(bbox.set_) {
  }

  CBBox3DT(const Point &point) {
    pmin_ = point;
    pmax_ = point;

    set_ = true;
  }

  CBBox3DT(const Point &point1, const Point &point2) {
    pmin_ = Point::min(point1, point2);
    pmax_ = Point::max(point1, point2);

    set_ = true;
  }

  CBBox3DT(T x1, T y1, T z1, T x2, T y2, T z2) {
    pmin_ = Point(std::min(x1,x2),std::min(y1,y2),std::min(z1,z2));
    pmax_ = Point(std::max(x1,x2),std::max(y1,y2),std::max(z1,z2));

    set_ = true;
  }

  void reset() { set_ = false; }

  bool isSet() const { return set_; }

  const Point &getMin() const { return pmin_; }
  const Point &getMax() const { return pmax_; }

  void setMin(const Point &pmin) { pmin_ = pmin; }
  void setMax(const Point &pmax) { pmax_ = pmax; }

  void setXMin(T x) { pmin_.x = x; }
  void setYMin(T y) { pmin_.y = y; }
  void setZMin(T z) { pmin_.z = z; }

  void setXMax(T x) { pmax_.x = x; }
  void setYMax(T y) { pmax_.y = y; }
  void setZMax(T z) { pmax_.z = z; }

  T getXMin() const { return pmin_.x; }
  T getYMin() const { return pmin_.y; }
  T getZMin() const { return pmin_.z; }

  T getXMax() const { return pmax_.x; }
  T getYMax() const { return pmax_.y; }
  T getZMax() const { return pmax_.z; }

  CBBox3DT operator+(const Point &rhs) const {
    CBBox3DT t(*this);

    t += rhs;

    return t;
  }

  CBBox3DT &operator+=(const Point &rhs) {
    add(rhs);

    return *this;
  }

  CBBox3DT operator+(const Vector &rhs) const {
    CBBox3DT t(*this);

    t += rhs;

    return t;
  }

  CBBox3DT &operator+=(const Vector &rhs) {
    add(rhs);

    return *this;
  }

  CBBox3DT operator+(const CBBox3DT &rhs) const {
    CBBox3DT t(*this);

    t += rhs;

    return t;
  }

  CBBox3DT &operator+=(const CBBox3DT &rhs) {
    add(rhs);

    return *this;
  }

  CBBox3DT &operator*=(T rhs) {
    scale(rhs);

    return *this;
  }

  void add(T x, T y, T z) {
    add(Point(x, y, z));
  }

  void add(const Point &point) {
    if (! set_) {
      pmin_ = point;
      pmax_ = point;

      set_ = true;
    }
    else {
      pmin_ = Point::min(pmin_, point);
      pmax_ = Point::max(pmax_, point);
    }
  }

  void add(const Vector &vector)  {
    add(vector.point());
  }

  void add(const CBBox3DT &bbox) {
    if (! bbox.set_) return;

    if (! set_) {
      pmin_ = bbox.pmin_;
      pmax_ = bbox.pmax_;

      set_ = true;
    }
    else {
      pmin_ = Point::min(pmin_, bbox.pmin_);
      pmax_ = Point::max(pmax_, bbox.pmax_);
    }
  }

  bool overlaps(const CBBox3DT &bbox) {
    return ((pmax_.x >= bbox.pmin_.x && pmin_.x <= bbox.pmax_.x) &&
            (pmax_.y >= bbox.pmin_.y && pmin_.y <= bbox.pmax_.y) &&
            (pmax_.z >= bbox.pmin_.z && pmin_.z <= bbox.pmax_.z));
  }

  bool inside(const CBBox3DT &bbox) {
    if ((pmin_.x > bbox.pmax_.x || pmax_.x < bbox.pmax_.x) ||
        (pmin_.y > bbox.pmax_.y || pmax_.y < bbox.pmax_.y) ||
        (pmin_.z > bbox.pmax_.z || pmax_.z < bbox.pmax_.z))
      return false;

    return true;
  }

  bool inside(const Point &point) {
    return ((point.x >= pmin_.x && point.x <= pmax_.x) &&
            (point.y >= pmin_.y && point.y <= pmax_.y) &&
            (point.z >= pmin_.z && point.z <= pmax_.z));
  }

  void expand(T delta) {
    pmin_ += delta;
    pmax_ += delta;
  }

  T volume() const {
    Vector diag = getSize();

    return diag.getX()*diag.getY()*diag.getZ();
  }

  CMathGen::AxisType3D maxAxis() const {
    Vector diag = getSize();

    if      (diag.getX() > diag.getY() && diag.getX() > diag.getZ())
      return CMathGen::X_AXIS_3D;
    else if (diag.getY() > diag.getZ())
      return CMathGen::Y_AXIS_3D;
    else
      return CMathGen::Z_AXIS_3D;
  }

  Point getCenter() const {
    return 0.5*(pmin_ + pmax_);
  }

  Vector getSize() const {
    return Vector(pmin_, pmax_);
  }

  T getXSize() const {
    return (pmax_.x - pmin_.x);
  }

  T getYSize() const {
    return (pmax_.y - pmin_.y);
  }

  T getZSize() const {
    return (pmax_.z - pmin_.z);
  }

  T getRadius() const {
    Vector radius = 0.5*getSize();

    return radius.length();
  }

  void scale(T factor) {
    Point center = getCenter();

    pmin_ = center + Vector(center, pmin_)*factor;
    pmax_ = center + Vector(center, pmax_)*factor;
  }

 public:
  bool lineInteracts(const Point &l1, const Point &l2) const {
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
  bool lineIntersect(const Point &l1, const Point &l2, Point &pi) const {
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
  int getIntersect(T d1, T d2, const Point &p1, const Point &p2, Point &pi) const {
    if ((d1*d2) >= 0.0) return 0;

    if (d1 == d2) return 0;

    pi = p1 + (p2 - p1)*(-d1/(d2 - d1));

    return 1;
  }

  bool inBox(const Point &pi, int axis) const {
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

  friend std::ostream &operator<<(std::ostream &os, const CBBox3DT &bbox) {
    bbox.print(os);

    return os;
  }
};

typedef CBBox3DT<double> CBBox3D;

#endif
