#ifndef CBBOX_2D_H
#define CBBOX_2D_H

#include <CMathGen.h>
#include <CPoint2D.h>
#include <CVector2D.h>
#include <CSize2D.h>

template<typename T>
class CBBox2DT {
 private:
  typedef CVector2DT<T> Vector;
  typedef CPoint2DT<T>  Point;
  typedef CBBox2DT<T>   BBox;
  typedef CSize2DT<T>   Size;

 public:
  CBBox2DT() :
   pmin_(), pmax_(), set_(false) {
  }

  CBBox2DT(const BBox &bbox) :
   pmin_(bbox.pmin_), pmax_(bbox.pmax_), set_(bbox.set_) {
  }

  CBBox2DT(const Point &point) :
   pmin_(point), pmax_(point), set_(true) {
  }

  CBBox2DT(const Point &pmin, const Point &pmax) :
   pmin_(pmin), pmax_(pmax), set_(true) {
    update();
  }

  CBBox2DT(T x1, T y1, T x2, T y2) :
   pmin_(x1, y1), pmax_(x2, y2), set_(true) {
    update();
  }

  CBBox2DT(const Point &o, const Size &s) :
   pmin_(o), pmax_(o + s), set_(true) {
    update();
  }

  void reset() { set_ = false; }

  bool isSet() const { return set_; }

  BBox operator+(const Point &rhs) const {
    BBox t(*this);

    t += rhs;

    return t;
  }

  BBox &operator+=(const Point &rhs) {
    add(rhs.x, rhs.y);

    return *this;
  }

  BBox operator+(const BBox &rhs) const {
    BBox t(*this);

    t += rhs;

    return t;
  }

  BBox &operator+=(const BBox &rhs) {
    add(rhs);

    return *this;
  }

  void add(const Point &point) {
    add(point.x, point.y);
  }

  void add(T x, T y) {
    if (! set_) {
      pmin_ = CPoint2D(x, y);
      pmax_ = pmin_;

      set_ = true;
    }
    else {
      pmin_.x = std::min(pmin_.x, x);
      pmin_.y = std::min(pmin_.y, y);
      pmax_.x = std::max(pmax_.x, x);
      pmax_.y = std::max(pmax_.y, y);
    }
  }

  void add(const BBox &bbox) {
    if (! bbox.set_) return;

    if (! set_) {
      pmin_ = bbox.pmin_;
      pmax_ = bbox.pmax_;

      set_ = true;
    }
    else {
      pmin_.x = std::min(pmin_.x, bbox.pmin_.x);
      pmin_.y = std::min(pmin_.y, bbox.pmin_.y);
      pmax_.x = std::max(pmax_.x, bbox.pmax_.x);
      pmax_.y = std::max(pmax_.y, bbox.pmax_.y);
    }
  }

  bool overlaps(const BBox &bbox) const {
    if (! set_ || ! bbox.set_) return false;

    return ((pmax_.x >= bbox.pmin_.x && pmin_.x <= bbox.pmax_.x) &&
            (pmax_.y >= bbox.pmin_.y && pmin_.y <= bbox.pmax_.y));
  }

  bool overlapsX(const BBox &bbox) const {
    if (! set_ || ! bbox.set_) return false;

    return (pmax_.x >= bbox.pmin_.x && pmin_.x <= bbox.pmax_.x);
  }

  bool overlapsY(const BBox &bbox) const {
    if (! set_ || ! bbox.set_) return false;

    return (pmax_.y >= bbox.pmin_.y && pmin_.y <= bbox.pmax_.y);
  }

  bool intersect(const BBox &bbox) const {
    if (! set_ || ! bbox.set_) return false;

    if ((pmax_.x < bbox.pmin_.x || pmin_.x > bbox.pmax_.x) ||
        (pmax_.y < bbox.pmin_.y || pmin_.y > bbox.pmax_.y))
      return false;

    return true;
  }

  bool intersect(const BBox &bbox, BBox &ibbox) const {
    if (! set_ || ! bbox.set_) return false;

    if ((pmax_.x < bbox.pmin_.x || pmin_.x > bbox.pmax_.x) ||
        (pmax_.y < bbox.pmin_.y || pmin_.y > bbox.pmax_.y))
      return false;

    ibbox.set_    = true;
    ibbox.pmin_.x = std::max(pmin_.x, bbox.pmin_.x);
    ibbox.pmin_.y = std::max(pmin_.y, bbox.pmin_.y);
    ibbox.pmax_.x = std::min(pmax_.x, bbox.pmax_.x);
    ibbox.pmax_.y = std::min(pmax_.y, bbox.pmax_.y);

    return true;
  }

  bool inside(T x, T y) const {
    if (! set_) return false;

    return ((x >= pmin_.x && x <= pmax_.x) &&
            (y >= pmin_.y && y <= pmax_.y));
  }

  bool inside(const Point &point) const {
    if (! set_) return false;

    return ((point.x >= pmin_.x && point.x <= pmax_.x) &&
            (point.y >= pmin_.y && point.y <= pmax_.y));
  }

  bool inside(const BBox &bbox) const {
    if (! set_) return false;

    return ((bbox.pmin_.x >= pmin_.x && bbox.pmax_.x <= pmax_.x) &&
            (bbox.pmin_.y >= pmin_.y && bbox.pmax_.y <= pmax_.y));
  }

  void expand(T delta) {
    if (! set_) return;

    pmin_ -= delta;
    pmax_ += delta;

    update();
  }

  void expand(T x1, T y1, T x2, T y2) {
    if (! set_) return;

    pmin_.x += x1;
    pmin_.y += y1;
    pmax_.x += x2;
    pmax_.y += y2;

    update();
  }

  T area() const {
    if (! set_) return 0.0;

    Vector diag = Vector(pmin_, pmax_);

    return fabs(diag.getX()*diag.getY());
  }

  CMathGen::AxisType2D maxAxis() const {
    if (! set_) return CMathGen::X_AXIS_2D;

    Vector diag = pmax_ - pmin_;

    if      (diag.getX() > diag.getY())
      return CMathGen::X_AXIS_2D;
    else
      return CMathGen::Y_AXIS_2D;
  }

  Point getMin() const { assert(set_); return pmin_; }
  Point getMax() const { assert(set_); return pmax_; }

  T getLeft  () const { return getXMin(); }
  T getBottom() const { return getYMin(); }
  T getRight () const { return getXMax(); }
  T getTop   () const { return getYMax(); }

  T getXMin() const { return getMin().x; }
  T getYMin() const { return getMin().y; }
  T getXMax() const { return getMax().x; }
  T getYMax() const { return getMax().y; }

  T getXMid() const { return (getXMin() + getXMax())/2; }
  T getYMid() const { return (getYMin() + getYMax())/2; }

  Point getCenter() const {
    return 0.5*(getMin() + getMax());
  }

  void setCenter(const Point &point) {
    T dx = point.x - getCenter().x;
    T dy = point.y - getCenter().y;

    moveBy(Point(dx, dy));
  }

  void setLL(const Point &point) {
    pmin_ = point;

    if (! set_) {
      pmax_ = point;
      set_  = true;
    }
  }

  void setUR(const Point &point) {
    pmax_ = point;

    if (! set_) {
      pmin_ = point;
      set_  = true;
    }
  }

  void setX(T x) { setXMin(x); }
  void setY(T y) { setYMin(y); }

  void setWidth(T width) {
    if (! set_)
      pmin_.x = 0;

    pmax_.x = pmin_.x + width;
    set_    = true;
  }

  void setHeight(T height) {
    if (! set_)
      pmin_.y = 0;

    pmax_.y = pmin_.y + height;
    set_    = true;
  }

  void setLeft  (T x) { setXMin(x); }
  void setBottom(T y) { setYMin(y); }
  void setRight (T x) { setXMax(x); }
  void setTop   (T y) { setYMax(y); }

  void setXMin(T x) {
    pmin_.x = x;
    set_    = true;
  }

  void setYMin(T y) {
    pmin_.y = y;
    set_    = true;
  }

  void setXMax(T x) {
    pmax_.x = x;
    set_    = true;
  }

  void setYMax(T y) {
    pmax_.y = y;
    set_    = true;
  }

  void setSize(const Size &size) {
    if (! set_) {
      pmin_.x = 0;
      pmin_.y = 0;
    }

    pmax_.x = pmin_.x + size.width;
    pmax_.y = pmin_.y + size.height;
    set_    = true;
  }

  Point getLL() const { return getMin(); }
  Point getLR() const { assert(set_); return Point(pmax_.x, pmin_.y); }
  Point getUL() const { assert(set_); return Point(pmin_.x, pmax_.y); }
  Point getUR() const { return getMax(); }

  Size getSize() const {
    return Size(getWidth(), getHeight());
  }

  T getRadius() const {
    Vector radius = 0.5*Vector(getMin(), getMax());

    return radius.length();
  }

  T getWidth() const {
    return fabs(getXMax() - getXMin());
  }

  T getHeight() const {
    return fabs(getYMax() - getYMin());
  }

  BBox &moveXTo(double x) {
    assert(set_);

    double dx = x - pmin_.x;

    pmin_.x += dx;
    pmax_.x += dx;

    update();

    return *this;
  }

  BBox &moveYTo(double y) {
    assert(set_);

    double dy = y - pmin_.y;

    pmin_.y += dy;
    pmax_.y += dy;

    update();

    return *this;
  }

  BBox &moveTo(const Point &p) {
    assert(set_);

    Point delta = p - pmin_;

    pmin_ += delta;
    pmax_ += delta;

    update();

    return *this;
  }

  BBox &moveBy(const Vector &delta) {
    assert(set_);

    pmin_ += delta;
    pmax_ += delta;

    update();

    return *this;
  }

  BBox &moveBy(const Point &delta) {
    assert(set_);

    pmin_ += delta;
    pmax_ += delta;

    update();

    return *this;
  }

  BBox &moveBy(const Point &dmin, const Point &dmax) {
    assert(set_);

    pmin_ += dmin;
    pmax_ += dmax;

    update();

    return *this;
  }

  void print(std::ostream &os) const {
    if (! set_)
      os << "( not set )";
    else
      os << "(" << pmin_ << ") (" << pmax_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const BBox &bbox) {
    bbox.print(os);

    return os;
  }

 private:
  void update() {
    assert(set_);

    if (pmin_.x > pmax_.x) std::swap(pmin_.x, pmax_.x);
    if (pmin_.y > pmax_.y) std::swap(pmin_.y, pmax_.y);
  }

 private:
  Point pmin_;
  Point pmax_;
  bool  set_;
};

typedef CBBox2DT<double> CBBox2D;

#endif
