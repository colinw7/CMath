#ifndef CIBBOX_2D_H
#define CIBBOX_2D_H

#include <CIPoint2D.h>
#include <CISize2D.h>
#include <CIVector2D.h>
#include <cmath>
#include <iostream>

template<typename T>
class CIBBox2DT {
 private:
  typedef CIPoint2DT<T>  Point;
  typedef CIBBox2DT<T>   BBox;
  typedef CISize2DT<T>   Size;
  typedef CIVector2DT<T> Vector;

 public:
  CIBBox2DT() : pmin_(), pmax_(), set_(false) { }

  CIBBox2DT(const CIBBox2DT &bbox) :
   pmin_(bbox.pmin_), pmax_(bbox.pmax_), set_(bbox.set_) {
  }

  CIBBox2DT(const Point &point) :
   pmin_(point), pmax_(point), set_(true) {
  }

  CIBBox2DT(const Point &point1, const Point &point2) :
   pmin_(), pmax_(), set_(true) {
    pmin_ = Point::min(point1, point2);
    pmax_ = Point::max(point1, point2);
  }

  CIBBox2DT(const Point &point1, const Size &size) :
   pmin_(), pmax_(), set_(true) {
    Point point2 = point1 + Point(size.getWidth(), size.getHeight());

    pmin_ = Point::min(point1, point2);
    pmax_ = Point::max(point1, point2);
  }

  CIBBox2DT(T x1, T y1, T x2, T y2) :
   pmin_(), pmax_(), set_(true) {
    pmin_ = Point::min(Point(x1, y1), Point(x2, y2));
    pmax_ = Point::max(Point(x1, y1), Point(x2, y2));
  }

  CIBBox2DT &operator=(const CIBBox2DT &bbox) {
    pmin_ = bbox.pmin_;
    pmax_ = bbox.pmax_;
    set_  = bbox.set_;

    return *this;
  }

  void reset() { set_ = false; }

  bool isSet() const { return set_; }

  void set(const CIBBox2DT &bbox) {
    pmin_ = bbox.pmin_;
    pmax_ = bbox.pmax_;
    set_  = bbox.set_;
  }

  void set(const Point &point1, const Point &point2) {
    pmin_ = Point::min(point1, point2);
    pmax_ = Point::max(point1, point2);

    set_ = true;
  }

  void set(T x1, T y1, T x2, T y2) {
    pmin_ = Point::min(Point(x1, y1), Point(x2, y2));
    pmax_ = Point::max(Point(x1, y1), Point(x2, y2));

    set_ = true;
  }

  void get(T *x1, T *y1, T *x2, T *y2) const {
    *x1 = getXMin();
    *y1 = getYMin();
    *x2 = getXMax();
    *y2 = getYMax();
  }

  const Point &getMin() const { return pmin_; }
  const Point &getMax() const { return pmax_; }

  Point getMin() { return set_ ? pmin_ : Point(0,0); }
  Point getMax() { return set_ ? pmax_ : Point(0,0); }

  T getLeft  () const { return set_ ? pmin_.x : 0; }
  T getBottom() const { return set_ ? pmin_.y : 0; }
  T getRight () const { return set_ ? pmax_.x : 0; }
  T getTop   () const { return set_ ? pmax_.y : 0; }

  T getXMin() const { return set_ ? pmin_.x : 0; }
  T getYMin() const { return set_ ? pmin_.y : 0; }
  T getXMax() const { return set_ ? pmax_.x : 0; }
  T getYMax() const { return set_ ? pmax_.y : 0; }

  T getXMid() const { return set_ ? (pmin_.x + pmax_.x)/2 : 0; }
  T getYMid() const { return set_ ? (pmin_.y + pmax_.y)/2 : 0; }

  void setMin(const Point &point) {
    pmin_ = point;

    if (! set_){
      pmax_ = point;
      set_  = true;
    }
  }

  void setMax(const Point &point) {
    pmax_ = point;

    if (! set_) {
      pmin_ = point;
      set_  = true;
    }
  }

  void setMin(T x, T y) {
    pmin_.x = x;
    pmin_.y = y;

    if (! set_) {
      pmax_ = pmin_;
      set_  = true;
    }
  }

  void setMax(T x, T y) {
    pmax_.x = x;
    pmax_.y = y;

    if (! set_) {
      pmin_ = pmax_;
      set_  = true;
    }
  }

  void setXMin(T x) {
    pmin_.x = x;

    if (! set_) {
      pmin_.y = 0;

      pmax_ = pmin_;
      set_  = true;
    }
  }

  void setYMin(T y) {
    pmin_.y = y;

    if (! set_) {
      pmin_.x = 0;

      pmax_ = pmin_;
      set_  = true;
    }
  }

  void setXMax(T x) {
    pmax_.x = x;

    if (! set_) {
      pmax_.y = 0;

      pmin_ = pmax_;
      set_  = true;
    }
  }

  void setYMax(T y) {
    pmax_.y = y;

    if (! set_) {
      pmax_.x = 0;

      pmin_ = pmax_;
      set_  = true;
    }
  }

  void setWidth(T width) {
    pmax_.x = pmin_.x + width;
  }

  void setHeight(T height) {
    pmax_.y = pmin_.y + height;
  }

  void setPosition(const Point &pos) {
    Size size = getSize();

    pmin_ = pos;

    setSize(size);
  }

  void setSize(const Size &size) {
    if (! set_) {
      pmin_.x = 0;
      pmin_.y = 0;

      set_ = true;
    }

    pmax_.x = pmin_.x + size.width;
    pmax_.y = pmin_.y + size.height;
  }

  Size getSize() const {
    return Size(getWidth(), getHeight());
  }

  T getRadius() const {
    Vector radius = Vector(getMin(), getMax())/2;

    return radius.length();
  }

  T getWidth() const {
    return std::abs((T) (getXMax() - getXMin()));
  }

  T getHeight() const {
    return std::abs((T) (getYMax() - getYMin()));
  }

  T getMinDim() const {
    return std::min(getWidth(), getHeight());
  }

  T getMaxDim() const {
    return std::max(getWidth(), getHeight());
  }

  void moveBy(const Vector &delta) {
    pmin_ = (Vector(pmin_) + delta).point();
    pmax_ = (Vector(pmax_) + delta).point();
  }

  bool isPoint() const {
    return (pmin_.x == pmax_.x && pmin_.y == pmax_.y);
  }

  CIBBox2DT operator+(const Point &rhs) const {
    CIBBox2DT t(*this);

    t += rhs;

    return t;
  }

  CIBBox2DT &operator+=(const Point &rhs) {
    add(rhs);

    return *this;
  }

  CIBBox2DT operator+(const CIBBox2DT &rhs) const {
    CIBBox2DT t(*this);

    t += rhs;

    return t;
  }

  CIBBox2DT &operator+=(const CIBBox2DT &rhs) {
    add(rhs);

    return *this;
  }

  void add(T x, T y) {
    add(Point(x, y));
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

  void add(const CIBBox2DT &bbox) {
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

  bool overlaps(const CIBBox2DT &bbox) const {
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

  bool intersect(const BBox &bbox, BBox &ibbox) const {
    if (! set_ || ! bbox.set_) return false;

    if ((pmax_.x < bbox.pmin_.x || pmin_.x > bbox.pmax_.x) ||
        (pmax_.y < bbox.pmin_.y || pmin_.y > bbox.pmax_.y))
      return false;

    ibbox.pmin_.x = std::max(pmin_.x, bbox.pmin_.x);
    ibbox.pmin_.y = std::max(pmin_.y, bbox.pmin_.y);
    ibbox.pmax_.x = std::min(pmax_.x, bbox.pmax_.x);
    ibbox.pmax_.y = std::min(pmax_.y, bbox.pmax_.y);

    return true;
  }

  bool insideX(int x) const {
    return (x >= pmin_.x && x <= pmax_.x);
  }

  bool insideY(int y) const {
    return (y >= pmin_.y && y <= pmax_.y);
  }

  bool inside(const Point &point) const {
    return ((point.x >= pmin_.x && point.x <= pmax_.x) &&
            (point.y >= pmin_.y && point.y <= pmax_.y));
  }

  void expand(T delta) {
    pmin_ += delta;
    pmax_ += delta;
  }

  Point getCenter() const {
    return (pmin_ + pmax_)/2;
  }

  Point getLL() const {
    return pmin_;
  }

  Point getLR() const {
    return Point(pmax_.x, pmin_.y);
  }

  Point getUL() const {
    return Point(pmin_.x, pmax_.y);
  }

  Point getUR() const {
    return pmax_;
  }

  double area() const {
    return getWidth()*getHeight();
  }

  void print(std::ostream &os) const {
    if (! set_)
      os << "( not set )";
    else
      os << "(" << pmin_ << ") (" << pmax_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CIBBox2DT &bbox) {
    bbox.print(os);

    return os;
  }

  friend bool operator==(const CIBBox2DT &lhs, const CIBBox2DT &rhs) {
    return (lhs.pmin_ == rhs.pmin_ && lhs.pmax_ == rhs.pmax_);
  }

  friend bool operator!=(const CIBBox2DT &lhs, const CIBBox2DT &rhs) {
    return ! (lhs == rhs);
  }

 private:
  Point pmin_;
  Point pmax_;
  bool  set_ { false };
};

typedef CIBBox2DT<int> CIBBox2D;

#endif
