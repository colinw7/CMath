#ifndef CIBBOX_2D_H
#define CIBBOX_2D_H

#include <iostream>
#include <CIPoint2D.h>
#include <CISize2D.h>
#include <CIVector2D.h>

class CIBBox2D {
 public:
  CIBBox2D() : pmin_(), pmax_(), set_(false) { }

  CIBBox2D(const CIBBox2D &bbox) :
   pmin_(bbox.pmin_), pmax_(bbox.pmax_), set_(bbox.set_) {
  }

  explicit CIBBox2D(const CIPoint2D &point) :
   pmin_(point), pmax_(point), set_(true) {
  }

  CIBBox2D(const CIPoint2D &point1, const CIPoint2D &point2) :
   pmin_(), pmax_(), set_(true) {
    pmin_ = CIPoint2D::min(point1, point2);
    pmax_ = CIPoint2D::max(point1, point2);
  }

  CIBBox2D(const CIPoint2D &point1, const CISize2D &size) :
   pmin_(), pmax_(), set_(true) {
    CIPoint2D point2 = point1 + CIPoint2D(size.getWidth(), size.getHeight());

    pmin_ = CIPoint2D::min(point1, point2);
    pmax_ = CIPoint2D::max(point1, point2);
  }

  CIBBox2D(int x1, int y1, int x2, int y2) :
   pmin_(), pmax_(), set_(true) {
    pmin_ = CIPoint2D::min(CIPoint2D(x1, y1), CIPoint2D(x2, y2));
    pmax_ = CIPoint2D::max(CIPoint2D(x1, y1), CIPoint2D(x2, y2));
  }

  CIBBox2D &operator=(const CIBBox2D &bbox) {
    pmin_ = bbox.pmin_;
    pmax_ = bbox.pmax_;
    set_  = bbox.set_;

    return *this;
  }

  void reset() { set_ = false; }

  bool isSet() const { return set_; }

  void set(const CIBBox2D &bbox) {
    pmin_ = bbox.pmin_;
    pmax_ = bbox.pmax_;
    set_  = bbox.set_;
  }

  void set(const CIPoint2D &point1, const CIPoint2D &point2) {
    pmin_ = CIPoint2D::min(point1, point2);
    pmax_ = CIPoint2D::max(point1, point2);

    set_ = true;
  }

  void set(int x1, int y1, int x2, int y2) {
    pmin_ = CIPoint2D::min(CIPoint2D(x1, y1), CIPoint2D(x2, y2));
    pmax_ = CIPoint2D::max(CIPoint2D(x1, y1), CIPoint2D(x2, y2));

    set_ = true;
  }

  void get(int *x1, int *y1, int *x2, int *y2) const {
    *x1 = getXMin();
    *y1 = getYMin();
    *x2 = getXMax();
    *y2 = getYMax();
  }

  const CIPoint2D &getMin() const { return pmin_; }
  const CIPoint2D &getMax() const { return pmax_; }

  CIPoint2D getMin() { return set_ ? pmin_ : CIPoint2D(0,0); }
  CIPoint2D getMax() { return set_ ? pmax_ : CIPoint2D(0,0); }

  int getLeft  () const { return set_ ? pmin_.x : 0; }
  int getBottom() const { return set_ ? pmin_.y : 0; }
  int getRight () const { return set_ ? pmax_.x : 0; }
  int getTop   () const { return set_ ? pmax_.y : 0; }

  int getXMin() const { return set_ ? pmin_.x : 0; }
  int getYMin() const { return set_ ? pmin_.y : 0; }
  int getXMax() const { return set_ ? pmax_.x : 0; }
  int getYMax() const { return set_ ? pmax_.y : 0; }

  int getXMid() const { return set_ ? (pmin_.x + pmax_.x)/2 : 0; }
  int getYMid() const { return set_ ? (pmin_.y + pmax_.y)/2 : 0; }

  void setMin(const CIPoint2D &point) {
    pmin_ = point;

    if (! set_){
      pmax_ = point;
      set_  = true;
    }
  }

  void setMax(const CIPoint2D &point) {
    pmax_ = point;

    if (! set_) {
      pmin_ = point;
      set_  = true;
    }
  }

  void setMin(int x, int y) {
    pmin_.x = x;
    pmin_.y = y;

    if (! set_) {
      pmax_ = pmin_;
      set_  = true;
    }
  }

  void setMax(int x, int y) {
    pmax_.x = x;
    pmax_.y = y;

    if (! set_) {
      pmin_ = pmax_;
      set_  = true;
    }
  }

  void setXMin(int x) {
    pmin_.x = x;

    if (! set_) {
      pmin_.y = 0;

      pmax_ = pmin_;
      set_  = true;
    }
  }

  void setYMin(int y) {
    pmin_.y = y;

    if (! set_) {
      pmin_.x = 0;

      pmax_ = pmin_;
      set_  = true;
    }
  }

  void setXMax(int x) {
    pmax_.x = x;

    if (! set_) {
      pmax_.y = 0;

      pmin_ = pmax_;
      set_  = true;
    }
  }

  void setYMax(int y) {
    pmax_.y = y;

    if (! set_) {
      pmax_.x = 0;

      pmin_ = pmax_;
      set_  = true;
    }
  }

  void setWidth(int width) {
    pmax_.x = pmin_.x + width;
  }

  void setHeight(int height) {
    pmax_.y = pmin_.y + height;
  }

  void setPosition(const CIPoint2D &pos) {
    CISize2D size = getSize();

    pmin_ = pos;

    setSize(size);
  }

  void setSize(const CISize2D &size) {
    if (! set_) {
      pmin_.x = 0;
      pmin_.y = 0;

      set_ = true;
    }

    pmax_.x = pmin_.x + size.width;
    pmax_.y = pmin_.y + size.height;
  }

  CISize2D getSize() const {
    return CISize2D(getWidth(), getHeight());
  }

  int getRadius() const {
    CIVector2D radius = CIVector2D(getMin(), getMax())/2;

    return radius.length();
  }

  int getWidth() const {
    return std::abs((int) (getXMax() - getXMin()));
  }

  int getHeight() const {
    return std::abs((int) (getYMax() - getYMin()));
  }

  int getMinDim() const {
    return std::min(getWidth(), getHeight());
  }

  int getMaxDim() const {
    return std::max(getWidth(), getHeight());
  }

  void moveBy(const CIVector2D &delta) {
    pmin_ = (CIVector2D(pmin_) + delta).point();
    pmax_ = (CIVector2D(pmax_) + delta).point();
  }

  bool isPoint() const {
    return (pmin_.x == pmax_.x && pmin_.y == pmax_.y);
  }

  CIBBox2D operator+(const CIPoint2D &rhs) const {
    CIBBox2D t(*this);

    t += rhs;

    return t;
  }

  CIBBox2D &operator+=(const CIPoint2D &rhs) {
    add(rhs);

    return *this;
  }

  CIBBox2D operator+(const CIBBox2D &rhs) const {
    CIBBox2D t(*this);

    t += rhs;

    return t;
  }

  CIBBox2D &operator+=(const CIBBox2D &rhs) {
    add(rhs);

    return *this;
  }

  void add(int x, int y) {
    add(CIPoint2D(x, y));
  }

  void add(const CIPoint2D &point) {
    if (! set_) {
      pmin_ = point;
      pmax_ = point;

      set_ = true;
    }
    else {
      pmin_ = CIPoint2D::min(pmin_, point);
      pmax_ = CIPoint2D::max(pmax_, point);
    }
  }

  void add(const CIBBox2D &bbox) {
    if (! bbox.set_) return;

    if (! set_) {
      pmin_ = bbox.pmin_;
      pmax_ = bbox.pmax_;

      set_ = true;
    }
    else {
      pmin_ = CIPoint2D::min(pmin_, bbox.pmin_);
      pmax_ = CIPoint2D::max(pmax_, bbox.pmax_);
    }
  }

  bool overlaps(const CIBBox2D &bbox) const {
    return ((pmax_.x >= bbox.pmin_.x && pmin_.x <= bbox.pmax_.x) &&
            (pmax_.y >= bbox.pmin_.y && pmin_.y <= bbox.pmax_.y));
  }

  bool overlapsX(const CIBBox2D &bbox) const {
    if (! set_ || ! bbox.set_) return false;

    return (pmax_.x >= bbox.pmin_.x && pmin_.x <= bbox.pmax_.x);
  }

  bool overlapsY(const CIBBox2D &bbox) const {
    if (! set_ || ! bbox.set_) return false;

    return (pmax_.y >= bbox.pmin_.y && pmin_.y <= bbox.pmax_.y);
  }

  bool intersect(const CIBBox2D &bbox, CIBBox2D &ibbox) const {
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
    if (! set_) return false;
    return (x >= pmin_.x && x <= pmax_.x);
  }

  bool insideY(int y) const {
    if (! set_) return false;
    return (y >= pmin_.y && y <= pmax_.y);
  }

  bool inside(const CIPoint2D &point) const {
    if (! set_) return false;
    return ((point.x >= pmin_.x && point.x <= pmax_.x) &&
            (point.y >= pmin_.y && point.y <= pmax_.y));
  }

  void expand(int delta) {
    pmin_ += delta;
    pmax_ += delta;
  }

  CIPoint2D getCenter() const {
    return (pmin_ + pmax_)/2;
  }

  CIPoint2D getLL() const {
    return pmin_;
  }

  CIPoint2D getLR() const {
    return CIPoint2D(pmax_.x, pmin_.y);
  }

  CIPoint2D getUL() const {
    return CIPoint2D(pmin_.x, pmax_.y);
  }

  CIPoint2D getUR() const {
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

  friend std::ostream &operator<<(std::ostream &os, const CIBBox2D &bbox) {
    bbox.print(os);

    return os;
  }

  friend bool operator==(const CIBBox2D &lhs, const CIBBox2D &rhs) {
    return (lhs.pmin_ == rhs.pmin_ && lhs.pmax_ == rhs.pmax_);
  }

  friend bool operator!=(const CIBBox2D &lhs, const CIBBox2D &rhs) {
    return ! (lhs == rhs);
  }

 private:
  CIPoint2D pmin_;
  CIPoint2D pmax_;
  bool      set_ { false };
};

#endif
