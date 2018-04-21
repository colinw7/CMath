#ifndef CSIZE_2D_H
#define CSIZE_2D_H

#include <CISize2D.h>
#include <CPoint2D.h>

class CSize2D {
 public:
  CSize2D() = default;

  CSize2D(double w, double h) :
   set_(true), width_(w), height_(h) {
  }

  CSize2D(const CSize2D &size) :
   set_(size.set_), width_(size.width_), height_(size.height_) {
  }

  explicit CSize2D(const CISize2D &size) :
   set_(true), width_(size.getWidth()), height_(size.getHeight()) {
  }

  bool isSet() const { return set_; }

  void set(double w, double h) {
    set_    = true;
    width_  = w;
    height_ = h;
  }

  void get(double *w, double *h) const {
    assert(set_);

    *w = width_;
    *h = height_;
  }

  double getWidth () const { assert(set_); return width_ ; }
  double getHeight() const { assert(set_); return height_; }

  void setWidth (double w) { set_ = true; width_  = w; }
  void setHeight(double h) { set_ = true; height_ = h; }

  double area() const { assert(set_); return width_*height_; }

  // m*size
  friend CSize2D operator*(double m, const CSize2D &size) {
    assert(size.set_);

    return CSize2D(m*size.width_, m*size.height_);
  }

  // size*m
  friend CSize2D operator*(const CSize2D &size, double m) {
    assert(size.set_);

    return CSize2D(m*size.width_, m*size.height_);
  }

  // size/m
  friend CSize2D operator/(const CSize2D &size, double m) {
    assert(size.set_);

    return CSize2D(size.width_/m, size.height_/m);
  }

  friend std::ostream &operator<<(std::ostream &os, const CSize2D &size) {
    if (size.set_)
      return os << "(" << size.width_ << "," << size.height_ << ")";
    else
      return os << "<not_set>";
  }

  // size + p
  friend CPoint2D operator+(const CSize2D &s, const CPoint2D &p) {
    assert(s.set_);

    return CPoint2D(p.x + s.width_, p.y + s.height_);
  }

  // p + size
  friend CPoint2D operator+(const CPoint2D &p, const CSize2D &s) {
    return (s + p);
  }

  friend bool operator==(const CSize2D &lhs, const CSize2D &rhs) {
    if (! lhs.set_ || ! rhs.set_) return false;

    return (lhs.width_ == rhs.width_ && lhs.height_ == rhs.height_);
  }

  friend bool operator!=(const CSize2D &lhs, const CSize2D &rhs) {
    return ! (lhs == rhs);
  }

 private:
  bool   set_    { false };
  double width_  { 0.0 };
  double height_ { 0.0 };
};

#endif
