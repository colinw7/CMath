#ifndef CILINE_2D_H
#define CILINE_2D_H

#include <CIPoint2D.h>
#include <CIVector2D.h>

class CILine2D {
 public:
  CILine2D() :
   p1_(), p2_(), v_() {
  }

  CILine2D(const CILine2D &line) :
    p1_(line.p1_), p2_(line.p2_), v_(line.v_) {
  }

  CILine2D(int x1, int y1, int x2, int y2) :
   p1_(x1, y1), p2_(x2, y2), v_(x2 - x1, y2 - y1) {
  }

  CILine2D(const CIPoint2D &p0, const CIPoint2D &p1) :
   p1_(p0), p2_(p1), v_(p0, p1) {
  }

  const CIPoint2D  &start () const { return p1_; }
  const CIPoint2D  &end   () const { return p2_; }
  const CIVector2D &vector() const { return v_ ; }

  void setStart(const CIPoint2D &start) { p1_ = start; }
  void setEnd  (const CIPoint2D &end  ) { p2_ = end  ; }

  void print(std::ostream &os) const {
    os << p1_ << " " << p2_;
  }

  friend std::ostream &operator<<(std::ostream &os, const CILine2D &line) {
    line.print(os);

    return os;
  }

 private:
  CIPoint2D  p1_, p2_;
  CIVector2D v_;
};

#endif
