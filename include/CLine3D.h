#ifndef CLINE_3D_H
#define CLINE_3D_H

// TODO: add transform, distance, closest points

#include <CPoint3D.h>
#include <CVector3D.h>
#include <CMathGen.h>

class CLine3D {
 private:
  typedef CMathGen::IntersectType IntersectType;

 public:
  // constructor/destructor
  CLine3D() :
   p0_(0,0,0), p1_(1,0,0), v_(1,0,0) {
  }

  CLine3D(double x1, double y1, double z1, double x2, double y2, double z2) :
   p0_(x1, y1, z1), p1_(x2, y2, z2), v_(x2 - x1, y2 - y1, z2 - z1) {
  }

  CLine3D(const CPoint3D &p0, const CPoint3D &p1) :
   p0_(p0), p1_(p1), v_(p0, p1) {
  }

  CLine3D(const CPoint3D &o, const CVector3D &d) :
   p0_(o), p1_(o + d), v_(d) {
  }

  CLine3D(const CLine3D &rhs) :
   p0_(rhs.p0_), p1_(rhs.p1_), v_(rhs.v_) {
  }

  CLine3D &operator=(const CLine3D &rhs) {
    p0_= rhs.p0_; p1_= rhs.p1_; v_= rhs.v_;

    return *this;
  }

 ~CLine3D() { }

  //------

  // accessors
  CPoint3D &start() { return p0_; }
  CPoint3D &end  () { return p1_; }

  const CPoint3D &start() const { return p0_; }
  const CPoint3D &end  () const { return p1_; }

  const CVector3D &vector() const { return v_; }

  CPoint3D point(double mu) const { return p0_ + mu*v_; }

  CPoint3D center() const { return (p0_ + p1_)/2.0; }

  void setStart(const CPoint3D &p0) {
    p0_ = p0;
    v_  = CVector3D(p0_, p1_);
  }

  void setEnd(const CPoint3D &p1) {
    p1_ = p1;
    v_  = CVector3D(p0_, p1_);
  }

  void setVector(const CVector3D &v) {
    v_  = v;
    p1_ = p0_ + v;
  }

  void flip() {
    std::swap(p0_, p1_);
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << p0_ << " " << p1_;
  }

  friend std::ostream &operator<<(std::ostream &os, const CLine3D &line) {
    line.print(os);

    return os;
  }

  //------

  // comparison
  int cmp(const CLine3D &v) const {
    if      (p0_ < v.p0_) return -1;
    else if (p0_ > v.p0_) return  1;
    else if (p1_ < v.p1_) return -1;
    else if (p1_ > v.p1_) return  1;
    else                  return  0;
  }

  friend bool operator==(const CLine3D &lhs, const CLine3D &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const CLine3D &lhs, const CLine3D &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const CLine3D &lhs, const CLine3D &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const CLine3D &lhs, const CLine3D &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const CLine3D &lhs, const CLine3D &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const CLine3D &lhs, const CLine3D &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  //------

  // intersect
  IntersectType intersect(const CLine3D &line, double *t1, double *t2) {
    double det = (v_.getX()*line.v_.getY() - v_.getY()*line.v_.getX());

    if (::fabs(det) == 0.0)
      return CMathGen::INTERSECT_NONE;

    *t1 = (line.v_.getX()*(p0_.y - line.p0_.y) -
           line.v_.getY()*(p0_.x - line.p0_.x))/det;
    *t2 = (     v_.getX()*(p0_.y - line.p0_.y) -
                v_.getY()*(p0_.x - line.p0_.x))/det;

    if (*t1 >= 0.0 && *t1 <= 1.0 &&
        *t2 >= 0.0 && *t2 <= 1.0)
      return CMathGen::INTERSECT_INSIDE;
    else
      return CMathGen::INTERSECT_OUTSIDE;
  }

  IntersectType intersect(const CLine3D &line, CPoint3D &point) {
    double t1, t2;

    IntersectType type = intersect(line, &t1, &t2);

    if (type != CMathGen::INTERSECT_NONE) {
      point.x = p0_.x + v_.getX()*t1;
      point.y = p0_.y + v_.getY()*t1;
    }

    return type;
  }

  //------

  // length
  double lengthSqr() const {
    return v_.lengthSqr();
  }

  double length() const {
    return ::sqrt(lengthSqr());
  }

  //------

  CPoint3D interp(double t) const {
    return p1_ + t*v_;
  }

 private:
  CPoint3D  p0_;
  CPoint3D  p1_;
  CVector3D v_;
};

#endif
