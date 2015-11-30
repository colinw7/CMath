#ifndef CLINE_3D_H
#define CLINE_3D_H

// TODO: add transform, distance, closest points

#include <CPoint3D.h>
#include <CVector3D.h>
#include <CMathGen.h>

template<typename T>
class CLine3DT {
 private:
  typedef CPoint3DT<T>            Point;
  typedef CVector3DT<T>           Vector;
  typedef CLine3DT<T>             Line;
  typedef CMathGen::IntersectType IntersectType;

 public:
  // constructor/destructor
  CLine3DT() :
   p0_(0,0,0), p1_(1,0,0), v_(1,0,0) {
  }

  CLine3DT(T x1, T y1, T z1, T x2, T y2, T z2) :
   p0_(x1, y1, z1), p1_(x2, y2, z2), v_(x2 - x1, y2 - y1, z2 - z1) {
  }

  CLine3DT(const Point &p0, const Point &p1) :
   p0_(p0), p1_(p1), v_(p0, p1) {
  }

  CLine3DT(const Point &o, const Vector &d) :
   p0_(o), p1_(o + d), v_(d) {
  }

  CLine3DT(const Line &rhs) :
   p0_(rhs.p0_), p1_(rhs.p1_), v_(rhs.v_) {
  }

  CLine3DT &operator=(const Line &rhs) {
    p0_= rhs.p0_; p1_= rhs.p1_; v_= rhs.v_;

    return *this;
  }

 ~CLine3DT() { }

  //------

  // accessors
  Point &start() { return p0_; }
  Point &end  () { return p1_; }

  const Point &start() const { return p0_; }
  const Point &end  () const { return p1_; }

  const Vector &vector() const { return v_; }

  Point point(double mu) const { return p0_ + mu*v_; }

  Point center() const { return (p0_ + p1_)/2.0; }

  void setStart(const Point &p0) {
    p0_ = p0;
    v_  = Vector(p0_, p1_);
  }

  void setEnd(const Point &p1) {
    p1_ = p1;
    v_  = Vector(p0_, p1_);
  }

  void setVector(const Vector &v) {
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

  friend std::ostream &operator<<(std::ostream &os, const Line &line) {
    line.print(os);

    return os;
  }

  //------

  // comparison
  int cmp(const Line &v) const {
    if      (p0_ < v.p0_) return -1;
    else if (p0_ > v.p0_) return  1;
    else if (p1_ < v.p1_) return -1;
    else if (p1_ > v.p1_) return  1;
    else                  return  0;
  }

  friend bool operator==(const Line &lhs, const Line &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const Line &lhs, const Line &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const Line &lhs, const Line &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const Line &lhs, const Line &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const Line &lhs, const Line &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const Line &lhs, const Line &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  //------

  // intersect
  IntersectType intersect(const Line &line, T *t1, T *t2) {
    T det = (v_.x_*line.v_.y_ - v_.y_*line.v_.x_);

    if (::fabs(det) == 0.0)
      return CMathGen::INTERSECT_NONE;

    *t1 = (line.v_.x_*(p0_.y - line.p0_.y) -
           line.v_.y_*(p0_.x - line.p0_.x))/det;
    *t2 = (     v_.x_*(p0_.y - line.p0_.y) -
                v_.y_*(p0_.x - line.p0_.x))/det;

    if (*t1 >= 0.0 && *t1 <= 1.0 &&
        *t2 >= 0.0 && *t2 <= 1.0)
      return CMathGen::INTERSECT_INSIDE;
    else
      return CMathGen::INTERSECT_OUTSIDE;
  }

  IntersectType intersect(const Line &line, Point &point) {
    T t1, t2;

    IntersectType type = intersect(line, &t1, &t2);

    if (type != CMathGen::INTERSECT_NONE) {
      point.x = p0_.x + v_.x_*t1;
      point.y = p0_.y + v_.y_*t1;
    }

    return type;
  }

  //------

  // length
  T lengthSqr() const {
    return v_.lengthSqr();
  }

  T length() const {
    return ::sqrt(lengthSqr());
  }

 private:
  Point  p0_;
  Point  p1_;
  Vector v_;
};

typedef CLine3DT<double> CLine3D;

#endif
