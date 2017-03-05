#ifndef CORTHNORMAL_BASIS_2D_H
#define CORTHNORMAL_BASIS_2D_H

#include <CVector2D.h>
#include <CMatrix2D.h>

class COrthonormalBasis2D {
 public:
  COrthonormalBasis2D() :
   u_(1,0), v_(0,1), m_(u_,v_) {
  }

  COrthonormalBasis2D(const CVector2D &u, const CVector2D &v) :
   u_(u), v_(v), m_(u_, v_) {
    if (! validate()) { std::cerr << "Invalid Basis" << std::endl; assert(false); }
  }

  const CVector2D &getU() const { return u_; }
  const CVector2D &getV() const { return v_; }

  const CMatrix2D &getMatrix() const { return m_; }

  void reset() {
    setUV(CVector2D(1,0),CVector2D(0,1));
  }

  void setUV(const CVector2D &u, const CVector2D &v) {
    if (validate(u, v)) {
      u_ = u; v_ = v;

      m_ = CMatrix2D(u_, v_);
    }
    else {
      std::cerr << "Invalid Basis" << std::endl; assert(false);
    }
  }

  void getUV(CVector2D &u, CVector2D &v) {
    u = u_; v = v_;
  }

  void initFromU(const CVector2D &u) {
    CVector2D u1 = u.unit();

    CVector2D v1(u1.getY(), -u1.getX());

    setUV(u1, v1);
  }

  void initFromV(const CVector2D &v) {
    CVector2D v1 = v.unit();

    CVector2D u1(v1.getY(), -v1.getX());

    setUV(u1, v1);
  }

  CVector2D fromBasis(const CVector2D &a) const {
    return CVector2D(a.dotProduct(u_), a.dotProduct(v_));
  }

  CPoint2D fromBasis(const CPoint2D &p) const {
    return fromBasis(CVector2D(p.x, p.y)).point();
  }

  CVector2D fromBasis(double x, double y) const {
    return fromBasis(CVector2D(x, y));
  }

  CVector2D toBasis(const CVector2D &a) const {
    return CVector2D(a.getX()*u_ + a.getY()*v_);
  }

  CPoint2D toBasis(const CPoint2D &p) const {
    return toBasis(CVector2D(p.x, p.y)).point();
  }

  CVector2D toBasis(double x, double y) const {
    return toBasis(CVector2D(x, y));
  }

  friend bool operator==(const COrthonormalBasis2D &lhs, const COrthonormalBasis2D &rhs) {
    return (lhs.u_ == rhs.u_ && lhs.v_ == rhs.v_);
  }

  // Right-Handed Rotations
  void rotate(double a) {
    CMatrix2D r;

    r.setRotation(a);

    m_ *= r;

    m_.getColumns(u_, v_);

    u_.normalize(); v_.normalize();

    m_.setColumns(u_, v_);
  }

  //-----------

  bool validate() const {
    return validate(u_, v_);
  }

  static bool validate(const CVector2D &u, const CVector2D &v) {
    if (u.isZero() || v.isZero())
      return false;

    double uv = u.dotProduct(v);

    if (! REAL_EQ(uv, 0))
      return false;

    return true;
  }

  //-----------

  void print(std::ostream &os) const {
    os << "(" << u_ << ", " << v_ << ") : [" << m_ << "]";
  }

  friend std::ostream &operator<<(std::ostream &os, const COrthonormalBasis2D &basis) {
    basis.print(os);

    return os;
  }

 private:
  /// Basis Vectors
  CVector2D u_, v_;

  /// CMatrix2D containing Basis Vectors as columns
  CMatrix2D m_;
};

#endif
