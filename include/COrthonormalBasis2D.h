#ifndef CORTHNORMAL_BASIS_2D_H
#define CORTHNORMAL_BASIS_2D_H

#include <CVector2D.h>
#include <CMatrix2D.h>

template<typename T>
class COrthonormalBasis2DT {
 private:
  typedef CVector2DT<T> Vector;
  typedef CPoint2DT<T>  Point;
  typedef CMatrix2DT<T> Matrix;

  /// Basis Vectors
  Vector u_, v_;

  /// Matrix containing Basis Vectors as columns
  Matrix m_;

 public:
  COrthonormalBasis2DT() :
   u_(1,0), v_(0,1), m_(u_,v_) {
  }

  COrthonormalBasis2DT(const Vector &u, const Vector &v) :
   u_(u), v_(v), m_(u_, v_) {
    if (! validate()) { std::cerr << "Invalid Basis" << std::endl; assert(false); }
  }

  const Vector &getU() const { return u_; }
  const Vector &getV() const { return v_; }

  const Matrix &getMatrix() const { return m_; }

  void reset() {
    setUV(Vector(1,0),Vector(0,1));
  }

  void setUV(const Vector &u, const Vector &v) {
    if (validate(u, v)) {
      u_ = u; v_ = v;

      m_ = Matrix(u_, v_);
    }
    else {
      std::cerr << "Invalid Basis" << std::endl; assert(false);
    }
  }

  void getUV(Vector &u, Vector &v) {
    u = u_; v = v_;
  }

  void initFromU(const Vector &u) {
    Vector u1 = u.unit();

    Vector v1(u1.getY(), -u1.getX());

    setUV(u1, v1);
  }

  void initFromV(const Vector &v) {
    Vector v1 = v.unit();

    Vector u1(v1.getY(), -v1.getX());

    setUV(u1, v1);
  }

  Vector fromBasis(const Vector &a) const {
    return Vector(a.dotProduct(u_), a.dotProduct(v_));
  }

  Point fromBasis(const Point &p) const {
    return fromBasis(Vector(p.x, p.y)).point();
  }

  Vector fromBasis(T x, T y) const {
    return fromBasis(Vector(x, y));
  }

  Vector toBasis(const Vector &a) const {
    return Vector(a.getX()*u_ + a.getY()*v_);
  }

  Point toBasis(const Point &p) const {
    return toBasis(Vector(p.x, p.y)).point();
  }

  Vector toBasis(T x, T y) const {
    return toBasis(Vector(x, y));
  }

  friend bool operator==(const COrthonormalBasis2DT &lhs,
                         const COrthonormalBasis2DT &rhs) {
    return (lhs.u_ == rhs.u_ && lhs.v_ == rhs.v_);
  }

  // Right-Handed Rotations
  void rotate(T a) {
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

  static bool validate(const Vector &u, const Vector &v) {
    if (u.isZero() || v.isZero())
      return false;

    T uv = u.dotProduct(v);

    if (! REAL_EQ(uv, 0))
      return false;

    return true;
  }

  //-----------

  void print(std::ostream &os) const {
    os << "(" << u_ << ", " << v_ << ") : [" << m_ << "]";
  }

  friend std::ostream &operator<<(std::ostream &os, const COrthonormalBasis2DT &basis) {
    basis.print(os);

    return os;
  }
};

typedef COrthonormalBasis2DT<double> COrthonormalBasis2D;

#endif
