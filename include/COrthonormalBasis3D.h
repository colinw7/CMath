#ifndef CORTHNORMAL_BASIS_3D_H
#define CORTHNORMAL_BASIS_3D_H

#include <CVector3D.h>
#include <CMatrix3D.h>
#include <CMathMacros.h>

template<typename T>
class COrthonormalBasis3DT {
 private:
  typedef CVector3DT<T> Vector;
  typedef CPoint3DT<T>  Point;
  typedef CMatrix3DT<T> Matrix;

  /// Basis Vectors
  Vector u_, v_, w_;

  /// Matrix containing Basis Vectors as columns
  Matrix m_;

 public:
  // create default
  COrthonormalBasis3DT() :
   u_(1,0,0), v_(0,1,0), w_(0,0,1), m_(u_,v_,w_) {
  }

  // create from three orthogonal vectors
  COrthonormalBasis3DT(const Vector &u, const Vector &v, const Vector &w) :
   u_(u), v_(v), w_(w), m_(u_, v_, w_) {
    if (! validate()) { std::cerr << "Invalid Basis" << std::endl; assert(false); }
  }

  // get vector components
  const Vector &getU() const { return u_; }
  const Vector &getV() const { return v_; }
  const Vector &getW() const { return w_; }

  // get matrix
  const Matrix &getMatrix() const { return m_; }

  // set to default
  void reset() {
    setUVW(Vector(1,0,0),Vector(0,1,0),Vector(0,0,1));
  }

  // set to three orthogonal vectors
  void setUVW(const Vector &u, const Vector &v, const Vector &w) {
    if (validate(u, v, w)) {
      u_ = u; v_ = v; w_ = w;

      m_ = Matrix(u_, v_, w_);
    }
    else {
      std::cerr << "Invalid Basis" << std::endl; assert(false);
    }
  }

  // get vector components
  void getUVW(Vector &u, Vector &v, Vector &w) {
    u = u_; v = v_; w = w_;
  }

  // initialize from u,v othogonal vectors (generate w)
  void initFromUV(const Vector &u, const Vector &v) {
    Vector u1 = u.unit();
    Vector w1 = u.crossProduct(v).unit();
    Vector v1 = w1.crossProduct(u1);

    setUVW(u1, v1, w1);
  }

  // initialize from v,w othogonal vectors (generate u)
  void initFromVW(const Vector &v, const Vector &w) {
    Vector v1 = v.unit();
    Vector u1 = v.crossProduct(w).unit();
    Vector w1 = u1.crossProduct(v1);

    setUVW(u1, v1, w1);
  }

  // initialize from w,u othogonal vectors (generate v)
  void initFromWU(const Vector &w, const Vector &u) {
    Vector w1 = w.unit();
    Vector v1 = w.crossProduct(u).unit();
    Vector u1 = v1.crossProduct(w1);

    setUVW(u1, v1, w1);
  }

  // initialize from w,v othogonal vectors (generate u)
  void initFromWV(const Vector &w, const Vector &v) {
    Vector w1 = w.unit();
    Vector u1 = v.crossProduct(w).unit();
    Vector v1 = w1.crossProduct(u1);

    setUVW(u1, v1, w1);
  }

  // initialize from v,u othogonal vectors (generate w)
  void initFromVU(const Vector &v, const Vector &u) {
    Vector v1 = v.unit();
    Vector w1 = u.crossProduct(v).unit();
    Vector u1 = v1.crossProduct(w1);

    setUVW(u1, v1, w1);
  }

  // initialize from u,w othogonal vectors (generate v)
  void initFromUW(const Vector &u, const Vector &w) {
    Vector u1 = u.unit();
    Vector v1 = w.crossProduct(u).unit();
    Vector w1 = u1.crossProduct(v1);

    setUVW(u1, v1, w1);
  }

  // initialize from u othogonal vector (generate v and w)
  void initFromU(const Vector &u) {
    Vector u1 = u.unit();
    Vector v1 = Vector(0,1,0);

    T vu = v1.dotProduct(u1);

    v1 = (v1 - vu*u1).unit();

    if (v1.isZero()) {
      v1 = Vector(0,0,1);

      vu = v1.dotProduct(u1);

      v1 = (v1 - vu*u1).unit();
    }

    Vector w1 = u1.crossProduct(v1);
    v1 = w1.crossProduct(u1);

    setUVW(u1, v1, w1);
  }

  // initialize from v othogonal vector (generate u and w)
  void initFromV(const Vector &v) {
    Vector v1 = v.unit();
    Vector w1 = Vector(0,0,1);

    T wv = w1.dotProduct(v1);

    w1 = (w1 - wv*v1).unit();

    if (w1.isZero()) {
      w1 = Vector(1,0,0);

      wv = w1.dotProduct(v1);

      w1 = (w1 - wv*v1).unit();
    }

    Vector u1 = v1.crossProduct(w1);
    w1 = u1.crossProduct(v1);

    setUVW(u1, v1, w1);
  }

  // initialize from w othogonal vector (generate u and v)
  void initFromW(const Vector &w) {
    Vector w1 = w.unit();
    Vector u1 = Vector(1,0,0);

    T uw = u1.dotProduct(w1);

    u1 = (u1 - uw*w1).unit();

    if (u1.isZero()) {
      u1 = Vector(0,1,0);

      uw = u1.dotProduct(w1);

      u1 = (u1 - uw*w1).unit();
    }

    Vector v1 = w1.crossProduct(u1);
    u1 = v1.crossProduct(w1);

    setUVW(u1, v1, w1);
  }

  // convert vector 'a' in world coordinates to our local coordinates
  Vector fromBasis(const Vector &a) const {
    return Vector(a.dotProduct(u_), a.dotProduct(v_), a.dotProduct(w_));
  }

  // convert point 'a' in world coordinates to our local coordinates
  Point fromBasis(const Point &p) const {
    return fromBasis(Vector(p.x, p.y, p.z)).point();
  }

  // convert point 'x, y, z' in world coordinates to our local coordinates
  Vector fromBasis(T x, T y, T z) const {
    return fromBasis(Vector(x, y, z));
  }

  // convert vector 'a' in local coordinates to work coordinates
  Vector toBasis(const Vector &a) const {
    return Vector(a.getX()*u_ + a.getY()*v_ + a.getZ()*w_);
  }

  // convert point 'a' in local coordinates to work coordinates
  Point toBasis(const Point &p) const {
    return toBasis(Vector(p.x, p.y, p.z)).point();
  }

  // convert point 'x, y, z' in local coordinates to work coordinates
  Vector toBasis(T x, T y, T z) const {
    return toBasis(Vector(x, y, z));
  }

  friend bool operator==(const COrthonormalBasis3DT &lhs, const COrthonormalBasis3DT &rhs) {
    return (lhs.u_ == rhs.u_ && lhs.v_ == rhs.v_ && lhs.w_ == rhs.w_);
  }

  // Right-Handed Rotation about X axis
  void rotateAboutX(T a) {
    CMatrix3D r;

    r.setRotation(CMathGen::X_AXIS_3D, a);

    m_ *= r;

    m_.getColumns(u_, v_, w_);

    u_.normalize(); v_.normalize(); w_.normalize();
  }

  // Right-Handed Rotation about Y axis
  void rotateAboutY(T a) {
    CMatrix3D r;

    r.setRotation(CMathGen::Y_AXIS_3D, a);

    m_ *= r;

    m_.getColumns(u_, v_, w_);

    u_.normalize(); v_.normalize(); w_.normalize();
  }

  // Right-Handed Rotation about Z axis
  void rotateAboutZ(T a) {
    CMatrix3D r;

    r.setRotation(CMathGen::Z_AXIS_3D, a);

    m_ *= r;

    m_.getColumns(u_, v_, w_);

    u_.normalize(); v_.normalize(); w_.normalize();
  }

  // rotate the basis about the unit axis u by theta (radians)
  void rotate(T theta, const Vector &u) {
    CMatrix3D r;

    r.setRotation(theta, u);

    m_ *= r;
  }

  // rotate, length of da is theta, unit direction of da is u
  void rotate(const Vector &v) {
    T theta = sqrt(v.dotProduct(v)); //angle to rotate by

    rotate(theta, v/theta); //unit vector is axis
  }

  //-----------

  // validate vectors are orthogonal
  bool validate() const {
    return validate(u_, v_, w_);
  }

  // validate vectors are orthogonal
  static bool validate(const Vector &u, const Vector &v, const Vector &w) {
    if (u.isZero() || v.isZero() || w.isZero())
      return false;

    T uv = u.dotProduct(v);
    T vw = v.dotProduct(w);
    T wu = w.dotProduct(u);

    if (! REAL_EQ(uv, 0) || ! REAL_EQ(vw, 0) || ! REAL_EQ(wu, 0))
      return false;

    Vector uvw = u.crossProduct(v).crossProduct(w);

    if (! REAL_EQ(uvw.length(), 0))
      return false;

    return true;
  }

  //-----------

  // display
  void print(std::ostream &os) const {
    os << "(" << u_ << ", " << v_ << ", " << w_ << ") : [" << m_ << "]";
  }

  // support operator <<
  friend std::ostream &operator<<(std::ostream &os, const COrthonormalBasis3DT &basis) {
    basis.print(os);

    return os;
  }
};

typedef COrthonormalBasis3DT<double> COrthonormalBasis3D;

#endif
