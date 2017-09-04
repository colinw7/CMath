#ifndef CORTHNORMAL_BASIS_3D_H
#define CORTHNORMAL_BASIS_3D_H

#include <CVector3D.h>
#include <CMatrix3D.h>
#include <CMathGen.h>
#include <CMathMacros.h>

class COrthonormalBasis3D {
 public:
  // create default
  COrthonormalBasis3D() :
   u_(1,0,0), v_(0,1,0), w_(0,0,1), m_(u_,v_,w_) {
  }

  // create from three orthogonal vectors
  COrthonormalBasis3D(const CVector3D &u, const CVector3D &v, const CVector3D &w) :
   u_(u), v_(v), w_(w), m_(u_, v_, w_) {
    if (! validate()) { std::cerr << "Invalid Basis" << std::endl; assert(false); }
  }

  // get vector components
  const CVector3D &getU() const { return u_; }
  const CVector3D &getV() const { return v_; }
  const CVector3D &getW() const { return w_; }

  // get matrix
  const CMatrix3D &getMatrix() const { return m_; }

  // set to default
  void reset() {
    setUVW(CVector3D(1,0,0),CVector3D(0,1,0),CVector3D(0,0,1));
  }

  // set to three orthogonal vectors
  void setUVW(const CVector3D &u, const CVector3D &v, const CVector3D &w) {
    if (validate(u, v, w)) {
      u_ = u; v_ = v; w_ = w;

      m_ = CMatrix3D(u_, v_, w_);
    }
    else {
      std::cerr << "Invalid Basis" << std::endl; assert(false);
    }
  }

  // get vector components
  void getUVW(CVector3D &u, CVector3D &v, CVector3D &w) {
    u = u_; v = v_; w = w_;
  }

  // initialize from u,v othogonal vectors (generate w)
  void initFromUV(const CVector3D &u, const CVector3D &v) {
    CVector3D u1 = u.unit();
    CVector3D w1 = u.crossProduct(v).unit();
    CVector3D v1 = w1.crossProduct(u1);

    setUVW(u1, v1, w1);
  }

  // initialize from v,w othogonal vectors (generate u)
  void initFromVW(const CVector3D &v, const CVector3D &w) {
    CVector3D v1 = v.unit();
    CVector3D u1 = v.crossProduct(w).unit();
    CVector3D w1 = u1.crossProduct(v1);

    setUVW(u1, v1, w1);
  }

  // initialize from w,u othogonal vectors (generate v)
  void initFromWU(const CVector3D &w, const CVector3D &u) {
    CVector3D w1 = w.unit();
    CVector3D v1 = w.crossProduct(u).unit();
    CVector3D u1 = v1.crossProduct(w1);

    setUVW(u1, v1, w1);
  }

  // initialize from w,v othogonal vectors (generate u)
  void initFromWV(const CVector3D &w, const CVector3D &v) {
    CVector3D w1 = w.unit();
    CVector3D u1 = v.crossProduct(w).unit();
    CVector3D v1 = w1.crossProduct(u1);

    setUVW(u1, v1, w1);
  }

  // initialize from v,u othogonal vectors (generate w)
  void initFromVU(const CVector3D &v, const CVector3D &u) {
    CVector3D v1 = v.unit();
    CVector3D w1 = u.crossProduct(v).unit();
    CVector3D u1 = v1.crossProduct(w1);

    setUVW(u1, v1, w1);
  }

  // initialize from u,w othogonal vectors (generate v)
  void initFromUW(const CVector3D &u, const CVector3D &w) {
    CVector3D u1 = u.unit();
    CVector3D v1 = w.crossProduct(u).unit();
    CVector3D w1 = u1.crossProduct(v1);

    setUVW(u1, v1, w1);
  }

  // initialize from u othogonal vector (generate v and w)
  void initFromU(const CVector3D &u) {
    CVector3D u1 = u.unit();
    CVector3D v1 = CVector3D(0,1,0);

    double vu = v1.dotProduct(u1);

    v1 = (v1 - vu*u1).unit();

    if (v1.isZero()) {
      v1 = CVector3D(0,0,1);

      vu = v1.dotProduct(u1);

      v1 = (v1 - vu*u1).unit();
    }

    CVector3D w1 = u1.crossProduct(v1);
    v1 = w1.crossProduct(u1);

    setUVW(u1, v1, w1);
  }

  // initialize from v othogonal vector (generate u and w)
  void initFromV(const CVector3D &v) {
    CVector3D v1 = v.unit();
    CVector3D w1 = CVector3D(0,0,1);

    double wv = w1.dotProduct(v1);

    w1 = (w1 - wv*v1).unit();

    if (w1.isZero()) {
      w1 = CVector3D(1,0,0);

      wv = w1.dotProduct(v1);

      w1 = (w1 - wv*v1).unit();
    }

    CVector3D u1 = v1.crossProduct(w1);
    w1 = u1.crossProduct(v1);

    setUVW(u1, v1, w1);
  }

  // initialize from w othogonal vector (generate u and v)
  void initFromW(const CVector3D &w) {
    CVector3D w1 = w.unit();
    CVector3D u1 = CVector3D(1,0,0);

    double uw = u1.dotProduct(w1);

    u1 = (u1 - uw*w1).unit();

    if (u1.isZero()) {
      u1 = CVector3D(0,1,0);

      uw = u1.dotProduct(w1);

      u1 = (u1 - uw*w1).unit();
    }

    CVector3D v1 = w1.crossProduct(u1);
    u1 = v1.crossProduct(w1);

    setUVW(u1, v1, w1);
  }

  // convert vector 'a' in world coordinates to our local coordinates
  CVector3D fromBasis(const CVector3D &a) const {
    return CVector3D(a.dotProduct(u_), a.dotProduct(v_), a.dotProduct(w_));
  }

  // convert point 'a' in world coordinates to our local coordinates
  CPoint3D fromBasis(const CPoint3D &p) const {
    return fromBasis(CVector3D(p.x, p.y, p.z)).point();
  }

  // convert point 'x, y, z' in world coordinates to our local coordinates
  CVector3D fromBasis(double x, double y, double z) const {
    return fromBasis(CVector3D(x, y, z));
  }

  // convert vector 'a' in local coordinates to work coordinates
  CVector3D toBasis(const CVector3D &a) const {
    return CVector3D(a.getX()*u_ + a.getY()*v_ + a.getZ()*w_);
  }

  // convert point 'a' in local coordinates to work coordinates
  CPoint3D toBasis(const CPoint3D &p) const {
    return toBasis(CVector3D(p.x, p.y, p.z)).point();
  }

  // convert point 'x, y, z' in local coordinates to work coordinates
  CVector3D toBasis(double x, double y, double z) const {
    return toBasis(CVector3D(x, y, z));
  }

  friend bool operator==(const COrthonormalBasis3D &lhs, const COrthonormalBasis3D &rhs) {
    return (lhs.u_ == rhs.u_ && lhs.v_ == rhs.v_ && lhs.w_ == rhs.w_);
  }

  // Right-Handed Rotation about X axis
  void rotateAboutX(double a) {
    CMatrix3D r;

    r.setRotation(CMathGen::X_AXIS_3D, a);

    m_ *= r;

    m_.getColumns(u_, v_, w_);

    u_.normalize(); v_.normalize(); w_.normalize();
  }

  // Right-Handed Rotation about Y axis
  void rotateAboutY(double a) {
    CMatrix3D r;

    r.setRotation(CMathGen::Y_AXIS_3D, a);

    m_ *= r;

    m_.getColumns(u_, v_, w_);

    u_.normalize(); v_.normalize(); w_.normalize();
  }

  // Right-Handed Rotation about Z axis
  void rotateAboutZ(double a) {
    CMatrix3D r;

    r.setRotation(CMathGen::Z_AXIS_3D, a);

    m_ *= r;

    m_.getColumns(u_, v_, w_);

    u_.normalize(); v_.normalize(); w_.normalize();
  }

  // rotate the basis about the unit axis u by theta (radians)
  void rotate(double theta, const CVector3D &u) {
    CMatrix3D r;

    r.setRotation(theta, u);

    m_ *= r;
  }

  // rotate, length of da is theta, unit direction of da is u
  void rotate(const CVector3D &v) {
    double theta = sqrt(v.dotProduct(v)); //angle to rotate by

    rotate(theta, v/theta); //unit vector is axis
  }

  //-----------

  // validate vectors are orthogonal
  bool validate() const {
    return validate(u_, v_, w_);
  }

  // validate vectors are orthogonal
  static bool validate(const CVector3D &u, const CVector3D &v, const CVector3D &w) {
    if (u.isZero() || v.isZero() || w.isZero())
      return false;

    double uv = u.dotProduct(v);
    double vw = v.dotProduct(w);
    double wu = w.dotProduct(u);

    if (! REAL_EQ(uv, 0) || ! REAL_EQ(vw, 0) || ! REAL_EQ(wu, 0))
      return false;

    CVector3D uvw = u.crossProduct(v).crossProduct(w);

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
  friend std::ostream &operator<<(std::ostream &os, const COrthonormalBasis3D &basis) {
    basis.print(os);

    return os;
  }

 private:
  /// Basis Vectors
  CVector3D u_, v_, w_;

  /// CMatrix3D containing Basis Vectors as columns
  CMatrix3D m_;
};

#endif
