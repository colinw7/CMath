#ifndef CCOORD_FRAME_3D_H
#define CCOORD_FRAME_3D_H

#include <CVector3D.h>
#include <COrthonormalBasis3D.h>

/// Coordinate Frame
class CCoordFrame3D {
 public:
  CCoordFrame3D() :
   origin_(0,0,0), basis_() {
  }

  explicit
  CCoordFrame3D(const CPoint3D &origin, const CVector3D &u=CVector3D(1,0,0),
                const CVector3D &v=CVector3D(0,1,0), const CVector3D &w=CVector3D(0,0,1)) :
   origin_(origin), basis_(u, v, w) {
  }

  explicit
  CCoordFrame3D(const CVector3D &origin, const CVector3D &u=CVector3D(1,0,0),
                const CVector3D &v=CVector3D(0,1,0), const CVector3D &w=CVector3D(0,0,1)) :
   origin_(origin), basis_(u, v, w) {
  }

  CCoordFrame3D(const CVector3D &origin, const COrthonormalBasis3D &basis) :
   origin_(origin), basis_(basis) {
  }

  CCoordFrame3D(const CCoordFrame3D &coord_frame) :
    origin_(coord_frame.origin_),
    basis_ (coord_frame.basis_ ) {
  }

  const CCoordFrame3D &operator=(const CCoordFrame3D &coord_frame) {
    origin_ = coord_frame.origin_;
    basis_  = coord_frame.basis_ ;

    return *this;
  }

  void init() {
    origin_ = CPoint3D(0,0,0);
    basis_  = COrthonormalBasis3D(CVector3D(1,0,0), CVector3D(0,1,0), CVector3D(0,0,1));
  }

  const CVector3D &getOrigin() const { return origin_; }

  CPoint3D getOriginPoint() const { return origin_.point(); }

  const COrthonormalBasis3D &getBasis() const { return basis_; }

  void setOrigin(const CPoint3D  &origin) { origin_ = origin; }
  void setOrigin(const CVector3D &origin) { origin_ = origin; }

  void setBasis(const CVector3D &u, const CVector3D &v, const CVector3D &w) {
    if (u.isZero() || v.isZero() || w.isZero()) {
      std::cerr << "Invalid Basis" << std::endl;
      return;
    }

    basis_.setUVW(u, v, w);
  }

  void getBasis(CVector3D &u, CVector3D &v, CVector3D &w) {
    basis_.getUVW(u, v, w);
  }

  const CMatrix3D &getMatrix() const {
    return basis_.getMatrix();
  }

  /// Move Origin by dv
  void move(const CVector3D &dv) {
    origin_ += dv;
  }

  /// Move Origin by dx
  void moveX(double dx) {
    origin_.setX(origin_.getX() + dx);
  }

  /// Move Origin by dy
  void moveY(double dy) {
    origin_.setY(origin_.getY() + dy);
  }

  /// Move Origin by dz
  void moveZ(double dz) {
    origin_.setZ(origin_.getZ() + dz);
  }

  /// Scale Origin by dx
  void scaleX(double dx) {
    origin_.setX(origin_.getX()*dx);
  }

  /// Scale Origin by dy
  void scaleY(double dy) {
    origin_.setY(origin_.getY()*dy);
  }

  /// Scale Origin by dz
  void scaleZ(double dz) {
    origin_.setZ(origin_.getZ()*dz);
  }

  void rotateAboutX(double a) {
    basis_.rotateAboutX(a);
  }

  void rotateAboutY(double a) {
    basis_.rotateAboutY(a);
  }

  void rotateAboutZ(double a) {
    basis_.rotateAboutZ(a);
  }

  void rotateAboutXYZ(double xa, double ya, double za) {
    rotateAboutX(xa);
    rotateAboutY(ya);
    rotateAboutZ(za);
  }

  // rotate, length of da is theta, unit direction of da is u
  void rotate(const CVector3D &v) {
    basis_.rotate(v);
  }

  /// Reset to (0,0,0) origin and (1,0,0),(0,1,0),(0,0,1) basis
  void reset() {
    origin_.zero();

    basis_.reset();
  }

  /// transform the coordinate vector and translate by this origin
  CPoint3D transformFrom(const CPoint3D &p) const {
    return basis_.fromBasis(p) + origin_;
  }

  /// transform the coordinate vector and translate by this origin
  CVector3D transformFrom(const CVector3D &p) const {
    return basis_.fromBasis(p) + origin_;
  }

  /// translate to this frame's origin, then project onto this basis
  CPoint3D transformTo(const CPoint3D &p) const {
    return basis_.toBasis(p - origin_);
  }

  /// translate to this frame's origin, then project onto this basis
  CVector3D transformTo(const CVector3D &p) const {
    return basis_.toBasis(p - origin_);
  }

  //-----------

  void print(std::ostream &os) const {
    os << "(" << basis_ << "), (" << origin_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const CCoordFrame3D &frame) {
    frame.print(os);

    return os;
  }

 private:
  CVector3D           origin_;
  COrthonormalBasis3D basis_;
};

#endif
