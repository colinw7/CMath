#ifndef CCOORD_FRAME_2D_H
#define CCOORD_FRAME_2D_H

#include <CVector2D.h>
#include <COrthonormalBasis2D.h>

/// Coordinate Frame
class CCoordFrame2D {
 public:
  CCoordFrame2D() :
   origin_(0,0), basis_() {
  }

  explicit
  CCoordFrame2D(const CPoint2D &origin, const CVector2D &u=CVector2D(1,0),
                const CVector2D &v=CVector2D(0,1)) :
   origin_(origin), basis_(u, v) {
  }

  explicit
  CCoordFrame2D(const CVector2D &origin, const CVector2D &u=CVector2D(1,0),
                const CVector2D &v=CVector2D(0,1)) :
   origin_(origin), basis_(u, v) {
  }

  CCoordFrame2D(const CVector2D &origin, const COrthonormalBasis2D &basis) :
   origin_(origin), basis_(basis) {
  }

  CCoordFrame2D(const CCoordFrame2D &coord_frame) :
   origin_(coord_frame.origin_), basis_ (coord_frame.basis_ ) {
  }

  const CCoordFrame2D &operator=(const CCoordFrame2D &coord_frame) {
    origin_ = coord_frame.origin_;
    basis_  = coord_frame.basis_ ;

    return *this;
  }

  const CVector2D &getOrigin() const { return origin_; }

  CPoint2D getOriginPoint() const { return origin_.point(); }

  const COrthonormalBasis2D &getBasis() const { return basis_; }

  void setOrigin(const CPoint2D  &origin) { origin_ = origin; }
  void setOrigin(const CVector2D &origin) { origin_ = origin; }

  void setBasis(const CVector2D &u, const CVector2D &v) {
    if (u.isZero() || v.isZero()) {
      std::cerr << "Invalid Basis" << std::endl;
      return;
    }

    basis_.setUV(u, v);
  }

  void getBasis(CVector2D &u, CVector2D &v) {
    basis_.getUV(u, v);
  }

  const CMatrix2D &getMatrix() const {
    return basis_.getMatrix();
  }

  // Move Origin by dv
  void move(const CVector2D &dv) {
    origin_ += dv;
  }

  // Move Origin by dx
  void moveX(double dx) {
    origin_.setX(origin_.getX() + dx);
  }

  // Move Origin by dy
  void moveY(double dy) {
    origin_.setY(origin_.getY() + dy);
  }

  // Scale Origin by dx
  void scaleX(double dx) {
    origin_.setX(origin_.getX()*dx);
  }

  // Scale Origin by dy
  void scaleY(double dy) {
    origin_.setY(origin_.getY()*dy);
  }

  // Rotate by angle a
  void rotate(double a) {
    basis_.rotate(a);
  }

  /// Reset to (0,0) origin and (1,0),(0,1) basis
  void reset() {
    origin_.zero();

    basis_.reset();
  }

  /// transform the coordinate vector and translate by this origin
  CPoint2D transformFrom(const CPoint2D &p) const {
    return basis_.fromBasis(p) + origin_;
  }

  /// transform the coordinate vector and translate by this origin
  CVector2D transformFrom(const CVector2D &p) const {
    return basis_.fromBasis(p) + origin_;
  }

  /// translate to this frame's origin, then project onto this basis
  CPoint2D transformTo(const CPoint2D &p) const {
    return basis_.toBasis(p - origin_);
  }

  /// translate to this frame's origin, then project onto this basis
  CVector2D transformTo(const CVector2D &p) const {
    return basis_.toBasis(p - origin_);
  }

  //-----------

  void print(std::ostream &os) const {
    os << "(" << basis_ << "), (" << origin_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CCoordFrame2D &frame) {
    frame.print(os);

    return os;
  }

 private:
  CVector2D           origin_;
  COrthonormalBasis2D basis_;
};

#endif
