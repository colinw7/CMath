#ifndef CCOORD_FRAME_2D_H
#define CCOORD_FRAME_2D_H

#include <CVector2D.h>
#include <COrthonormalBasis2D.h>

/// Coordinate Frame
template<typename T>
class CCoordFrame2DT {
 private:
  typedef CPoint2DT<T>            Point;
  typedef CVector2DT<T>           Vector;
  typedef CMatrix2DT<T>           Matrix;
  typedef COrthonormalBasis2DT<T> OrthonormalBasis;

  Vector           origin_;
  OrthonormalBasis basis_;

 public:
  CCoordFrame2DT() :
   origin_(0,0), basis_() {
  }

  explicit
  CCoordFrame2DT(const Point &origin, const Vector &u=Vector(1,0), const Vector &v=Vector(0,1)) :
   origin_(origin), basis_(u, v) {
  }

  explicit
  CCoordFrame2DT(const Vector &origin, const Vector &u=Vector(1,0), const Vector &v=Vector(0,1)) :
   origin_(origin), basis_(u, v) {
  }

  CCoordFrame2DT(const Vector &origin, const OrthonormalBasis &basis) :
   origin_(origin), basis_(basis) {
  }

  CCoordFrame2DT(const CCoordFrame2DT &coord_frame) :
    origin_(coord_frame.origin_), basis_ (coord_frame.basis_ ) {
  }

  const CCoordFrame2DT &operator=(const CCoordFrame2DT &coord_frame) {
    origin_ = coord_frame.origin_;
    basis_  = coord_frame.basis_ ;

    return *this;
  }

  const Vector &getOrigin() const { return origin_; }

  Point getOriginPoint() const { return origin_.point(); }

  const OrthonormalBasis &getBasis() const { return basis_; }

  void setOrigin(const Point  &origin) { origin_ = origin; }
  void setOrigin(const Vector &origin) { origin_ = origin; }

  void setBasis(const Vector &u, const Vector &v) {
    if (u.isZero() || v.isZero()) {
      std::cerr << "Invalid Basis" << std::endl;
      return;
    }

    basis_.setUV(u, v);
  }

  void getBasis(Vector &u, Vector &v) {
    basis_.getUV(u, v);
  }

  const Matrix &getMatrix() const {
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

  friend std::ostream &operator<<(std::ostream &os, const CCoordFrame2DT &frame) {
    frame.print(os);

    return os;
  }
};

typedef class CCoordFrame2DT<double> CCoordFrame2D;

#endif
