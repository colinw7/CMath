#ifndef CREF_TRANSFORM_3D_H
#define CREF_TRANSFORM_3D_H

#include <CRefPtr.h>
#include <CMatrix3D.h>
#include <CPoint3D.h>
#include <CVector3D.h>
#include <CMathGen.h>

class CRefTransform3D {
 private:
  typedef CMatrix3D       Matrix;
  typedef CVector3D       Vector;
  typedef CPoint3D        Point;
  typedef CNormal3D       Normal;
  typedef CRefPtr<Matrix> MatrixPtr;

 public:
  CRefTransform3D() {
    m_ = new Matrix(Matrix::CMATRIX_3D_IDENTITY);

    im_ = m_;

    inverse_set_ = true;
  }

  CRefTransform3D(const CRefTransform3D &t) {
    m_           = t.m_;
    inverse_set_ = t.inverse_set_;

    if (inverse_set_)
      im_ = t.im_;
  }

  CRefTransform3D(const Matrix &m) {
    m_ = new Matrix(m);

    inverse_set_ = false;
  }

  CRefTransform3D(const Matrix &m, const Matrix &im) {
    m_  = new Matrix(m);
    im_ = new Matrix(im);

    inverse_set_ = true;
  }

  CRefTransform3D(const MatrixPtr &m) {
    m_ = m;

    inverse_set_ = false;
  }

  CRefTransform3D(const MatrixPtr &m, const MatrixPtr &im) {
    m_  = m;
    im_ = im;

    inverse_set_ = true;
  }

  const Matrix &getMatrix() {
    return m_;
  }

  const Matrix &getIMatrix() {
    calcInverse();

    return im_;
  }

  CRefTransform3D getInverse() {
    return CRefTransform3D(im_, m_);
  }

  CRefTransform3D translation(const Vector &v) {
    Matrix m, im;

    m .setTranslation( v.getX(),  v.getY(),  v.getZ());
    im.setTranslation(-v.getX(), -v.getY(), -v.getZ());

    return CRefTransform3D(m, im);
  }

  CRefTransform3D scale(double x, double y, double z) {
    Matrix m, im;

    m.setScale(x, y, z);

    im.setScale(1.0/x, 1.0/y, 1.0/z);

    return CRefTransform3D(m, im);
  }

  CRefTransform3D rotation(CMathGen::AxisType3D axis, double angle,
                           CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    Matrix m, im;

    m.setRotation(axis, angle, handedness);

    im = m.transposed();

    return CRefTransform3D(m, im);
  }

  CRefTransform3D genRotation(const Vector &axis, double angle,
                              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    Matrix m, im;

    m.setGenRotation(axis, angle, handedness);

    im = m.transposed();

    return CRefTransform3D(m, im);
  }

  CRefTransform3D lookAt(const Point &eye, const Point &center,
                          const Vector &up) {
    Matrix m, im;

    m.lookAt(eye, center, up);

    im = m.transposed();

    return CRefTransform3D(m, im);
  }

  bool swapsHandedness() const {
    return m_.upperDeterminant() < 0.0;
  }

  CRefTransform3D operator*(const CRefTransform3D &rhs) const {
    if (inverse_set_ && rhs.inverse_set_)
      return CRefTransform3D(m_*rhs.m_, im_*rhs.im_)
    else
      return CRefTransform3D(m_*rhs.m_);
  }

  Point operator()(const Point &p) {
    return m_*p;
  }

  void operator()(const Point &ip, Point &op) {
    op = m_*ip;
  }

  Vector operator()(const Vector &v) {
    return m_*v;
  }

  void operator()(const Vector &iv, Vector &ov) {
    ov = m_*iv;
  }

  Normal operator()(const Normal &n) {
    double ix, iy, iz;

    n.getXYZ(&ix, &iy, &iz);

    calcInverse();

    double ox = im_.m00_*ix + im_.m10_*iy + im_.m20_*iz;
    double oy = im_.m01_*ix + im_.m11_*iy + im_.m21_*iz;
    double oz = im_.m02_*ix + im_.m12_*iy + im_.m22_*iz;

    return Normal(ox, oy, oz);
  }

  void operator()(const Normal &in, Normal &on) {
    double ix, iy, iz;

    in.getXYZ(&ix, &iy, &iz);

    calcInverse();

    double ox = im_.m00_*ix + im_.m10_*iy + im_.m20_*iz;
    double oy = im_.m01_*ix + im_.m11_*iy + im_.m21_*iz;
    double oz = im_.m02_*ix + im_.m12_*iy + im_.m22_*iz;

    on.setXYZ(ox, oy, oz);
  }

 private:
  void calcInverse() {
    if (! inverse_set_) {
      im_ = new Matrix(m_.inverse());

      inverse_set_ = true;
    }
  }

 private:
  MatrixPtr m_;
  MatrixPtr im_;
  bool      inverse_set_ { false };
};

#endif
