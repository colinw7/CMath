#ifndef CREF_TRANSFORM_3D_H
#define CREF_TRANSFORM_3D_H

#include <CRefPtr.h>
#include <CMatrix3D.h>
#include <CPoint3D.h>
#include <CVector3D.h>

template<typename T>
class CRefTransform3DT {
 private:
  typedef CMatrix3DT<T>   Matrix;
  typedef CVector3DT<T>   Vector;
  typedef CPoint3DT<T>    Point;
  typedef CNormal3DT<T>   Normal;
  typedef CRefPtr<Matrix> MatrixPtr;

  MatrixPtr m_;
  MatrixPtr im_;
  bool      inverse_set_;

 public:
  CRefTransform3DT() {
    m_ = new Matrix(Matrix::CMATRIX_3D_IDENTITY);

    im_ = m_;

    inverse_set_ = true;
  }

  CRefTransform3DT(const CRefTransform3DT &t) {
    m_           = t.m_;
    inverse_set_ = t.inverse_set_;

    if (inverse_set_)
      im_ = t.im_;
  }

  CRefTransform3DT(const Matrix &m) {
    m_ = new Matrix(m);

    inverse_set_ = false;
  }

  CRefTransform3DT(const Matrix &m, const Matrix &im) {
    m_  = new Matrix(m);
    im_ = new Matrix(im);

    inverse_set_ = true;
  }

  CRefTransform3DT(const MatrixPtr &m) {
    m_ = m;

    inverse_set_ = false;
  }

  CRefTransform3DT(const MatrixPtr &m, const MatrixPtr &im) {
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

  CRefTransform3DT getInverse() {
    return CRefTransform3DT(im_, m_);
  }

  CRefTransform3DT translation(const Vector &v) {
    Matrix m, im;

    m .setTranslation( v.getX(),  v.getY(),  v.getZ());
    im.setTranslation(-v.getX(), -v.getY(), -v.getZ());

    return CRefTransform3DT(m, im);
  }

  CRefTransform3DT scale(T x, T y, T z) {
    Matrix m, im;

    m.setScale(x, y, z);

    im.setScale(1.0/x, 1.0/y, 1.0/z);

    return CRefTransform3DT(m, im);
  }

  CRefTransform3DT rotation(CMathGen::AxisType3D axis, T angle,
                           CMathGen::Handedness handedness =
                            CMathGen::RIGHT_HANDEDNESS) {
    Matrix m, im;

    m.setRotation(axis, angle, handedness);

    im = m.transposed();

    return CRefTransform3DT(m, im);
  }

  CRefTransform3DT genRotation(const Vector &axis, T angle,
                              CMathGen::Handedness handedness =
                               CMathGen::RIGHT_HANDEDNESS) {
    Matrix m, im;

    m.setGenRotation(axis, angle, handedness);

    im = m.transposed();

    return CRefTransform3DT(m, im);
  }

  CRefTransform3DT lookAt(const Point &eye, const Point &center,
                          const Vector &up) {
    Matrix m, im;

    m.lookAt(eye, center, up);

    im = m.transposed();

    return CRefTransform3DT(m, im);
  }

  bool swapsHandedness() const {
    return m_.upperDeterminant() < 0.0;
  }

  CRefTransform3DT operator*(const CRefTransform3DT &rhs) const {
    if (inverse_set_ && rhs.inverse_set_)
      return CRefTransform3DT(m_*rhs.m_, im_*rhs.im_)
    else
      return CRefTransform3DT(m_*rhs.m_);
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
    T ix, iy, iz;

    n.getXYZ(&ix, &iy, &iz);

    calcInverse();

    T ox = im_.m00_*ix + im_.m10_*iy + im_.m20_*iz;
    T oy = im_.m01_*ix + im_.m11_*iy + im_.m21_*iz;
    T oz = im_.m02_*ix + im_.m12_*iy + im_.m22_*iz;

    return Normal(ox, oy, oz);
  }

  void operator()(const Normal &in, Normal &on) {
    T ix, iy, iz;

    in.getXYZ(&ix, &iy, &iz);

    calcInverse();

    T ox = im_.m00_*ix + im_.m10_*iy + im_.m20_*iz;
    T oy = im_.m01_*ix + im_.m11_*iy + im_.m21_*iz;
    T oz = im_.m02_*ix + im_.m12_*iy + im_.m22_*iz;

    on.setXYZ(ox, oy, oz);
  }

 private:
  void calcInverse() {
    if (! inverse_set_) {
      im_ = new Matrix(m_.inverse());

      inverse_set_ = true;
    }
  }
};

typedef CRefTransform3DT<double> CRefTransform3D;

#endif
