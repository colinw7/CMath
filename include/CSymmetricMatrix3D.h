#ifndef CSYMMETRIC_MATRIX_3D_H
#define CSYMMETRIC_MATRIX_3D_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CVector3D.h>

template<typename T>
class CSymmetricMatrix3DT {
 private:
  union {
    T m1_[16];
    T m2_[4][4];

    struct {
      T m00_, m01_, m02_, m03_;
      T m10_, m11_, m12_, m13_;
      T m20_, m21_, m22_, m23_;
      T m30_, m31_, m32_, m33_;
    };
  };

 public:
  CSymmetricMatrix3DT() {
  }

  CSymmetricMatrix3DT(T a, T b, T c, T d, T e, T f) {
    setValues(a, b, c, d, e, f, 0.0, 0.0, 0.0);
  }

  CSymmetricMatrix3DT(T a, T b, T c, T d, T e, T f, T tx, T ty, T tz) {
    setValues(a, b, c, d, e, f, tx, ty, tz);
  }

  CSymmetricMatrix3DT(const CSymmetricMatrix3DT &a) {
    memcpy(m1_, a.m1_, sizeof(m1_));
  }

 ~CSymmetricMatrix3DT() { }

  void setIdentity() {
    setValues(1.0, 0.0, 0.0,
                   1.0, 0.0,
                        1.0 ,
              0.0, 0.0, 0.0);
  }

  void setTranslation(T tx, T ty, T tz) {
    setValues(1.0, 0.0, 0.0,
                   1.0, 0.0,
                        1.0,
               tx,  ty, tz );
  }

  void setScale(T sx, T sy, T sz) {
    setValues(sx , 0.0, 0.0,
                    sy, 0.0,
                         sz,
              0.0, 0.0, 0.0);
  }

  void multiplyPoint(T xi, T yi, T zi, T *xo, T *yo, T *zo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi + m03_;
    *yo = m10_*xi + m11_*yi + m12_*zi + m13_;
    *zo = m20_*xi + m21_*yi + m22_*zi + m23_;
  }

  void multiplyPoint(const CPoint3DT<T> &ipoint,
                     CPoint3DT<T> &opoint) const {
    opoint.x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z + m03_;
    opoint.y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z + m13_;
    opoint.z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z + m23_;
  }

  void multiplyVector(const CVector3DT<T> &ivector,
                      CVector3DT<T> &ovector) const {
    ovector.setX(m00_*ivector.getX() +
                 m01_*ivector.getY() +
                 m02_*ivector.getZ() +
                 m03_);
    ovector.setY(m10_*ivector.getX() +
                 m11_*ivector.getY() +
                 m12_*ivector.getZ() +
                 m13_);
    ovector.setZ(m20_*ivector.getX() +
                 m21_*ivector.getY() +
                 m22_*ivector.getZ() +
                 m23_);
  }

  void setValues(T a, T b, T c, T d, T e, T f) {
    m00_ =   a; m01_ =   b; m02_ =   c; m03_ = 0.0;
    m10_ =   b; m11_ =   d; m12_ =   e; m13_ = 0.0;
    m20_ =   c; m21_ =   e; m22_ =   f; m23_ = 0.0;
    m30_ = 0.0; m31_ = 0.0; m32_ = 0.0; m33_ = 1.0;
  }

  void setValues(T a, T b, T c, T d, T e, T f, T tx, T ty, T tz) {
    m00_ = a  ; m01_ = b  ; m02_ = c  ; m03_ = tx;
    m10_ = b  ; m11_ = d  ; m12_ = e  ; m13_ = ty;
    m20_ = c  ; m21_ = e  ; m22_ = f  ; m23_ = tz;
    m30_ = 0.0; m31_ = 0.0; m32_ = 0.0; m33_ = 1.0;
  }

  void getValues(T *a, T *b, T *c, T *d, T *e, T *f,
                 T *tx, T *ty, T *tz) const {
    if (a ) *a  = m00_;
    if (b ) *b  = m01_;
    if (c ) *c  = m02_;
    if (d ) *d  = m11_;
    if (e ) *e  = m12_;
    if (f ) *f  = m13_;
    if (tx) *tx = m03_;
    if (ty) *ty = m13_;
    if (tz) *tz = m23_;
  }

  void getValues(T *a, T *b, T *c, T *d, T *e, T *f) const {
    if (a) *a  = m00_;
    if (b) *b  = m01_;
    if (c) *c  = m02_;
    if (d) *d  = m11_;
    if (e) *e  = m12_;
    if (f) *f  = m13_;
  }

  void getValues(T v[9]) const {
    v[0] = m00_; v[1] = m01_; v[2] = m02_;
    v[3] = m10_; v[4] = m11_; v[5] = m12_;
    v[6] = m20_; v[7] = m21_; v[8] = m22_;
  }

  CSymmetricMatrix3DT &operator=(const CSymmetricMatrix3DT &a) {
    memcpy(m1_, a.m1_, 16*sizeof(T));

    return *this;
  }

  CSymmetricMatrix3DT &operator*=(T s) {
    CSymmetricMatrix3DT a = *this;

    m00_ = a.m00_*s; m01_ = a.m00_*s; m02_ = a.m00_*s; m03_ = a.m00_*s;
    m10_ = a.m10_*s; m11_ = a.m10_*s; m12_ = a.m10_*s; m13_ = a.m10_*s;
    m20_ = a.m20_*s; m21_ = a.m20_*s; m22_ = a.m20_*s; m23_ = a.m20_*s;
    m30_ = 0.0     ; m31_ = 0.0     ; m32_ = 0.0     ; m33_ = 1.0     ;

    return *this;
  }

  CSymmetricMatrix3DT operator*(T s) {
    CSymmetricMatrix3DT c = *this;

    c *= s;

    return c;
  }

  void setValue(unsigned int i, T value) { m1_[i] = value; }

  T getValue(unsigned int i) { return m1_[i]; }

  void zero() { memset(m1_, 0, sizeof(m1_)); }

  T operator[](unsigned int i) { return m1_[i]; }

  const T &operator[](unsigned int i) const { return m1_[i]; }

  friend std::ostream &operator<<(std::ostream &os, const CSymmetricMatrix3DT &matrix) {
    os << "(" << matrix.m00_ << "," <<
                 matrix.m01_ << "," <<
                 matrix.m02_ << "," <<
                 matrix.m03_ << ")" << std::endl;
    os << "(" << matrix.m10_ << "," <<
                 matrix.m11_ << "," <<
                 matrix.m12_ << "," <<
                 matrix.m13_ << ")" << std::endl;
    os << "(" << matrix.m20_ << "," <<
                 matrix.m21_ << "," <<
                 matrix.m22_ << "," <<
                 matrix.m23_ << ")" << std::endl;
    os << "(" << matrix.m30_ << "," <<
                 matrix.m31_ << "," <<
                 matrix.m32_ << "," <<
                 matrix.m33_ << ")" << std::endl;

    return os;
  }

 private:
  static T calcDeterminant(T a, T b, T c, T d, T e, T f);
};

typedef CSymmetricMatrix3DT<double> CSymmetricMatrix3D;

#endif
