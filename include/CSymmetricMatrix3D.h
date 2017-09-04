#ifndef CSYMMETRIC_MATRIX_3D_H
#define CSYMMETRIC_MATRIX_3D_H

#include <CPoint3D.h>
#include <CVector3D.h>

class CSymmetricMatrix3D {
 public:
  CSymmetricMatrix3D() {
  }

  CSymmetricMatrix3D(double a, double b, double c, double d, double e, double f) {
    setValues(a, b, c, d, e, f, 0.0, 0.0, 0.0);
  }

  CSymmetricMatrix3D(double a, double b, double c, double d, double e, double f,
                     double tx, double ty, double tz) {
    setValues(a, b, c, d, e, f, tx, ty, tz);
  }

  CSymmetricMatrix3D(const CSymmetricMatrix3D &a) {
    memcpy(&m00_, &a.m00_, 16*sizeof(double));
  }

 ~CSymmetricMatrix3D() { }

  void setIdentity() {
    setValues(1.0, 0.0, 0.0,
                   1.0, 0.0,
                        1.0 ,
              0.0, 0.0, 0.0);
  }

  void setTranslation(double tx, double ty, double tz) {
    setValues(1.0, 0.0, 0.0,
                   1.0, 0.0,
                        1.0,
               tx,  ty, tz );
  }

  void setScale(double sx, double sy, double sz) {
    setValues(sx , 0.0, 0.0,
                    sy, 0.0,
                         sz,
              0.0, 0.0, 0.0);
  }

  void multiplyPoint(double xi, double yi, double zi, double *xo, double *yo, double *zo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi + m03_;
    *yo = m10_*xi + m11_*yi + m12_*zi + m13_;
    *zo = m20_*xi + m21_*yi + m22_*zi + m23_;
  }

  void multiplyPoint(const CPoint3D &ipoint, CPoint3D &opoint) const {
    opoint.x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z + m03_;
    opoint.y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z + m13_;
    opoint.z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z + m23_;
  }

  void multiplyVector(const CVector3D &ivector, CVector3D &ovector) const {
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

  void setValues(double a, double b, double c, double d, double e, double f) {
    m00_ =   a; m01_ =   b; m02_ =   c; m03_ = 0.0;
    m10_ =   b; m11_ =   d; m12_ =   e; m13_ = 0.0;
    m20_ =   c; m21_ =   e; m22_ =   f; m23_ = 0.0;
    m30_ = 0.0; m31_ = 0.0; m32_ = 0.0; m33_ = 1.0;
  }

  void setValues(double a, double b, double c, double d, double e, double f,
                 double tx, double ty, double tz) {
    m00_ = a  ; m01_ = b  ; m02_ = c  ; m03_ = tx;
    m10_ = b  ; m11_ = d  ; m12_ = e  ; m13_ = ty;
    m20_ = c  ; m21_ = e  ; m22_ = f  ; m23_ = tz;
    m30_ = 0.0; m31_ = 0.0; m32_ = 0.0; m33_ = 1.0;
  }

  void getValues(double *a, double *b, double *c, double *d, double *e, double *f,
                 double *tx, double *ty, double *tz) const {
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

  void getValues(double *a, double *b, double *c, double *d, double *e, double *f) const {
    if (a) *a  = m00_;
    if (b) *b  = m01_;
    if (c) *c  = m02_;
    if (d) *d  = m11_;
    if (e) *e  = m12_;
    if (f) *f  = m13_;
  }

  void getValues(double v[9]) const {
    v[0] = m00_; v[1] = m01_; v[2] = m02_;
    v[3] = m10_; v[4] = m11_; v[5] = m12_;
    v[6] = m20_; v[7] = m21_; v[8] = m22_;
  }

  CSymmetricMatrix3D &operator=(const CSymmetricMatrix3D &a) {
    memcpy(&m00_, &a.m00_, 16*sizeof(double));

    return *this;
  }

  CSymmetricMatrix3D &operator*=(double s) {
    CSymmetricMatrix3D a = *this;

    m00_ = a.m00_*s; m01_ = a.m00_*s; m02_ = a.m00_*s; m03_ = a.m00_*s;
    m10_ = a.m10_*s; m11_ = a.m10_*s; m12_ = a.m10_*s; m13_ = a.m10_*s;
    m20_ = a.m20_*s; m21_ = a.m20_*s; m22_ = a.m20_*s; m23_ = a.m20_*s;
    m30_ = 0.0     ; m31_ = 0.0     ; m32_ = 0.0     ; m33_ = 1.0     ;

    return *this;
  }

  CSymmetricMatrix3D operator*(double s) {
    CSymmetricMatrix3D c = *this;

    c *= s;

    return c;
  }

  void setValue(unsigned int i, double value) { (&m00_)[i] = value; }

  double getValue(unsigned int i) { return (&m00_)[i]; }

  void zero() { memset(&m00_, 0, 16*sizeof(double)); }

  double operator[](unsigned int i) { return (&m00_)[i]; }

  const double &operator[](unsigned int i) const { return (&m00_)[i]; }

  friend std::ostream &operator<<(std::ostream &os, const CSymmetricMatrix3D &matrix) {
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
  static double calcDeterminant(double a, double b, double c, double d, double e, double f);

 private:
  double m00_, m01_, m02_, m03_;
  double m10_, m11_, m12_, m13_;
  double m20_, m21_, m22_, m23_;
  double m30_, m31_, m32_, m33_;
};

#endif
