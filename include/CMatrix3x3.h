#ifndef CMATRIX_3X3_H
#define CMATRIX_3X3_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CVector3D.h>
#include <iostream>
#include <cstring>

/* / m00 m01 m02 \ */
/* | m10 m11 m12 | */
/* \ m20 m21 m22 / */

class CMatrix3x3 {
 public:
  enum Type {
    CMATRIX_3x3_IDENTITY
  };

 private:
  typedef CPoint3D  Point;
  typedef CVector3D Vector;

 public:
  CMatrix3x3() { }

  explicit CMatrix3x3(typename CMatrix3x3::Type type) {
    if (type == CMATRIX_3x3_IDENTITY)
      setIdentity();
    else
      assert(false && "Bad Matrix Type");
  }

  CMatrix3x3(double m00, double m01, double m02, double m10, double m11, double m12,
             double m20, double m21, double m22) :
    m00_(m00), m01_(m01), m02_(m02),
    m10_(m10), m11_(m11), m12_(m12),
    m20_(m20), m21_(m21), m22_(m22) {
  }

  CMatrix3x3(const double *m, int n) {
    if      (n == 9)
      memcpy(&m00_, m, 9*sizeof(double));
    else
      assert(false && "Invalid size");
  }

  CMatrix3x3(CMathGen::AxisType3D axis, double angle,
             CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    setRotation(axis, angle, handedness);
  }

  CMatrix3x3(const CMatrix3x3 &a) {
    memcpy(&m00_, &a.m00_, 9*sizeof(double));
  }

 ~CMatrix3x3() { }

  //----------

  void setIdentity() {
    m00_ = 1.0, m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = 1.0, m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = 1.0;
  }

  //----------

  void setTranslation(double tx, double ty, double tz) {
    setIdentity();

    m02_ = tx;
    m10_ = ty;
    m20_ = tz;
  }

  //----------

  void setScale(double s) {
    m00_ = s  , m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = s  , m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = s  ;
  }

  void setScale(double sx, double sy, double sz) {
    m00_ = sx , m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = sy , m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = sz ;
  }

  //----------

  void setRotation(CMathGen::AxisType3D axis, double angle,
                   CMathGen::Handedness handedness =
                    CMathGen::RIGHT_HANDEDNESS) {
    if (handedness == CMathGen::RIGHT_HANDEDNESS)
      setInnerRotationRHS(axis, angle);
    else
      setInnerRotationLHS(axis, angle);
  }

  void setXYZRotation(double x_angle, double y_angle, double z_angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3x3 xmatrix(CMathGen::X_AXIS_3D, x_angle, handedness);
    CMatrix3x3 ymatrix(CMathGen::Y_AXIS_3D, y_angle, handedness);
    CMatrix3x3 zmatrix(CMathGen::Z_AXIS_3D, z_angle, handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setXYZRotation(const Vector &angles,
                      CMathGen::Handedness handedness =
                       CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3x3 xmatrix(CMathGen::X_AXIS_3D, angles.getX(), handedness);
    CMatrix3x3 ymatrix(CMathGen::Y_AXIS_3D, angles.getY(), handedness);
    CMatrix3x3 zmatrix(CMathGen::Z_AXIS_3D, angles.getZ(), handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setGenRotation(double x1, double y1, double z1,
                      double x2, double y2, double z2, double angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3x3 matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7;

    matrix1.setTranslation(-x1, -y1, -z1);
    matrix2.setTranslation( x1,  y1,  z1);

    double theta = CMathGen::atan2(x2 - x1, y2 - y1);

    matrix3.setRotation(CMathGen::Z_AXIS_3D,  theta, handedness);
    matrix4.setRotation(CMathGen::Z_AXIS_3D, -theta, handedness);

    double v = ::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));

    double beta = CMathGen::atan2(z2 - z1, v);

    matrix5.setRotation(CMathGen::Y_AXIS_3D,  beta, handedness);
    matrix6.setRotation(CMathGen::Y_AXIS_3D, -beta, handedness);

    matrix7.setRotation(CMathGen::Z_AXIS_3D, angle, handedness);

    *this = matrix2*(matrix4*(matrix6*(matrix7*(matrix5*(matrix3*matrix1)))));
  }

  void setGenRotation(const Vector &axis, double angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    Vector a = axis.normalized();

    double c = ::cos(angle);
    double s = ::sin(angle);

    double c1 = 1.0 - c;

    double axx = a.getX()*a.getX();
    double axy = a.getX()*a.getY();
    double axz = a.getX()*a.getZ();
    double ayy = a.getY()*a.getY();
    double ayz = a.getY()*a.getZ();
    double azz = a.getZ()*a.getZ();

    double axs = a.getX()*s;
    double ays = a.getY()*s;
    double azs = a.getZ()*s;

    if (handedness == CMathGen::RIGHT_HANDEDNESS) {
      m00_ = axx*c1 + c;
      m01_ = axy*c1 + azs;
      m02_ = axz*c1 - ays;

      m10_ = axy*c1 - azs;
      m11_ = ayy*c1 + c;
      m12_ = ayz*c1 + axs;

      m20_ = axz*c1 + ays;
      m21_ = ayz*c1 - axs;
      m22_ = azz*c1 + c;
    }
    else {
      m00_ = axx*c1 + c;
      m01_ = axy*c1 - azs;
      m02_ = axz*c1 + ays;

      m10_ = axy*c1 + azs;
      m11_ = ayy*c1 + c;
      m12_ = ayz*c1 - axs;

      m20_ = axz*c1 - ays;
      m21_ = ayz*c1 + axs;
      m22_ = azz*c1 + c;
    }
  }

  void getValues(double *a, double *b, double *c, double *d, double *e, double *f,
                 double *g, double *h, double *i) const {
    if (a) *a = m00_; if (b) *b = m01_; if (c) *c = m02_;
    if (d) *d = m10_; if (e) *e = m11_; if (f) *f = m12_;
    if (g) *g = m20_; if (h) *h = m21_; if (i) *i = m22_;
  }

  void getValues(double *v, int n) const {
    if      (n == 9) {
      v[0] = m00_; v[1] = m01_; v[2] = m02_;
      v[3] = m10_; v[4] = m11_; v[5] = m12_;
      v[6] = m20_; v[7] = m21_; v[8] = m22_;
    }
    else
      assert(false && "Invalid size");
  }

  void setColumn(int c, double x, double y, double z) {
    assert(c < 3);

    double *m = &(&m00_)[c];

    m[0] = x, m[3] = y; m[6] = z;
  }

  void setColumn(int c, const Point &point) {
    setColumn(c, point.x, point.y, point.z);
  }

  void setColumn(int c, const Vector &vector) {
    setColumn(c, vector.getX(), vector.getY(), vector.getZ());
  }

  void getColumn(int c, double *x, double *y, double *z) {
    assert(c < 3);

    double *m = &(&m00_)[c];

    if (x) *x = m[0];
    if (y) *y = m[3];
    if (z) *z = m[6];
  }

  void setRow(int r, double x, double y, double z) {
    assert(r < 3);

    double *m = &(&m00_)[r*3];

    m[0] = x, m[1] = y; m[2] = z;
  }

  void setRow(int r, const Point &point) {
    setRow(r, point.x, point.y, point.z);
  }

  void setRow(int r, const Vector &vector) {
    setRow(r, vector.getX(), vector.getY(), vector.getZ());
  }

  void getRow(int r, double *x, double *y, double *z) {
    assert(r < 3);

    double *m = &(&m00_)[r*3];

    if (x) *x = m[0];
    if (y) *y = m[1];
    if (z) *z = m[2];
  }

  void multiplyPoint(double  xi, double  yi, double  zi,
                     double *xo, double *yo, double *zo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi;
    *yo = m10_*xi + m11_*yi + m12_*zi;
    *zo = m20_*xi + m21_*yi + m22_*zi;
  }

  void multiplyPoint(const Point &point, Point &opoint) const {
    opoint.x = m00_*point.x + m01_*point.y + m02_*point.z;
    opoint.y = m10_*point.x + m11_*point.y + m12_*point.z;
    opoint.z = m20_*point.x + m21_*point.y + m22_*point.z;
  }

  void multiplyVector(const Vector &ivector, Vector &ovector) const {
    double ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    double ox = m00_*ix + m01_*iy + m02_*iz;
    double oy = m10_*ix + m11_*iy + m12_*iz;
    double oz = m20_*ix + m21_*iy + m22_*iz;

    ovector.setXYZ(ox, oy, oz);
  }

  void preMultiplyPoint(double xi, double yi, double zi,
                        double *xo, double *yo, double *zo) const {
    *xo = m00_*xi + m10_*yi + m20_*zi;
    *yo = m01_*xi + m11_*yi + m21_*zi;
    *zo = m02_*xi + m12_*yi + m22_*zi;
  }

  void preMultiplyPoint(const Point &ipoint, Point &opoint) const {
    opoint.x = m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z;
    opoint.z = m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z;
  }

  void preMultiplyVector(const Vector &ivector, Vector &ovector) const {
    double ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    double ox = m00_*ix + m10_*iy + m20_*iz;
    double oy = m01_*ix + m11_*iy + m21_*iz;
    double oz = m02_*ix + m12_*iy + m22_*iz;

    ovector.setXYZ(ox, oy, oz);
  }

  void transpose() {
    std::swap(m10_, m01_);
    std::swap(m20_, m02_);
    std::swap(m21_, m12_);
  }

  CMatrix3x3 transposed() const {
    return CMatrix3x3(m00_, m10_, m20_,
                       m01_, m11_, m21_,
                       m02_, m12_, m22_);
  }

  bool invert(CMatrix3x3 &imatrix) const {
    double d = determinant();

    if (::fabs(d) == 0.0)
      return false;

    double id = 1.0/d;

    imatrix.m00_ =  id*calcDeterminant(m11_, m12_, m21_, m22_);
    imatrix.m10_ = -id*calcDeterminant(m10_, m12_, m20_, m22_);
    imatrix.m20_ =  id*calcDeterminant(m10_, m11_, m20_, m21_);

    imatrix.m01_ = -id*calcDeterminant(m01_, m02_, m21_, m22_);
    imatrix.m11_ =  id*calcDeterminant(m00_, m02_, m20_, m22_);
    imatrix.m21_ = -id*calcDeterminant(m00_, m01_, m20_, m21_);

    imatrix.m02_ =  id*calcDeterminant(m01_, m02_, m11_, m12_);
    imatrix.m12_ = -id*calcDeterminant(m00_, m02_, m10_, m12_);
    imatrix.m22_ =  id*calcDeterminant(m00_, m01_, m10_, m11_);

    return true;
  }

  CMatrix3x3 inverse() const {
    CMatrix3x3 imatrix;

    if (! invert(imatrix))
      assert(false && "Divide by 0.0");

    return imatrix;
  }

  double determinant() const {
    return (m00_*calcDeterminant(m11_, m12_, m21_, m22_) -
            m01_*calcDeterminant(m10_, m12_, m20_, m22_) +
            m02_*calcDeterminant(m10_, m11_, m20_, m21_));
  }

  void normalize() {
    double d = determinant();

    double id = 1.0/d;

    for (int i = 0; i < 9; ++i)
      (&m00_)[i] *= id;
  }

  static CMatrix3x3 *newIdentityMatrix() {
    CMatrix3x3 *m = new CMatrix3x3();

    m->setIdentity();

    return m;
  }

  static bool solveAXeqB(const CMatrix3x3 &a, Point &x, const Point &b) {
    double det_a = a.determinant();

    if (::abs(det_a) < 0.0)
      return false;

    double idet_a = 1.0/det_a;

    CMatrix3x3 t(a);

    t.setColumn(0, b.x, b.y, b.z);

    double det_t = t.determinant();

    x.x = det_t*idet_a;

    t = a;

    t.setColumn(1, b.x, b.y, b.z);

    det_t = t.determinant();

    x.y = det_t*idet_a;

    t = a;

    t.setColumn(2, b.x, b.y, b.z);

    det_t = t.determinant();

    x.z = det_t*idet_a;

    return true;
  }

  void zero() { memset(&m00_, 0, 9*sizeof(double)); }

  CMatrix3x3 &operator=(const CMatrix3x3 &a) {
    memcpy(&m00_, &a.m00_, 9*sizeof(double));

    return *this;
  }

  CMatrix3x3 &operator+=(const CMatrix3x3 &b) {
    m00_ += b.m00_; m01_ += b.m01_; m02_ += b.m02_;
    m10_ += b.m10_; m11_ += b.m11_; m12_ += b.m12_;
    m20_ += b.m20_; m21_ += b.m21_; m22_ += b.m22_;

    return *this;
  }

  CMatrix3x3 operator+(const CMatrix3x3 &b) {
    CMatrix3x3 c = *this;

    c += b;

    return c;
  }

  CMatrix3x3 &operator-=(const CMatrix3x3 &b) {
    m00_ -= b.m00_; m01_ -= b.m01_, m02_ -= b.m02_;
    m10_ -= b.m10_; m11_ -= b.m11_; m12_ -= b.m12_;
    m20_ -= b.m20_; m21_ -= b.m21_; m22_ -= b.m22_;

    return *this;
  }

  CMatrix3x3 operator-(const CMatrix3x3 &b) {
    CMatrix3x3 c = *this;

    c -= b;

    return c;
  }

  CMatrix3x3 &operator*=(const CMatrix3x3 &b) {
    CMatrix3x3 a = *this;

    m00_ = a.m00_*b.m00_ + a.m01_*b.m10_ + a.m02_*b.m20_;
    m01_ = a.m00_*b.m01_ + a.m01_*b.m11_ + a.m02_*b.m21_;
    m02_ = a.m00_*b.m02_ + a.m01_*b.m12_ + a.m02_*b.m22_;

    m10_ = a.m10_*b.m00_ + a.m11_*b.m10_ + a.m12_*b.m20_;
    m11_ = a.m10_*b.m01_ + a.m11_*b.m11_ + a.m12_*b.m21_;
    m12_ = a.m10_*b.m02_ + a.m11_*b.m12_ + a.m12_*b.m22_;

    m20_ = a.m20_*b.m00_ + a.m21_*b.m10_ + a.m22_*b.m20_;
    m21_ = a.m20_*b.m01_ + a.m21_*b.m11_ + a.m22_*b.m21_;
    m22_ = a.m20_*b.m02_ + a.m21_*b.m12_ + a.m22_*b.m22_;

    return *this;
  }

  CMatrix3x3 operator*(const CMatrix3x3 &b) const {
    CMatrix3x3 c = *this;

    c *= b;

    return c;
  }

  CMatrix3x3 &operator*=(double s) {
    CMatrix3x3 a = *this;

    m00_ = a.m00_*s; m01_ = a.m01_*s; m02_ = a.m02_*s;
    m10_ = a.m10_*s; m11_ = a.m11_*s; m12_ = a.m12_*s;
    m20_ = a.m20_*s; m21_ = a.m21_*s; m22_ = a.m22_*s;

    return *this;
  }

  CMatrix3x3 operator*(double s) {
    CMatrix3x3 c = *this;

    c *= s;

    return c;
  }

  friend Point operator*(const CMatrix3x3 &m, const Point &p) {
    Point p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend Point operator*(const Point &p, const CMatrix3x3 &m) {
    Point p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend Vector operator*(const CMatrix3x3 &m, const Vector &v) {
    Vector v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend Vector operator*(const Vector &v, const CMatrix3x3 &m) {
    Vector v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  CMatrix3x3 &operator/=(const CMatrix3x3 &b) {
    CMatrix3x3 bi;

    if (! b.invert(bi)) {
      assert(false && "Divide by 0.0");
      return *this;
    }

    return (*this) *= bi;
  }

  CMatrix3x3 operator/(const CMatrix3x3 &b) {
    CMatrix3x3 c = *this;

    c /= b;

    return c;
  }

  void setValue(unsigned int i, double value) {
    (&m00_)[i] = value;
  }

  void setValue(unsigned int i, unsigned int j, double value) {
    assert(i < 3 && j < 3);

    double &m = (&m00_)[3*j + i];

    m = value;
  }

  double getValue(unsigned int i) const {
    return (&m00_)[i];
  }

  double getValue(unsigned int i, unsigned int j) const {
    assert(i < 3 && j < 3);

    const double &m = (&m00_)[3*j + i];

    return m;
  }

  double operator[](unsigned int i) { return (&m00_)[i]; }

  const double &operator[](unsigned int i) const { return (&m00_)[i]; }

  //------

  void print(std::ostream &os) const {
    os << "(" << m00_ << "," << m01_ << "," << m02_ << ")" << std::endl;
    os << "(" << m10_ << "," << m11_ << "," << m12_ << ")" << std::endl;
    os << "(" << m20_ << "," << m21_ << "," << m22_ << ")" << std::endl;
  }

  friend std::ostream &operator<<(std::ostream &os, const CMatrix3x3 &matrix) {
    matrix.print(os);

    return os;
  }

  //------

 private:
  void setInnerRotationRHS(CMathGen::AxisType3D axis, double angle) {
    double c = ::cos(angle);
    double s = ::sin(angle);

    if      (axis == CMathGen::X_AXIS_3D) {
      m00_ = 1.0; m01_ = 0.0; m02_ = 0.0;
      m10_ = 0.0; m11_ =   c; m12_ =   s;
      m20_ = 0.0; m21_ =  -s; m22_ =   c;
    }
    else if (axis == CMathGen::Y_AXIS_3D) {
      m00_ =   c; m01_ = 0.0; m02_ =  -s;
      m10_ = 0.0; m11_ = 1.0; m12_ = 0.0;
      m20_ =   s; m21_ = 0.0; m22_ =   c;
    }
    else {
      m00_ =   c; m01_ =   s; m02_ = 0.0;
      m10_ =  -s; m11_ =   c; m12_ = 0.0;
      m20_ = 0.0; m21_ = 0.0; m22_ = 1.0;
    }
  }

  void setInnerRotationLHS(CMathGen::AxisType3D axis, double angle) {
    double c = ::cos(angle);
    double s = ::sin(angle);

    if      (axis == CMathGen::X_AXIS_3D) {
      m00_ = 1.0; m01_ = 0.0; m02_ = 0.0;
      m10_ = 1.0; m11_ =   c; m12_ =  -s;
      m20_ = 0.0; m21_ =   s; m22_ =   c;
    }
    else if (axis == CMathGen::Y_AXIS_3D) {
      m00_ =   c; m01_ = 0.0; m02_ =   s;
      m10_ = 0.0; m11_ = 1.0; m12_ = 0.0;
      m20_ =  -s; m21_ = 0.0; m22_ =   c;
    }
    else {
      m00_ =   c; m01_ =  -s; m02_ = 0.0;
      m10_ =   s; m11_ =   c; m12_ = 0.0;
      m20_ = 0.0; m21_ = 0.0; m22_ = 1.0;
    }
  }

 private:
  static double calcDeterminant(double m00, double m01, double m10, double m11) {
    return m00*m11 - m01*m10;
  }

 private:
  double m00_ { 0 }, m01_ { 0 }, m02_ { 0 };
  double m10_ { 0 }, m11_ { 0 }, m12_ { 0 };
  double m20_ { 0 }, m21_ { 0 }, m22_ { 0 };
};

#endif
