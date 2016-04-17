#ifndef CMATRIX_3X3_H
#define CMATRIX_3X3_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CVector3D.h>
#include <CThrow.h>

/* / m00 m01 m02 \ */
/* | m10 m11 m12 | */
/* \ m20 m21 m22 / */

template<typename T>
class CMatrix3x3T {
 public:
  enum Type {
    CMATRIX_3x3_IDENTITY
  };

 private:
  typedef CPoint3DT<T>  Point;
  typedef CVector3DT<T> Vector;

 private:
  union {
    T m1_[9];
    T m2_[3][3];

    struct {
      T m00_, m01_, m02_;
      T m10_, m11_, m12_;
      T m20_, m21_, m22_;
    };

    struct {
      T a_, b_, c_;
      T d_, e_, f_;
      T g_, h_, i_;
    };
  };

 public:
  CMatrix3x3T() { }

  explicit CMatrix3x3T(typename CMatrix3x3T::Type type) {
    if (type == CMATRIX_3x3_IDENTITY)
      setIdentity();
    else
      CTHROW("Bad Matrix Type");
  }

  CMatrix3x3T(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) :
    m00_(m00), m01_(m01), m02_(m02),
    m10_(m10), m11_(m11), m12_(m12),
    m20_(m20), m21_(m21), m22_(m22) {
  }

  CMatrix3x3T(const T *m, int n) {
    if      (n == 9)
      setValues(m[ 0], m[ 1], m[ 2],
                m[ 3], m[ 4], m[ 5],
                m[ 6], m[ 7], m[ 8]);
    else
      CTHROW("Invalid size");
  }

  CMatrix3x3T(CMathGen::AxisType3D axis, T angle,
             CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    setRotation(axis, angle, handedness);
  }

  CMatrix3x3T(const CMatrix3x3T &a) {
    memcpy(m1_, a.m1_, sizeof(m1_));
  }

 ~CMatrix3x3T() { }

  void setIdentity() {
    m00_ = 1.0, m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = 1.0, m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = 1.0;
  }

  void setScale(T s) {
    m00_ = s  , m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = s  , m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = s  ;
  }

  void setScale(T sx, T sy, T sz) {
    m00_ = sx , m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = sy , m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = sz ;
  }

  void setRotation(CMathGen::AxisType3D axis, T angle,
                   CMathGen::Handedness handedness =
                    CMathGen::RIGHT_HANDEDNESS) {
    if (handedness == CMathGen::RIGHT_HANDEDNESS)
      setInnerRotationRHS(axis, angle);
    else
      setInnerRotationLHS(axis, angle);
  }

  void setXYZRotation(T x_angle, T y_angle, T z_angle,
                      CMathGen::Handedness handedness =
                       CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3x3T xmatrix(CMathGen::X_AXIS_3D, x_angle, handedness);
    CMatrix3x3T ymatrix(CMathGen::Y_AXIS_3D, y_angle, handedness);
    CMatrix3x3T zmatrix(CMathGen::Z_AXIS_3D, z_angle, handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setXYZRotation(const Vector &angles,
                      CMathGen::Handedness handedness =
                       CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3x3T xmatrix(CMathGen::X_AXIS_3D, angles.getX(), handedness);
    CMatrix3x3T ymatrix(CMathGen::Y_AXIS_3D, angles.getY(), handedness);
    CMatrix3x3T zmatrix(CMathGen::Z_AXIS_3D, angles.getZ(), handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setGenRotation(T x1, T y1, T z1, T x2, T y2, T z2, T angle,
                      CMathGen::Handedness handedness =
                       CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3x3T matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7;

    matrix1.setTranslation(-x1, -y1, -z1);
    matrix2.setTranslation( x1,  y1,  z1);

    T theta = CMathGen::atan2(x2 - x1, y2 - y1);

    matrix3.setRotation(CMathGen::Z_AXIS_3D,  theta, handedness);
    matrix4.setRotation(CMathGen::Z_AXIS_3D, -theta, handedness);

    T v = ::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));

    T beta = CMathGen::atan2(z2 - z1, v);

    matrix5.setRotation(CMathGen::Y_AXIS_3D,  beta, handedness);
    matrix6.setRotation(CMathGen::Y_AXIS_3D, -beta, handedness);

    matrix7.setRotation(CMathGen::Z_AXIS_3D, angle, handedness);

    *this = matrix2*(matrix4*(matrix6*(matrix7*(matrix5*(matrix3*matrix1)))));
  }

  void setGenRotation(const Vector &axis, T angle,
                      CMathGen::Handedness handedness =
                       CMathGen::RIGHT_HANDEDNESS) {
    Vector a = axis.normalized();

    T c = ::cos(angle);
    T s = ::sin(angle);

    T c1 = 1.0 - c;

    T axx = a.x*a.x;
    T axy = a.x*a.y;
    T axz = a.x*a.z;
    T ayy = a.y*a.y;
    T ayz = a.y*a.z;
    T azz = a.z*a.z;

    T axs = a.x*s;
    T ays = a.y*s;
    T azs = a.z*s;

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

  void setLookAt(const Point &eye, const Point &center, const Vector &up) {
    Vector dir = (center - eye).normalize();

    CVector3D right = dir.crossProduct(up.normalize());
    CVector3D newUp = right.crossProduct(dir);

    dir = -dir;

    setColumn(0, right);
    setColumn(1, newUp);
    setColumn(2, dir  );

    setOuterTranslate(eye, eye, eye);
  }

  void setEye(T x1, T y1, T z1, T x2, T y2, T z2,
              CMathGen::Handedness handedness =
               CMathGen::RIGHT_HANDEDNESS) {
    T angle1, angle2, angle3;

    calcEye(x1, y1, z1, x2, y2, z2, &angle1, &angle2, &angle3);

    CMatrix3x3T matrix1, matrix2, matrix3, matrix4;

    matrix1.setTranslation(-x1, -y1, -z1);

    matrix2.setRotation(CMathGen::Z_AXIS_3D,  angle1, handedness);
    matrix3.setRotation(CMathGen::Y_AXIS_3D,  angle2, handedness);
    matrix4.setRotation(CMathGen::Z_AXIS_3D, -angle3, handedness);

    *this = matrix4*(matrix3*(matrix2*matrix1));
  }

  static void calcEye(T x1, T y1, T z1, T x2, T y2, T z2,
                      T *angle1, T *angle2, T *angle3) {
    T dx = x2 - x1;
    T dy = y2 - y1;
    T dz = z2 - z1;

    *angle1 = CMathGen::atan2(-dx, -dy);

    T v = ::sqrt(dx*dx + dy*dy);

    *angle2 = CMathGen::atan2(-dz, v);

    T w = ::sqrt(v*v + dz*dz);

    *angle3 = CMathGen::atan2(-dx*w, dy*dz);
  }

  void setValues(T a, T b, T c, T d, T e, T f, T g, T h, T i) {
    m00_ = a, m01_ = b, m02_ = c;
    m10_ = d, m11_ = e, m12_ = f;
    m20_ = g, m21_ = h, m22_ = i;
  }

  void getValues(T *a, T *b, T *c, T *d, T *e, T *f, T *g, T *h, T *i) const {
    if (a) *a = m00_; if (b) *b = m01_; if (c) *c = m02_;
    if (d) *d = m10_; if (e) *e = m11_; if (f) *f = m12_;
    if (g) *g = m20_; if (h) *h = m21_; if (i) *i = m22_;
  }

  void getValues(T *v, int n) const {
    if      (n == 9) {
      v[0] = m00_; v[1] = m01_; v[2] = m02_;
      v[3] = m10_; v[4] = m11_; v[5] = m12_;
      v[6] = m20_; v[7] = m21_; v[8] = m22_;
    }
    else
      CTHROW("Invalid size");
  }

  void setColumn(int c, T x, T y, T z) {
    m2_[0][c] = x, m2_[1][c] = y; m2_[2][c] = z;
  }

  void setColumn(int c, const Point &point) {
    m2_[0][c] = point.x, m2_[1][c] = point.y, m2_[2][c] = point.z;
  }

  void setColumn(int c, const Vector &vector) {
    vector.getXYZ(&m2_[0][c], &m2_[1][c], &m2_[2][c]);
  }

  void getColumn(int c, T *x, T *y, T *z) {
    if (x) *x = m2_[0][c];
    if (y) *y = m2_[1][c];
    if (z) *z = m2_[2][c];
  }

  void setRow(int r, T x, T y, T z) {
    m2_[r][0] = x, m2_[r][1] = y; m2_[r][2] = z;
  }

  void setRow(int r, const Point &point) {
    m2_[r][0] = point.x, m2_[r][1] = point.y, m2_[r][2] = point.z;
  }

  void setRow(int r, const Vector &vector) {
    vector.getXYZ(&m2_[r][0], &m2_[r][1], &m2_[r][2]);
  }

  void getRow(int r, T *x, T *y, T *z) {
    if (x) *x = m2_[r][0];
    if (y) *y = m2_[r][1];
    if (z) *z = m2_[r][2];
  }

  void multiplyPoint(T  xi, T  yi, T  zi, T *xo, T *yo, T *zo) const {
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
    T ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    T ox = m00_*ix + m01_*iy + m02_*iz;
    T oy = m10_*ix + m11_*iy + m12_*iz;
    T oz = m20_*ix + m21_*iy + m22_*iz;

    ovector.setXYZ(ox, oy, oz);
  }

  void preMultiplyPoint(T xi, T yi, T zi, T *xo, T *yo, T *zo) const {
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
    T ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    T ox = m00_*ix + m10_*iy + m20_*iz;
    T oy = m01_*ix + m11_*iy + m21_*iz;
    T oz = m02_*ix + m12_*iy + m22_*iz;

    ovector.setXYZ(ox, oy, oz);
  }

  void transpose() {
    swap(m10_, m01_);
    swap(m20_, m02_);
    swap(m21_, m12_);
  }

  CMatrix3x3T transposed() const {
    return CMatrix3x3T(m00_, m10_, m20_,
                       m01_, m11_, m21_,
                       m02_, m12_, m22_);
  }

  bool invert(CMatrix3x3T &imatrix) const {
    T d = determinant();

    if (::fabs(d) == 0.0)
      return false;

    T id = 1.0/d;

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

  CMatrix3x3T inverse() const {
    CMatrix3x3T imatrix;

    if (! invert(imatrix))
      CTHROW("Divide by 0.0");

    return imatrix;
  }

  T determinant() const {
    return (m00_*calcDeterminant(m11_, m12_, m21_, m22_) -
            m01_*calcDeterminant(m10_, m12_, m20_, m22_) +
            m02_*calcDeterminant(m10_, m11_, m20_, m21_));
  }

  void normalize() {
    T d = determinant();

    T id = 1.0/d;

    for (int i = 0; i < 9; ++i)
      m1_[i] *= id;
  }

  static CMatrix3x3T *newIdentityMatrix() {
    CMatrix3x3T *m = new CMatrix3x3T();

    m->setIdentity();

    return m;
  }

  static bool solveAXeqB(const CMatrix3x3T &a, Point &x, const Point &b) {
    T det_a = a.determinant();

    if (::abs(det_a) < 0.0)
      return false;

    T idet_a = 1.0/det_a;

    CMatrix3x3T t(a);

    t.setColumn(0, b.x, b.y, b.z);

    T det_t = t.determinant();

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

  void zero() { memset(m1_, 0, sizeof(m1_)); }

  CMatrix3x3T &operator=(const CMatrix3x3T &a) {
    memcpy(m1_, a.m1_, sizeof(m1_));

    return *this;
  }

  CMatrix3x3T &operator+=(const CMatrix3x3T &b) {
    m00_ += b.m00_; m01_ += b.m01_; m02_ += b.m02_;
    m10_ += b.m10_; m11_ += b.m11_; m12_ += b.m12_;
    m20_ += b.m20_; m21_ += b.m21_; m22_ += b.m22_;

    return *this;
  }

  CMatrix3x3T operator+(const CMatrix3x3T &b) {
    CMatrix3x3T c = *this;

    c += b;

    return c;
  }

  CMatrix3x3T &operator-=(const CMatrix3x3T &b) {
    m00_ -= b.m00_; m01_ -= b.m01_, m02_ -= b.m02_;
    m10_ -= b.m10_; m11_ -= b.m11_; m12_ -= b.m12_;
    m20_ -= b.m20_; m21_ -= b.m21_; m22_ -= b.m22_;

    return *this;
  }

  CMatrix3x3T operator-(const CMatrix3x3T &b) {
    CMatrix3x3T c = *this;

    c -= b;

    return c;
  }

  CMatrix3x3T &operator*=(const CMatrix3x3T &b) {
    CMatrix3x3T a = *this;

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

  CMatrix3x3T operator*(const CMatrix3x3T &b) const {
    CMatrix3x3T c = *this;

    c *= b;

    return c;
  }

  CMatrix3x3T &operator*=(T s) {
    CMatrix3x3T a = *this;

    m00_ = a.m00_*s; m01_ = a.m01_*s; m02_ = a.m02_*s;
    m10_ = a.m10_*s; m11_ = a.m11_*s; m12_ = a.m12_*s;
    m20_ = a.m20_*s; m21_ = a.m21_*s; m22_ = a.m22_*s;

    return *this;
  }

  CMatrix3x3T operator*(T s) {
    CMatrix3x3T c = *this;

    c *= s;

    return c;
  }

  friend Point operator*(const CMatrix3x3T &m, const Point &p) {
    Point p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend Point operator*(const Point &p, const CMatrix3x3T &m) {
    Point p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend Vector operator*(const CMatrix3x3T &m, const Vector &v) {
    Vector v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend Vector operator*(const Vector &v, const CMatrix3x3T &m) {
    Vector v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  CMatrix3x3T &operator/=(const CMatrix3x3T &b) {
    CMatrix3x3T bi;

    if (! b.invert(bi)) {
      CTHROW("Divide by 0.0");
      return *this;
    }

    return (*this) *= bi;
  }

  CMatrix3x3T operator/(const CMatrix3x3T &b) {
    CMatrix3x3T c = *this;

    c /= b;

    return c;
  }

  void setValue(unsigned int i, T value) {
    m1_[i] = value;
  }

  void setValue(unsigned int i, unsigned int j, T value) {
    m2_[i][j] = value;
  }

  T getValue(unsigned int i) const {
    return m1_[i];
  }

  T getValue(unsigned int i, unsigned int j) const {
    return m2_[i][j];
  }

  T operator[](unsigned int i) { return m1_[i]; }

  const T &operator[](unsigned int i) const { return m1_[i]; }

  void print(ostream &os) const {
    os << "(" << m00_ << "," << m01_ << "," << m02_ << ")" << endl;
    os << "(" << m10_ << "," << m11_ << "," << m12_ << ")" << endl;
    os << "(" << m20_ << "," << m21_ << "," << m22_ << ")" << endl;
  }

  friend ostream &operator<<(ostream &os, const CMatrix3x3T &matrix) {
    matrix.print(os);

    return os;
  }

 private:
  static T calcDeterminant(T m00, T m01, T m10, T m11) {
    return m00*m11 - m01*m10;
  }
};

typedef CMatrix3x3T<double> CMatrix3x3;

#endif
