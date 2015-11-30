#ifndef CMATRIX_3DH_H
#define CMATRIX_3DH_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CVector3D.h>
#include <CMatrix3D.h>
#include <CThrow.h>
#include <cstring>

/* Homogeneous 3D Matrix */

/* / m00 m01 m02 m03 \ */
/* | m10 m11 m12 m13 | */
/* | m20 m21 m22 m23 | */
/* \ m30 m31 m32 m33 / */

template<typename T>
class CMatrix3DHT {
 public:
  enum Type {
    CMATRIX_3DH_IDENTITY
  };

 private:
  typedef CPoint3DT<T>   Point;
  typedef CVector3DT<T>  Vector;
  typedef CMatrix3DHT<T> MatrixH;
  typedef CMatrix3DT<T>  Matrix3;

 private:
  T m00_, m01_, m02_, m03_;
  T m10_, m11_, m12_, m13_;
  T m20_, m21_, m22_, m23_;
  T m30_, m31_, m32_, m33_;

 public:
  // constructor/destructor
  CMatrix3DHT() :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
  }

 ~CMatrix3DHT() { }

  explicit CMatrix3DHT(Type type) :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
    if (type == CMATRIX_3DH_IDENTITY)
      setIdentity();
    else
      CTHROW("Bad Matrix Type");
  }

  CMatrix3DHT(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) :
   m00_(m00), m01_(m01), m02_(m02), m03_(0),
   m10_(m10), m11_(m11), m12_(m12), m13_(0),
   m20_(m20), m21_(m21), m22_(m22), m23_(0),
   m30_(0  ), m31_(0  ), m32_(0  ), m33_(0) {
    setOuterIdentity();
  }

  CMatrix3DHT(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22, T tx, T ty, T tz) :
   m00_(m00), m01_(m01), m02_(m02), m03_(tx),
   m10_(m10), m11_(m11), m12_(m12), m13_(ty),
   m20_(m20), m21_(m21), m22_(m22), m23_(tz),
   m30_(0  ), m31_(0  ), m32_(0  ), m33_(0 ) {
    setBottomIdentity();
  }

  CMatrix3DHT(T m00, T m01, T m02, T m03, T m10, T m11, T m12, T m13,
              T m20, T m21, T m22, T m23, T m30, T m31, T m32, T m33) :
   m00_(m00), m01_(m01), m02_(m02), m03_(m03),
   m10_(m10), m11_(m11), m12_(m12), m13_(m13),
   m20_(m20), m21_(m21), m22_(m22), m23_(m23),
   m30_(m30), m31_(m31), m32_(m32), m33_(m33) {
  }

  CMatrix3DHT(const T *m, uint n) :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
    if      (n == 9) {
      setValues(m[ 0], m[ 1], m[ 2],
                m[ 3], m[ 4], m[ 5],
                m[ 6], m[ 7], m[ 8]);

      setOuterIdentity();
    }
    else if (n == 12) {
      setValues(m[ 0], m[ 1], m[ 2],
                m[ 3], m[ 4], m[ 5],
                m[ 6], m[ 7], m[ 8],
                m[ 9], m[10], m[11]);

      setBottomIdentity();
    }
    else if (n == 16)
      setValues(m[ 0], m[ 1], m[ 2], m[ 3],
                m[ 4], m[ 5], m[ 6], m[ 7],
                m[ 8], m[ 9], m[10], m[11],
                m[12], m[13], m[14], m[15]);
    else
      CTHROW("Invalid size");
  }

  CMatrix3DHT(CMathGen::AxisType3D axis, T angle,
              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
    setRotation(axis, angle, handedness);
  }

  //------

  // copy operations
  CMatrix3DHT(const MatrixH &a) :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
    memcpy(&m00_, &a.m00_, 16*sizeof(T));
  }

  MatrixH &operator=(const MatrixH &a) {
    memcpy(&m00_, &a.m00_, 16*sizeof(T));

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "((" << m00_ << "," << m01_ << "," <<
                  m02_ << "," << m03_ << "),";
    os << " (" << m10_ << "," << m11_ << "," <<
                  m12_ << "," << m13_ << "),";
    os << " (" << m20_ << "," << m21_ << "," <<
                  m22_ << "," << m23_ << "),";
    os << " (" << m30_ << "," << m31_ << "," <<
                  m32_ << "," << m33_ << "))";
  }

  friend std::ostream &operator<<(std::ostream &os, const MatrixH &matrix) {
    matrix.print(os);

    return os;
  }

  //------

  // comparison
  int cmp(const MatrixH &v) const {
    if      (m00_ < v.m00_) return -1;
    else if (m00_ > v.m00_) return  1;
    else if (m01_ < v.m01_) return -1;
    else if (m01_ > v.m01_) return  1;
    else if (m02_ < v.m02_) return -1;
    else if (m02_ > v.m02_) return  1;
    else if (m03_ < v.m03_) return -1;
    else if (m03_ > v.m03_) return  1;

    else if (m10_ < v.m10_) return -1;
    else if (m10_ > v.m10_) return  1;
    else if (m11_ < v.m11_) return -1;
    else if (m11_ > v.m11_) return  1;
    else if (m12_ < v.m12_) return -1;
    else if (m12_ > v.m12_) return  1;
    else if (m13_ < v.m13_) return -1;
    else if (m13_ > v.m13_) return  1;

    else if (m20_ < v.m20_) return -1;
    else if (m20_ > v.m20_) return  1;
    else if (m21_ < v.m21_) return -1;
    else if (m21_ > v.m21_) return  1;
    else if (m22_ < v.m22_) return -1;
    else if (m22_ > v.m22_) return  1;
    else if (m23_ < v.m23_) return -1;
    else if (m23_ > v.m23_) return  1;

    else if (m30_ < v.m30_) return -1;
    else if (m30_ > v.m30_) return  1;
    else if (m31_ < v.m31_) return -1;
    else if (m31_ > v.m31_) return  1;
    else if (m32_ < v.m32_) return -1;
    else if (m32_ > v.m32_) return  1;
    else if (m33_ < v.m33_) return -1;
    else if (m33_ > v.m33_) return  1;

    else                    return  0;
  }

  friend bool operator==(const MatrixH &lhs, const MatrixH &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const MatrixH &lhs, const MatrixH &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const MatrixH &lhs, const MatrixH &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const MatrixH &lhs, const MatrixH &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const MatrixH &lhs, const MatrixH &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const MatrixH &lhs, const MatrixH &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  // TODO: isZero(), isIdentity()

  //------

  bool isIdentity() {
    return (m00_ == 1 && m01_ == 0 && m02_ == 0 && m03_ == 0 &&
            m10_ == 0 && m11_ == 1 && m12_ == 0 && m13_ == 0 &&
            m20_ == 0 && m21_ == 0 && m22_ == 1 && m23_ == 0 &&
            m30_ == 0 && m31_ == 0 && m32_ == 0 && m33_ == 1);
  }

  //------

  Matrix3 getMatrix() {
    return Matrix3(m00_, m01_, m02_,
                   m10_, m11_, m12_,
                   m20_, m21_, m22_,
                   m03_, m13_, m23_);
  }

  //------

  void setIdentity() {
    setInnerIdentity();

    setOuterIdentity();
  }

  static CMatrix3DHT translation(T tx, T ty, T tz) {
    CMatrix3DHT m;

    m.setTranslation(tx, ty, tz);

    return m;
  }

  void setTranslation(T tx, T ty, T tz) {
    setInnerIdentity();

    setOuterTranslate(tx, ty, tz);
  }

  void setTranslation(const Point &point) {
    setInnerIdentity();

    setOuterTranslate(point.x, point.y, point.z);
  }

  void setTranslation(const Vector &vector) {
    setInnerIdentity();

    setOuterTranslate(vector.getX(), vector.getY(), vector.getZ());
  }

  static CMatrix3DHT scale(T sx, T sy, T sz) {
    CMatrix3DHT m;

    m.setScale(sx, sy, sz);

    return m;
  }

  void setScale(T s) {
    setInnerScale(s, s, s);

    setOuterIdentity();
  }

  void setScale(T sx, T sy, T sz) {
    setInnerScale(sx, sy, sz);

    setOuterIdentity();
  }

  void setScaleTranslation(T s, T tx, T ty, T tz) {
    setInnerScale(s, s, s);

    setOuterTranslate(tx, ty, tz);
  }

  void setScaleTranslation(T sx, T sy, T sz, T tx, T ty, T tz) {
    setInnerScale(sx, sy, sz);

    setOuterTranslate(tx, ty, tz);
  }

  void setRotation(CMathGen::AxisType3D axis, T angle,
                   CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    if (handedness == CMathGen::RIGHT_HANDEDNESS)
      setInnerRotationRHS(axis, angle);
    else
      setInnerRotationLHS(axis, angle);

    setOuterIdentity();
  }

  void setRotationTranslation(CMathGen::AxisType3D axis, T angle, T tx, T ty, T tz,
                              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    if (handedness == CMathGen::RIGHT_HANDEDNESS)
      setInnerRotationRHS(axis, angle);
    else
      setInnerRotationLHS(axis, angle);

    setOuterTranslate(tx, ty, tz);
  }

  void setXYZRotation(T x_angle, T y_angle, T z_angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    MatrixH xmatrix(CMathGen::X_AXIS_3D, x_angle, handedness);
    MatrixH ymatrix(CMathGen::Y_AXIS_3D, y_angle, handedness);
    MatrixH zmatrix(CMathGen::Z_AXIS_3D, z_angle, handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setXYZRotation(const Vector &angles,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    MatrixH xmatrix(CMathGen::X_AXIS_3D, angles.getX(), handedness);
    MatrixH ymatrix(CMathGen::Y_AXIS_3D, angles.getY(), handedness);
    MatrixH zmatrix(CMathGen::Z_AXIS_3D, angles.getZ(), handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setGenRotation(T x1, T y1, T z1, T x2, T y2, T z2, T angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    MatrixH matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7;

    matrix1.setTranslation(-x1, -y1, -z1);
    matrix2.setTranslation( x1,  y1,  z1);

    T theta = atan2(y2 - y1, x2 - x1);

    matrix3.setRotation(CMathGen::Z_AXIS_3D,  theta, handedness);
    matrix4.setRotation(CMathGen::Z_AXIS_3D, -theta, handedness);

    T v = ::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));

    T beta = atan2(v, z2 - z1);

    matrix5.setRotation(CMathGen::Y_AXIS_3D,  beta, handedness);
    matrix6.setRotation(CMathGen::Y_AXIS_3D, -beta, handedness);

    matrix7.setRotation(CMathGen::Z_AXIS_3D, angle, handedness);

    *this = matrix2*(matrix4*(matrix6*(matrix7*(matrix5*(matrix3*matrix1)))));
  }

  void setGenRotation(const Vector &axis, T angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    setOuterIdentity();

    Vector a = axis.normalized();

    T c = ::cos(angle);
    T s = ::sin(angle);

    T c1 = 1.0 - c;

    T axx = a.getX()*a.getX();
    T axy = a.getX()*a.getY();
    T axz = a.getX()*a.getZ();
    T ayy = a.getY()*a.getY();
    T ayz = a.getY()*a.getZ();
    T azz = a.getZ()*a.getZ();

    T axs = a.getX()*s; T ays = a.getY()*s; T azs = a.getZ()*s;

    if (handedness == CMathGen::RIGHT_HANDEDNESS) {
      m00_ = axx*c1 + c  ; m01_ = axy*c1 + azs; m02_ = axz*c1 - ays;
      m10_ = axy*c1 - azs; m11_ = ayy*c1 + c  ; m12_ = ayz*c1 + axs;
      m20_ = axz*c1 + ays; m21_ = ayz*c1 - axs; m22_ = azz*c1 + c  ;
    }
    else {
      m00_ = axx*c1 + c  ; m01_ = axy*c1 - azs; m02_ = axz*c1 + ays;
      m10_ = axy*c1 + azs; m11_ = ayy*c1 + c  ; m12_ = ayz*c1 - axs;
      m20_ = axz*c1 - ays; m21_ = ayz*c1 + axs; m22_ = azz*c1 + c  ;
    }
  }

  void setLookAt(const Point &eye, const Point &center) {
    Vector dir(eye, center);

    setLookAt(eye, dir);
  }

  void setLookAt(const Point &eye, const Point &center, const Vector &up) {
    Vector dir(eye, center);

    setLookAt(eye, dir, up);
  }

  void setLookAt(const Point &eye, const Vector &dir) {
    Vector dir1 = dir.normalized();

    Vector up(0,0,1);

    Vector up1 = up - (up.dotProduct(dir1))*dir1;

    setLookAt(eye, dir1, up1);
  }

  void setLookAt(const Point &eye, const Vector &dir, const Vector &up) {
    Vector dir1 = dir.normalized();
    Vector up1  = up .normalized();

    Vector right = dir1 .crossProduct(up1 );
    Vector newUp = right.crossProduct(dir1);

    //dir1 = -dir1;

    setColumn(0, right);
    setColumn(1, newUp);
    setColumn(2, dir1 );

    setOuterTranslate(eye.x, eye.y, eye.z);
  }

  void setEye(T x1, T y1, T z1, T x2, T y2, T z2,
              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    T angle1, angle2, angle3;

    calcEye(x1, y1, z1, x2, y2, z2, &angle1, &angle2, &angle3);

    MatrixH matrix1, matrix2, matrix3, matrix4;

    matrix1.setTranslation(-x1, -y1, -z1);

    matrix2.setRotation(CMathGen::Z_AXIS_3D,  angle1, handedness);
    matrix3.setRotation(CMathGen::Y_AXIS_3D,  angle2, handedness);
    matrix4.setRotation(CMathGen::Z_AXIS_3D, -angle3, handedness);

    *this = matrix4*(matrix3*(matrix2*matrix1));
  }

  static void calcEye(T x1, T y1, T z1, T x2, T y2, T z2, T *angle1, T *angle2, T *angle3) {
    T dx = x2 - x1;
    T dy = y2 - y1;
    T dz = z2 - z1;

    *angle1 = atan2(-dy, -dx);

    T v = ::sqrt(dx*dx + dy*dy);

    *angle2 = atan2(v, -dz);

    T w = ::sqrt(v*v + dz*dz);

    *angle3 = atan2(dy*dz, -dx*w);
  }

  void setValues(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) {
    m00_ = m00, m01_ = m01, m02_ = m02;
    m10_ = m10, m11_ = m11, m12_ = m12;
    m20_ = m20, m21_ = m21, m22_ = m22;

    setOuterIdentity();
  }

  void setValues(T m00, T m01, T m02, T m10, T m11, T m12,
                 T m20, T m21, T m22, T tx , T ty , T tz ) {
    m00_ = m00, m01_ = m01, m02_ = m02, m03_ = tx;
    m10_ = m10, m11_ = m11, m12_ = m12, m13_ = ty;
    m20_ = m20, m21_ = m21, m22_ = m22, m23_ = tz;

    setBottomIdentity();
  }

  void setValues(T m00, T m01, T m02, T m03, T m10, T m11, T m12, T m13,
                 T m20, T m21, T m22, T m23, T m30, T m31, T m32, T m33) {
    m00_ = m00; m01_ = m01; m02_ = m02; m03_ = m03;
    m10_ = m10; m11_ = m11; m12_ = m12; m13_ = m13;
    m20_ = m20; m21_ = m21; m22_ = m22; m23_ = m23;
    m30_ = m30; m31_ = m31; m32_ = m32; m33_ = m33;
  }

  void getValues(T *m00, T *m01, T *m02,
                 T *m10, T *m11, T *m12,
                 T *m20, T *m21, T *m22) const {
    if (m00) *m00 = m00_; if (m01) *m01 = m01_; if (m02) *m02 = m02_;
    if (m10) *m10 = m10_; if (m11) *m11 = m11_; if (m12) *m12 = m12_;
    if (m20) *m20 = m20_; if (m21) *m21 = m21_; if (m22) *m22 = m22_;
  }

  void getValues(T *m00, T *m01, T *m02,
                 T *m10, T *m11, T *m12,
                 T *m20, T *m21, T *m22,
                 T *tx , T *ty , T *tz ) const {
    if (m00) *m00 = m00_; if (m01) *m01 = m01_; if (m02) *m02 = m02_;
    if (m10) *m10 = m10_; if (m11) *m11 = m11_; if (m12) *m12 = m12_;
    if (m20) *m20 = m20_; if (m21) *m21 = m21_; if (m22) *m22 = m22_;

    if (tx) *tx = m03_ ; if (ty) *ty = m13_ ; if (tz) *tz = m23_ ;
  }

  void getValues(T *m00, T *m01, T *m02, T *m03,
                 T *m10, T *m11, T *m12, T *m13,
                 T *m20, T *m21, T *m22, T *m23,
                 T *m30, T *m31, T *m32, T *m33) const {
    if (m00) *m00 = m00_; if (m01) *m01 = m01_;
    if (m02) *m02 = m02_; if (m03) *m03 = m03_;
    if (m10) *m10 = m10_; if (m11) *m11 = m11_;
    if (m12) *m12 = m12_; if (m13) *m13 = m13_;
    if (m20) *m20 = m20_; if (m21) *m21 = m21_;
    if (m22) *m22 = m22_; if (m23) *m23 = m23_;
    if (m30) *m30 = m30_; if (m31) *m31 = m31_;
    if (m32) *m32 = m32_; if (m33) *m33 = m33_;
  }

  void getValues(T *v, int n) const {
    if      (n == 9) {
      v[0] = m00_; v[1] = m01_; v[2] = m02_;
      v[3] = m10_; v[4] = m11_; v[5] = m12_;
      v[6] = m20_; v[7] = m21_; v[8] = m22_;
    }
    else if (n == 12) {
      v[ 0] = m00_; v[ 1] = m01_; v[ 2] = m02_;
      v[ 3] = m10_; v[ 4] = m11_; v[ 5] = m12_;
      v[ 6] = m20_; v[ 7] = m21_; v[ 8] = m22_;
      v[ 9] = m03_; v[10] = m13_; v[11] = m23_;
    }
    else if (n == 16) {
      v[ 0] = m00_; v[ 1] = m01_; v[ 2] = m02_; v[ 3] = m03_;
      v[ 4] = m10_; v[ 5] = m11_; v[ 6] = m12_; v[ 7] = m13_;
      v[ 8] = m20_; v[ 9] = m21_; v[10] = m22_; v[11] = m23_;
      v[12] = m30_; v[13] = m31_; v[14] = m32_; v[15] = m33_;
    }
    else
      CTHROW("Invalid size");
  }

  //---------

  const T *getData() const { return &m00_; }

  void setColumn(uint c, T x, T y, T z) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    m[0] = x, m[4] = y, m[8] = z;
  }

  void setColumn(uint c, T x, T y, T z, T w) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    m[0] = x, m[4] = y, m[8] = z, m[12] = w;
  }

  void setColumn(uint c, const Vector &vector) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    vector.getXYZ(&m[0], &m[4], &m[8]);
  }

  void getColumn(uint c, T *x, T *y, T *z) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    if (x) *x = m[0];
    if (y) *y = m[4];
    if (z) *z = m[8];
  }

  void getColumn(uint c, T *x, T *y, T *z, T *w) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    if (x) *x = m[ 0];
    if (y) *y = m[ 4];
    if (z) *z = m[ 8];
    if (w) *w = m[12];
  }

  //---------

  void setRow(uint r, T x, T y, T z) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    m[0] = x, m[1] = y, m[2] = z;
  }

  void setRow(uint r, T x, T y, T z, T w) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    m[0] = x, m[1] = y, m[2] = z, m[3] = w;
  }

  void setRow(uint r, const Vector &vector) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    vector.getXYZ(&m[0], &m[1], &m[2]);
  }

  void getRow(uint r, T *x, T *y, T *z) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    if (x) *x = m[0];
    if (y) *y = m[1];
    if (z) *z = m[2];
  }

  void getRow(uint r, T *x, T *y, T *z, T *w) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    if (x) *x = m[0];
    if (y) *y = m[1];
    if (z) *z = m[2];
    if (w) *w = m[3];
  }

  //---------

  void multiplyPoint(T xi, T yi, T zi, T *xo, T *yo, T *zo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi + m03_;
    *yo = m10_*xi + m11_*yi + m12_*zi + m13_;
    *zo = m20_*xi + m21_*yi + m22_*zi + m23_;

    T w = m30_*xi + m31_*yi + m32_*zi + m33_;

    if (w > 1E-6) {
      T iw = 1.0/w;

      *xo *= iw;
      *yo *= iw;
      *zo *= iw;
    }
  }

  void multiplyPoint(T xi, T yi, T zi, T wi, T *xo, T *yo, T *zo, T *wo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi + m03_*wi;
    *yo = m10_*xi + m11_*yi + m12_*zi + m13_*wi;
    *zo = m20_*xi + m21_*yi + m22_*zi + m23_*wi;
    *wo = m30_*xi + m31_*yi + m32_*zi + m33_*wi;
  }

  void multiplyPoint(const Point &ipoint, Point &opoint) const {
    T x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z + m03_;
    T y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z + m13_;
    T z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z + m23_;

    T w = m30_*ipoint.x + m31_*ipoint.y + m32_*ipoint.z + m33_;

    T iw = (fabs(w) > 1E-6 ? 1.0/w : 1.0);

    opoint.x = x*iw;
    opoint.y = y*iw;
    opoint.z = z*iw;
  }

  Point multiplyPoint(const Point &ipoint) const {
    Point opoint;

    T x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z + m03_;
    T y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z + m13_;
    T z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z + m23_;

    T w = m30_*ipoint.x + m31_*ipoint.y + m32_*ipoint.z + m33_;

    T iw = (fabs(w) > 1E-6 ? 1.0/w : 1.0);

    opoint.x = x*iw;
    opoint.y = y*iw;
    opoint.z = z*iw;

    return opoint;
  }

  void multiplyVector(const Vector &ivector, Vector &ovector) const {
    T ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    T ox = m00_*ix + m01_*iy + m02_*iz + m03_;
    T oy = m10_*ix + m11_*iy + m12_*iz + m13_;
    T oz = m20_*ix + m21_*iy + m22_*iz + m23_;

    ovector.setXYZ(ox, oy, oz);
  }

  Vector multiplyVector(const Vector &ivector) const {
    Vector ovector;

    T ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    T ox = m00_*ix + m01_*iy + m02_*iz + m03_;
    T oy = m10_*ix + m11_*iy + m12_*iz + m13_;
    T oz = m20_*ix + m21_*iy + m22_*iz + m23_;

    ovector.setXYZ(ox, oy, oz);

    return ovector;
  }

  void preMultiplyPoint(T xi, T yi, T zi, T *xo, T *yo, T *zo) const {
    *xo = m00_*xi + m10_*yi + m20_*zi;
    *yo = m01_*xi + m11_*yi + m21_*zi;
    *zo = m02_*xi + m12_*yi + m22_*zi;

    T w = m03_*xi + m13_*yi + m23_*zi + m33_;

    T iw = (fabs(w) > 1E-6 ? 1.0/w : 1.0);

    *xo *= iw;
    *yo *= iw;
    *xo *= iw;
  }

  void preMultiplyPoint(T xi, T yi, T zi, T wi, T *xo, T *yo, T *zo, T *wo) const {
    *xo = m00_*xi + m10_*yi + m20_*zi + m30_*wi;
    *yo = m01_*xi + m11_*yi + m21_*zi + m31_*wi;
    *zo = m02_*xi + m12_*yi + m22_*zi + m32_*wi;
    *wo = m03_*xi + m13_*yi + m23_*zi + m33_*wi;
  }

  void preMultiplyPoint(const Point &ipoint, Point &opoint) const {
    opoint.x = m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z;
    opoint.z = m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z;

    T w = m03_*ipoint.x + m13_*ipoint.y + m23_*ipoint.z + m33_;

    T iw = (fabs(w) > 1E-6 ? 1.0/w : 1.0);

    opoint *= iw;
  }

  void preMultiplyVector(const Vector &ivector, Vector &ovector) const {
    T ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    T ox = m00_*ix + m10_*iy + m20_*iz;
    T oy = m01_*ix + m11_*iy + m21_*iz;
    T oz = m02_*ix + m12_*iy + m22_*iz;

    ovector.setXYZ(ox, oy, oz);
  }

  MatrixH &translate(T x, T y, T z) {
    m03_ += x;
    m13_ += y;
    m23_ += z;

    return *this;
  }

  //------

  void transpose() {
    std::swap(m10_, m01_);
    std::swap(m20_, m02_);
    std::swap(m21_, m12_);
    std::swap(m30_, m03_);
    std::swap(m31_, m13_);
    std::swap(m32_, m23_);
  }

  MatrixH transposed() const {
    MatrixH matrix = *this;

    matrix.transpose();

    return matrix;
  }

  //------

  bool invert(MatrixH &imatrix) const {
    T d = determinant();

    if (::fabs(d) == 0.0)
      return false;

    T id = 1.0/d;

    imatrix.m00_ =  id*det3x3(m11_, m12_, m13_,
                              m21_, m22_, m23_,
                              m31_, m32_, m33_);
    imatrix.m10_ = -id*det3x3(m10_, m12_, m13_,
                              m20_, m22_, m23_,
                              m30_, m32_, m33_);
    imatrix.m20_ =  id*det3x3(m10_, m11_, m13_,
                              m20_, m21_, m23_,
                              m30_, m31_, m33_);
    imatrix.m30_ = -id*det3x3(m10_, m11_, m12_,
                              m20_, m21_, m22_,
                              m30_, m31_, m32_);

    imatrix.m01_ = -id*det3x3(m01_, m02_, m03_,
                              m21_, m22_, m23_,
                              m31_, m32_, m33_);
    imatrix.m11_ =  id*det3x3(m00_, m02_, m03_,
                              m20_, m22_, m23_,
                              m30_, m32_, m33_);
    imatrix.m21_ = -id*det3x3(m00_, m01_, m03_,
                              m20_, m21_, m23_,
                              m30_, m31_, m33_);
    imatrix.m31_ =  id*det3x3(m00_, m01_, m02_,
                              m20_, m21_, m22_,
                              m30_, m31_, m32_);

    imatrix.m02_ =  id*det3x3(m01_, m02_, m03_,
                              m11_, m12_, m13_,
                              m31_, m32_, m33_);
    imatrix.m12_ = -id*det3x3(m00_, m02_, m03_,
                              m10_, m12_, m13_,
                              m30_, m32_, m33_);
    imatrix.m22_ =  id*det3x3(m00_, m01_, m03_,
                              m10_, m11_, m13_,
                              m30_, m31_, m33_);
    imatrix.m32_ = -id*det3x3(m00_, m01_, m02_,
                              m10_, m11_, m12_,
                              m30_, m31_, m32_);

    imatrix.m03_ = -id*det3x3(m01_, m02_, m03_,
                              m11_, m12_, m13_,
                              m21_, m22_, m23_);
    imatrix.m13_ =  id*det3x3(m00_, m02_, m03_,
                              m10_, m12_, m13_,
                              m20_, m22_, m23_);
    imatrix.m23_ = -id*det3x3(m00_, m01_, m03_,
                              m10_, m11_, m13_,
                              m20_, m21_, m23_);
    imatrix.m33_ =  id*det3x3(m00_, m01_, m02_,
                              m10_, m11_, m12_,
                              m20_, m21_, m22_);

    return true;
  }

  MatrixH inverse() const {
    MatrixH imatrix;

    if (! invert(imatrix))
      CTHROW("Divide by zero");

    return imatrix;
  }

  T determinant() const {
    return
      (m00_*det3x3(m11_, m12_, m13_, m21_, m22_, m23_, m31_, m32_, m33_) -
       m01_*det3x3(m10_, m12_, m13_, m20_, m22_, m23_, m30_, m32_, m33_) +
       m02_*det3x3(m10_, m11_, m13_, m20_, m21_, m23_, m30_, m31_, m33_) -
       m03_*det3x3(m10_, m11_, m12_, m20_, m21_, m22_, m30_, m31_, m32_));
  }

  void normalize() {
    T d = determinant();

    T id = 1.0/d;

    for (int i = 0; i < 9; ++i)
      (&m00_)[i] *= id;
  }

  //------

  bool affineInvert(MatrixH &imatrix) const {
    T d = affineDeterminant();

    if (::fabs(d) == 0.0)
      return false;

    T id = 1.0/d;

    imatrix.m00_ =  id*det2x2(m11_, m12_, m21_, m22_);
    imatrix.m10_ = -id*det2x2(m10_, m12_, m20_, m22_);
    imatrix.m20_ =  id*det2x2(m10_, m11_, m20_, m21_);

    imatrix.m01_ = -id*det2x2(m01_, m02_, m21_, m22_);
    imatrix.m11_ =  id*det2x2(m00_, m02_, m20_, m22_);
    imatrix.m21_ = -id*det2x2(m00_, m01_, m20_, m21_);

    imatrix.m02_ =  id*det2x2(m01_, m02_, m11_, m12_);
    imatrix.m12_ = -id*det2x2(m00_, m02_, m10_, m12_);
    imatrix.m22_ =  id*det2x2(m00_, m01_, m10_, m11_);

    imatrix.m03_ = -id*det3x3(m01_, m02_, m03_,
                              m11_, m12_, m13_,
                              m21_, m22_, m23_);
    imatrix.m13_ =  id*det3x3(m00_, m02_, m03_,
                              m10_, m12_, m13_,
                              m20_, m22_, m23_);
    imatrix.m23_ = -id*det3x3(m00_, m01_, m03_,
                              m10_, m11_, m13_,
                              m20_, m21_, m23_);

    imatrix.m30_ = 0;
    imatrix.m31_ = 0;
    imatrix.m32_ = 0;
    imatrix.m33_ = 1;

    return true;
  }

  MatrixH affineInverse() const {
    MatrixH imatrix;

    if (! affineInvert(imatrix))
      CTHROW("Divide by zero");

    return imatrix;
  }

  T affineDeterminant() const {
    return (m00_*det2x2(m11_, m12_, m21_, m22_) -
            m01_*det2x2(m10_, m12_, m20_, m22_) +
            m02_*det2x2(m10_, m11_, m20_, m21_));
  }

  void affineNormalize() {
    T d = affineDeterminant();

    T id = 1.0/d;

    for (int i = 0; i < 9; ++i)
      (&m00_)[i] *= id;
  }

  //------

  void setTransform(T xmin1, T ymin1, T zmin1, T xmax1, T ymax1, T zmax1,
                    T xmin2, T ymin2, T zmin2, T xmax2, T ymax2, T zmax2) {
    T sx = (xmax2 - xmin2)/(xmax1 - xmin1);
    T sy = (ymax2 - ymin2)/(ymax1 - ymin1);
    T sz = (zmax2 - zmin2)/(zmax1 - zmin1);

    T tx = -xmin1*sx + xmin2;
    T ty = -ymin1*sy + ymin2;
    T tz = -zmin1*sz + zmin2;

    setInnerScale(sx, sy, sz);

    setOuterTranslate(tx, ty, tz);
  }

  //------

  static MatrixH *newIdentityMatrix() {
    MatrixH *m = new MatrixH();

    m->setIdentity();

    return m;
  }

  //------

  static bool solveAXeqB(const MatrixH &a, Point &x, const Point &b) {
    T det_a = a.determinant();

    if (::fabs(det_a) == 0.0)
      return false;

    T idet_a = 1.0/det_a;

    MatrixH t(a);

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

  //------

  void zero() { memset(&m00_, 0, 16*sizeof(T)); }

  //------

  MatrixH &operator+=(const MatrixH &b) {
    m00_ += b.m00_; m01_ += b.m01_; m02_ += b.m02_; m03_ += b.m03_;
    m10_ += b.m10_; m11_ += b.m11_; m12_ += b.m12_; m13_ += b.m13_;
    m20_ += b.m20_; m21_ += b.m21_; m22_ += b.m22_; m23_ += b.m23_;
    m30_ += b.m30_; m31_ += b.m31_; m32_ += b.m32_; m33_ += b.m33_;

    return *this;
  }

  MatrixH operator+(const MatrixH &b) {
    MatrixH c = *this;

    c += b;

    return c;
  }

  MatrixH &operator-=(const MatrixH &b) {
    m00_ -= b.m00_; m01_ -= b.m01_; m02_ -= b.m02_; m03_ -= b.m03_;
    m10_ -= b.m10_; m11_ -= b.m11_; m12_ -= b.m12_; m13_ -= b.m13_;
    m20_ -= b.m20_; m21_ -= b.m21_; m22_ -= b.m22_; m23_ -= b.m23_;
    m30_ -= b.m30_; m31_ -= b.m31_; m32_ -= b.m32_; m33_ -= b.m33_;

    return *this;
  }

  MatrixH operator-(const MatrixH &b) {
    MatrixH c = *this;

    c -= b;

    return c;
  }

  MatrixH &operator*=(const MatrixH &b) {
    MatrixH a;

    memcpy(&a.m00_, &m00_, 16*sizeof(T));

    m00_ = a.m00_*b.m00_ + a.m01_*b.m10_ + a.m02_*b.m20_ + a.m03_*b.m30_;
    m01_ = a.m00_*b.m01_ + a.m01_*b.m11_ + a.m02_*b.m21_ + a.m03_*b.m31_;
    m02_ = a.m00_*b.m02_ + a.m01_*b.m12_ + a.m02_*b.m22_ + a.m03_*b.m32_;
    m03_ = a.m00_*b.m03_ + a.m01_*b.m13_ + a.m02_*b.m23_ + a.m03_*b.m33_;

    m10_ = a.m10_*b.m00_ + a.m11_*b.m10_ + a.m12_*b.m20_ + a.m13_*b.m30_;
    m11_ = a.m10_*b.m01_ + a.m11_*b.m11_ + a.m12_*b.m21_ + a.m13_*b.m31_;
    m12_ = a.m10_*b.m02_ + a.m11_*b.m12_ + a.m12_*b.m22_ + a.m13_*b.m32_;
    m13_ = a.m10_*b.m03_ + a.m11_*b.m13_ + a.m12_*b.m23_ + a.m13_*b.m33_;

    m20_ = a.m20_*b.m00_ + a.m21_*b.m10_ + a.m22_*b.m20_ + a.m23_*b.m30_;
    m21_ = a.m20_*b.m01_ + a.m21_*b.m11_ + a.m22_*b.m21_ + a.m23_*b.m31_;
    m22_ = a.m20_*b.m02_ + a.m21_*b.m12_ + a.m22_*b.m22_ + a.m23_*b.m32_;
    m23_ = a.m20_*b.m03_ + a.m21_*b.m13_ + a.m22_*b.m23_ + a.m23_*b.m33_;

    m30_ = a.m30_*b.m00_ + a.m31_*b.m10_ + a.m32_*b.m20_ + a.m33_*b.m30_;
    m31_ = a.m30_*b.m01_ + a.m31_*b.m11_ + a.m32_*b.m21_ + a.m33_*b.m31_;
    m32_ = a.m30_*b.m02_ + a.m31_*b.m12_ + a.m32_*b.m22_ + a.m33_*b.m32_;
    m33_ = a.m30_*b.m03_ + a.m31_*b.m13_ + a.m32_*b.m23_ + a.m33_*b.m33_;

    return *this;
  }

  MatrixH operator*(const MatrixH &b) const {
    MatrixH c = *this;

    c *= b;

    return c;
  }

  MatrixH &operator*=(T s) {
    m00_ *= s; m01_ *= s; m02_ *= s; m03_ *= s;
    m10_ *= s; m11_ *= s; m12_ *= s; m13_ *= s;
    m20_ *= s; m21_ *= s; m22_ *= s; m23_ *= s;
    m30_ *= s; m31_ *= s; m32_ *= s; m33_ *= s;

    return *this;
  }

  MatrixH operator*(T s) {
    MatrixH c = *this;

    c *= s;

    return c;
  }

  friend Point operator*(const MatrixH &m, const Point &p) {
    Point p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend Point operator*(const Point &p, const MatrixH &m) {
    Point p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend Vector operator*(const MatrixH &m, const Vector &v) {
    Vector v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend Vector operator*(const Vector &v, const MatrixH &m) {
    Vector v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  MatrixH &operator/=(const MatrixH &b) {
    MatrixH bi;

    if (! b.invert(bi)) {
      CTHROW("Divide by zero");
      return *this;
    }

    return (*this) *= bi;
  }

  MatrixH operator/(const MatrixH &b) {
    MatrixH c = *this;

    c /= b;

    return c;
  }

  //------

  void setValue(uint i, T value) {
    (&m00_)[i] = value;
  }

  void setValue(uint i, uint j, T value) {
    assert(i < 4 && j < 4);

    T *m  = (&m00_)[4*j + i];

    *m = value;
  }

  T getValue(uint i) const {
    assert(i < 16);

    return (&m00_)[i];
  }

  T getValue(uint i, uint j) const {
    assert(i < 4 && j < 4);

    T *m = (&m00_)[4*j + i];

    return *m;
  }

  T operator[](uint i) const {
    assert(i < 16);

    return (&m00_)[i];
  }

  T &operator[](uint i) {
    assert(i < 16);

    return (&m00_)[i];
  }

  T operator()(uint i, uint j) const {
    assert(i < 4 && j < 4);

    T *m = (&m00_)[4*j + i];

    return *m;
  }

  T &operator()(uint i, uint j) {
    assert(i < 4 && j < 4);

    T *m = (&m00_)[4*j + i];

    return *m;
  }

  //------

  /*! build perspective projection matrix
      (near = -ve, far = -ve) look down negative Z-axis
      (near = +ve, far = +ve) look down position Z-axis
  */
  void buildPerspective(T fov, T aspect, T near, T far) {
    // can't have near/far on other sides of origin
    assert(near*far > 0);

    if (near < 0) {
      double tf2 = tan(CMathGen::DegToRad(fov)*0.5);

      double itf2 = 1.0/tf2;

      double npf = near + far;
      double ntf = near * far;

      double d = near - far; double id = 1.0/d;

      T a  = itf2/aspect;
      T e  = itf2;
      T i  = npf*id;
      T tz = 2.0*ntf*id;

      m00_ = a  , m01_ = 0.0, m02_ =  0.0; m03_ = 0.0;
      m10_ = 0.0, m11_ = e  , m12_ =  0.0; m13_ = 0.0;
      m20_ = 0.0, m21_ = 0.0, m22_ =  i  ; m23_ = tz ;
      m30_ = 0.0, m31_ = 0.0, m32_ = -1.0; m33_ = 0.0;
    }
    else {
      double tf2 = tan(CMathGen::DegToRad(fov)*0.5);

      double itf2 = 1.0/tf2;

      double ntf = near * far;

      double d = far - near; double id = 1.0/d;

      T a  = itf2/aspect;
      T e  = itf2;
      T i  = far*id;
      T tz = -ntf*id;

      m00_ = a  , m01_ = 0.0, m02_ = 0.0; m03_ = 0.0;
      m10_ = 0.0, m11_ = e  , m12_ = 0.0; m13_ = 0.0;
      m20_ = 0.0, m21_ = 0.0, m22_ = i  ; m23_ = tz ;
      m30_ = 0.0, m31_ = 0.0, m32_ = 1.0; m33_ = 0.0;
    }
  }

  void buildOrtho(T left, T right, T bottom, T top, T near, T far) {
    T w = right - left  ; T iw = 1.0/w;
    T h = top   - bottom; T ih = 1.0/h;
    T d = near  - far   ; T id = 1.0/d; // far - near ?

    T rpl = right + left  ;
    T tpb = top   + bottom;
    T npf = near  + far   ;

    T a = 2.0*iw;
    T e = 2.0*ih;
    T i = 2.0*id;

    T tx = -rpl*iw;
    T ty = -tpb*ih;
    T tz =  npf*id;

    setInnerScale(a, e, i);

    setOuterTranslate(tx, ty, tz);
  }

  void buildFrustrum(T left, T right, T bottom, T top, T near, T far) {
    T w = right - left  ; T iw = 1.0/w;
    T h = top   - bottom; T ih = 1.0/h;
    T d = near  - far   ; T id = 1.0/d; // far - near ?

    T npf = near + far;
    T ntf = near * far;

    T a  = 2.0*near*iw;
    T e  = 2.0*near*ih;
    T i  = npf*id;

    T c  = (right + left  )*iw;
    T f  = (top   + bottom)*ih;
    T tz = 2.0*ntf*id;

    m00_ = a  , m01_ = 0.0, m02_ =  c  ; m03_ = 0.0;
    m10_ = 0.0, m11_ = e  , m12_ =  f  ; m13_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ =  i  ; m23_ = tz ;
    m30_ = 0.0, m31_ = 0.0, m32_ = -1.0; m33_ = 0.0;
  }

  //------

 private:
  void setInnerIdentity() {
    m00_ = 1.0, m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = 1.0, m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = 1.0;
  }

  void setInnerScale(T sx, T sy, T sz) {
    m00_ = sx , m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = sy , m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = sz ;
  }

  void setInnerRotationRHS(CMathGen::AxisType3D axis, T angle) {
    T c = ::cos(angle);
    T s = ::sin(angle);

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

  void setInnerRotationLHS(CMathGen::AxisType3D axis, T angle) {
    T c = ::cos(angle);
    T s = ::sin(angle);

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

  void setOuterIdentity() {
    m03_ = 0.0; m13_ = 0.0; m23_ = 0.0;

    setBottomIdentity();
  }

  void setOuterTranslate(T tx, T ty, T tz) {
    m03_ = tx; m13_ = ty; m23_ = tz;

    setBottomIdentity();
  }

  void setBottomIdentity() {
    m30_ = 0.0, m31_ = 0.0, m32_ = 0.0, m33_ = 1.0;
  }

  static T det3x3(T a, T b, T c, T d, T e, T f, T g, T h, T i) {
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
  }

  static T det2x2(T a, T b, T c, T d) {
    return a*c - b*d;
  }
};

typedef CMatrix3DHT<double> CMatrix3DH;
typedef CMatrix3DHT<float>  CMatrix3DHF;

#endif
