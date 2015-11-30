#ifndef CMATRIX_3D_H
#define CMATRIX_3D_H

#include <cassert>
#include <cstring>

#include <CMatrixType.h>
#include <CMathGen.h>
#include <CPoint3D.h>
#include <CVector3D.h>
#include <CThrow.h>

#ifdef GAUSSIAN_MATRIX
#include <CGaussianMatrix.h>
#endif

/* / m00 m01 m02 m03 \ */
/* | m10 m11 m12 m13 | */
/* | m20 m21 m22 m23 | */
/* \ m30 m31 m32 m33 / */

template<typename T>
class CMatrix3DT {
 private:
  typedef CPoint3DT<T>   Point;
  typedef CVector3DT<T>  Vector;
  typedef CMatrix3DT<T>  Matrix;

 private:
  T m00_, m01_, m02_, m03_;
  T m10_, m11_, m12_, m13_;
  T m20_, m21_, m22_, m23_;
  T m30_, m31_, m32_, m33_;

 public:
  // constructor/destructor
  CMatrix3DT() :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
  }

 ~CMatrix3DT() { }

  explicit CMatrix3DT(CMatrixType type) :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
    if (type == CMATRIX_TYPE_IDENTITY)
      setIdentity();
    else
      CTHROW("Bad Matrix Type");
  }

  CMatrix3DT(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) :
   m00_(m00), m01_(m01), m02_(m02), m03_(0),
   m10_(m10), m11_(m11), m12_(m12), m13_(0),
   m20_(m20), m21_(m21), m22_(m22), m23_(0),
   m30_(0  ), m31_(0  ), m32_(0  ), m33_(0) {
    setOuterIdentity();
  }

  CMatrix3DT(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22, T tx, T ty, T tz) :
   m00_(m00), m01_(m01), m02_(m02), m03_(tx),
   m10_(m10), m11_(m11), m12_(m12), m13_(ty),
   m20_(m20), m21_(m21), m22_(m22), m23_(tz),
   m30_(0  ), m31_(0  ), m32_(0  ), m33_(0 ) {
    setBottomIdentity();
  }

  CMatrix3DT(T m00, T m01, T m02, T m03, T m10, T m11, T m12, T m13,
             T m20, T m21, T m22, T m23, T m30, T m31, T m32, T m33) :
   m00_(m00), m01_(m01), m02_(m02), m03_(m03),
   m10_(m10), m11_(m11), m12_(m12), m13_(m13),
   m20_(m20), m21_(m21), m22_(m22), m23_(m23),
   m30_(m30), m31_(m31), m32_(m32), m33_(m33) {
  }

  CMatrix3DT(const T *m, uint n) :
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

  CMatrix3DT(const Vector &v0, const Vector &v1, const Vector &v2) :
   m00_(v0.getX()), m01_(v1.getX()), m02_(v2.getX()), m03_(0),
   m10_(v0.getY()), m11_(v1.getY()), m12_(v2.getY()), m13_(0),
   m20_(v0.getZ()), m21_(v1.getZ()), m22_(v2.getZ()), m23_(0),
   m30_(0        ), m31_(0        ), m32_(0        ), m33_(0) {
    setOuterIdentity();
  }

  CMatrix3DT(CMathGen::AxisType3D axis, T angle,
             CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    setRotation(axis, angle, handedness);
  }

  CMatrix3DT *dup() const {
    return new CMatrix3DT(*this);
  }

  //------

  // copy operations
  CMatrix3DT(const Matrix &a) :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
    memcpy(&m00_, &a.m00_, 16*sizeof(T));
  }

  Matrix &operator=(const Matrix &a) {
    memcpy(&m00_, &a.m00_, 16*sizeof(T));

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "((" << m00_ << "," << m01_ << "," << m02_ << "," << m03_ << "),";
    os << " (" << m10_ << "," << m11_ << "," << m12_ << "," << m13_ << "),";
    os << " (" << m20_ << "," << m21_ << "," << m22_ << "," << m23_ << "),";
    os << " (" << m30_ << "," << m31_ << "," << m32_ << "," << m33_ << "))";
  }

  friend std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
    matrix.print(os);

    return os;
  }

  //------

  // comparison
  int cmp(const Matrix &v) const {
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

  friend bool operator==(const Matrix &lhs, const Matrix &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const Matrix &lhs, const Matrix &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const Matrix &lhs, const Matrix &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const Matrix &lhs, const Matrix &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const Matrix &lhs, const Matrix &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const Matrix &lhs, const Matrix &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  // TODO: isZero(), isIdentity()

  //------

  void setIdentity() {
    setInnerIdentity();

    setOuterIdentity();
  }

  //----------

  static CMatrix3DT translation(T tx, T ty, T tz) {
    CMatrix3DT m;

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

  Matrix &translate(T x, T y, T z) {
    m03_ += x;
    m13_ += y;
    m23_ += z;

    return *this;
  }

  void getTranslate(T *tx, T *ty, T *tz) {
    *tx = m03_;
    *ty = m13_;
    *tz = m23_;
  }

  //----------

  static CMatrix3DT scale(T s) {
    return scale(s, s, s);
  }

  static CMatrix3DT scale(T sx, T sy, T sz) {
    CMatrix3DT m;

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

  void getScale(T *sx, T *sy, T *sz) {
    *sx = m00_;
    *sy = m11_;
    *sz = m22_;
  }

  //----------

  void setScaleTranslation(T s, T tx, T ty, T tz) {
    setInnerScale(s, s, s);

    setOuterTranslate(tx, ty, tz);
  }

  void setScaleTranslation(T sx, T sy, T sz, T tx, T ty, T tz) {
    setInnerScale(sx, sy, sz);

    setOuterTranslate(tx, ty, tz);
  }

  //----------

  static CMatrix3DT rotation(T rx, T ry, T rz) {
    CMatrix3DT m;

    m.setXYZRotation(rx, ry, rz);

    return m;
  }

  static CMatrix3DT rotation(CMathGen::AxisType3D axis, T angle,
                             CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3DT m;

    m.setRotation(axis, angle, handedness);

    return m;
  }

  void setRotation(CMathGen::AxisType3D axis, T angle,
                   CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    if (handedness == CMathGen::RIGHT_HANDEDNESS)
      setInnerRotationRHS(axis, angle);
    else
      setInnerRotationLHS(axis, angle);

    setOuterIdentity();
  }

  void setRotation(T theta, const Vector &u) {
    T theta2 = 0.5*theta;

    // rotate around the line
    T w = cos(theta2);

    //TODO: shouldn't have to normalize u
    Vector v = u*sin(theta2);

    //assign matrix
    T x = v.getX(); T y = v.getY(); T z = v.getZ();

    T x2 = 2.0*x; T y2 = 2.0*y; T z2 = 2.0*z;

    T wx2 = x2*w; T wy2 = y2*w; T wz2 = z2*w;
    T xx2 = x2*x; T xy2 = y2*x; T xz2 = z2*x;
    T yy2 = y2*y; T yz2 = z2*y; T zz2 = z2*z;

    T a = 1.0 - (yy2 + zz2); T b =        xy2 - wz2 ; T c =        xz2 + wy2 ;
    T d =        xy2 + wz2 ; T e = 1.0 - (xx2 + zz2); T f =        yz2 - wx2 ;
    T g =        xz2 - wy2 ; T h =        yz2 + wx2 ; T i = 1.0 - (xx2 + yy2);

    setValues(a, b, c, d, e, f, g, h, i, 0.0, 0.0, 0.0);
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
    Matrix xmatrix(CMathGen::X_AXIS_3D, x_angle, handedness);
    Matrix ymatrix(CMathGen::Y_AXIS_3D, y_angle, handedness);
    Matrix zmatrix(CMathGen::Z_AXIS_3D, z_angle, handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setXYZRotation(const Vector &angles,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    Matrix xmatrix(CMathGen::X_AXIS_3D, angles.getX(), handedness);
    Matrix ymatrix(CMathGen::Y_AXIS_3D, angles.getY(), handedness);
    Matrix zmatrix(CMathGen::Z_AXIS_3D, angles.getZ(), handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setGenRotation(T x1, T y1, T z1, T x2, T y2, T z2, T angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    Matrix matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7;

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

    T axx = a.x*a.x; T axy = a.x*a.y; T axz = a.x*a.z;
                     T ayy = a.y*a.y; T ayz = a.y*a.z;
                                      T azz = a.z*a.z;

    T axs = a.x*s; T ays = a.y*s; T azs = a.z*s;

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

  //----------

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

    dir1 = -dir1;

    setColumn(0, right);
    setColumn(1, newUp);
    setColumn(2, dir1 );

    setOuterTranslate(-eye.x, -eye.y, -eye.z);
  }

  //----------

  void setEye(T x1, T y1, T z1, T x2, T y2, T z2,
              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    T angle1, angle2, angle3;

    calcEye(x1, y1, z1, x2, y2, z2, &angle1, &angle2, &angle3);

    Matrix matrix1, matrix2, matrix3, matrix4;

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

  //----------

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
    m00_ = m00, m01_ = m01, m02_ = m02, m03_ = m03;
    m10_ = m10, m11_ = m11, m12_ = m12, m13_ = m13;
    m20_ = m20, m21_ = m21, m22_ = m22, m23_ = m23;
    m30_ = m30, m31_ = m31, m32_ = m32, m33_ = m33;
  }

  void setValues(const Vector &v0, const Vector &v1, const Vector &v2) {
    m00_ = v0.getX(); m01_ = v1.getX(); m02_ = v2.getX();
    m10_ = v0.getY(); m11_ = v1.getY(); m12_ = v2.getY();
    m20_ = v0.getZ(); m21_ = v1.getZ(); m22_ = v2.getZ();

    setOuterIdentity();
  }

  void getValues(T *m00, T *m01, T *m02, T *m10, T *m11, T *m12, T *m20, T *m21, T *m22) const {
    if (m00) *m00 = m00_; if (m01) *m01 = m01_; if (m02) *m02 = m02_;
    if (m10) *m10 = m10_; if (m11) *m11 = m11_; if (m12) *m12 = m12_;
    if (m20) *m20 = m20_; if (m21) *m21 = m21_; if (m22) *m22 = m22_;
  }

  void getValues(T *m00, T *m01, T *m02, T *m10, T *m11, T *m12,
                 T *m20, T *m21, T *m22, T *tx , T *ty , T *tz ) const {
    if (m00) *m00 = m00_; if (m01) *m01 = m01_; if (m02) *m02 = m02_;
    if (m10) *m10 = m10_; if (m11) *m11 = m11_; if (m12) *m12 = m12_;
    if (m20) *m20 = m20_; if (m21) *m21 = m21_; if (m22) *m22 = m22_;

    if (tx ) *tx  = m03_; if (ty ) *ty  = m13_; if (tz ) *tz  = m23_;
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

  void getValuesI(T *v, int n) const {
    if      (n == 9) {
      v[0] = m00_; v[1] = m10_; v[2] = m20_;
      v[3] = m01_; v[4] = m11_; v[5] = m21_;
      v[6] = m02_; v[7] = m12_; v[8] = m22_;
    }
    else if (n == 12) {
      v[ 0] = m00_; v[ 1] = m10_; v[ 2] = m20_;
      v[ 3] = m01_; v[ 4] = m11_; v[ 5] = m21_;
      v[ 6] = m02_; v[ 7] = m12_; v[ 8] = m22_;
      v[ 9] = m03_; v[10] = m13_; v[11] = m23_;
    }
    else if (n == 16) {
      v[ 0] = m00_; v[ 1] = m10_; v[ 2] = m20_; v[ 3] = m30_;
      v[ 4] = m01_; v[ 5] = m11_; v[ 6] = m21_; v[ 7] = m31_;
      v[ 8] = m02_; v[ 9] = m12_; v[10] = m22_; v[11] = m32_;
      v[12] = m03_; v[13] = m13_; v[14] = m23_; v[15] = m33_;
    }
    else
      CTHROW("Invalid size");
  }

  //---------

  const T *getData() const {
    return &m00_;
  }

  void setColumn(int c, T x, T y, T z) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    m[0] = x, m[4] = y, m[8] = z;
  }

  void setColumn(int c, const Point &point) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    m[0] = point.x, m[4] = point.y, m[8] = point.z;
  }

  void setColumn(int c, const Vector &vector) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    vector.getXYZ(&m[0], &m[4], &m[8]);
  }

  void getColumn(int c, T *x, T *y, T *z) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    if (x) *x = m[0];
    if (y) *y = m[4];
    if (z) *z = m[8];
  }

  void getColumn(int c, Point &point) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    point.x = m[0];
    point.y = m[4];
    point.z = m[8];
  }

  void getColumn(int c, Vector &vector) {
    assert(c < 4);

    T *m = &(&m00_)[c];

    vector = Vector(m[0], m[4], m[8]);
  }

  void getColumns(Vector &u, Vector &v, Vector &w) {
    getColumn(0, u);
    getColumn(1, v);
    getColumn(2, w);
  }

  //---------

  void setRow(int r, T x, T y, T z) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    m[0] = x, m[1] = y, m[2] = z;
  }

  void setRow(int r, const Point &point) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    m[0] = point.x, m[1] = point.y, m[2] = point.z;
  }

  void setRow(int r, const Vector &vector) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    vector.getXYZ(&m[0], &m[1], &m[2]);
  }

  void getRow(int r, T *x, T *y, T *z) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    if (x) *x = m[0];
    if (y) *y = m[1];
    if (z) *z = m[2];
  }

  void getRow(int r, Point &point) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    point.x = m[0];
    point.y = m[1];
    point.z = m[2];
  }

  void getRow(int r, Vector &vector) {
    assert(r < 4);

    T *m = &(&m00_)[r*4];

    vector = Vector(m[0], m[1], m[2]);
  }

  //---------

  // Point can be expressed as a 1x4 matrix (x,y,z,1)
  void multiplyPoint(T xi, T yi, T zi, T *xo, T *yo, T *zo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi + m03_;
    *yo = m10_*xi + m11_*yi + m12_*zi + m13_;
    *zo = m20_*xi + m21_*yi + m22_*zi + m23_;
  }

  void multiplyPoint(const Point &ipoint, Point &opoint) const {
    opoint.x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z + m03_;
    opoint.y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z + m13_;
    opoint.z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z + m23_;
  }

  Point multiplyPoint(const Point &ipoint) const {
    return Point(m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z + m03_,
                 m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z + m13_,
                 m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z + m23_);
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

  Point preMultiplyPoint(const Point &ipoint) const {
    return Point(m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z,
                 m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z,
                 m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z);
  }

  // Vector can be expressed as a 1x4 matrix (x,y,z,0)
  void multiplyVector(T xi, T yi, T zi, T *xo, T *yo, T *zo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi;
    *yo = m10_*xi + m11_*yi + m12_*zi;
    *zo = m20_*xi + m21_*yi + m22_*zi;
  }

  void multiplyVector(const Vector &ivector, Vector &ovector) const {
    Point ipoint, opoint;

    ivector.getXYZ(&ipoint.x, &ipoint.y, &ipoint.z);

    opoint.x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z;
    opoint.y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z;
    opoint.z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z;

    ovector = opoint;
  }

  Vector multiplyVector(const Vector &ivector) const {
    Point ipoint, opoint;

    ivector.getXYZ(&ipoint.x, &ipoint.y, &ipoint.z);

    opoint.x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z;
    opoint.y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z;
    opoint.z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z;

    return Vector(opoint);
  }

  void preMultiplyVector(T xi, T yi, T zi, T *xo, T *yo, T *zo) const {
    *xo = m00_*xi + m10_*yi + m20_*zi;
    *yo = m01_*xi + m11_*yi + m21_*zi;
    *zo = m02_*xi + m12_*yi + m22_*zi;
  }

  void preMultiplyVector(const Vector &ivector, Vector &ovector) const {
    Point ipoint, opoint;

    ivector.getXYZ(&ipoint.x, &ipoint.y, &ipoint.z);

    opoint.x = m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z;
    opoint.z = m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z;

    ovector = opoint;
  }

  Vector preMultiplyVector(const Vector &ivector) const {
    Point ipoint, opoint;

    ivector.getXYZ(&ipoint.x, &ipoint.y, &ipoint.z);

    opoint.x = m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z;
    opoint.z = m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z;

    return Vector(opoint);
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

  Matrix transposed() const {
    Matrix matrix = *this;

    matrix.transpose();

    return matrix;
  }

  //------

  bool invert(Matrix &imatrix) const {
    T d = determinant();

    if (::fabs(d) == 0.0)
      return false;

    T id = 1.0/d;

    imatrix.m00_ = det2x2(m11_, m12_, m21_, m22_)*id;
    imatrix.m10_ = det2x2(m12_, m10_, m22_, m20_)*id;
    imatrix.m20_ = det2x2(m10_, m11_, m20_, m21_)*id;

    imatrix.m01_ = det2x2(m21_, m22_, m01_, m02_)*id;
    imatrix.m11_ = det2x2(m22_, m20_, m02_, m00_)*id;
    imatrix.m21_ = det2x2(m20_, m21_, m00_, m01_)*id;

    imatrix.m02_ = det2x2(m01_, m02_, m11_, m12_)*id;
    imatrix.m12_ = det2x2(m02_, m00_, m12_, m10_)*id;
    imatrix.m22_ = det2x2(m00_, m01_, m10_, m11_)*id;

    T adjoint;

    adjoint = det3x3(m01_, m02_, m03_, m11_, m12_, m13_, m21_, m22_, m23_);

    imatrix.m03_ = -adjoint*id;

    adjoint = det3x3(m02_, m03_, m00_, m12_, m13_, m10_, m22_, m23_, m20_);

    imatrix.m13_ =  adjoint*id;

    adjoint = det3x3(m03_, m00_, m01_, m13_, m10_, m11_, m23_, m20_, m21_);

    imatrix.m23_ = -adjoint*id;

    imatrix.setBottomIdentity();

    return true;
  }

  Matrix inverse() const {
    Matrix imatrix;

    if (! invert(imatrix))
      CTHROW("Divide by 0.0");

    return imatrix;
  }

  //------

  T determinant() const {
    return (m00_*det2x2(m11_, m12_, m21_, m22_) -
            m01_*det2x2(m10_, m12_, m20_, m22_) +
            m02_*det2x2(m10_, m11_, m20_, m21_));
  }

  T upperDeterminant() const {
    return (m00_*det2x2(m11_, m12_, m21_, m22_) -
            m01_*det2x2(m10_, m12_, m20_, m22_) +
            m02_*det2x2(m10_, m11_, m20_, m21_));
  }

  //------

  void normalize() {
    T d = determinant();

    T id = 1.0/d;

    for (int i = 0; i < 9; ++i)
      (&m00_)[i] *= id;
  }

  Matrix normalized() {
    Matrix nmatrix = *this;

    nmatrix.normalize();

    return nmatrix;
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

  static Matrix *newIdentityMatrix() {
    Matrix *m = new Matrix();

    m->setIdentity();

    return m;
  }

  //------

  static bool solveAXeqB(const Matrix &a, Point &x, const Point &b) {
    T det_a = a.determinant();

    if (::fabs(det_a) < 0.0)
      return false;

    T idet_a = 1.0/det_a;

    Matrix t(a);

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

#ifdef GAUSSIAN_MATRIX
  static bool gaussianSolve(const Matrix &m, const Vector &v, Vector &v1) {
    T d[12];

    memset(d, 0, sizeof(d));

    m.getValues(&d[ 0], &d[ 1], &d[ 2], &d[ 4], &d[ 5], &d[ 6], &d[ 8], &d[ 9], &d[10]);

    v.getXYZ(&d[3], &d[7], &d[11]);

    CGaussianMatrix<T,3,4> matrix(d);

    if (! matrix.solve(d))
      return false;

    v1.setXYZ(d[0], d[1], d[2]);

    return true;
  }
#endif

  //------

  void zero() { memset(&m00_, 0, 16*sizeof(T)); }

  //------

  Matrix &operator+=(const Matrix &b) {
    m00_ += b.m00_; m01_ += b.m01_; m02_ += b.m02_; m03_ += b.m03_;
    m10_ += b.m10_; m11_ += b.m11_; m12_ += b.m12_; m13_ += b.m13_;
    m20_ += b.m20_; m21_ += b.m21_; m22_ += b.m22_; m23_ += b.m23_;

    return *this;
  }

  Matrix operator+(const Matrix &b) const {
    Matrix c = *this;

    c += b;

    return c;
  }

  Matrix &operator-=(const Matrix &b) {
    m00_ -= b.m00_; m01_ -= b.m01_, m02_ -= b.m02_, m03_ -= b.m03_;
    m10_ -= b.m10_; m11_ -= b.m11_; m12_ -= b.m12_; m13_ -= b.m13_;
    m20_ -= b.m20_; m21_ -= b.m21_; m22_ -= b.m22_; m23_ -= b.m23_;

    return *this;
  }

  Matrix operator-(const Matrix &b) const {
    Matrix c = *this;

    c -= b;

    return c;
  }

  Matrix &operator*=(const Matrix &b) {
    Matrix a;

    memcpy(&a.m00_, &m00_, 16*sizeof(T));

    m00_ = a.m00_*b.m00_ + a.m01_*b.m10_ + a.m02_*b.m20_;
    m01_ = a.m00_*b.m01_ + a.m01_*b.m11_ + a.m02_*b.m21_;
    m02_ = a.m00_*b.m02_ + a.m01_*b.m12_ + a.m02_*b.m22_;
    m03_ = a.m00_*b.m03_ + a.m01_*b.m13_ + a.m02_*b.m23_ + a.m03_;

    m10_ = a.m10_*b.m00_ + a.m11_*b.m10_ + a.m12_*b.m20_;
    m11_ = a.m10_*b.m01_ + a.m11_*b.m11_ + a.m12_*b.m21_;
    m12_ = a.m10_*b.m02_ + a.m11_*b.m12_ + a.m12_*b.m22_;
    m13_ = a.m10_*b.m03_ + a.m11_*b.m13_ + a.m12_*b.m23_ + a.m13_;

    m20_ = a.m20_*b.m00_ + a.m21_*b.m10_ + a.m22_*b.m20_;
    m21_ = a.m20_*b.m01_ + a.m21_*b.m11_ + a.m22_*b.m21_;
    m22_ = a.m20_*b.m02_ + a.m21_*b.m12_ + a.m22_*b.m22_;
    m23_ = a.m20_*b.m03_ + a.m21_*b.m13_ + a.m22_*b.m23_ + a.m23_;

    return *this;
  }

  Matrix operator*(const Matrix &b) const {
    Matrix c = *this;

    c *= b;

    return c;
  }

  Matrix &operator*=(T s) {
    m00_ *= s; m01_ *= s; m02_ *= s;
    m10_ *= s; m11_ *= s; m12_ *= s;
    m20_ *= s; m21_ *= s; m22_ *= s;

    return *this;
  }

  Matrix operator*(T s) const {
    Matrix c = *this;

    c *= s;

    return c;
  }

  friend Point operator*(const Matrix &m, const Point &p) {
    Point p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend Point operator*(const Point &p, const Matrix &m) {
    Point p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend Vector operator*(const Matrix &m, const Vector &v) {
    Vector v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend Vector operator*(const Vector &v, const Matrix &m) {
    Vector v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  Matrix &operator/=(const Matrix &b) {
    Matrix bi;

    if (! b.invert(bi)) {
      CTHROW("Divide by 0.0");
      return *this;
    }

    return (*this) *= bi;
  }

  Matrix operator/(const Matrix &b) const {
    Matrix c = *this;

    c /= b;

    return c;
  }

  //------

  void setValue(uint i, T value) {
    (&m00_)[i] = value;
  }

  void setValue(uint i, uint j, T value) {
    assert(i < 4 && j < 4);

    T &m = (&m00_)[4*j + i];

    m = value;
  }

  const T &getValue(uint i) const {
    assert(i < 16);

    return (&m00_)[i];
  }

  const T &getValue(uint i, uint j) const {
    assert(i < 4 && j < 4);

    const T &m = (&m00_)[4*j + i];

    return m;
  }

  const T &operator[](uint i) const {
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

  static T det3x3(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) {
    return m00*(m11*m22 - m12*m21) -
           m01*(m10*m22 - m12*m20) +
           m02*(m10*m21 - m11*m20);
  }

  static T det2x2(T m00, T m01, T m10, T m11) {
    return m00*m11 - m01*m10;
  }
};

typedef CMatrix3DT<double> CMatrix3D;
typedef CMatrix3DT<float>  CMatrix3DF;

#endif
