#ifndef CGLMatrix3D_H
#define CGLMatrix3D_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CGLVector3D.h>
#include <cstring>

/* Homogeneous 3D Matrix for OpenGL */

/* / m00 m01 m02 m03 \ */
/* | m10 m11 m12 m13 | */
/* | m20 m21 m22 m23 | */
/* \ m30 m31 m32 m33 / */

class CGLMatrix3D {
 public:
  enum class Type {
    IDENTITY
  };

 public:
  static CGLMatrix3D identity() {
    CGLMatrix3D m; m.setIdentity(); return m;
  }

  //---

  // constructor/destructor
  CGLMatrix3D() { }

 ~CGLMatrix3D() { }

  explicit CGLMatrix3D(Type type) {
    if (type == Type::IDENTITY)
      setIdentity();
    else
      assert(false && "Bad Matrix Type");
  }

  CGLMatrix3D(float m00, float m01, float m02,
              float m10, float m11, float m12,
              float m20, float m21, float m22) :
   m00_(m00), m10_(m10), m20_(m20),
   m01_(m01), m11_(m11), m21_(m21),
   m02_(m02), m12_(m12), m22_(m22) {
    setOuterIdentity();
  }

  CGLMatrix3D(float m00, float m01, float m02,
              float m10, float m11, float m12,
              float m20, float m21, float m22, float tx, float ty, float tz) :
   m00_(m00), m10_(m10), m20_(m20),
   m01_(m01), m11_(m11), m21_(m21),
   m02_(m02), m12_(m12), m22_(m22),
   m03_(tx), m13_(ty), m23_(tz) {
    setBottomIdentity();
  }

  CGLMatrix3D(float m00, float m01, float m02, float m03,
              float m10, float m11, float m12, float m13,
              float m20, float m21, float m22, float m23,
              float m30, float m31, float m32, float m33) :
   m00_(m00), m10_(m10), m20_(m20), m30_(m30),
   m01_(m01), m11_(m11), m21_(m21), m31_(m31),
   m02_(m02), m12_(m12), m22_(m22), m32_(m32),
   m03_(m03), m13_(m13), m23_(m23), m33_(m33) {
  }

  CGLMatrix3D(const float *m, uint n) {
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
      assert(false && "Invalid size");
  }

  CGLMatrix3D(CMathGen::AxisType3D axis, float angle,
              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    setRotation(axis, angle, handedness);
  }

  //------

  // copy operations
  CGLMatrix3D(const CGLMatrix3D &a) {
    memcpy(&m00_, &a.m00_, 16*sizeof(float));
  }

  CGLMatrix3D &operator=(const CGLMatrix3D &a) {
    memcpy(&m00_, &a.m00_, 16*sizeof(float));

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

  friend std::ostream &operator<<(std::ostream &os, const CGLMatrix3D &matrix) {
    matrix.print(os);

    return os;
  }

  //------

  // comparison
  int cmp(const CGLMatrix3D &v) const {
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

  friend bool operator==(const CGLMatrix3D &lhs, const CGLMatrix3D &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const CGLMatrix3D &lhs, const CGLMatrix3D &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const CGLMatrix3D &lhs, const CGLMatrix3D &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const CGLMatrix3D &lhs, const CGLMatrix3D &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const CGLMatrix3D &lhs, const CGLMatrix3D &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const CGLMatrix3D &lhs, const CGLMatrix3D &rhs) {
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

  void setIdentity() {
    setInnerIdentity();

    setOuterIdentity();
  }

  //------

  static CGLMatrix3D translation(float tx, float ty, float tz) {
    CGLMatrix3D m;

    m.setTranslation(tx, ty, tz);

    return m;
  }

  static CGLMatrix3D translation(const CGLVector3D &vector) {
    CGLMatrix3D m;

    m.setTranslation(vector);

    return m;
  }

  CGLMatrix3D &translated(float tx, float ty, float tz) {
    CGLMatrix3D m;

    m.setTranslation(tx, ty, tz);

    *this *= m;

    return *this;
  }

  CGLMatrix3D &setTranslation(float tx, float ty, float tz) {
    setInnerIdentity();

    setOuterTranslate(tx, ty, tz);

    return *this;
  }

  CGLMatrix3D &setTranslation(const CPoint3D &point) {
    setInnerIdentity();

    setOuterTranslate(float(point.x), float(point.y), float(point.z));

    return *this;
  }

  CGLMatrix3D &setTranslation(const CGLVector3D &vector) {
    setInnerIdentity();

    setOuterTranslate(vector.getX(), vector.getY(), vector.getZ());

    return *this;
  }

  //---

  static CGLMatrix3D scale(float sx, float sy, float sz) {
    CGLMatrix3D m;

    m.setScale(sx, sy, sz);

    return m;
  }

  CGLMatrix3D &scaled(float sx, float sy, float sz) {
    CGLMatrix3D m;

    m.setScale(sx, sy, sz);

    *this *= m;

    return *this;
  }

  void setScale(float s) {
    setInnerScale(s, s, s);

    setOuterIdentity();
  }

  void setScale(float sx, float sy, float sz) {
    setInnerScale(sx, sy, sz);

    setOuterIdentity();
  }

  void setScaleTranslation(float s, float tx, float ty, float tz) {
    setInnerScale(s, s, s);

    setOuterTranslate(tx, ty, tz);
  }

  void setScaleTranslation(float sx, float sy, float sz, float tx, float ty, float tz) {
    setInnerScale(sx, sy, sz);

    setOuterTranslate(tx, ty, tz);
  }

  //---

  static CGLMatrix3D rotation(float theta, const CGLVector3D &u) {
    CGLMatrix3D m;

    m.setRotation(theta, u);

    return m;
  }

  CGLMatrix3D &rotated(float theta, const CGLVector3D &u) {
    CGLMatrix3D m;

    m.setRotation(theta, u);

    *this *= m;

    return *this;
  }

  CGLMatrix3D &setRotation(CMathGen::AxisType3D axis, float angle,
                           CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    if (handedness == CMathGen::RIGHT_HANDEDNESS)
      setInnerRotationRHS(axis, angle);
    else
      setInnerRotationLHS(axis, angle);

    setOuterIdentity();

    return *this;
  }

  CGLMatrix3D setRotation(float theta, const CGLVector3D &u) {
    float theta2 = 0.5f*theta;

    // rotate around the line
    float w = std::cos(theta2);

    //TODO: shouldn't have to normalize u
    CGLVector3D v = u.normalized()*std::sin(theta2);

    // assign matrix
    float x = v.getX(); float y = v.getY(); float z = v.getZ();

    float x2 = 2.0f*x; float y2 = 2.0f*y; float z2 = 2.0f*z;

    float wx2 = x2*w; float wy2 = y2*w; float wz2 = z2*w;
    float xx2 = x2*x; float xy2 = y2*x; float xz2 = z2*x;
    float yy2 = y2*y; float yz2 = z2*y; float zz2 = z2*z;

    float a = 1.0f - (yy2 + zz2); float b =         xy2 - wz2 ; float c =         xz2 + wy2 ;
    float d =         xy2 + wz2 ; float e = 1.0f - (xx2 + zz2); float f =         yz2 - wx2 ;
    float g =         xz2 - wy2 ; float h =         yz2 + wx2 ; float i = 1.0f - (xx2 + yy2);

    setValues(a, b, c, d, e, f, g, h, i, 0.0, 0.0, 0.0);

    return *this;
  }

  CGLMatrix3D &setRotationTranslation(CMathGen::AxisType3D axis, float angle,
                                      float tx, float ty, float tz,
                                      CMathGen::Handedness handedness=CMathGen::RIGHT_HANDEDNESS) {
    if (handedness == CMathGen::RIGHT_HANDEDNESS)
      setInnerRotationRHS(axis, angle);
    else
      setInnerRotationLHS(axis, angle);

    setOuterTranslate(tx, ty, tz);

    return *this;
  }

  CGLMatrix3D &setXYZRotation(float x_angle, float y_angle, float z_angle,
                              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    CGLMatrix3D xmatrix(CMathGen::X_AXIS_3D, x_angle, handedness);
    CGLMatrix3D ymatrix(CMathGen::Y_AXIS_3D, y_angle, handedness);
    CGLMatrix3D zmatrix(CMathGen::Z_AXIS_3D, z_angle, handedness);

    *this = xmatrix*ymatrix*zmatrix;

    return *this;
  }

  CGLMatrix3D &setXYZRotation(const CGLVector3D &angles,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    CGLMatrix3D xmatrix(CMathGen::X_AXIS_3D, angles.getX(), handedness);
    CGLMatrix3D ymatrix(CMathGen::Y_AXIS_3D, angles.getY(), handedness);
    CGLMatrix3D zmatrix(CMathGen::Z_AXIS_3D, angles.getZ(), handedness);

    *this = xmatrix*ymatrix*zmatrix;

    return *this;
  }

  CGLMatrix3D &setGenRotation(float x1, float y1, float z1, float x2, float y2, float z2,
                              float angle,
                              CMathGen::Handedness handedness=CMathGen::RIGHT_HANDEDNESS) {
    CGLMatrix3D matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7;

    matrix1.setTranslation(-x1, -y1, -z1);
    matrix2.setTranslation( x1,  y1,  z1);

    float theta = CMathGen::atan2(x2 - x1, y2 - y1);

    matrix3.setRotation(CMathGen::Z_AXIS_3D,  theta, handedness);
    matrix4.setRotation(CMathGen::Z_AXIS_3D, -theta, handedness);

    float v = float(std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)));

    float beta = CMathGen::atan2(z2 - z1, v);

    matrix5.setRotation(CMathGen::Y_AXIS_3D,  beta, handedness);
    matrix6.setRotation(CMathGen::Y_AXIS_3D, -beta, handedness);

    matrix7.setRotation(CMathGen::Z_AXIS_3D, angle, handedness);

    *this = matrix2*(matrix4*(matrix6*(matrix7*(matrix5*(matrix3*matrix1)))));

    return *this;
  }

  CGLMatrix3D &setGenRotation(const CGLVector3D &axis, float angle,
                              CMathGen::Handedness handedness=CMathGen::RIGHT_HANDEDNESS) {
    setOuterIdentity();

    CGLVector3D a = axis.normalized();

    float c = std::cos(angle);
    float s = std::sin(angle);

    float c1 = 1.0f - c;

    float axx = a.getX()*a.getX();
    float axy = a.getX()*a.getY();
    float axz = a.getX()*a.getZ();
    float ayy = a.getY()*a.getY();
    float ayz = a.getY()*a.getZ();
    float azz = a.getZ()*a.getZ();

    float axs = a.getX()*s; float ays = a.getY()*s; float azs = a.getZ()*s;

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

    return *this;
  }

  //---

  // NOTE: does not match setLookAt code below (OpenGL specific ?)
  static CGLMatrix3D lookAt(const CGLVector3D &eye, const CGLVector3D &center,
                            const CGLVector3D &up) {
    CGLVector3D f(center.x() - eye.x(), center.y() - eye.y(), center.z() - eye.z());

    auto up1 = up;

    f  .normalize();
    up1.normalize();

    CGLVector3D s = f.crossProduct(up1);
    CGLVector3D u = s.crossProduct(f);

    f = -f;

    CGLMatrix3D m;

    m.setColumn(0, s);
    m.setColumn(1, u);
    m.setColumn(2, f);

    auto ex = s.dotProduct(eye);
    auto ey = u.dotProduct(eye);
    auto ez = f.dotProduct(eye);

    m.setOuterTranslate(-ex, -ey, -ez);

    return m;
  }

  void setLookAt(const CPoint3D &eye, const CPoint3D &center) {
    CGLVector3D dir(eye, center);

    setLookAt(eye, dir);
  }

  void setLookAt(const CPoint3D &eye, const CPoint3D &center, const CGLVector3D &up) {
    CGLVector3D dir(eye, center);

    setLookAt(eye, dir, up);
  }

  void setLookAt(const CPoint3D &eye, const CGLVector3D &dir) {
    CGLVector3D dir1 = dir.normalized();

    CGLVector3D up(0, 0, 1);

    CGLVector3D up1 = up - (up.dotProduct(dir1))*dir1;

    setLookAt(eye, dir1, up1);
  }

  void setLookAt(const CPoint3D &eye, const CGLVector3D &dir, const CGLVector3D &up) {
    CGLVector3D dir1 = dir.normalized();
    CGLVector3D up1  = up .normalized();

    CGLVector3D right = dir1 .crossProduct(up1 );
    CGLVector3D newUp = right.crossProduct(dir1);

    //dir1 = -dir1;

    setColumn(0, right);
    setColumn(1, newUp);
    setColumn(2, -dir1);

    setOuterTranslate(float(eye.x), float(eye.y), float(eye.z));
  }

  void setEye(float x1, float y1, float z1, float x2, float y2, float z2,
              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    float angle1, angle2, angle3;

    calcEye(x1, y1, z1, x2, y2, z2, &angle1, &angle2, &angle3);

    CGLMatrix3D matrix1, matrix2, matrix3, matrix4;

    matrix1.setTranslation(-x1, -y1, -z1);

    matrix2.setRotation(CMathGen::Z_AXIS_3D,  angle1, handedness);
    matrix3.setRotation(CMathGen::Y_AXIS_3D,  angle2, handedness);
    matrix4.setRotation(CMathGen::Z_AXIS_3D, -angle3, handedness);

    *this = matrix4*(matrix3*(matrix2*matrix1));
  }

  static void calcEye(float x1, float y1, float z1, float x2, float y2, float z2,
                      float *angle1, float *angle2, float *angle3) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;

    *angle1 = CMathGen::atan2(-dx, -dy);

    float v = float(std::sqrt(dx*dx + dy*dy));

    *angle2 = CMathGen::atan2(-dz, v);

    float w = float(std::sqrt(v*v + dz*dz));

    *angle3 = CMathGen::atan2(-dx*w, dy*dz);
  }

  void setValues(float m00, float m01, float m02,
                 float m10, float m11, float m12,
                 float m20, float m21, float m22) {
    m00_ = m00, m01_ = m01, m02_ = m02;
    m10_ = m10, m11_ = m11, m12_ = m12;
    m20_ = m20, m21_ = m21, m22_ = m22;

    setOuterIdentity();
  }

  void setValues(float m00, float m01, float m02, float m10, float m11, float m12,
                 float m20, float m21, float m22, float tx , float ty , float tz ) {
    m00_ = m00, m01_ = m01, m02_ = m02, m03_ = tx;
    m10_ = m10, m11_ = m11, m12_ = m12, m13_ = ty;
    m20_ = m20, m21_ = m21, m22_ = m22, m23_ = tz;

    setBottomIdentity();
  }

  void setValues(float m00, float m01, float m02, float m03,
                 float m10, float m11, float m12, float m13,
                 float m20, float m21, float m22, float m23,
                 float m30, float m31, float m32, float m33) {
    m00_ = m00; m01_ = m01; m02_ = m02; m03_ = m03;
    m10_ = m10; m11_ = m11; m12_ = m12; m13_ = m13;
    m20_ = m20; m21_ = m21; m22_ = m22; m23_ = m23;
    m30_ = m30; m31_ = m31; m32_ = m32; m33_ = m33;
  }

  void getValues(float *m00, float *m01, float *m02,
                 float *m10, float *m11, float *m12,
                 float *m20, float *m21, float *m22) const {
    if (m00) { *m00 = m00_; } if (m01) { *m01 = m01_; } if (m02) { *m02 = m02_; }
    if (m10) { *m10 = m10_; } if (m11) { *m11 = m11_; } if (m12) { *m12 = m12_; }
    if (m20) { *m20 = m20_; } if (m21) { *m21 = m21_; } if (m22) { *m22 = m22_; }
  }

  void getValues(float *m00, float *m01, float *m02,
                 float *m10, float *m11, float *m12,
                 float *m20, float *m21, float *m22,
                 float *tx , float *ty , float *tz ) const {
    if (m00) { *m00 = m00_; } if (m01) { *m01 = m01_; } if (m02) { *m02 = m02_; }
    if (m10) { *m10 = m10_; } if (m11) { *m11 = m11_; } if (m12) { *m12 = m12_; }
    if (m20) { *m20 = m20_; } if (m21) { *m21 = m21_; } if (m22) { *m22 = m22_; }
    if (tx ) { *tx  = m03_; } if (ty ) { *ty  = m13_; } if (tz ) { *tz  = m23_; }
  }

  void getValues(float *m00, float *m01, float *m02, float *m03,
                 float *m10, float *m11, float *m12, float *m13,
                 float *m20, float *m21, float *m22, float *m23,
                 float *m30, float *m31, float *m32, float *m33) const {
    if (m00) { *m00 = m00_; } if (m01) { *m01 = m01_; }
    if (m02) { *m02 = m02_; } if (m03) { *m03 = m03_; }
    if (m10) { *m10 = m10_; } if (m11) { *m11 = m11_; }
    if (m12) { *m12 = m12_; } if (m13) { *m13 = m13_; }
    if (m20) { *m20 = m20_; } if (m21) { *m21 = m21_; }
    if (m22) { *m22 = m22_; } if (m23) { *m23 = m23_; }
    if (m30) { *m30 = m30_; } if (m31) { *m31 = m31_; }
    if (m32) { *m32 = m32_; } if (m33) { *m33 = m33_; }
  }

  void getValues(float *v, int n) const {
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
      assert(false && "Invalid size");
  }

  //---------

  const float *getData() const { return &m00_; }

  void setColumn(uint c, float x, float y, float z) {
    assert(c < 4);

    float *m = &(&m00_)[c];

    m[0] = x, m[4] = y, m[8] = z;
  }

  void setColumn(uint c, float x, float y, float z, float w) {
    assert(c < 4);

    float *m = &(&m00_)[c];

    m[0] = x, m[4] = y, m[8] = z, m[12] = w;
  }

  void setColumn(uint c, const CGLVector3D &vector) {
    assert(c < 4);

    float *m = &(&m00_)[c];

    float x, y, z;

    vector.getXYZ(&x, &y, &z);

    m[0] = float(x);
    m[4] = float(y);
    m[8] = float(z);
  }

  void getColumn(uint c, float *x, float *y, float *z) {
    assert(c < 4);

    float *m = &(&m00_)[c];

    if (x) *x = m[0];
    if (y) *y = m[4];
    if (z) *z = m[8];
  }

  void getColumn(uint c, float *x, float *y, float *z, float *w) {
    assert(c < 4);

    float *m = &(&m00_)[c];

    if (x) *x = m[ 0];
    if (y) *y = m[ 4];
    if (z) *z = m[ 8];
    if (w) *w = m[12];
  }

  //---------

  void setRow(uint r, float x, float y, float z) {
    assert(r < 4);

    float *m = &(&m00_)[r*4];

    m[0] = x, m[1] = y, m[2] = z;
  }

  void setRow(uint r, float x, float y, float z, float w) {
    assert(r < 4);

    float *m = &(&m00_)[r*4];

    m[0] = x, m[1] = y, m[2] = z, m[3] = w;
  }

  void setRow(uint r, const CGLVector3D &vector) {
    assert(r < 4);

    float *m = &(&m00_)[r*4];

    float x, y, z;

    vector.getXYZ(&x, &y, &z);

    m[0] = float(x);
    m[1] = float(y);
    m[3] = float(z);
  }

  void getRow(uint r, float *x, float *y, float *z) {
    assert(r < 4);

    float *m = &(&m00_)[r*4];

    if (x) *x = m[0];
    if (y) *y = m[1];
    if (z) *z = m[2];
  }

  void getRow(uint r, float *x, float *y, float *z, float *w) {
    assert(r < 4);

    float *m = &(&m00_)[r*4];

    if (x) *x = m[0];
    if (y) *y = m[1];
    if (z) *z = m[2];
    if (w) *w = m[3];
  }

  //---------

  void multiplyPoint(float xi, float yi, float zi, float *xo, float *yo, float *zo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi + m03_;
    *yo = m10_*xi + m11_*yi + m12_*zi + m13_;
    *zo = m20_*xi + m21_*yi + m22_*zi + m23_;

    float w = m30_*xi + m31_*yi + m32_*zi + m33_;

    if (std::abs(w) > 1E-6) {
      float iw = 1.0f/w;

      *xo *= iw;
      *yo *= iw;
      *zo *= iw;
    }
  }

  void multiplyPoint(float xi, float yi, float zi, float wi,
                     float *xo, float *yo, float *zo, float *wo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi + m03_*wi;
    *yo = m10_*xi + m11_*yi + m12_*zi + m13_*wi;
    *zo = m20_*xi + m21_*yi + m22_*zi + m23_*wi;
    *wo = m30_*xi + m31_*yi + m32_*zi + m33_*wi;
  }

  void multiplyPoint(const CPoint3D &ipoint, CPoint3D &opoint) const {
    float x = m00_*float(ipoint.x) + m01_*float(ipoint.y) + m02_*float(ipoint.z) + m03_;
    float y = m10_*float(ipoint.x) + m11_*float(ipoint.y) + m12_*float(ipoint.z) + m13_;
    float z = m20_*float(ipoint.x) + m21_*float(ipoint.y) + m22_*float(ipoint.z) + m23_;
    float w = m30_*float(ipoint.x) + m31_*float(ipoint.y) + m32_*float(ipoint.z) + m33_;

    float iw = (std::abs(w) > 1E-6 ? 1.0f/w : 1.0f);

    opoint.x = x*iw;
    opoint.y = y*iw;
    opoint.z = z*iw;
  }

  CPoint3D multiplyPoint(const CPoint3D &ipoint) const {
    CPoint3D opoint;

    float x = m00_*float(ipoint.x) + m01_*float(ipoint.y) + m02_*float(ipoint.z) + m03_;
    float y = m10_*float(ipoint.x) + m11_*float(ipoint.y) + m12_*float(ipoint.z) + m13_;
    float z = m20_*float(ipoint.x) + m21_*float(ipoint.y) + m22_*float(ipoint.z) + m23_;
    float w = m30_*float(ipoint.x) + m31_*float(ipoint.y) + m32_*float(ipoint.z) + m33_;

    float iw = (std::abs(w) > 1E-6 ? 1.0f/w : 1.0f);

    opoint.x = x*iw;
    opoint.y = y*iw;
    opoint.z = z*iw;

    return opoint;
  }

  void multiplyVector(const CGLVector3D &ivector, CGLVector3D &ovector) const {
    float ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    float ox = m00_*ix + m01_*iy + m02_*iz + m03_;
    float oy = m10_*ix + m11_*iy + m12_*iz + m13_;
    float oz = m20_*ix + m21_*iy + m22_*iz + m23_;

    ovector.setXYZ(ox, oy, oz);
  }

  CGLVector3D multiplyVector(const CGLVector3D &ivector) const {
    CGLVector3D ovector;

    float ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    float ox = m00_*ix + m01_*iy + m02_*iz + m03_;
    float oy = m10_*ix + m11_*iy + m12_*iz + m13_;
    float oz = m20_*ix + m21_*iy + m22_*iz + m23_;

    ovector.setXYZ(ox, oy, oz);

    return ovector;
  }

  void preMultiplyPoint(float xi, float yi, float zi, float *xo, float *yo, float *zo) const {
    *xo = m00_*xi + m10_*yi + m20_*zi;
    *yo = m01_*xi + m11_*yi + m21_*zi;
    *zo = m02_*xi + m12_*yi + m22_*zi;

    float w = m03_*xi + m13_*yi + m23_*zi + m33_;

    float iw = (std::abs(w) > 1E-6 ? 1.0f/w : 1.0f);

    *xo *= iw;
    *yo *= iw;
    *xo *= iw;
  }

  void preMultiplyPoint(float xi, float yi, float zi, float wi,
                        float *xo, float *yo, float *zo, float *wo) const {
    *xo = m00_*xi + m10_*yi + m20_*zi + m30_*wi;
    *yo = m01_*xi + m11_*yi + m21_*zi + m31_*wi;
    *zo = m02_*xi + m12_*yi + m22_*zi + m32_*wi;
    *wo = m03_*xi + m13_*yi + m23_*zi + m33_*wi;
  }

  void preMultiplyPoint(const CPoint3D &ipoint, CPoint3D &opoint) const {
    opoint.x = m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z;
    opoint.z = m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z;

    float w = float(m03_*ipoint.x + m13_*ipoint.y + m23_*ipoint.z + m33_);

    float iw = (std::abs(w) > 1E-6 ? 1.0f/w : 1.0f);

    opoint *= iw;
  }

  void preMultiplyVector(const CGLVector3D &ivector, CGLVector3D &ovector) const {
    float ix, iy, iz;

    ivector.getXYZ(&ix, &iy, &iz);

    float ox = m00_*ix + m10_*iy + m20_*iz;
    float oy = m01_*ix + m11_*iy + m21_*iz;
    float oz = m02_*ix + m12_*iy + m22_*iz;

    ovector.setXYZ(ox, oy, oz);
  }

  CGLMatrix3D &translate(float x, float y, float z) {
    m03_ += x;
    m13_ += y;
    m23_ += z;

    return *this;
  }

  //------

  CGLMatrix3D &transpose() {
    std::swap(m10_, m01_);
    std::swap(m20_, m02_);
    std::swap(m21_, m12_);
    std::swap(m30_, m03_);
    std::swap(m31_, m13_);
    std::swap(m32_, m23_);

    return *this;
  }

  CGLMatrix3D transposed() const {
    auto matrix = *this;

    matrix.transpose();

    return matrix;
  }

  //------

  bool invert(CGLMatrix3D &imatrix) const {
    float d = determinant();

    if (::fabs(d) == 0.0)
      return false;

    float id = 1.0f/d;

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

  CGLMatrix3D inverse() const {
    CGLMatrix3D imatrix;

    if (! invert(imatrix))
      assert(false && "Divide by zero");

    return imatrix;
  }

  float determinant() const {
    return
      (m00_*det3x3(m11_, m12_, m13_, m21_, m22_, m23_, m31_, m32_, m33_) -
       m01_*det3x3(m10_, m12_, m13_, m20_, m22_, m23_, m30_, m32_, m33_) +
       m02_*det3x3(m10_, m11_, m13_, m20_, m21_, m23_, m30_, m31_, m33_) -
       m03_*det3x3(m10_, m11_, m12_, m20_, m21_, m22_, m30_, m31_, m32_));
  }

  void normalize() {
    float d = determinant();

    float id = 1.0f/d;

    for (int i = 0; i < 9; ++i)
      (&m00_)[i] *= id;
  }

  //------

  bool affineInvert(CGLMatrix3D &imatrix) const {
    float d = affineDeterminant();

    if (::fabs(d) == 0.0)
      return false;

    float id = 1.0f/d;

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

  CGLMatrix3D affineInverse() const {
    CGLMatrix3D imatrix;

    if (! affineInvert(imatrix))
      assert(false && "Divide by zero");

    return imatrix;
  }

  float affineDeterminant() const {
    return (m00_*det2x2(m11_, m12_, m21_, m22_) -
            m01_*det2x2(m10_, m12_, m20_, m22_) +
            m02_*det2x2(m10_, m11_, m20_, m21_));
  }

  void affineNormalize() {
    float d = affineDeterminant();

    float id = 1.0f/d;

    for (int i = 0; i < 9; ++i)
      (&m00_)[i] *= id;
  }

  //------

  CGLMatrix3D &rotate(float theta, const CGLVector3D &u) {
    CGLMatrix3D m; m.setRotation(theta, u);

    *this *= m;

    return *this;
  }

  //------

  void setTransform(float xmin1, float ymin1, float zmin1,
                    float xmax1, float ymax1, float zmax1,
                    float xmin2, float ymin2, float zmin2,
                    float xmax2, float ymax2, float zmax2) {
    float sx = (xmax2 - xmin2)/(xmax1 - xmin1);
    float sy = (ymax2 - ymin2)/(ymax1 - ymin1);
    float sz = (zmax2 - zmin2)/(zmax1 - zmin1);

    float tx = -xmin1*sx + xmin2;
    float ty = -ymin1*sy + ymin2;
    float tz = -zmin1*sz + zmin2;

    setInnerScale(sx, sy, sz);

    setOuterTranslate(tx, ty, tz);
  }

  //------

  static CGLMatrix3D *newIdentityMatrix() {
    CGLMatrix3D *m = new CGLMatrix3D();

    m->setIdentity();

    return m;
  }

  //------

  static bool solveAXeqB(const CGLMatrix3D &a, CPoint3D &x, const CPoint3D &b) {
    float det_a = a.determinant();

    if (::fabs(det_a) == 0.0)
      return false;

    float idet_a = 1.0f/det_a;

    CGLMatrix3D t(a);

    t.setColumn(0, float(b.x), float(b.y), float(b.z));

    float det_t = t.determinant();

    x.x = det_t*idet_a;

    t = a;

    t.setColumn(1, float(b.x), float(b.y), float(b.z));

    det_t = t.determinant();

    x.y = det_t*idet_a;

    t = a;

    t.setColumn(2, float(b.x), float(b.y), float(b.z));

    det_t = t.determinant();

    x.z = det_t*idet_a;

    return true;
  }

  //------

  void zero() { memset(&m00_, 0, 16*sizeof(float)); }

  //------

  CGLMatrix3D &operator+=(const CGLMatrix3D &b) {
    m00_ += b.m00_; m01_ += b.m01_; m02_ += b.m02_; m03_ += b.m03_;
    m10_ += b.m10_; m11_ += b.m11_; m12_ += b.m12_; m13_ += b.m13_;
    m20_ += b.m20_; m21_ += b.m21_; m22_ += b.m22_; m23_ += b.m23_;
    m30_ += b.m30_; m31_ += b.m31_; m32_ += b.m32_; m33_ += b.m33_;

    return *this;
  }

  CGLMatrix3D operator+(const CGLMatrix3D &b) {
    CGLMatrix3D c = *this;

    c += b;

    return c;
  }

  CGLMatrix3D &operator-=(const CGLMatrix3D &b) {
    m00_ -= b.m00_; m01_ -= b.m01_; m02_ -= b.m02_; m03_ -= b.m03_;
    m10_ -= b.m10_; m11_ -= b.m11_; m12_ -= b.m12_; m13_ -= b.m13_;
    m20_ -= b.m20_; m21_ -= b.m21_; m22_ -= b.m22_; m23_ -= b.m23_;
    m30_ -= b.m30_; m31_ -= b.m31_; m32_ -= b.m32_; m33_ -= b.m33_;

    return *this;
  }

  CGLMatrix3D operator-(const CGLMatrix3D &b) {
    CGLMatrix3D c = *this;

    c -= b;

    return c;
  }

  CGLMatrix3D &operator*=(const CGLMatrix3D &b) {
    CGLMatrix3D a;

    memcpy(&a.m00_, &m00_, 16*sizeof(float));

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

  CGLMatrix3D operator*(const CGLMatrix3D &b) const {
    CGLMatrix3D c = *this;

    c *= b;

    return c;
  }

  CGLMatrix3D &operator*=(float s) {
    m00_ *= s; m01_ *= s; m02_ *= s; m03_ *= s;
    m10_ *= s; m11_ *= s; m12_ *= s; m13_ *= s;
    m20_ *= s; m21_ *= s; m22_ *= s; m23_ *= s;
    m30_ *= s; m31_ *= s; m32_ *= s; m33_ *= s;

    return *this;
  }

  CGLMatrix3D operator*(float s) {
    CGLMatrix3D c = *this;

    c *= s;

    return c;
  }

  friend CPoint3D operator*(const CGLMatrix3D &m, const CPoint3D &p) {
    CPoint3D p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend CPoint3D operator*(const CPoint3D &p, const CGLMatrix3D &m) {
    CPoint3D p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend CGLVector3D operator*(const CGLMatrix3D &m, const CGLVector3D &v) {
    CGLVector3D v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend CGLVector3D operator*(const CGLVector3D &v, const CGLMatrix3D &m) {
    CGLVector3D v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  CGLMatrix3D &operator/=(const CGLMatrix3D &b) {
    CGLMatrix3D bi;

    if (! b.invert(bi)) {
      assert(false && "Divide by zero");
      return *this;
    }

    return (*this) *= bi;
  }

  CGLMatrix3D operator/(const CGLMatrix3D &b) {
    CGLMatrix3D c = *this;

    c /= b;

    return c;
  }

  //------

  void setValue(uint i, float value) {
    (&m00_)[i] = value;
  }

  void setValue(uint i, uint j, float value) {
    assert(i < 4 && j < 4);

    float *m = &(&m00_)[4*j + i];

    *m = value;
  }

  float getValue(uint i) const {
    assert(i < 16);

    return (&m00_)[i];
  }

  float getValue(uint i, uint j) const {
    assert(i < 4 && j < 4);

    const float &m = (&m00_)[4*j + i];

    return m;
  }

  float operator[](uint i) const {
    assert(i < 16);

    return (&m00_)[i];
  }

  float &operator[](uint i) {
    assert(i < 16);

    return (&m00_)[i];
  }

  float operator()(uint i, uint j) const {
    assert(i < 4 && j < 4);

    const float &m = (&m00_)[4*j + i];

    return m;
  }

  float &operator()(uint i, uint j) {
    assert(i < 4 && j < 4);

    float &m = (&m00_)[4*j + i];

    return m;
  }

  //------

  static CGLMatrix3D perspective(float fov, float aspect, float near, float far) {
    CGLMatrix3D m;

    m.buildPerspective(fov, aspect, near, far);

    return m;
  }

  /*! build perspective projection matrix
      (near = -ve, far = -ve) look down negative Z-axis
      (near = +ve, far = +ve) look down position Z-axis
  */
  CGLMatrix3D &buildPerspective(float fov, float aspect, float near, float far) {
    // can't have near/far on other sides of origin
    if (near*far <= 0) return *this;

    if (near < 0) {
      float tf2  = float(std::tan(0.5*CMathGen::DegToRad(fov)));
      float itf2 = 1.0f/tf2;

      float npf = near + far;
      float ntf = near * far;

      float d  = near - far;
      float id = 1.0f/d;

      float a  = itf2/aspect;
      float e  = itf2;
      float i  = npf*id;
      float tz = 2.0f*ntf*id;

      m00_ = a  , m01_ = 0.0, m02_ = 0.0 ; m03_ = 0.0;
      m10_ = 0.0, m11_ = e  , m12_ = 0.0 ; m13_ = 0.0;
      m20_ = 0.0, m21_ = 0.0, m22_ = i   ; m23_ = tz ;
      m30_ = 0.0, m31_ = 0.0, m32_ = 1.0f; m33_ = 0.0;
    }
    else {
      float tf2  = float(std::tan(0.5*CMathGen::DegToRad(fov)));
      float itf2 = 1.0f/tf2;

      float ntf = near*far;

      float d  = far - near;
      float id = 1.0f/d;

      float a  = itf2/aspect;
      float e  = itf2;
      float i  = far*id;
      float tz = -2.0f*ntf*id;

      m00_ = a  , m01_ = 0.0, m02_ =  0.0 ; m03_ = 0.0;
      m10_ = 0.0, m11_ = e  , m12_ =  0.0 ; m13_ = 0.0;
      m20_ = 0.0, m21_ = 0.0, m22_ = -i   ; m23_ = tz ;
      m30_ = 0.0, m31_ = 0.0, m32_ = -1.0f; m33_ = 0.0;
    }

    return *this;
  }

  void buildOrtho(float left, float right, float bottom, float top, float near, float far) {
    float w = right - left  ; float iw = 1.0f/w;
    float h = top   - bottom; float ih = 1.0f/h;
    float d = near  - far   ; float id = 1.0f/d; // far - near ?

    float rpl = right + left  ;
    float tpb = top   + bottom;
    float npf = near  + far   ;

    float a = 2.0f*iw;
    float e = 2.0f*ih;
    float i = 2.0f*id;

    float tx = -rpl*iw;
    float ty = -tpb*ih;
    float tz =  npf*id;

    setInnerScale(a, e, i);

    setOuterTranslate(tx, ty, tz);
  }

  void buildFrustrum(float left, float right, float bottom, float top,
                     float near, float far) {
    float w = right - left  ; float iw = 1.0f/w;
    float h = top   - bottom; float ih = 1.0f/h;
    float d = near  - far   ; float id = 1.0f/d; // far - near ?

    float npf = near + far;
    float ntf = near * far;

    float a  = 2.0f*near*iw;
    float e  = 2.0f*near*ih;
    float i  = npf*id;

    float c  = (right + left  )*iw;
    float f  = (top   + bottom)*ih;
    float tz = 2.0f*ntf*id;

    m00_ = a  , m01_ = 0.0, m02_ =  c   ; m03_ = 0.0;
    m10_ = 0.0, m11_ = e  , m12_ =  f   ; m13_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ =  i   ; m23_ = tz ;
    m30_ = 0.0, m31_ = 0.0, m32_ = -1.0f; m33_ = 0.0;
  }

  //------

 private:
  void setInnerIdentity() {
    m00_ = 1.0f, m01_ = 0.0f, m02_ = 0.0f;
    m10_ = 0.0f, m11_ = 1.0f, m12_ = 0.0f;
    m20_ = 0.0f, m21_ = 0.0f, m22_ = 1.0f;
  }

  void setInnerScale(float sx, float sy, float sz) {
    m00_ = sx , m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = sy , m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = sz ;
  }

  void setInnerRotationRHS(CMathGen::AxisType3D axis, float angle) {
    float c = std::cos(angle);
    float s = std::sin(angle);

    if      (axis == CMathGen::X_AXIS_3D) {
      m00_ = 1.0f; m01_ = 0.0; m02_ = 0.0;
      m10_ = 0.0f; m11_ =   c; m12_ =   s;
      m20_ = 0.0f; m21_ =  -s; m22_ =   c;
    }
    else if (axis == CMathGen::Y_AXIS_3D) {
      m00_ =   c; m01_ = 0.0f; m02_ =  -s;
      m10_ = 0.0; m11_ = 1.0f; m12_ = 0.0;
      m20_ =   s; m21_ = 0.0f; m22_ =   c;
    }
    else {
      m00_ =   c; m01_ =   s; m02_ = 0.0f;
      m10_ =  -s; m11_ =   c; m12_ = 0.0f;
      m20_ = 0.0; m21_ = 0.0; m22_ = 1.0f;
    }
  }

  void setInnerRotationLHS(CMathGen::AxisType3D axis, float angle) {
    float c = std::cos(angle);
    float s = std::sin(angle);

    if      (axis == CMathGen::X_AXIS_3D) {
      m00_ = 1.0f; m01_ = 0.0; m02_ = 0.0;
      m10_ = 0.0f; m11_ =   c; m12_ =  -s;
      m20_ = 0.0f; m21_ =   s; m22_ =   c;
    }
    else if (axis == CMathGen::Y_AXIS_3D) {
      m00_ =   c; m01_ = 0.0f; m02_ =   s;
      m10_ = 0.0; m11_ = 1.0f; m12_ = 0.0;
      m20_ =  -s; m21_ = 0.0f; m22_ =   c;
    }
    else {
      m00_ =   c; m01_ =  -s; m02_ = 0.0f;
      m10_ =   s; m11_ =   c; m12_ = 0.0f;
      m20_ = 0.0; m21_ = 0.0; m22_ = 1.0f;
    }
  }

  void setOuterIdentity() {
    m03_ = 0.0; m13_ = 0.0; m23_ = 0.0;

    setBottomIdentity();
  }

  void setOuterTranslate(float tx, float ty, float tz) {
    m03_ = tx; m13_ = ty; m23_ = tz;

    setBottomIdentity();
  }

  void setBottomIdentity() {
    m30_ = 0.0, m31_ = 0.0, m32_ = 0.0, m33_ = 1.0f;
  }

  static float det3x3(float a, float b, float c, float d, float e,
                      float f, float g, float h, float i) {
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
  }

  static float det2x2(float a, float b, float c, float d) {
    return a*c - b*d;
  }

 private:
  // transposed data order
  float m00_ { 0.0 }, m10_ { 0.0 }, m20_ { 0.0 }, m30_ { 0.0 };
  float m01_ { 0.0 }, m11_ { 0.0 }, m21_ { 0.0 }, m31_ { 0.0 };
  float m02_ { 0.0 }, m12_ { 0.0 }, m22_ { 0.0 }, m32_ { 0.0 };
  float m03_ { 0.0 }, m13_ { 0.0 }, m23_ { 0.0 }, m33_ { 0.0 };
};

#endif
