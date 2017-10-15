#ifndef CMATRIX_3D_H
#define CMATRIX_3D_H

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CVector3D.h>
#include <CGaussianMatrix.h>

#include <cassert>
#include <cstring>

/* / m00 m01 m02 m03 \ */
/* | m10 m11 m12 m13 | */
/* | m20 m21 m22 m23 | */
/* \ m30 m31 m32 m33 / */

class CMatrix3D {
 public:
  enum class Type {
    IDENTITY
  };

 public:
  // constructor/destructor
  CMatrix3D() :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
  }

 ~CMatrix3D() { }

  explicit CMatrix3D(Type type) :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
    if (type == Type::IDENTITY)
      setIdentity();
    else
      assert(false && "Bad CMatrix3D Type");
  }

  CMatrix3D(double m00, double m01, double m02,
            double m10, double m11, double m12,
            double m20, double m21, double m22) :
   m00_(m00), m01_(m01), m02_(m02), m03_(0),
   m10_(m10), m11_(m11), m12_(m12), m13_(0),
   m20_(m20), m21_(m21), m22_(m22), m23_(0),
   m30_(0  ), m31_(0  ), m32_(0  ), m33_(0) {
    setOuterIdentity();
  }

  CMatrix3D(double m00, double m01, double m02,
            double m10, double m11, double m12,
            double m20, double m21, double m22, double tx, double ty, double tz) :
   m00_(m00), m01_(m01), m02_(m02), m03_(tx),
   m10_(m10), m11_(m11), m12_(m12), m13_(ty),
   m20_(m20), m21_(m21), m22_(m22), m23_(tz),
   m30_(0  ), m31_(0  ), m32_(0  ), m33_(0 ) {
    setBottomIdentity();
  }

  CMatrix3D(double m00, double m01, double m02, double m03,
            double m10, double m11, double m12, double m13,
            double m20, double m21, double m22, double m23,
            double m30, double m31, double m32, double m33) :
   m00_(m00), m01_(m01), m02_(m02), m03_(m03),
   m10_(m10), m11_(m11), m12_(m12), m13_(m13),
   m20_(m20), m21_(m21), m22_(m22), m23_(m23),
   m30_(m30), m31_(m31), m32_(m32), m33_(m33) {
  }

  CMatrix3D(const double *m, uint n) :
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
      assert(false && "Invalid size");
  }

  CMatrix3D(const CVector3D &v0, const CVector3D &v1, const CVector3D &v2) :
   m00_(v0.getX()), m01_(v1.getX()), m02_(v2.getX()), m03_(0),
   m10_(v0.getY()), m11_(v1.getY()), m12_(v2.getY()), m13_(0),
   m20_(v0.getZ()), m21_(v1.getZ()), m22_(v2.getZ()), m23_(0),
   m30_(0        ), m31_(0        ), m32_(0        ), m33_(0) {
    setOuterIdentity();
  }

  CMatrix3D(CMathGen::AxisType3D axis, double angle,
             CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    setRotation(axis, angle, handedness);
  }

  CMatrix3D *dup() const {
    return new CMatrix3D(*this);
  }

  //------

  // copy operations
  CMatrix3D(const CMatrix3D &a) :
   m00_(0), m01_(0), m02_(0), m03_(0),
   m10_(0), m11_(0), m12_(0), m13_(0),
   m20_(0), m21_(0), m22_(0), m23_(0),
   m30_(0), m31_(0), m32_(0), m33_(0) {
    memcpy(&m00_, &a.m00_, 16*sizeof(double));
  }

  CMatrix3D &operator=(const CMatrix3D &a) {
    memcpy(&m00_, &a.m00_, 16*sizeof(double));

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

  friend std::ostream &operator<<(std::ostream &os, const CMatrix3D &matrix) {
    matrix.print(os);

    return os;
  }

  //------

  // comparison
  int cmp(const CMatrix3D &v) const {
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

  friend bool operator==(const CMatrix3D &lhs, const CMatrix3D &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const CMatrix3D &lhs, const CMatrix3D &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const CMatrix3D &lhs, const CMatrix3D &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const CMatrix3D &lhs, const CMatrix3D &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const CMatrix3D &lhs, const CMatrix3D &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const CMatrix3D &lhs, const CMatrix3D &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  // TODO: isZero(), isIdentity()

  //------

  void setIdentity() {
    setInnerIdentity();

    setOuterIdentity();
  }

  //----------

  static CMatrix3D translation(double tx, double ty, double tz) {
    CMatrix3D m;

    m.setTranslation(tx, ty, tz);

    return m;
  }

  void setTranslation(double tx, double ty, double tz) {
    setInnerIdentity();

    setOuterTranslate(tx, ty, tz);
  }

  void setTranslation(const CPoint3D &point) {
    setInnerIdentity();

    setOuterTranslate(point.x, point.y, point.z);
  }

  void setTranslation(const CVector3D &vector) {
    setInnerIdentity();

    setOuterTranslate(vector.getX(), vector.getY(), vector.getZ());
  }

  CMatrix3D &translate(double x, double y, double z) {
    m03_ += x;
    m13_ += y;
    m23_ += z;

    return *this;
  }

  void getTranslate(double *tx, double *ty, double *tz) {
    *tx = m03_;
    *ty = m13_;
    *tz = m23_;
  }

  //----------

  static CMatrix3D scale(double s) {
    return scale(s, s, s);
  }

  static CMatrix3D scale(double sx, double sy, double sz) {
    CMatrix3D m;

    m.setScale(sx, sy, sz);

    return m;
  }

  void setScale(double s) {
    setInnerScale(s, s, s);

    setOuterIdentity();
  }

  void setScale(double sx, double sy, double sz) {
    setInnerScale(sx, sy, sz);

    setOuterIdentity();
  }

  void getScale(double *sx, double *sy, double *sz) {
    *sx = m00_;
    *sy = m11_;
    *sz = m22_;
  }

  //----------

  void setScaleTranslation(double s, double tx, double ty, double tz) {
    setInnerScale(s, s, s);

    setOuterTranslate(tx, ty, tz);
  }

  void setScaleTranslation(double sx, double sy, double sz, double tx, double ty, double tz) {
    setInnerScale(sx, sy, sz);

    setOuterTranslate(tx, ty, tz);
  }

  //----------

  static CMatrix3D rotation(double rx, double ry, double rz) {
    CMatrix3D m;

    m.setXYZRotation(rx, ry, rz);

    return m;
  }

  static CMatrix3D rotation(CMathGen::AxisType3D axis, double angle,
                             CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3D m;

    m.setRotation(axis, angle, handedness);

    return m;
  }

  void setRotation(CMathGen::AxisType3D axis, double angle,
                   CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    if (handedness == CMathGen::RIGHT_HANDEDNESS)
      setInnerRotationRHS(axis, angle);
    else
      setInnerRotationLHS(axis, angle);

    setOuterIdentity();
  }

  void setRotation(double theta, const CVector3D &u) {
    double theta2 = 0.5*theta;

    // rotate around the line
    double w = cos(theta2);

    //TODO: shouldn't have to normalize u
    CVector3D v = u*sin(theta2);

    //assign matrix
    double x = v.getX(); double y = v.getY(); double z = v.getZ();

    double x2 = 2.0*x; double y2 = 2.0*y; double z2 = 2.0*z;

    double wx2 = x2*w; double wy2 = y2*w; double wz2 = z2*w;
    double xx2 = x2*x; double xy2 = y2*x; double xz2 = z2*x;
    double yy2 = y2*y; double yz2 = z2*y; double zz2 = z2*z;

    double a = 1.0 - (yy2 + zz2); double b =        xy2 - wz2 ; double c =        xz2 + wy2 ;
    double d =        xy2 + wz2 ; double e = 1.0 - (xx2 + zz2); double f =        yz2 - wx2 ;
    double g =        xz2 - wy2 ; double h =        yz2 + wx2 ; double i = 1.0 - (xx2 + yy2);

    setValues(a, b, c, d, e, f, g, h, i, 0.0, 0.0, 0.0);
  }

  void setRotationTranslation(CMathGen::AxisType3D axis,
                              double angle, double tx, double ty, double tz,
                              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    if (handedness == CMathGen::RIGHT_HANDEDNESS)
      setInnerRotationRHS(axis, angle);
    else
      setInnerRotationLHS(axis, angle);

    setOuterTranslate(tx, ty, tz);
  }

  void setXYZRotation(double x_angle, double y_angle, double z_angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3D xmatrix(CMathGen::X_AXIS_3D, x_angle, handedness);
    CMatrix3D ymatrix(CMathGen::Y_AXIS_3D, y_angle, handedness);
    CMatrix3D zmatrix(CMathGen::Z_AXIS_3D, z_angle, handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setXYZRotation(const CVector3D &angles,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3D xmatrix(CMathGen::X_AXIS_3D, angles.getX(), handedness);
    CMatrix3D ymatrix(CMathGen::Y_AXIS_3D, angles.getY(), handedness);
    CMatrix3D zmatrix(CMathGen::Z_AXIS_3D, angles.getZ(), handedness);

    *this = xmatrix*ymatrix*zmatrix;
  }

  void setGenRotation(double x1, double y1, double z1,
                      double x2, double y2, double z2, double angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    CMatrix3D matrix1, matrix2, matrix3, matrix4, matrix5, matrix6, matrix7;

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

  void setGenRotation(const CVector3D &axis, double angle,
                      CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    setOuterIdentity();

    CVector3D a = axis.normalized();

    double c = ::cos(angle);
    double s = ::sin(angle);

    double c1 = 1.0 - c;

    double axx = a.getX()*a.getX(); double axy = a.getX()*a.getY(); double axz = a.getX()*a.getZ();
                                    double ayy = a.getY()*a.getY(); double ayz = a.getY()*a.getZ();
                                                                    double azz = a.getZ()*a.getZ();

    double axs = a.getX()*s; double ays = a.getY()*s; double azs = a.getZ()*s;

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

  void setLookAt(const CPoint3D &eye, const CPoint3D &center) {
    CVector3D dir(eye, center);

    setLookAt(eye, dir);
  }

  void setLookAt(const CPoint3D &eye, const CPoint3D &center, const CVector3D &up) {
    CVector3D dir(eye, center);

    setLookAt(eye, dir, up);
  }

  void setLookAt(const CPoint3D &eye, const CVector3D &dir) {
    CVector3D dir1 = dir.normalized();

    CVector3D up(0,0,1);

    CVector3D up1 = up - (up.dotProduct(dir1))*dir1;

    setLookAt(eye, dir1, up1);
  }

  void setLookAt(const CPoint3D &eye, const CVector3D &dir, const CVector3D &up) {
    CVector3D dir1 = dir.normalized();
    CVector3D up1  = up .normalized();

    CVector3D right = dir1 .crossProduct(up1 );
    CVector3D newUp = right.crossProduct(dir1);

    dir1 = -dir1;

    setColumn(0, right);
    setColumn(1, newUp);
    setColumn(2, dir1 );

    setOuterTranslate(-eye.x, -eye.y, -eye.z);
  }

  //----------

  void setEye(double x1, double y1, double z1, double x2, double y2, double z2,
              CMathGen::Handedness handedness = CMathGen::RIGHT_HANDEDNESS) {
    double angle1, angle2, angle3;

    calcEye(x1, y1, z1, x2, y2, z2, &angle1, &angle2, &angle3);

    CMatrix3D matrix1, matrix2, matrix3, matrix4;

    matrix1.setTranslation(-x1, -y1, -z1);

    matrix2.setRotation(CMathGen::Z_AXIS_3D,  angle1, handedness);
    matrix3.setRotation(CMathGen::Y_AXIS_3D,  angle2, handedness);
    matrix4.setRotation(CMathGen::Z_AXIS_3D, -angle3, handedness);

    *this = matrix4*(matrix3*(matrix2*matrix1));
  }

  static void calcEye(double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double *angle1, double *angle2, double *angle3) {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;

    *angle1 = CMathGen::atan2(-dx, -dy);

    double v = ::sqrt(dx*dx + dy*dy);

    *angle2 = CMathGen::atan2(-dz, v);

    double w = ::sqrt(v*v + dz*dz);

    *angle3 = CMathGen::atan2(-dx*w, dy*dz);
  }

  //----------

  void setValues(double m00, double m01, double m02,
                 double m10, double m11, double m12,
                 double m20, double m21, double m22) {
    m00_ = m00, m01_ = m01, m02_ = m02;
    m10_ = m10, m11_ = m11, m12_ = m12;
    m20_ = m20, m21_ = m21, m22_ = m22;

    setOuterIdentity();
  }

  void setValues(double m00, double m01, double m02, double m10, double m11, double m12,
                 double m20, double m21, double m22, double tx , double ty , double tz ) {
    m00_ = m00, m01_ = m01, m02_ = m02, m03_ = tx;
    m10_ = m10, m11_ = m11, m12_ = m12, m13_ = ty;
    m20_ = m20, m21_ = m21, m22_ = m22, m23_ = tz;

    setBottomIdentity();
  }

  void setValues(double m00, double m01, double m02, double m03,
                 double m10, double m11, double m12, double m13,
                 double m20, double m21, double m22, double m23,
                 double m30, double m31, double m32, double m33) {
    m00_ = m00, m01_ = m01, m02_ = m02, m03_ = m03;
    m10_ = m10, m11_ = m11, m12_ = m12, m13_ = m13;
    m20_ = m20, m21_ = m21, m22_ = m22, m23_ = m23;
    m30_ = m30, m31_ = m31, m32_ = m32, m33_ = m33;
  }

  void setValues(const CVector3D &v0, const CVector3D &v1, const CVector3D &v2) {
    m00_ = v0.getX(); m01_ = v1.getX(); m02_ = v2.getX();
    m10_ = v0.getY(); m11_ = v1.getY(); m12_ = v2.getY();
    m20_ = v0.getZ(); m21_ = v1.getZ(); m22_ = v2.getZ();

    setOuterIdentity();
  }

  void getValues(double *m00, double *m01, double *m02,
                 double *m10, double *m11, double *m12,
                 double *m20, double *m21, double *m22) const {
    if (m00) *m00 = m00_; if (m01) *m01 = m01_; if (m02) *m02 = m02_;
    if (m10) *m10 = m10_; if (m11) *m11 = m11_; if (m12) *m12 = m12_;
    if (m20) *m20 = m20_; if (m21) *m21 = m21_; if (m22) *m22 = m22_;
  }

  void getValues(double *m00, double *m01, double *m02,
                 double *m10, double *m11, double *m12,
                 double *m20, double *m21, double *m22,
                 double *tx , double *ty , double *tz ) const {
    if (m00) *m00 = m00_; if (m01) *m01 = m01_; if (m02) *m02 = m02_;
    if (m10) *m10 = m10_; if (m11) *m11 = m11_; if (m12) *m12 = m12_;
    if (m20) *m20 = m20_; if (m21) *m21 = m21_; if (m22) *m22 = m22_;

    if (tx ) *tx  = m03_; if (ty ) *ty  = m13_; if (tz ) *tz  = m23_;
  }

  void getValues(double *v, int n) const {
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

  void getValuesI(double *v, int n) const {
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
      assert(false && "Invalid size");
  }

  //---------

  const double *getData() const {
    return &m00_;
  }

  void setColumn(int c, double x, double y, double z) {
    assert(c < 4);

    double *m = &(&m00_)[c];

    m[0] = x, m[4] = y, m[8] = z;
  }

  void setColumn(int c, const CPoint3D &point) {
    setColumn(c, point.x, point.y, point.z);
  }

  void setColumn(int c, const CVector3D &vector) {
    assert(c < 4);

    double *m = &(&m00_)[c];

    vector.getXYZ(&m[0], &m[4], &m[8]);
  }

  void getColumn(int c, double *x, double *y, double *z) {
    assert(c < 4);

    double *m = &(&m00_)[c];

    if (x) *x = m[0];
    if (y) *y = m[4];
    if (z) *z = m[8];
  }

  void getColumn(int c, CPoint3D &point) {
    assert(c < 4);

    double *m = &(&m00_)[c];

    point.x = m[0];
    point.y = m[4];
    point.z = m[8];
  }

  void getColumn(int c, CVector3D &vector) {
    assert(c < 4);

    double *m = &(&m00_)[c];

    vector = CVector3D(m[0], m[4], m[8]);
  }

  void getColumns(CVector3D &u, CVector3D &v, CVector3D &w) {
    getColumn(0, u);
    getColumn(1, v);
    getColumn(2, w);
  }

  //---------

  void setRow(int r, double x, double y, double z) {
    assert(r < 4);

    double *m = &(&m00_)[r*4];

    m[0] = x, m[1] = y, m[2] = z;
  }

  void setRow(int r, const CPoint3D &point) {
    assert(r < 4);

    double *m = &(&m00_)[r*4];

    m[0] = point.x, m[1] = point.y, m[2] = point.z;
  }

  void setRow(int r, const CVector3D &vector) {
    assert(r < 4);

    double *m = &(&m00_)[r*4];

    vector.getXYZ(&m[0], &m[1], &m[2]);
  }

  void getRow(int r, double *x, double *y, double *z) {
    assert(r < 4);

    double *m = &(&m00_)[r*4];

    if (x) *x = m[0];
    if (y) *y = m[1];
    if (z) *z = m[2];
  }

  void getRow(int r, CPoint3D &point) {
    assert(r < 4);

    double *m = &(&m00_)[r*4];

    point.x = m[0];
    point.y = m[1];
    point.z = m[2];
  }

  void getRow(int r, CVector3D &vector) {
    assert(r < 4);

    double *m = &(&m00_)[r*4];

    vector = CVector3D(m[0], m[1], m[2]);
  }

  //---------

  // CPoint3D can be expressed as a 1x4 matrix (x,y,z,1)
  void multiplyPoint(double xi, double yi, double zi,
                    double *xo, double *yo, double *zo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi + m03_;
    *yo = m10_*xi + m11_*yi + m12_*zi + m13_;
    *zo = m20_*xi + m21_*yi + m22_*zi + m23_;
  }

  void multiplyPoint(const CPoint3D &ipoint, CPoint3D &opoint) const {
    opoint.x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z + m03_;
    opoint.y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z + m13_;
    opoint.z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z + m23_;
  }

  CPoint3D multiplyPoint(const CPoint3D &ipoint) const {
    return CPoint3D(m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z + m03_,
                 m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z + m13_,
                 m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z + m23_);
  }

  void preMultiplyPoint(double xi, double yi, double zi,
                        double *xo, double *yo, double *zo) const {
    *xo = m00_*xi + m10_*yi + m20_*zi;
    *yo = m01_*xi + m11_*yi + m21_*zi;
    *zo = m02_*xi + m12_*yi + m22_*zi;
  }

  void preMultiplyPoint(const CPoint3D &ipoint, CPoint3D &opoint) const {
    opoint.x = m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z;
    opoint.z = m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z;
  }

  CPoint3D preMultiplyPoint(const CPoint3D &ipoint) const {
    return CPoint3D(m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z,
                 m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z,
                 m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z);
  }

  // CVector3D can be expressed as a 1x4 matrix (x,y,z,0)
  void multiplyVector(double xi, double yi, double zi,
                      double *xo, double *yo, double *zo) const {
    *xo = m00_*xi + m01_*yi + m02_*zi;
    *yo = m10_*xi + m11_*yi + m12_*zi;
    *zo = m20_*xi + m21_*yi + m22_*zi;
  }

  void multiplyVector(const CVector3D &ivector, CVector3D &ovector) const {
    CPoint3D ipoint, opoint;

    ivector.getXYZ(&ipoint.x, &ipoint.y, &ipoint.z);

    opoint.x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z;
    opoint.y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z;
    opoint.z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z;

    ovector = opoint;
  }

  CVector3D multiplyVector(const CVector3D &ivector) const {
    CPoint3D ipoint, opoint;

    ivector.getXYZ(&ipoint.x, &ipoint.y, &ipoint.z);

    opoint.x = m00_*ipoint.x + m01_*ipoint.y + m02_*ipoint.z;
    opoint.y = m10_*ipoint.x + m11_*ipoint.y + m12_*ipoint.z;
    opoint.z = m20_*ipoint.x + m21_*ipoint.y + m22_*ipoint.z;

    return CVector3D(opoint);
  }

  void preMultiplyVector(double xi, double yi, double zi,
                         double *xo, double *yo, double *zo) const {
    *xo = m00_*xi + m10_*yi + m20_*zi;
    *yo = m01_*xi + m11_*yi + m21_*zi;
    *zo = m02_*xi + m12_*yi + m22_*zi;
  }

  void preMultiplyVector(const CVector3D &ivector, CVector3D &ovector) const {
    CPoint3D ipoint, opoint;

    ivector.getXYZ(&ipoint.x, &ipoint.y, &ipoint.z);

    opoint.x = m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z;
    opoint.z = m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z;

    ovector = opoint;
  }

  CVector3D preMultiplyVector(const CVector3D &ivector) const {
    CPoint3D ipoint, opoint;

    ivector.getXYZ(&ipoint.x, &ipoint.y, &ipoint.z);

    opoint.x = m00_*ipoint.x + m10_*ipoint.y + m20_*ipoint.z;
    opoint.y = m01_*ipoint.x + m11_*ipoint.y + m21_*ipoint.z;
    opoint.z = m02_*ipoint.x + m12_*ipoint.y + m22_*ipoint.z;

    return CVector3D(opoint);
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

  CMatrix3D transposed() const {
    CMatrix3D matrix = *this;

    matrix.transpose();

    return matrix;
  }

  //------

  bool invert(CMatrix3D &imatrix) const {
    double d = determinant();

    if (::fabs(d) == 0.0)
      return false;

    double id = 1.0/d;

    imatrix.m00_ = det2x2(m11_, m12_, m21_, m22_)*id;
    imatrix.m10_ = det2x2(m12_, m10_, m22_, m20_)*id;
    imatrix.m20_ = det2x2(m10_, m11_, m20_, m21_)*id;

    imatrix.m01_ = det2x2(m21_, m22_, m01_, m02_)*id;
    imatrix.m11_ = det2x2(m22_, m20_, m02_, m00_)*id;
    imatrix.m21_ = det2x2(m20_, m21_, m00_, m01_)*id;

    imatrix.m02_ = det2x2(m01_, m02_, m11_, m12_)*id;
    imatrix.m12_ = det2x2(m02_, m00_, m12_, m10_)*id;
    imatrix.m22_ = det2x2(m00_, m01_, m10_, m11_)*id;

    double adjoint;

    adjoint = det3x3(m01_, m02_, m03_, m11_, m12_, m13_, m21_, m22_, m23_);

    imatrix.m03_ = -adjoint*id;

    adjoint = det3x3(m02_, m03_, m00_, m12_, m13_, m10_, m22_, m23_, m20_);

    imatrix.m13_ =  adjoint*id;

    adjoint = det3x3(m03_, m00_, m01_, m13_, m10_, m11_, m23_, m20_, m21_);

    imatrix.m23_ = -adjoint*id;

    imatrix.setBottomIdentity();

    return true;
  }

  CMatrix3D inverse() const {
    CMatrix3D imatrix;

    if (! invert(imatrix))
      assert(false && "Divide by 0.0");

    return imatrix;
  }

  //------

  double determinant() const {
    return (m00_*det2x2(m11_, m12_, m21_, m22_) -
            m01_*det2x2(m10_, m12_, m20_, m22_) +
            m02_*det2x2(m10_, m11_, m20_, m21_));
  }

  double upperDeterminant() const {
    return (m00_*det2x2(m11_, m12_, m21_, m22_) -
            m01_*det2x2(m10_, m12_, m20_, m22_) +
            m02_*det2x2(m10_, m11_, m20_, m21_));
  }

  //------

  void normalize() {
    double d = determinant();

    double id = 1.0/d;

    for (int i = 0; i < 9; ++i)
      (&m00_)[i] *= id;
  }

  CMatrix3D normalized() {
    CMatrix3D nmatrix = *this;

    nmatrix.normalize();

    return nmatrix;
  }

  //------

  void setTransform(double xmin1, double ymin1, double zmin1,
                    double xmax1, double ymax1, double zmax1,
                    double xmin2, double ymin2, double zmin2,
                    double xmax2, double ymax2, double zmax2) {
    double sx = (xmax2 - xmin2)/(xmax1 - xmin1);
    double sy = (ymax2 - ymin2)/(ymax1 - ymin1);
    double sz = (zmax2 - zmin2)/(zmax1 - zmin1);

    double tx = -xmin1*sx + xmin2;
    double ty = -ymin1*sy + ymin2;
    double tz = -zmin1*sz + zmin2;

    setInnerScale(sx, sy, sz);

    setOuterTranslate(tx, ty, tz);
  }

  //------

  static CMatrix3D *newIdentityMatrix() {
    CMatrix3D *m = new CMatrix3D();

    m->setIdentity();

    return m;
  }

  //------

  static bool solveAXeqB(const CMatrix3D &a, CPoint3D &x, const CPoint3D &b) {
    double det_a = a.determinant();

    if (::fabs(det_a) < 0.0)
      return false;

    double idet_a = 1.0/det_a;

    CMatrix3D t(a);

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

  //------

  static bool gaussianSolve(const CMatrix3D &m, const CVector3D &v, CVector3D &v1) {
    double d[12];

    memset(d, 0, sizeof(d));

    m.getValues(&d[ 0], &d[ 1], &d[ 2], &d[ 4], &d[ 5], &d[ 6], &d[ 8], &d[ 9], &d[10]);

    v.getXYZ(&d[3], &d[7], &d[11]);

    CGaussianMatrix<3,4> matrix(d);

    if (! matrix.solve(d))
      return false;

    v1.setXYZ(d[0], d[1], d[2]);

    return true;
  }

  //------

  void zero() { memset(&m00_, 0, 16*sizeof(double)); }

  //------

  CMatrix3D &operator+=(const CMatrix3D &b) {
    m00_ += b.m00_; m01_ += b.m01_; m02_ += b.m02_; m03_ += b.m03_;
    m10_ += b.m10_; m11_ += b.m11_; m12_ += b.m12_; m13_ += b.m13_;
    m20_ += b.m20_; m21_ += b.m21_; m22_ += b.m22_; m23_ += b.m23_;

    return *this;
  }

  CMatrix3D operator+(const CMatrix3D &b) const {
    CMatrix3D c = *this;

    c += b;

    return c;
  }

  CMatrix3D &operator-=(const CMatrix3D &b) {
    m00_ -= b.m00_; m01_ -= b.m01_, m02_ -= b.m02_, m03_ -= b.m03_;
    m10_ -= b.m10_; m11_ -= b.m11_; m12_ -= b.m12_; m13_ -= b.m13_;
    m20_ -= b.m20_; m21_ -= b.m21_; m22_ -= b.m22_; m23_ -= b.m23_;

    return *this;
  }

  CMatrix3D operator-(const CMatrix3D &b) const {
    CMatrix3D c = *this;

    c -= b;

    return c;
  }

  CMatrix3D &operator*=(const CMatrix3D &b) {
    CMatrix3D a;

    memcpy(&a.m00_, &m00_, 16*sizeof(double));

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

  CMatrix3D operator*(const CMatrix3D &b) const {
    CMatrix3D c = *this;

    c *= b;

    return c;
  }

  CMatrix3D &operator*=(double s) {
    m00_ *= s; m01_ *= s; m02_ *= s;
    m10_ *= s; m11_ *= s; m12_ *= s;
    m20_ *= s; m21_ *= s; m22_ *= s;

    return *this;
  }

  CMatrix3D operator*(double s) const {
    CMatrix3D c = *this;

    c *= s;

    return c;
  }

  friend CPoint3D operator*(const CMatrix3D &m, const CPoint3D &p) {
    CPoint3D p1;

    m.multiplyPoint(p, p1);

    return p1;
  }

  friend CPoint3D operator*(const CPoint3D &p, const CMatrix3D &m) {
    CPoint3D p1;

    m.preMultiplyPoint(p, p1);

    return p1;
  }

  friend CVector3D operator*(const CMatrix3D &m, const CVector3D &v) {
    CVector3D v1;

    m.multiplyVector(v, v1);

    return v1;
  }

  friend CVector3D operator*(const CVector3D &v, const CMatrix3D &m) {
    CVector3D v1;

    m.preMultiplyVector(v, v1);

    return v1;
  }

  CMatrix3D &operator/=(const CMatrix3D &b) {
    CMatrix3D bi;

    if (! b.invert(bi)) {
      assert(false && "Divide by 0.0");
      return *this;
    }

    return (*this) *= bi;
  }

  CMatrix3D operator/(const CMatrix3D &b) const {
    CMatrix3D c = *this;

    c /= b;

    return c;
  }

  //------

  void setValue(uint i, double value) {
    (&m00_)[i] = value;
  }

  void setValue(uint i, uint j, double value) {
    assert(i < 4 && j < 4);

    double &m = (&m00_)[4*j + i];

    m = value;
  }

  const double &getValue(uint i) const {
    assert(i < 16);

    return (&m00_)[i];
  }

  const double &getValue(uint i, uint j) const {
    assert(i < 4 && j < 4);

    const double &m = (&m00_)[4*j + i];

    return m;
  }

  const double &operator[](uint i) const {
    assert(i < 16);

    return (&m00_)[i];
  }

  double &operator[](uint i) {
    assert(i < 16);

    return (&m00_)[i];
  }

  double operator()(uint i, uint j) const {
    assert(i < 4 && j < 4);

    const double &m = (&m00_)[4*j + i];

    return m;
  }

  double &operator()(uint i, uint j) {
    assert(i < 4 && j < 4);

    double &m = (&m00_)[4*j + i];

    return m;
  }

  //------

 private:
  void setInnerIdentity() {
    m00_ = 1.0, m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = 1.0, m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = 1.0;
  }

  void setInnerScale(double sx, double sy, double sz) {
    m00_ = sx , m01_ = 0.0, m02_ = 0.0;
    m10_ = 0.0, m11_ = sy , m12_ = 0.0;
    m20_ = 0.0, m21_ = 0.0, m22_ = sz ;
  }

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

  void setOuterIdentity() {
    m03_ = 0.0; m13_ = 0.0; m23_ = 0.0;

    setBottomIdentity();
  }

  void setOuterTranslate(double tx, double ty, double tz) {
    m03_ = tx; m13_ = ty; m23_ = tz;

    setBottomIdentity();
  }

  void setBottomIdentity() {
    m30_ = 0.0, m31_ = 0.0, m32_ = 0.0, m33_ = 1.0;
  }

  static double det3x3(double m00, double m01, double m02,
                       double m10, double m11, double m12,
                       double m20, double m21, double m22) {
    return m00*(m11*m22 - m12*m21) -
           m01*(m10*m22 - m12*m20) +
           m02*(m10*m21 - m11*m20);
  }

  static double det2x2(double m00, double m01, double m10, double m11) {
    return m00*m11 - m01*m10;
  }

 private:
  double m00_ { 0 }, m01_ { 0 }, m02_ { 0 }, m03_ { 0 };
  double m10_ { 0 }, m11_ { 0 }, m12_ { 0 }, m13_ { 0 };
  double m20_ { 0 }, m21_ { 0 }, m22_ { 0 }, m23_ { 0 };
  double m30_ { 0 }, m31_ { 0 }, m32_ { 0 }, m33_ { 0 };
};

#endif
