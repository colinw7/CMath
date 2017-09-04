#ifndef CQUATERNION_H
#define CQUATERNION_H

#include <CMatrix3D.h>
#include <CPoint3D.h>
#include <CVector3D.h>

class CQuaternion {
 public:
  CQuaternion() { }

  explicit CQuaternion(double w, double x=0.0, double y=0.0, double z=0.0) :
   w_(w), v_(x, y, z) {
  }

  CQuaternion(double w, const CVector3D &v) :
   w_(w), v_(v) {
  }

  explicit CQuaternion(const CVector3D &v) :
   w_(0.0), v_(v) {
  }

  CQuaternion(const CQuaternion &q) :
   w_(q.w_), v_(q.v_) {
  }

 ~CQuaternion() { }

  double getW() const { return w_; }
  double getX() const { return v_.getX(); }
  double getY() const { return v_.getY(); }
  double getZ() const { return v_.getZ(); }

  void setW(double w) { w_ = w; }
  void setX(double x) { v_.setX(x); }
  void setY(double y) { v_.setY(y); }
  void setZ(double z) { v_.setZ(z); }

  CQuaternion &operator=(const CQuaternion &q) {
    w_ = q.w_; v_ = q.v_;

    return *this;
  }

  //----------

  // operator +=, -=, *=, /=

  CQuaternion &operator+=(const CQuaternion &rhs) {
    w_ += rhs.w_; v_ += rhs.v_;

    return *this;
  }

  CQuaternion &operator+=(double w) {
    w_ += w;

    return *this;
  }

  CQuaternion &operator-=(const CQuaternion &rhs) {
    w_ -= rhs.w_; v_ -= rhs.v_;

    return *this;
  }

  CQuaternion &operator-=(double w) {
    w_ -= w;

    return *this;
  }

  CQuaternion &operator*=(const CQuaternion &rhs) {
    double t = w_*rhs.w_ - v_.dotProduct(rhs.v_);

    v_ = w_*rhs.v_ + v_*rhs.w_ + v_.crossProduct(rhs.v_);
    w_ = t;

    return *this;
  }

  CQuaternion &operator*=(double s) {
    w_ *= s; v_ *= s;

    return *this;
  }

  CQuaternion &operator/=(const CQuaternion &rhs) {
    return *this *= rhs.inverted();
  }

  CQuaternion &operator/=(double rhs) {
    double irhs = 1.0/rhs;

    return *this *= irhs;
  }

  //----------

  // operator +, -, *, / (implemented in terms of +=, -=, *=, /=)

  CQuaternion operator+(const CQuaternion &rhs) const {
    CQuaternion t = *this; t += rhs; return t;
  }

  CQuaternion operator+(double rhs) {
    CQuaternion t = *this; t += rhs; return t;
  }

  CQuaternion operator-(const CQuaternion &rhs) const {
    CQuaternion t = *this; t -= rhs; return t;
  }

  CQuaternion operator-(double rhs) {
    CQuaternion t = *this; t -= rhs; return t;
  }

  CQuaternion operator*(const CQuaternion &rhs) const {
    CQuaternion t = *this; t *= rhs; return t;
  }

  CQuaternion operator*(double rhs) {
    CQuaternion t = *this; t *= rhs; return t;
  }

  CQuaternion operator/(const CQuaternion &rhs) const {
    CQuaternion t = *this; t /= rhs; return t;
  }

  CQuaternion operator/(double rhs) {
    CQuaternion t = *this; t /= rhs; return t;
  }

  //----------

  CQuaternion operator+() const {
    return *this;
  }

  CQuaternion operator-() const {
    return CQuaternion(-w_, -v_);
  }

  //----------

  friend CQuaternion operator+(double w, const CQuaternion &q) {
    return CQuaternion(w + q.w_, q.v_);
  }

  friend CQuaternion operator-(double w, const CQuaternion &q) {
    return CQuaternion(w - q.w_, -q.v_);
  }

  friend CQuaternion operator*(double s, const CQuaternion &q) {
    return CQuaternion(s*q.w_, s*q.v_);
  }

  friend CQuaternion operator/(double s, const CQuaternion &q) {
    return s*q.inverted();
  }

  //-----

  double operator^(const CQuaternion &rhs) const {
    return dotProduct(rhs);
  }

  CQuaternion operator~() const {
    return conjugated();
  }

  //-----

  int cmp(const CQuaternion &rhs) const {
    double dw = w_ - rhs.w_;

    if (dw != 0) return dw;

    return v_.cmp(rhs.v_);
  }

  //-----

  CQuaternion sqr() const {
    return (*this)*(*this);
  }

  //-----

  double length() const {
    return std::sqrt(lengthSqr());
  }

  double lengthSqr() const {
    return dotProductSelf();
  }

  double normal() const {
    return length();
  }

  double normalSqr() const {
    return lengthSqr();
  }

  const CQuaternion &normalize() {
    double l = length();

    if (l <= 0.0) {
      assert(false && "Divide by zero");
      return *this;
    }

    double li = 1.0/l;

    w_ *= li; v_ *= li;

    return *this;
  }

  CQuaternion normalized() const {
    CQuaternion t = *this;

    return t.normalize();
  }

  const CQuaternion &invert() {
    double l = dotProductSelf();

    if (l <= 0.0) {
      assert(false && "Divide by zero");
      return *this;
    }

    double li = 1.0/l;

    w_ *= li; v_ *= -li;

    return *this;
  }

  CQuaternion inverted() const {
    CQuaternion t = *this;

    return t.invert();
  }

  const CQuaternion &conjugate() {
    v_ = -v_;

    return *this;
  }

  CQuaternion conjugated() const {
    CQuaternion t = *this;

    return t.conjugate();
  }

  const CQuaternion &zero() {
    w_ = 0.0; v_.zero();

    return *this;
  }

  double dotProduct(const CQuaternion &rhs) const {
    return w_*rhs.w_ + v_.dotProduct(rhs.v_);
  }

  double dotProductSelf() const {
    return w_*w_ + v_.dotProductSelf();
  }

  void fromRotationMatrix(const CMatrix3D &matrix) {
    double a, b, c, d, e, f, g, h, i, tx, ty, tz;

    matrix.getValues(&a, &b, &c, &d, &e, &f, &g, &h, &i, &tx, &ty, &tz);

    double trace = a + e + i + 1.0;

    if (trace > 0.0) {
      double s = 0.5/std::sqrt(trace);

      w_ = trace*s;
      v_ = CVector3D(h - f, c - g, d - b)*s;
    }
    else {
      if      (a > e && a > i) {
        trace = a - e - i + 1.0;

        double s = 0.5/std::sqrt(trace);

        w_ = (h - f)*s;
        v_ = CVector3D(trace, d + b, g + c)*s;
      }
      else if (e > i) {
        trace = e - i - a + 1.0;

        double s = 0.5/std::sqrt(trace);

        w_ = (c - g)*s;
        v_ = CVector3D(b + d, trace, h + f)*s;
      }
      else {
        trace = i - a - e + 1.0;

        double s = 0.5/std::sqrt(trace);

        w_ = (d - b)*s;
        v_ = CVector3D(c + g, f + h, trace)*s;
      }
    }
  }

  void toRotationMatrix(CMatrix3D &matrix) const {
    CPoint3D v = v_.point();

    double a = 1.0 - 2.0*(v.y*v.y + v.z*v.z);
    double b =       2.0*(v.x*v.y -  w_*v.z);
    double c =       2.0*(v.x*v.z +  w_*v.y);
    double d =       2.0*(v.x*v.y +  w_*v.z);
    double e = 1.0 - 2.0*(v.x*v.x + v.z*v.z);
    double f =       2.0*(v.y*v.z -  w_*v.x);
    double g =       2.0*(v.x*v.z -  w_*v.y);
    double h =       2.0*(v.y*v.z +  w_*v.x);
    double i = 1.0 - 2.0*(v.x*v.x + v.y*v.y);

    matrix.setValues(a, b, c, d, e, f, g, h, i, 0.0, 0.0, 0.0);
  }

  void fromAngleAxis(double angle, const CPoint3D &axis) {
    fromAngleAxis(angle, CVector3D(axis));
  }

  void fromAngleAxis(double angle, const CVector3D &axis) {
    double angle2 = 0.5*angle;

    double s = std::sin(angle2);
    double c = std::cos(angle2);

    w_ = c;
    v_ = s*axis.normalized();
  }

  void toAngleAxis(double &angle, CPoint3D &axis) const {
    CVector3D v;

    toAngleAxis(angle, v);

    axis = v.point();
  }

  void toAngleAxis(double &angle, CVector3D &axis) const {
    angle = 2.0*acos(w_);

    axis = v_.normalized();
  }

  void fromAxes(const CPoint3D *axis) {
    CMatrix3D matrix(axis[0].x, axis[0].y, axis[0].z,
                  axis[1].x, axis[1].y, axis[1].z,
                  axis[2].x, axis[2].y, axis[2].z,
                  0.0, 0.0, 0.0);

    fromRotationMatrix(matrix);
  }

  void toAxes(CPoint3D *axis) const {
    CMatrix3D matrix;

    toRotationMatrix(matrix);

    double tx, ty, tz;

    matrix.getValues(&axis[0].x, &axis[0].y, &axis[0].z,
                     &axis[1].x, &axis[1].y, &axis[1].z,
                     &axis[2].x, &axis[2].y, &axis[2].z,
                     &tx, &ty, &tz);
  }

  void fromEulerZYX(double theta_z, double theta_y, double theta_x) {
    double cos_z_2 = 0.5*std::cos(theta_z);
    double cos_y_2 = 0.5*std::cos(theta_y);
    double cos_x_2 = 0.5*std::cos(theta_x);

    double sin_z_2 = 0.5*std::sin(theta_z);
    double sin_y_2 = 0.5*std::sin(theta_y);
    double sin_x_2 = 0.5*std::sin(theta_x);

    w_ = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;

    double x = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
    double y = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
    double z = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;

    v_ = CVector3D(x, y, z);
  }

  static CQuaternion lerp(double t, const CQuaternion &p, const CQuaternion &q) {
    CQuaternion r = p + t*(q - p);

    return r.normalize();
  }

  static CQuaternion slerp(double t, const CQuaternion &p, const CQuaternion &q) {
    double c = p.dotProduct(q);

    // Special case: q1 and q2 are the same, so just return one of them.
    // This also catches the case where cosomega is very slightly > 1.0
    if (c >= 1.0)
      return p;

    CQuaternion q1(q);

    if (c < 0.0) {
      c  = -c;
      q1 = -q;
    }

    double tol = double(95)/100;

    if (c < tol) {
      double angle = acos(c);

      if (::abs(angle) < 1E-6)
        return p;

      double s = std::sin(angle);

      double si = 1.0/s;

      double coeff0 = std::sin((1.0 - t)*angle)*si;
      double coeff1 = std::sin(t*angle)*si;

      return coeff0*p + coeff1*q1;
    }
    else
      return lerp(t, p, q1);
  }

  static CQuaternion slerpNoInvert(double t, const CQuaternion &p,
                                    const CQuaternion &q) {
    double c = p.dotProduct(q);

    double tol = double(95)/100;

    if (c > -tol && c < tol) {
      double angle = acos(c);

      if (::abs(angle) < 1E-6)
        return p;

      double s = std::sin(angle);

      double si = 1.0/s;

      double coeff0 = std::sin((1.0 - t)*angle)*si;
      double coeff1 = std::sin(t*angle)*si;

      return coeff0*p + coeff1*q;
    }
    else
      return lerp(t, p, q);
  }

  static CQuaternion slerpExtraSpins(double t, const CQuaternion &p,
                                      const CQuaternion &q, int nspins) {
    double c = p.dotProduct(q);

    double angle = acos(c);

    if (::abs(angle) < 1E-6)
      return p;

    double s = std::sin(angle);

    double phase = M_PI*nspins*t;

    double si = 1.0/s;

    double coeff0 = std::sin((1.0 - t)*angle - phase)*si;
    double coeff1 = std::sin(t*angle + phase)*si;

    return coeff0*p + coeff1*q;
  }

  static CQuaternion squad(double t, const CQuaternion &p,
                            const CQuaternion &a, const CQuaternion &b,
                            const CQuaternion &q) {
    CQuaternion slerpp = slerp(t, p, q);
    CQuaternion slerpq = slerp(t, a, b);

    double tslerp = 2.0*t*(1.0 - t);

    return slerp(tslerp, slerpp, slerpq);
  }

  static CQuaternion rotationArc(const CVector3D &v0, const CVector3D &v1) {
    CQuaternion q;

    CVector3D c = v0.crossProduct(v1);
    double      d = v0.dotProduct(v1);
    double      s = std::sqrt((1.0 + d)*2.0);

    double s1 = 1.0/s;

    q.v_ = c*s1;
    q.w_ = s / 2.0;

    return q;
  }

  static void toMatrix(const CQuaternion &q, const CVector3D &v, CMatrix3D &matrix) {
    CMatrix3D m;

    q.toRotationMatrix(m);

    matrix.translate(v.getX(), v.getY(), v.getZ());
  }

  CQuaternion log() const {
    if (::abs(w_) < 1.0) {
      double angle = acos(w_);

      double s = std::sin(angle);

      if (::abs(s) >= 1E-6) {
        double coeff = angle/s;

        return CQuaternion(0.0, coeff*v_);
      }
    }

    return CQuaternion(0.0, 0.0, 0.0, 0.0);
  }

  CQuaternion exp() const {
    double angle = v_.length();

    double s = std::sin(angle);

    if (::abs(s) >= 1E-6) {
      double coeff = s/angle;

      return CQuaternion(std::cos(angle), coeff*v_);
    }

    return CQuaternion(0.0, 0.0, 0.0, 0.0);
  }

  void print(ostream &os) const {
    os << "[" << w_ << ", " << v_ << "]";
  }

 private:
  double    w_;
  CVector3D v_;
};

//------

inline ostream &operator<<(ostream &os, const CQuaternion &q) {
  q.print(os);

  return os;
}

inline bool operator==(const CQuaternion &lhs, const CQuaternion &rhs) {
  return (lhs.cmp(rhs) == 0);
}

inline bool operator!=(const CQuaternion &lhs, const CQuaternion &rhs) {
  return ! operator==(lhs, rhs);
}

inline bool operator<(const CQuaternion &lhs, const CQuaternion &rhs) {
  return (lhs.cmp(rhs) < 0);
}

inline bool operator>(const CQuaternion &lhs, const CQuaternion &rhs) {
  return (lhs.cmp(rhs) > 0);
}

inline bool operator>=(const CQuaternion &lhs, const CQuaternion &rhs) {
  return ! operator<(lhs, rhs);
}

inline bool operator<=(const CQuaternion &lhs, const CQuaternion &rhs) {
  return ! operator>(lhs, rhs);
}

inline double dotProduct(const CQuaternion &lhs, const CQuaternion &rhs) {
  return lhs.dotProduct(rhs);
}

inline CQuaternion mul3(const CQuaternion &q1, const CQuaternion &q2, const CQuaternion &q3) {
  return (q1*q2)*q3;
}

#endif
