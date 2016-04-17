#ifndef CQUATERNION_H
#define CQUATERNION_H

#include <CThrow.h>
#include <CMatrix3D.h>
#include <CPoint3D.h>
#include <CVector3D.h>

template<typename T>
class CQuaternionT {
 private:
  typedef CVector3DT<T> Vector;
  typedef CPoint3DT<T>  Point;
  typedef CMatrix3DT<T> Matrix;

  T      w_;
  Vector v_;

 public:
  CQuaternionT() { }

  explicit CQuaternionT(T w, T x=0.0, T y=0.0, T z=0.0) :
   w_(w), v_(x, y, z) {
  }

  CQuaternionT(T w, const Vector &v) :
   w_(w), v_(v) {
  }

  explicit CQuaternionT(const Vector &v) :
   w_(0.0), v_(v) {
  }

  CQuaternionT(const CQuaternionT &q) :
   w_(q.w_), v_(q.v_) {
  }

 ~CQuaternionT() { }

  T getW() const { return w_; }
  T getX() const { return v_.getX(); }
  T getY() const { return v_.getY(); }
  T getZ() const { return v_.getZ(); }

  void setW(T w) { w_ = w; }
  void setX(T x) { v_.setX(x); }
  void setY(T y) { v_.setY(y); }
  void setZ(T z) { v_.setZ(z); }

  CQuaternionT &operator=(const CQuaternionT &q) {
    w_ = q.w_; v_ = q.v_;

    return *this;
  }

  //----------

  // operator +=, -=, *=, /=

  CQuaternionT &operator+=(const CQuaternionT &rhs) {
    w_ += rhs.w_; v_ += rhs.v_;

    return *this;
  }

  CQuaternionT &operator+=(T w) {
    w_ += w;

    return *this;
  }

  CQuaternionT &operator-=(const CQuaternionT &rhs) {
    w_ -= rhs.w_; v_ -= rhs.v_;

    return *this;
  }

  CQuaternionT &operator-=(T w) {
    w_ -= w;

    return *this;
  }

  CQuaternionT &operator*=(const CQuaternionT &rhs) {
    T t = w_*rhs.w_ - v_.dotProduct(rhs.v_);

    v_ = w_*rhs.v_ + v_*rhs.w_ + v_.crossProduct(rhs.v_);
    w_ = t;

    return *this;
  }

  CQuaternionT &operator*=(T s) {
    w_ *= s; v_ *= s;

    return *this;
  }

  CQuaternionT &operator/=(const CQuaternionT &rhs) {
    return *this *= rhs.inverted();
  }

  CQuaternionT &operator/=(T rhs) {
    T irhs = 1.0/rhs;

    return *this *= irhs;
  }

  //----------

  // operator +, -, *, / (implemented in terms of +=, -=, *=, /=)

  CQuaternionT operator+(const CQuaternionT &rhs) const {
    CQuaternionT t = *this; t += rhs; return t;
  }

  CQuaternionT operator+(T rhs) {
    CQuaternionT t = *this; t += rhs; return t;
  }

  CQuaternionT operator-(const CQuaternionT &rhs) const {
    CQuaternionT t = *this; t -= rhs; return t;
  }

  CQuaternionT operator-(T rhs) {
    CQuaternionT t = *this; t -= rhs; return t;
  }

  CQuaternionT operator*(const CQuaternionT &rhs) const {
    CQuaternionT t = *this; t *= rhs; return t;
  }

  CQuaternionT operator*(T rhs) {
    CQuaternionT t = *this; t *= rhs; return t;
  }

  CQuaternionT operator/(const CQuaternionT &rhs) const {
    CQuaternionT t = *this; t /= rhs; return t;
  }

  CQuaternionT operator/(T rhs) {
    CQuaternionT t = *this; t /= rhs; return t;
  }

  //----------

  CQuaternionT operator+() const {
    return *this;
  }

  CQuaternionT operator-() const {
    return CQuaternionT(-w_, -v_);
  }

  //----------

  friend CQuaternionT operator+(T w, const CQuaternionT &q) {
    return CQuaternionT(w + q.w_, q.v_);
  }

  friend CQuaternionT operator-(T w, const CQuaternionT &q) {
    return CQuaternionT(w - q.w_, -q.v_);
  }

  friend CQuaternionT operator*(T s, const CQuaternionT &q) {
    return CQuaternionT(s*q.w_, s*q.v_);
  }

  friend CQuaternionT operator/(T s, const CQuaternionT &q) {
    return s*q.inverted();
  }

  //-----

  T operator^(const CQuaternionT &rhs) const {
    return dotProduct(rhs);
  }

  CQuaternionT operator~() const {
    return conjugated();
  }

  //-----

  int cmp(const CQuaternionT &rhs) const {
    T dw = w_ - rhs.w_;

    if (dw != 0) return dw;

    return v_.cmp(rhs.v_);
  }

  //-----

  CQuaternionT sqr() const {
    return (*this)*(*this);
  }

  //-----

  T length() const {
    return std::sqrt(lengthSqr());
  }

  T lengthSqr() const {
    return dotProductSelf();
  }

  T normal() const {
    return length();
  }

  T normalSqr() const {
    return lengthSqr();
  }

  const CQuaternionT &normalize() {
    T l = length();

    if (l <= 0.0) {
      CTHROW("Divide by zero");
      return *this;
    }

    T li = 1.0/l;

    w_ *= li; v_ *= li;

    return *this;
  }

  CQuaternionT normalized() const {
    CQuaternionT t = *this;

    return t.normalize();
  }

  const CQuaternionT &invert() {
    T l = dotProductSelf();

    if (l <= 0.0) {
      CTHROW("Divide by zero");
      return *this;
    }

    T li = 1.0/l;

    w_ *= li; v_ *= -li;

    return *this;
  }

  CQuaternionT inverted() const {
    CQuaternionT t = *this;

    return t.invert();
  }

  const CQuaternionT &conjugate() {
    v_ = -v_;

    return *this;
  }

  CQuaternionT conjugated() const {
    CQuaternionT t = *this;

    return t.conjugate();
  }

  const CQuaternionT &zero() {
    w_ = 0.0; v_.zero();

    return *this;
  }

  T dotProduct(const CQuaternionT &rhs) const {
    return w_*rhs.w_ + v_.dotProduct(rhs.v_);
  }

  T dotProductSelf() const {
    return w_*w_ + v_.dotProductSelf();
  }

  void fromRotationMatrix(const Matrix &matrix) {
    T a, b, c, d, e, f, g, h, i, tx, ty, tz;

    matrix.getValues(&a, &b, &c, &d, &e, &f, &g, &h, &i, &tx, &ty, &tz);

    T trace = a + e + i + 1.0;

    if (trace > 0.0) {
      T s = 0.5/std::sqrt(trace);

      w_ = trace*s;
      v_ = Vector(h - f, c - g, d - b)*s;
    }
    else {
      if      (a > e && a > i) {
        trace = a - e - i + 1.0;

        T s = 0.5/std::sqrt(trace);

        w_ = (h - f)*s;
        v_ = Vector(trace, d + b, g + c)*s;
      }
      else if (e > i) {
        trace = e - i - a + 1.0;

        T s = 0.5/std::sqrt(trace);

        w_ = (c - g)*s;
        v_ = Vector(b + d, trace, h + f)*s;
      }
      else {
        trace = i - a - e + 1.0;

        T s = 0.5/std::sqrt(trace);

        w_ = (d - b)*s;
        v_ = Vector(c + g, f + h, trace)*s;
      }
    }
  }

  void toRotationMatrix(Matrix &matrix) const {
    Point v = v_.point();

    T a = 1.0 - 2.0*(v.y*v.y + v.z*v.z);
    T b =       2.0*(v.x*v.y -  w_*v.z);
    T c =       2.0*(v.x*v.z +  w_*v.y);
    T d =       2.0*(v.x*v.y +  w_*v.z);
    T e = 1.0 - 2.0*(v.x*v.x + v.z*v.z);
    T f =       2.0*(v.y*v.z -  w_*v.x);
    T g =       2.0*(v.x*v.z -  w_*v.y);
    T h =       2.0*(v.y*v.z +  w_*v.x);
    T i = 1.0 - 2.0*(v.x*v.x + v.y*v.y);

    matrix.setValues(a, b, c, d, e, f, g, h, i, 0.0, 0.0, 0.0);
  }

  void fromAngleAxis(T angle, const Point &axis) {
    fromAngleAxis(angle, Vector(axis));
  }

  void fromAngleAxis(T angle, const Vector &axis) {
    T angle2 = 0.5*angle;

    T s = std::sin(angle2);
    T c = std::cos(angle2);

    w_ = c;
    v_ = s*axis.normalized();
  }

  void toAngleAxis(T &angle, Point &axis) const {
    Vector v;

    toAngleAxis(angle, v);

    axis = v.point();
  }

  void toAngleAxis(T &angle, Vector &axis) const {
    angle = 2.0*acos(w_);

    axis = v_.normalized();
  }

  void fromAxes(const Point *axis) {
    Matrix matrix(axis[0].x, axis[0].y, axis[0].z,
                  axis[1].x, axis[1].y, axis[1].z,
                  axis[2].x, axis[2].y, axis[2].z,
                  0.0, 0.0, 0.0);

    fromRotationMatrix(matrix);
  }

  void toAxes(Point *axis) const {
    Matrix matrix;

    toRotationMatrix(matrix);

    T tx, ty, tz;

    matrix.getValues(&axis[0].x, &axis[0].y, &axis[0].z,
                     &axis[1].x, &axis[1].y, &axis[1].z,
                     &axis[2].x, &axis[2].y, &axis[2].z,
                     &tx, &ty, &tz);
  }

  void fromEulerZYX(T theta_z, T theta_y, T theta_x) {
    T cos_z_2 = 0.5*std::cos(theta_z);
    T cos_y_2 = 0.5*std::cos(theta_y);
    T cos_x_2 = 0.5*std::cos(theta_x);

    T sin_z_2 = 0.5*std::sin(theta_z);
    T sin_y_2 = 0.5*std::sin(theta_y);
    T sin_x_2 = 0.5*std::sin(theta_x);

    w_ = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;

    T x = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
    T y = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
    T z = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;

    v_ = Vector(x, y, z);
  }

  static CQuaternionT lerp(T t, const CQuaternionT &p, const CQuaternionT &q) {
    CQuaternionT r = p + t*(q - p);

    return r.normalize();
  }

  static CQuaternionT slerp(T t, const CQuaternionT &p, const CQuaternionT &q) {
    T c = p.dotProduct(q);

    // Special case: q1 and q2 are the same, so just return one of them.
    // This also catches the case where cosomega is very slightly > 1.0
    if (c >= 1.0)
      return p;

    CQuaternionT q1(q);

    if (c < 0.0) {
      c  = -c;
      q1 = -q;
    }

    T tol = T(95)/100;

    if (c < tol) {
      T angle = acos(c);

      if (::abs(angle) < 1E-6)
        return p;

      T s = std::sin(angle);

      T si = 1.0/s;

      T coeff0 = std::sin((1.0 - t)*angle)*si;
      T coeff1 = std::sin(t*angle)*si;

      return coeff0*p + coeff1*q1;
    }
    else
      return lerp(t, p, q1);
  }

  static CQuaternionT slerpNoInvert(T t, const CQuaternionT &p,
                                    const CQuaternionT &q) {
    T c = p.dotProduct(q);

    T tol = T(95)/100;

    if (c > -tol && c < tol) {
      T angle = acos(c);

      if (::abs(angle) < 1E-6)
        return p;

      T s = std::sin(angle);

      T si = 1.0/s;

      T coeff0 = std::sin((1.0 - t)*angle)*si;
      T coeff1 = std::sin(t*angle)*si;

      return coeff0*p + coeff1*q;
    }
    else
      return lerp(t, p, q);
  }

  static CQuaternionT slerpExtraSpins(T t, const CQuaternionT &p,
                                      const CQuaternionT &q, int nspins) {
    T c = p.dotProduct(q);

    T angle = acos(c);

    if (::abs(angle) < 1E-6)
      return p;

    T s = std::sin(angle);

    T phase = M_PI*nspins*t;

    T si = 1.0/s;

    T coeff0 = std::sin((1.0 - t)*angle - phase)*si;
    T coeff1 = std::sin(t*angle + phase)*si;

    return coeff0*p + coeff1*q;
  }

  static CQuaternionT squad(T t, const CQuaternionT &p,
                            const CQuaternionT &a, const CQuaternionT &b,
                            const CQuaternionT &q) {
    CQuaternionT slerpp = slerp(t, p, q);
    CQuaternionT slerpq = slerp(t, a, b);

    T tslerp = 2.0*t*(1.0 - t);

    return slerp(tslerp, slerpp, slerpq);
  }

  static CQuaternionT rotationArc(const Vector &v0, const Vector &v1) {
    CQuaternionT q;

    Vector c = v0.crossProduct(v1);
    T      d = v0.dotProduct(v1);
    T      s = std::sqrt((1.0 + d)*2.0);

    T s1 = 1.0/s;

    q.v_ = c*s1;
    q.w_ = s / 2.0;

    return q;
  }

  static void toMatrix(const CQuaternionT &q, const Vector &v, Matrix &matrix) {
    Matrix m;

    q.toRotationMatrix(m);

    matrix.translate(v.getX(), v.getY(), v.getZ());
  }

  CQuaternionT log() const {
    if (::abs(w_) < 1.0) {
      T angle = acos(w_);

      T s = std::sin(angle);

      if (::abs(s) >= 1E-6) {
        T coeff = angle/s;

        return CQuaternionT(0.0, coeff*v_);
      }
    }

    return CQuaternionT(0.0, 0.0, 0.0, 0.0);
  }

  CQuaternionT exp() const {
    T angle = v_.length();

    T s = std::sin(angle);

    if (::abs(s) >= 1E-6) {
      T coeff = s/angle;

      return CQuaternionT(std::cos(angle), coeff*v_);
    }

    return CQuaternionT(0.0, 0.0, 0.0, 0.0);
  }

  void print(ostream &os) const {
    os << "[" << w_ << ", " << v_ << "]";
  }
};

typedef CQuaternionT<double> CQuaternion;

//------

template<typename T>
inline ostream &operator<<(ostream &os, const CQuaternionT<T> &q) {
  q.print(os);

  return os;
}

template<typename T>
inline bool operator==(const CQuaternionT<T> &lhs, const CQuaternionT<T> &rhs) {
  return (lhs.cmp(rhs) == 0);
}

template<typename T>
inline bool operator!=(const CQuaternionT<T> &lhs, const CQuaternionT<T> &rhs) {
  return ! operator==(lhs, rhs);
}

template<typename T>
inline bool operator<(const CQuaternionT<T> &lhs, const CQuaternionT<T> &rhs) {
  return (lhs.cmp(rhs) < 0);
}

template<typename T>
inline bool operator>(const CQuaternionT<T> &lhs, const CQuaternionT<T> &rhs) {
  return (lhs.cmp(rhs) > 0);
}

template<typename T>
inline bool operator>=(const CQuaternionT<T> &lhs, const CQuaternionT<T> &rhs) {
  return ! operator<(lhs, rhs);
}

template<typename T>
inline bool operator<=(const CQuaternionT<T> &lhs, const CQuaternionT<T> &rhs) {
  return ! operator>(lhs, rhs);
}

template<typename T>
inline T dotProduct(const CQuaternionT<T> &lhs, const CQuaternionT<T> &rhs) {
  return lhs.dotProduct(rhs);
}

template<typename T>
inline CQuaternionT<T> mul3(const CQuaternionT<T> &q1, const CQuaternionT<T> &q2,
                            const CQuaternionT<T> &q3) {
  return (q1*q2)*q3;
}

#endif
