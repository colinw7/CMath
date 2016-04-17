#ifndef CUNIT_QUATERNION_H
#define CUNIT_QUATERNION_H

template<typename T>
class CUnitQuaternionT {
 private:
  T w_, x_, y_, z_;

 public:
  CUnitQuaternionT(T w, T x, T y, T z) :
   w_(w), x_(x), y_(y), z_(z) {
    normalize();
  }

  CUnitQuaternionT(const CUnitQuaternionT &q) :
   w_(q.w_), x_(q.x_), y_(q.y_), z_(q.z_) {
  }

  CUnitQuaternionT(const CQuaternionT<T> &q) :
   w_(q.w_), x_(q.x_), y_(q.y_), z_(q.z_) {
    normalize();
  }

 ~CUnitQuaternionT() { }

  CQuaternionT<T> quaternion() const {
    return CQuaternionT<T>(w_, x_, y_, z_);
  }

  CUnitQuaternionT &operator=(const CUnitQuaternionT &q) {
    w_ = q.w_; x_ = q.x_; y_ = q.y_; z_ = q.z_;

    return *this;
  }

  CUnitQuaternionT &operator=(const CQuaternionT<T> &q) {
    w_ = q.w_; x_ = q.x_; y_ = q.y_; z_ = q.z_;

    normalize();

    return *this;
  }

  CUnitQuaternionT &operator+=(const CUnitQuaternionT &rhs) {
    w_ += rhs.w_;
    x_ += rhs.x_;
    y_ += rhs.y_;
    z_ += rhs.z_;

    normalize();

    return *this;
  }

  CUnitQuaternionT operator+(const CUnitQuaternionT &rhs) const {
    CUnitQuaternionT lhs = *this;

    lhs += rhs;

    return lhs;
  }

  CUnitQuaternionT &operator-=(const CUnitQuaternionT &rhs) {
    w_ -= rhs.w_;
    x_ -= rhs.x_;
    y_ -= rhs.y_;
    z_ -= rhs.z_;

    normalize();

    return *this;
  }

  CUnitQuaternionT operator-(const CUnitQuaternionT &rhs) const {
    CUnitQuaternionT lhs = *this;

    lhs -= rhs;

    return lhs;
  }

  CUnitQuaternionT &operator*=(const CUnitQuaternionT &rhs) {
    CUnitQuaternionT lhs = *this;

    w_ = lhs.w_*rhs.w_ - lhs.x_*rhs.x_ - lhs.y_*rhs.y_ - lhs.z_*rhs.z_;
    x_ = lhs.w_*rhs.x_ + lhs.x_*rhs.w_ + lhs.y_*rhs.z_ - lhs.z_*rhs.y_;
    y_ = lhs.w_*rhs.y_ - lhs.x_*rhs.z_ + lhs.y_*rhs.w_ + lhs.z_*rhs.x_;
    z_ = lhs.w_*rhs.z_ + lhs.x_*rhs.y_ - lhs.y_*rhs.x_ + lhs.z_*rhs.w_;

    normalize();

    return *this;
  }

  CUnitQuaternionT operator*(const CUnitQuaternionT &rhs) const {
    CUnitQuaternionT lhs = *this;

    lhs *= rhs;

    return lhs;
  }

  CUnitQuaternionT operator-() const {
    return CUnitQuaternionT(-w_, -x_, -y_, -z_);
  }

  const CUnitQuaternionT &invert() {
    x_ = -x_; y_ = -y_; z_ = -z_;

    return *this;
  }

  CUnitQuaternionT inverted() const {
    return CUnitQuaternionT(w_, -x_, -y_, -z_);
  }

  CUnitQuaternionT exp() const {
    T angle = ::sqrt(x_*x_ + y_*y_ + z_*z_);

    T s = ::sin(angle);

    CUnitQuaternionT q(0.0, 0.0, 0.0, 0.0);

    q.w_ = ::cos(angle);

    if (::fabs(s) >= 1E-6) {
      T coeff = s/angle;

      q.x_ = coeff*x_;
      q.y_ = coeff*y_;
      q.z_ = coeff*z_;
    }
    else {
      q.x_ = x_;
      q.y_ = y_;
      q.z_ = z_;
    }

    q.normalize();

    return q;
  }

  CUnitQuaternionT log() const {
    CUnitQuaternionT q(0.0, 0.0, 0.0, 0.0);

    q.w_ = 0.0;

    if (::fabs(w_) < 1.0) {
      T angle = acos(w_);

      T s = ::sin(angle);

      if (::fabs(s) >= 1E-6) {
        T coeff = angle/s;

        q.x_ = coeff*x_;
        q.y_ = coeff*y_;
        q.z_ = coeff*z_;

        return q;
      }
    }

    q.x_ = x_;
    q.y_ = y_;
    q.z_ = z_;

    q.normalize();

    return q;
  }

  static void intermediate(const CUnitQuaternionT &q0,
                           const CUnitQuaternionT &q1,
                           const CUnitQuaternionT &q2,
                           CUnitQuaternionT &a,
                           CUnitQuaternionT &b) {
    CUnitQuaternionT q0i = q0.inverted();
    CUnitQuaternionT q1i = q1.inverted();

    CQuaternionT<T> p0 = q0i.quaternion()*q1.quaternion();
    CQuaternionT<T> p1 = q1i.quaternion()*q2.quaternion();

    CQuaternionT<T> arg = quarter_*(p0.log() - p1.log());

    CQuaternionT<T> argi = -arg;

    a = q1*arg .exp();
    b = q1*argi.exp();
  }

 private:
  void normalize() {
    T l = ::sqrt(w_*w_ + x_*x_ + y_*y_ + z_*z_);

    if (l <= 0.0)
      CTHROW("Divide by zero");
    else {
      T li = 1.0/l;

      w_ *= li; x_ *= li; y_ *= li; z_ *= li;
    }
  }
};

typedef CUnitQuaternionT<double> CUnitQuaternion;

#endif
