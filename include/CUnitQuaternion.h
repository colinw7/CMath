#ifndef CUNIT_QUATERNION_H
#define CUNIT_QUATERNION_H

class CUnitQuaternion {
 public:
  CUnitQuaternion(double w, double x, double y, double z) :
   w_(w), x_(x), y_(y), z_(z) {
    normalize();
  }

  CUnitQuaternion(const CUnitQuaternion &q) :
   w_(q.w_), x_(q.x_), y_(q.y_), z_(q.z_) {
  }

  CUnitQuaternion(const CQuaternion &q) :
   w_(q.w_), x_(q.x_), y_(q.y_), z_(q.z_) {
    normalize();
  }

 ~CUnitQuaternion() { }

  CQuaternion quaternion() const {
    return CQuaternion(w_, x_, y_, z_);
  }

  CUnitQuaternion &operator=(const CUnitQuaternion &q) {
    w_ = q.w_; x_ = q.x_; y_ = q.y_; z_ = q.z_;

    return *this;
  }

  CUnitQuaternion &operator=(const CQuaternion &q) {
    w_ = q.w_; x_ = q.x_; y_ = q.y_; z_ = q.z_;

    normalize();

    return *this;
  }

  CUnitQuaternion &operator+=(const CUnitQuaternion &rhs) {
    w_ += rhs.w_;
    x_ += rhs.x_;
    y_ += rhs.y_;
    z_ += rhs.z_;

    normalize();

    return *this;
  }

  CUnitQuaternion operator+(const CUnitQuaternion &rhs) const {
    CUnitQuaternion lhs = *this;

    lhs += rhs;

    return lhs;
  }

  CUnitQuaternion &operator-=(const CUnitQuaternion &rhs) {
    w_ -= rhs.w_;
    x_ -= rhs.x_;
    y_ -= rhs.y_;
    z_ -= rhs.z_;

    normalize();

    return *this;
  }

  CUnitQuaternion operator-(const CUnitQuaternion &rhs) const {
    CUnitQuaternion lhs = *this;

    lhs -= rhs;

    return lhs;
  }

  CUnitQuaternion &operator*=(const CUnitQuaternion &rhs) {
    CUnitQuaternion lhs = *this;

    w_ = lhs.w_*rhs.w_ - lhs.x_*rhs.x_ - lhs.y_*rhs.y_ - lhs.z_*rhs.z_;
    x_ = lhs.w_*rhs.x_ + lhs.x_*rhs.w_ + lhs.y_*rhs.z_ - lhs.z_*rhs.y_;
    y_ = lhs.w_*rhs.y_ - lhs.x_*rhs.z_ + lhs.y_*rhs.w_ + lhs.z_*rhs.x_;
    z_ = lhs.w_*rhs.z_ + lhs.x_*rhs.y_ - lhs.y_*rhs.x_ + lhs.z_*rhs.w_;

    normalize();

    return *this;
  }

  CUnitQuaternion operator*(const CUnitQuaternion &rhs) const {
    CUnitQuaternion lhs = *this;

    lhs *= rhs;

    return lhs;
  }

  CUnitQuaternion operator-() const {
    return CUnitQuaternion(-w_, -x_, -y_, -z_);
  }

  const CUnitQuaternion &invert() {
    x_ = -x_; y_ = -y_; z_ = -z_;

    return *this;
  }

  CUnitQuaternion inverted() const {
    return CUnitQuaternion(w_, -x_, -y_, -z_);
  }

  CUnitQuaternion exp() const {
    double angle = ::sqrt(x_*x_ + y_*y_ + z_*z_);

    double s = ::sin(angle);

    CUnitQuaternion q(0.0, 0.0, 0.0, 0.0);

    q.w_ = ::cos(angle);

    if (::fabs(s) >= 1E-6) {
      double coeff = s/angle;

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

  CUnitQuaternion log() const {
    CUnitQuaternion q(0.0, 0.0, 0.0, 0.0);

    q.w_ = 0.0;

    if (::fabs(w_) < 1.0) {
      double angle = acos(w_);

      double s = ::sin(angle);

      if (::fabs(s) >= 1E-6) {
        double coeff = angle/s;

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

  static void intermediate(const CUnitQuaternion &q0,
                           const CUnitQuaternion &q1,
                           const CUnitQuaternion &q2,
                           CUnitQuaternion &a,
                           CUnitQuaternion &b) {
    CUnitQuaternion q0i = q0.inverted();
    CUnitQuaternion q1i = q1.inverted();

    CQuaternion p0 = q0i.quaternion()*q1.quaternion();
    CQuaternion p1 = q1i.quaternion()*q2.quaternion();

    CQuaternion arg = quarter_*(p0.log() - p1.log());

    CQuaternion argi = -arg;

    a = q1*arg .exp();
    b = q1*argi.exp();
  }

 private:
  void normalize() {
    double l = ::sqrt(w_*w_ + x_*x_ + y_*y_ + z_*z_);

    if (l <= 0.0)
      assert(false && "Divide by zero");
    else {
      double li = 1.0/l;

      w_ *= li; x_ *= li; y_ *= li; z_ *= li;
    }
  }

 private:
  double w_ { 0 }, x_ { 0 }, y_ { 0 }, z_ { 0 };
};

#endif
