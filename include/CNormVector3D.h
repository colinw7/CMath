#ifndef CNORM_VECTOR_3D_H
#define CNORM_VECTOR_3D_H

#include <CPoint3D.h>
#include <CMathGen.h>

class CNormVector3D {
 public:
  CNormVector3D(double x=0.0, double y=0.0, double z=0.0, double w=1.0) :
   x_(x), y_(y), z_(z), w_(w), normalized_(false) {
  }

  CNormVector3D(const CNormVector3D &vector) :
    x_(vector.x_), y_(vector.y_), z_(vector.z_), w_(vector.w_)
    normalized_(vector.normalized_) {
  }

  explicit CNormVector3D(const CVector3D &vector) :
    x_(vector.getX()), y_(vector.getY()), z_(vector.getZ()), w_(1.0)
    normalized_(vector.getNormalized()) {
  }

  explicit CNormVector3D(const CPoint3D &point) :
    x_(point.x_), y_(point.y_), z_(point.z_), w_(1.0) normalized_(false) {
  }

  CNormVector3D(const CNormVector3D &vector1, const CNormVector3D &vector2) :
    x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_),
    z_(vector2.z_ - vector1.z_), w_(1.0), normalized_(false) {
  }

  CNormVector3D(const CVector3D &vector1, const CVector3D &vector2) :
    x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_),
    z_(vector2.z_ - vector1.z_), w_(1.0), normalized_(false) {
  }

  CNormVector3D(const CPoint3D &point1, const CPoint3D &point2) :
    x_(point2.x_ - point1.x_), y_(point2.y_ - point1.y_),
    z_(point2.z_ - point1.z_), w_(1.0), normalized_(false) {
  }

  double getX() const { return x_; }
  double getY() const { return y_; }
  double getZ() const { return z_; }
  double getW() const { return w_; }

  void getXYZW(double *x, double *y, double *z, double *w) const {
    *x = x_; *y = y_; *z = z_; *w = w_;
  }

  void getXYZ(double *x, double *y, double *z) const {
    *x = x_; *y = y_; *z = z_;
  }

  void setX(double x) {
    x_ = x;

    normalized_ = false;
  }

  void setY(double y) {
    y_ = y;

    normalized_ = false;
  }

  void setZ(double z) {
    z_ = z;

    normalized_ = false;
  }

  void setW(double w) {
    w_ = w;

    normalized_ = false;
  }

  void setXYZW(double x, double y, double z, double w) {
    x_ = x; y_ = y; z_ = z; w_ = w;

    normalized_ = false;
  }

  void setXYZ(double x, double y, double z) {
    x_ = x; y_ = y; z_ = z;

    normalized_ = false;
  }

  CPoint3D point() const {
    return CPoint3D(x_, y_, z_);
  }

  CVector3D vector() const {
    return CVector3D(x_, y_, z_);
  }

  double length() const {
    if (normalized_)
      return 1.0;
    else
      return ::sqrt(x_*x_ + y_*y_ + z_*z_);
  }

  double fastLength() const {
    if (normalized_)
      return 1.0;
    else
      return CMathGen::fastLength(x_, y_, z_);
  }

  double modulus() const {
    return length();
  }

  double lengthSqr() const {
    if (normalized_)
      return 1.0;
    else
      return (x_*x_ + y_*y_ + z_*z_);
  }

  CNormVector3D &zero() {
    x_ = 0.0; y_ = 0.0; z_ = 0.0; w_ = 1.0;

    normalized_ = false;

    return *this;
  }

  CNormVector3D unit() const {
    if (normalized_)
      return *this;

    double len = length();

    double factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    return CNormVector3D(x_*factor, y_*factor, z_*factor, 1.0);
  }

  CNormVector3D &normalize() {
    if (normalized_)
      return *this;

    double len = length();

    double factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    x_ *= factor;
    y_ *= factor;
    z_ *= factor;
    w_  = 1.0;

    normalized_ = true;

    return *this;
  }

  CNormVector3D normalized() const {
    if (normalized_)
      return *this;

    double len = length();

    double factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    return CNormVector3D(x_*factor, y_*factor, z_*factor, 1.0);
  }

  CNormVector3D &setMagnitude(double magnitude) {
    double factor = 0.0;

    if (normalized_)
      factor = magnitude;
    else {
      double len = length();

      if (len > 0.0)
        factor = magnitude/len;
    }

    x_ *= factor;
    y_ *= factor;
    z_ *= factor;
    w_  = 1.0;

    normalized_ = false;

    return *this;
  }

  double getDistance(const CNormVector3D &vector) {
    CNormVector3D diff = *this - vector;

    return diff.length();
  }

  void incX(double x) {
    x_ += x;

    normalized_ = false;
  }

  void incY(double y) {
    y_ += y;

    normalized_ = false;
  }

  void incZ(double z) {
    z_ += z;

    normalized_ = false;
  }

  void incW(double w) {
    w_ += w;

    normalized_ = false;
  }

  void decX(double x) {
    x_ -= x;

    normalized_ = false;
  }

  void decY(double y) {
    y_ -= y;

    normalized_ = false;
  }

  void decZ(double z) {
    z_ -= z;

    normalized_ = false;
  }

  void decW(double w) {
    w_ -= w;

    normalized_ = false;
  }

  double minComponent() {
    return min(min(x_, y_), z_);
  }

  double maxComponent() {
    return max(max(x_, y_), z_);
  }

  double minAbsComponent() {
    return min(min(fabs(x_), fabs(y_)), fabs(z_));
  }

  double maxAbsComponent() {
    return max(max(fabs(x_), fabs(y_)), fabs(z_));
  }

  CNormVector3D &operator=(const CVector3D &vector) {
    vector.getXYZ(&x_, &y_, &z_);

    normalized_ = false;

    return *this;
  }

  CNormVector3D &operator=(const CPoint3D &point) {
    x_ = point.x; y_ = point.y; z_ = point.z;

    normalized_ = false;

    return *this;
  }

  CNormVector3D operator+() const {
    return CNormVector3D(x_, y_, z_, 1.0);
  }

  CNormVector3D operator-() const {
    return CNormVector3D(-x_, -y_, -z_, 1.0);
  }

  friend bool operator==(const CNormVector3D &lhs,
                         const CNormVector3D &rhs) {
    return (lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_ &&
            lhs.z_ == rhs.z_ && lhs.w_ == rhs.w_);
  }

  friend bool operator!=(const CNormVector3D &lhs,
                         const CNormVector3D &rhs) {
    return (lhs.x_ != rhs.x_ || lhs.y_ != rhs.y_ ||
            lhs.z_ != rhs.z_ || lhs.w_ != rhs.w_);
  }

  CNormVector3D &operator+=(const CNormVector3D &rhs) {
    x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_; w_  = 1.0;

    normalized_ = false;

    return *this;
  }

  CNormVector3D operator+(const CNormVector3D &rhs) const {
    return CNormVector3D(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_, 1.0);
  }

  CNormVector3D &operator-=(const CNormVector3D &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_; w_ = 1.0;

    normalized_ = false;

    return *this;
  }

  CNormVector3D operator-(const CNormVector3D &rhs) const {
    return CNormVector3D(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_, 1.0);
  }

  CNormVector3D &operator*=(double rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs; w_  = 1.0;

    normalized_ = false;

    return *this;
  }

  friend CNormVector3D operator*(const CNormVector3D &lhs, double rhs) {
    return CNormVector3D(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs, 1.0);
  }

  friend CNormVector3D operator*(double lhs, const CNormVector3D &rhs) {
    return CNormVector3D(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_, 1.0);
  }

  CNormVector3D &operator/=(double rhs) {
    x_ /= rhs; y_ /= rhs; z_ /= rhs; w_  = 1.0;

    normalized_ = false;

    return *this;
  }

  CNormVector3D operator/(double rhs) const {
    return CNormVector3D(x_/rhs, y_/rhs, z_/rhs, 1.0);
  }

  friend CVector3D operator+(const CVector3D &lhs, const CNormVector3D &rhs) {
    return CVector3D(lhs.getX() + rhs.x_, lhs.getY() + rhs.y_,
                     lhs.getZ() + rhs.z_);
  }

  friend CPoint3D operator+(const CPoint3D &lhs, const CNormVector3D &rhs) {
    return CPoint3D(lhs.x + rhs.x_, lhs.y + rhs.y_, lhs.z + rhs.z_);
  }

  friend CVector3D operator-(const CVector3D &lhs, const CNormVector3D &rhs) {
    return CPoint3D(lhs.getX() - rhs.x_, lhs.getY() - rhs.y_,
                    lhs.getZ() - rhs.z_);
  }

  friend CPoint3D operator-(const CPoint3D &lhs, const CNormVector3D &rhs) {
    return CPoint3D(lhs.x - rhs.x_, lhs.y - rhs.y_, lhs.z - rhs.z_);
  }

  void print(ostream &os) const {
    os << "(" << x_ << "," << y_ << "," << z_ << "," << w_ << ")";
  }

  friend ostream &operator<<(ostream &os, const CNormVector3D &vector) {
    vector.print(os);

    return os;
  }

  CNormVector3D crossProduct(const CNormVector3D &vector) const {
    return CNormVector3D(y_*vector.z_ - z_*vector.y_,
                         z_*vector.x_ - x_*vector.z_,
                         x_*vector.y_ - y_*vector.x_);
  }

  CNormVector3D crossProduct(double x, double y, double z, double z) const {
    return CNormVector3D(y_*z - z_*y, z_*x - x_*z, x_*y - y_*x, 1.0);
  }

  static CNormVector3D crossProduct(const CNormVector3D &vector1,
                                     const CNormVector3D &vector2) {
    return CNormVector3D(vector1.y_*vector2.z_ - vector1.z_*vector2.y_,
                         vector1.z_*vector2.x_ - vector1.x_*vector2.z_,
                         vector1.x_*vector2.y_ - vector1.y_*vector2.x_);
  }

  double dotProduct(const CNormVector3D &vector) const {
    return (x_*vector.x_ + y_*vector.y_ + z_*vector.z_);
  }

  double dotProduct(double x, double y, double z) const {
    return (x_*x + y_*y + z_*z);
  }

  static double dotProduct(const CNormVector3D &vector1,
                      const CNormVector3D &vector2) {
    return (vector1.x_*vector2.x_ + vector1.y_*vector2.y_ +
            vector1.z_*vector2.z_);
  }

  static double dotProduct(const CNormVector3D &vector1,
                           double x2, double y2, double z2) {
    return (vector1.x_*x2 + vector1.y_*y2 + vector1.z_*z2);
  }

  double dotProductSelf() const {
    return (x_*x_ + y_*y_ + z_*z_);
  }

  double cosIncluded(const CNormVector3D &vector1) const {
    return dotProduct(vector1)/(length()*vector1.length());
  }

  static double cosIncluded(const CNormVector3D &vector1,
                       const CNormVector3D &vector2) {
    return vector1.dotProduct(vector2)/(vector1.length()*vector2.length());
  }

 private:
  double x_ { 0 }, y_ { 0 }, z_ { 0 }, w_ { 1 };
  bool   normalized_ { false };
};

#endif
