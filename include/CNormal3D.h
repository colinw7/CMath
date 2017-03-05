#ifndef CNORMAL_3D_H
#define CNORMAL_3D_H

#include <CMathGen.h>
#include <CVector3D.h>

class CNormal3D {
 public:
  explicit CNormal3D(double x=0.0, double y=0.0, double z=0.0) :
   x_(x), y_(y), z_(z), normalized_(false) {
  }

  CNormal3D(const CNormal3D &normal) :
    x_(normal.x_), y_(normal.y_), z_(normal.z_),
    normalized_(normal.normalized_) {
  }

  explicit CNormal3D(const CVector3D &vector) :
    x_(vector.getX()), y_(vector.getY()), z_(vector.getZ()),
    normalized_(vector.getNormalized()) {
  }

  //------

  double getX() const { return x_; }
  double getY() const { return y_; }
  double getZ() const { return z_; }

  void getXYZ(double *x, double *y, double *z) const {
    *x = x_; *y = y_; *z = z_;
  }

  double operator[](int i) const { assert(i < 3); return (&x_)[i]; }

  // Reference routine would break encapsulation

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

  void setXYZ(double x, double y, double z) {
    x_ = x; y_ = y; z_ = z;

    normalized_ = false;
  }

  bool getNormalized() const {
    return normalized_;
  }

  //------

  double length() const {
    if (normalized_)
      return 1.0 ;
    else
      return ::sqrt(x_*x_ + y_*y_ + z_*z_);
  }

  double fastLength() const {
    if (normalized_)
      return 1.0 ;
    else
      return CMathGen::fastDistance(x_, y_, z_);
  }

  double modulus() const {
    return length();
  }

  double lengthSqr() const {
    if (normalized_)
      return 1.0 ;
    else
      return (x_*x_ + y_*y_ + z_*z_);
  }

  //------

  CNormal3D &zero() {
    x_ = 0.0; y_ = 0.0; z_ = 0.0;

    normalized_ = false;

    return *this;
  }

  CNormal3D unit() const {
    if (normalized_)
      return *this;

    double len = length();

    double factor = 0.0;

    if (len > 0.0)
      factor = 1.0 /len;

    return CNormal3D(x_*factor, y_*factor, z_*factor, true);
  }

  //------

  CNormal3D &normalize() {
    if (normalized_)
      return *this;

    double len = length();

    double factor = 0.0;

    if (len > 0.0)
      factor = 1.0 /len;

    x_ *= factor;
    y_ *= factor;
    z_ *= factor;

    normalized_ = true;

    return *this;
  }

  CNormal3D normalized() const {
    if (normalized_)
      return *this;

    double len = length();

    double factor = 0.0;

    if (len > 0.0)
      factor = 1.0 /len;

    return CNormal3D(x_*factor, y_*factor, z_*factor, true);
  }

  //------

  CNormal3D &setMagnitude(double magnitude) {
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

    normalized_ = false;

    return *this;
  }

  //------

  double getDistance(const CNormal3D &normal) {
    CNormal3D diff = *this - normal;

    return diff.length();
  }

  //------

  void incX(double x = 1.0 ) {
    x_ += x;

    normalized_ = false;
  }

  void incY(double y = 1.0 ) {
    y_ += y;

    normalized_ = false;
  }

  void incZ(double z = 1.0 ) {
    z_ += z;

    normalized_ = false;
  }

  void decX(double x = 1.0 ) {
    x_ -= x;

    normalized_ = false;
  }

  void decY(double y = 1.0 ) {
    y_ -= y;

    normalized_ = false;
  }

  void decZ(double z = 1.0 ) {
    z_ -= z;

    normalized_ = false;
  }

  //------

  double minComponent() {
    return std::min(std::min(x_, y_), z_);
  }

  double maxComponent() {
    return std::max(std::max(x_, y_), z_);
  }

  double minAbsComponent() {
    return std::min(std::min(std::fabs(x_), std::fabs(y_)), std::fabs(z_));
  }

  double maxAbsComponent() {
    return std::max(std::max(std::fabs(x_), std::fabs(y_)), std::fabs(z_));
  }

  //------

  CNormal3D &operator=(const CNormal3D &normal) {
    x_ = normal.x_; y_ = normal.y_; z_ = normal.z_;

    normalized_ = normal.normalized_;

    return *this;
  }

  //------

  CNormal3D operator+() const {
    return CNormal3D(x_, y_, z_);
  }

  CNormal3D operator-() const {
    return CNormal3D(-x_, -y_, -z_);
  }

  //------

  friend bool operator==(const CNormal3D &lhs, const CNormal3D &rhs) {
    return (lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_ && lhs.z_ == rhs.z_);
  }

  friend bool operator!=(const CNormal3D &lhs, const CNormal3D &rhs) {
    return (lhs.x_ != rhs.x_ || lhs.y_ != rhs.y_ || lhs.z_ != rhs.z_);
  }

  //------

  CNormal3D &operator+=(const CNormal3D &rhs) {
    x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_;

    normalized_ = false;

    return *this;
  }

  CNormal3D operator+(const CNormal3D &rhs) const {
    return CNormal3D(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
  }

  //------

  CNormal3D &operator-=(const CNormal3D &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_;

    normalized_ = false;

    return *this;
  }

  CNormal3D operator-(const CNormal3D &rhs) const {
    return CNormal3D(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
  }

  friend CVector3D operator-(const CNormal3D &lhs, const CVector3D &rhs) {
    return CVector3D(lhs.x_ - rhs.getX(), lhs.y_ - rhs.getY(), lhs.z_ - rhs.getZ());
  }

  friend CVector3D operator-(const CVector3D &lhs, const CNormal3D &rhs) {
    return CVector3D(lhs.getX() - rhs.x_, lhs.getY() - rhs.y_, lhs.getZ() - rhs.z_);
  }

  //------

  CNormal3D &operator*=(double rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs;

    normalized_ = false;

    return *this;
  }

  friend CNormal3D operator*(const CNormal3D &lhs, double rhs) {
    return CNormal3D(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs);
  }

  friend CNormal3D operator*(double lhs, const CNormal3D &rhs) {
    return CNormal3D(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_);
  }

  //------

  CNormal3D &operator/=(double rhs) {
    double irhs = 1.0 /rhs;

    x_ *= irhs; y_ *= irhs; z_ *= irhs;

    normalized_ = false;

    return *this;
  }

  CNormal3D operator/(double rhs) const {
    double irhs = 1.0 /rhs;

    return CNormal3D(x_*irhs, y_*irhs, z_*irhs);
  }

  //------

  void print(std::ostream &os=std::cout) const {
    os << "(" << x_ << "," << y_ << "," << z_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CNormal3D &normal) {
    normal.print(os);

    return os;
  }

  //------

  double dotProduct(const CNormal3D &normal) const {
    return (x_*normal.x_ + y_*normal.y_ + z_*normal.z_);
  }

  double dotProduct(const CVector3D &vector) const {
    return (x_*vector.getX() + y_*vector.getY() + z_*vector.getZ());
  }

  friend double dotProduct(const CVector3D &vector, const CNormal3D &normal) {
    return (vector.getX()*normal.x_ + vector.getY()*normal.y_ + vector.getZ()*normal.z_);
  }

  double dotProduct(double x, double y, double z) const {
    return (x_*x + y_*y + z_*z);
  }

  static double dotProduct(const CNormal3D &normal1, const CNormal3D &normal2) {
    return (normal1.x_*normal2.x_ +
            normal1.y_*normal2.y_ +
            normal1.z_*normal2.z_);
  }

  static double dotProduct(const CNormal3D &normal1, const CVector3D &vector2) {
    return (normal1.x_*vector2.getX() +
            normal1.y_*vector2.getY() +
            normal1.z_*vector2.getZ());
  }

  static double dotProduct(const CVector3D &vector1, const CNormal3D &normal2) {
    return (vector1.getX()*normal2.x_ +
            vector1.getY()*normal2.y_ +
            vector1.getZ()*normal2.z_);
  }

  static double dotProduct(const CNormal3D &normal1, double x2, double y2, double z2) {
    return (normal1.x_*x2 + normal1.y_*y2 + normal1.z_*z2);
  }

  double dotProductSelf() const {
    return (x_*x_ + y_*y_ + z_*z_);
  }

  static double absDotProduct(const CNormal3D &normal1, const CNormal3D &normal2) {
    return ::fabs(normal1.x_*normal2.x_ + normal1.y_*normal2.y_ + normal1.z_*normal2.z_);
  }

  static double absDotProduct(const CNormal3D &normal1, const CVector3D &vector2) {
    return ::fabs(normal1.x_*vector2.getX() +
                  normal1.y_*vector2.getY() +
                  normal1.z_*vector2.getZ());
  }

  static double absDotProduct(const CVector3D &vector1, const CNormal3D &normal2) {
    return ::fabs(vector1.getX()*normal2.x_ +
                  vector1.getY()*normal2.y_ +
                  vector1.getZ()*normal2.z_);
  }

  //------

  void flip() {
    x_ = -x_; y_ = -y_; z_ = -z_;
  }

  //------

  double cosIncluded(const CNormal3D &normal1) const {
    double dot = dotProduct(normal1);

    if (! normalized_)
      dot /= length();

    if (! normal1.normalized_)
      dot /= normal1.length();

    return dot;
  }

  static double cosIncluded(const CNormal3D &normal1, const CNormal3D &normal2) {
    double dot = normal1.dotProduct(normal2);

    if (! normal1.normalized_)
      dot /= normal1.length();

    if (! normal2.normalized_)
      dot /= normal2.length();

    return dot;
  }

  //------

 private:
  CNormal3D(double x, double y, double z, bool normalized) :
   x_(x), y_(y), z_(z), normalized_(normalized) {
  }

 private:
  double x_ { 0 }, y_ { 0 }, z_ { 0 };
  bool   normalized_ { false };
};

#endif
