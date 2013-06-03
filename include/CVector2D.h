#ifndef CVECTOR_2D_H
#define CVECTOR_2D_H

#include <CPoint2D.h>
#include <CMathGen.h>

#include <cassert>
#include <cmath>

template<typename T>
class CPoint2DT;

template<typename T>
class CNormal2DT;

template<typename T>
class CVector2DT {
 private:
  typedef CVector2DT<T> Vector;
  typedef CPoint2DT<T>  Point;
  typedef CNormal2DT<T> Normal;

 private:
  T    x_, y_;
  bool normalized_;

 public:
  enum Type {
    ZERO,
    UNIT,
    UNIT_X,
    UNIT_Y
  };

 public:
  // constructor/destructor
  CVector2DT() :
   x_(0), y_(0), normalized_(false) {
  }

  CVector2DT(T x, T y) :
   x_(x), y_(y), normalized_(false) {
  }

 ~CVector2DT() { }

  //------

  // copy operations
  CVector2DT(const Vector &vector) :
   x_(vector.x_), y_(vector.y_),
   normalized_(vector.normalized_) {
  }

  Vector &operator=(const Vector &vector) {
    x_ = vector.x_; y_ = vector.y_;

    normalized_ = false;

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const Vector &vector) {
    vector.print(os);

    return os;
  }

  //------

  // accessors

  // get
  T x() const { return x_; }
  T y() const { return y_; }

  T getX() const { return x_; }
  T getY() const { return y_; }

  void getXY(T *x, T *y) const {
    *x = x_; *y = y_;
  }

  const T *getValues() const { return &x_; }

  // Reference routine would break encapsulation

  T operator[](uint i) const {
    assert(i < 2);

    return (&x_)[i];
  }

  // set
  void setX(T x) {
    x_ = x;

    normalized_ = false;
  }

  void setY(T y) {
    y_ = y;

    normalized_ = false;
  }

  void setXY(T x, T y) {
    x_ = x; y_ = y;

    normalized_ = false;
  }

  void iset(uint i, T v) {
    assert(i < 2);

    (&x_)[i] = v;

    normalized_ = false;
  }

  // more get accessors
  bool getNormalized() const {
    return normalized_;
  }

  //------

  T length() const {
    if (normalized_)
      return 1.0;

    return ::sqrt(lengthSqr());
  }

  T modulus() const {
    return length();
  }

  T lengthSqr() const {
    if (normalized_)
      return 1.0;

    return (x_*x_ + y_*y_);
  }

  //------

  // angle

  double angle() const {
    return atan2(y_, x_);
  }

  //------

  // comparison
  int cmp(const Vector &v) const {
    if      (x_ < v.x_) return -1;
    else if (x_ > v.x_) return 1;
    else if (y_ < v.y_) return -1;
    else if (y_ > v.y_) return 1;
    else                return 0;
  }

  friend bool operator==(const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) < 0;
  }

  friend bool operator<=(const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) > 0;
  }

  friend bool operator>=(const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  bool eq(const Vector &rhs) const {
    return fabs(x_ - rhs.x_) < 1E-6 &&
           fabs(y_ - rhs.y_) < 1E-6;
  }

  //------

  explicit CVector2DT(const Point &point) :
   x_(point.x), y_(point.y), normalized_(false) {
  }

  CVector2DT(const Point &point1, const Point &point2) :
   x_(point2.x - point1.x), y_(point2.y - point1.y),
   normalized_(false) {
  }

  CVector2DT(const Vector &vector1, const Vector &vector2) :
   x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_),
   normalized_(false) {
  }

  //------

  CVector2DT(Type type) {
    if      (type == ZERO  ) { x_ = 0.0; y_ = 0.0; }
    else if (type == UNIT  ) { x_ = 1.0; y_ = 1.0; }
    else if (type == UNIT_X) { x_ = 1.0; y_ = 0.0; }
    else if (type == UNIT_Y) { x_ = 0.0; y_ = 1.0; }
  }

  //------

  Point point() const {
    return Point(x_, y_);
  }

  //------

  Vector &zero() {
    x_ = 0.0; y_ = 0.0;

    normalized_ = false;

    return *this;
  }

  bool isZero() const {
    return REAL_EQ(lengthSqr(), 0);
  }

  bool isUnit() const {
    return REAL_EQ(lengthSqr(), 1);
  }

  //------

  Vector &normalize() {
    if (normalized_)
      return *this;

    T len = length();

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    x_ *= factor;
    y_ *= factor;

    normalized_ = true;

    return *this;
  }

  Vector normalized() const {
    if (normalized_)
      return *this;

    T len = length();

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    return Vector(x_*factor, y_*factor, true);
  }

  Vector unit() const {
    return normalized();
  }

  //------

  Vector perpendicular() const {
    return Vector(y_, -x_, normalized_);
  }

  Vector unitPerpendicular() const {
    return perpendicular().normalize();
  }

  T dotPerpendicular(const Vector &v) const {
    return x_*v.y_ - y_*v.x_;
  }

  //------

  Vector & setMagnitude(T magnitude) {
    T factor = 0.0;

    if (normalized_)
      factor = magnitude;
    else {
      T len = length();

      if (len > 0.0)
        factor = magnitude/len;
    }

    x_ *= factor; y_ *= factor;

    normalized_ = false;

    return *this;
  }

  //------

  T getDistance(const Vector &vector) {
    Vector diff = *this - vector;

    return diff.length();
  }

  T getDistanceSqr(const Vector &vector) {
    Vector diff = *this - vector;

    return diff.lengthSqr();
  }

  //------

  void incX(T x = 1.0) {
    x_ += x;

    normalized_ = false;
  }

  void incY(T y = 1.0) {
    y_ += y;

    normalized_ = false;
  }

  void decX(T x = 1.0) {
    x_ -= x;

    normalized_ = false;
  }

  void decY(T y = 1.0) {
    y_ -= y;

    normalized_ = false;
  }

  //------

  T minComponent() {
    return std::min(x_, y_);
  }

  T maxComponent() {
    return std::max(x_, y_);
  }

  T minAbsComponent() {
    return std::min(::fabs(x_), ::fabs(y_));
  }

  T maxAbsComponent() {
    return std::max(::fabs(x_), ::fabs(y_));
  }

  //------

  static Vector min(const Vector &lhs, const Vector &rhs) {
    return Vector(min(lhs.x_, rhs.x_), min(lhs.y_, rhs.y_));
  }

  static Vector max(const Vector &lhs, const Vector &rhs) {
    return Vector(max(lhs.x_, rhs.x_), max(lhs.y_, rhs.y_));
  }

  //------

  Vector &operator=(const Point &point) {
    x_ = point.x; y_ = point.y;

    normalized_ = false;

    return *this;
  }

  //------

  // operators

  // unary +/-
  Vector operator+() const {
    return Vector(x_, y_, normalized_);
  }

  Vector operator-() const {
    return Vector(-x_, -y_, normalized_);
  }

  // addition
  Vector &operator+=(const Vector &rhs) {
    x_ += rhs.x_; y_ += rhs.y_;

    normalized_ = false;

    return *this;
  }

  Vector operator+(const Vector &rhs) const {
    return Vector(x_ + rhs.x_, y_ + rhs.y_);
  }

  // subtraction
  Vector &operator-=(const Vector &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_;

    normalized_ = false;

    return *this;
  }

  Vector operator-(const Vector &rhs) const {
    return Vector(x_ - rhs.x_, y_ - rhs.y_);
  }

  // scalar multiplication/division
  Vector &operator*=(T rhs) {
    x_ *= rhs; y_ *= rhs;

    normalized_ = false;

    return *this;
  }

  Vector operator*(const Vector &rhs) {
    Vector t(*this);

    t *= rhs;

    return t;
  }

  friend Vector operator*(const Vector &lhs, T rhs) {
    return Vector(lhs.x_*rhs, lhs.y_*rhs);
  }

  friend Vector operator*(T lhs, const Vector &rhs) {
    return Vector(lhs*rhs.x_, lhs*rhs.y_);
  }

  Vector &operator/=(T rhs) {
    T irhs = 1.0/rhs;

    x_ *= irhs; y_ *= irhs;

    normalized_ = false;

    return *this;
  }

  Vector operator/(T rhs) const {
    T irhs = 1.0/rhs;

    return Vector(x_*irhs, y_*irhs);
  }

  //------

  // dot product
  T dotProduct(const Vector &vector) const {
    return (x_*vector.x_ + y_*vector.y_);
  }

  T dotProduct(const Point &point) const {
    return dotProduct(Vector(point.x, point.y));
  }

  T dotProduct(T x, T y) const {
    return dotProduct(Vector(x, y));
  }

  static T dotProduct(const Vector &vector1, const Vector &vector2) {
    return vector1.dotProduct(vector2);
  }

  static T dotProduct(const Vector &vector1, T x2, T y2) {
    return vector1.dotProduct(Vector(x2, y2));
  }

  T dotProductSelf() const {
    return dotProduct(*this);
  }

  static T absDotProduct(const Vector &vector1, const Vector &vector2) {
    return ::fabs(dotProduct(vector1, vector2));
  }

  //------

  // cross product
  Vector crossProduct(const Vector &vector) const {
    return Vector(y_*vector.x_ - x_*vector.y_,
                  x_*vector.y_ - y_*vector.x_,
                  normalized_ && vector.normalized_);
  }

  Vector crossProduct(const Point &point) const {
    return crossProduct(Vector(point.x, point.y));
  }

  Vector crossProduct(T x, T y) const {
    return crossProduct(Vector(x, y));
  }

  static Vector crossProduct(const Vector &vector1, const Vector &vector2) {
    return vector1.crossProduct(vector2);
  }

  Vector unitCrossProduct(const Vector &vector) const {
    return crossProduct(vector).normalize();
  }

  // Note: area of parallelogram is
  // CVector v1(x2 - x1, y2 - y1), v2(x4 - x1, y4 - y1);
  // area = v1.crossProduct(v2).length();

  //------

  friend Point operator+(const Point &lhs, const Vector &rhs) {
    return Point(lhs.x + rhs.x_, lhs.y + rhs.y_);
  }

  friend Point operator+=(Point &lhs, const Vector &rhs) {
    lhs.x += rhs.x_; lhs.y += rhs.y_;

    return lhs;
  }

  friend Point operator-(const Point &lhs, const Vector &rhs) {
    return Point(lhs.x - rhs.x_, lhs.y - rhs.y_);
  }

  friend Point operator-=(Point &lhs, const Vector &rhs) {
    lhs.x -= rhs.x_; lhs.y -= rhs.y_;

    return lhs;
  }

/*
  friend Vector operator-(const Point &lhs, const Point &rhs) {
    return Vector(lhs.x - rhs.x, lhs.y - rhs.y);
  }
*/

  //------

  Vector normal(const Vector &vector2) {
    return unitCrossProduct(vector2);
  }

  static Vector normal(const Vector &vector1, const Vector &vector2) {
    return vector1.normal(vector2);
  }

  //------

  T cosIncluded(const Vector &vector1) const {
    T dot = dotProduct(vector1);

    if (! normalized_)
      dot /= length();

    if (! vector1.normalized_)
      dot /= vector1.length();

    return dot;
  }

  static T cosIncluded(const Vector &vector1, const Vector &vector2) {
    return vector1.cosIncluded(vector2);
  }

  //------

  void getBarycentrics(const Vector &vector1, const Vector &vector2,
                       const Vector &vector3, T barycentrics[3]) const {
    // compute the vectors relative to V2 of the triangle
    Vector diff[3] = {
      vector1 - vector3,
      vector2 - vector3,
      *this   - vector3
    };

    // If the vertices have large magnitude, the linear system of equations
    // for computing barycentric coordinates can be ill-conditioned.  To avoid
    // this, uniformly scale the triangle edges to be of order 1.  The scaling
    // of all differences does not change the barycentric coordinates.
    T maxval = 0.0;

    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        maxval = max(maxval, ::fabs(diff[i][j]));

    // scale down only large data
    if (maxval > 1.0) {
      T imaxval = 1.0/maxval;

      for (int i = 0; i < 3; i++)
        diff[i] *= imaxval;
    }

    T det = diff[0].dotPerp(diff[1]);
    if (::fabs(det) > 1E-6 ) {
      T idet = 1.0/det;

      barycentrics[0] = diff[2].DotPerp(diff[1])*idet;
      barycentrics[1] = diff[0].DotPerp(diff[2])*idet;
      barycentrics[2] = 1.0 - barycentrics[0] - barycentrics[1];
    }
    else {
      // The triangle is a sliver. Determine the longest edge and
      // compute barycentric coordinates with respect to that edge.
      Vector v12 = vector1 - vector2;

      T max_length = v12.lengthSqr();

      int max_ind = 2;

      T fSqrLength = diff[1].lengthSqr();

      if (fSqrLength > max_length ) {
        max_ind    = 1;
        max_length = fSqrLength;
      }

      fSqrLength = diff[0].lengthSqr();

      if (fSqrLength > max_length) {
        max_ind    = 2;
        max_length = fSqrLength;
      }

      if (max_length > 1E-6) {
        T imax_length = 1.0/max_length;

        if      (max_ind == 0) {
          // P-V2 = t(V0-V2)
          barycentrics[0] = diff[2].Dot(diff[0])*imax_length;
          barycentrics[1] = 0.0;
          barycentrics[2] = 1.0 - barycentrics[0];
        }
        else if (max_ind == 1) {
          // P-V2 = t(V1-V2)
          barycentrics[0] = 0.0;
          barycentrics[1] = diff[2].Dot(diff[1])*imax_length;
          barycentrics[2] = 1.0 - barycentrics[1];
        }
        else {
          // P-V1 = t(V0-V1)
          diff[2] = *this - vector2;
          barycentrics[0] = diff[2].Dot(v12)*imax_length;
          barycentrics[1] = 1.0 - barycentrics[0];
          barycentrics[2] = 0.0;
        }
      }
      else {
        // triangle is a nearly a point, just return equal weights
        barycentrics[0] = 0.333333333;
        barycentrics[1] = barycentrics[0];
        barycentrics[2] = barycentrics[0];
      }
    }
  }

  //------

  void orthonormalize(Vector &u, Vector &v) {
    // If the input vectors are v0 and v1, then the Gram-Schmidt
    // orthonormalization produces vectors u0 and u1 as follows,
    //
    //   u0 = v0/|v0|
    //   u1 = (v1-(u0*v1)u0)/|v1-(u0*v1)u0|
    //
    // where |A| indicates length of vector A and A*B indicates dot
    // product of vectors A and B.

    // compute u0
    u.normalize();

    // compute u1
    T d = u.dotProduc(v);

    v -= u*d;

    v.normalize();
  }

  //------

  void generateOrthonormalBasis(Vector &u, Vector &v) {
    v.normalize();

    u = v.perpendicular();
  }

  //------

 private:
  CVector2DT(T x, T y, bool normalized) :
   x_(x), y_(y), normalized_(normalized) {
  }
};

typedef CVector2DT<double> CVector2D;
typedef CVector2DT<float>  CVector2DF;

#endif
