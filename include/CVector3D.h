#ifndef CVECTOR_3D_H
#define CVECTOR_3D_H

// TODO: use CMath::sqrt() and CMath::isqrt() to allow replacement

#include <CMathGen.h>
#include <CPoint3D.h>

template<typename T>
class CPoint3DT;

template<typename T>
class CVector3DT {
 private:
  typedef CPoint3DT<T>  Point;
  typedef CVector3DT<T> Vector;

 private:
  T    x_, y_, z_;
  bool normalized_;

 public:
  enum Type {
    ZERO,
    UNIT,
    UNIT_X,
    UNIT_Y,
    UNIT_Z
  };

 public:
  // constructor/destructor
  CVector3DT() :
   x_(0), y_(0), z_(0), normalized_(false) {
  }

  CVector3DT(T x, T y, T z) :
   x_(x), y_(y), z_(z), normalized_(false) {
  }

 ~CVector3DT() { }

  //------

  // copy operations
  CVector3DT(const Vector &vector) :
   x_(vector.x_), y_(vector.y_), z_(vector.z_), normalized_(vector.normalized_) {
  }

  Vector &operator=(const Vector &v) {
    x_ = v.x_; y_ = v.y_; z_ = v.z_;

    normalized_ = v.normalized_;

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << "," << z_ << ")";
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
  T z() const { return z_; }

  T getX() const { return x_; }
  T getY() const { return y_; }
  T getZ() const { return z_; }

  void getXYZ(T *x, T *y, T *z) const {
    *x = x_; *y = y_; *z = z_;
  }

  const T *getValues() const { return &x_; }

  // Reference routine would break encapsulation

  T operator[](uint i) const {
    assert(i < 3);

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

  void setZ(T z) {
    z_ = z;

    normalized_ = false;
  }

  void setXYZ(T x, T y, T z) {
    x_ = x; y_ = y; z_ = z;

    normalized_ = false;
  }

  void iset(uint i, T v) {
    assert(i < 3);

    (&x_)[i] = v;

    normalized_ = false;
  }

  // more get accessors
  bool getNormalized() const {
    return normalized_;
  }

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

    return (x_*x_ + y_*y_ + z_*z_);
  }

  //------

  // comparison
  // TODO: tolerance ? use eq()
  int cmp(const Vector &v) const {
    if      (x_ < v.x_) return -1;
    else if (x_ > v.x_) return  1;
    else if (y_ < v.y_) return -1;
    else if (y_ > v.y_) return  1;
    else if (z_ < v.z_) return -1;
    else if (z_ > v.z_) return  1;
    else                return  0;
  }

  friend bool operator==(const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const Vector &lhs, const Vector &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  bool eq(const Vector &rhs) const {
    return fabs(x_ - rhs.x_) < 1E-6 &&
           fabs(y_ - rhs.y_) < 1E-6 &&
           fabs(z_ - rhs.z_) < 1E-6;
  }

  //------

  explicit CVector3DT(const Point &point) :
    x_(point.x), y_(point.y), z_(point.z), normalized_(false) {
  }

  // v = point2 - point1
  CVector3DT(const Point &point1, const Point &point2) :
   x_(point2.x - point1.x), y_(point2.y - point1.y),
   z_(point2.z - point1.z), normalized_(false) {
  }

  // v = vector2 - vector1
  CVector3DT(const Vector &vector1, const Vector &vector2) :
   x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_),
   z_(vector2.z_ - vector1.z_), normalized_(false) {
  }

  //------

  // parametric creation
  CVector3DT(Type type, T s=1.0) :
   normalized_(false) {
    if      (type == ZERO  ) { x_ = 0.0; y_ = 0.0; z_ = 0.0; }
    else if (type == UNIT  ) { x_ =   s; y_ =   s; z_ =   s; }
    else if (type == UNIT_X) { x_ =   s; y_ = 0.0; z_ = 0.0; }
    else if (type == UNIT_Y) { x_ = 0.0; y_ =   s; z_ = 0.0; }
    else if (type == UNIT_Z) { x_ = 0.0; y_ = 0.0; z_ =   s; }
  }

  //------

  // to point
  Point point() const {
    return Point(x_, y_, z_);
  }

  //------

  Vector &zero() {
    x_ = 0.0; y_ = 0.0; z_ = 0.0;

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

    // assert on len == 0 ?

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    x_ *= factor;
    y_ *= factor;
    z_ *= factor;

    normalized_ = true;

    return *this;
  }

  Vector normalized() const {
    if (normalized_)
      return *this;

    T len = length();

    // assert on len == 0 ?

    T factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    return Vector(x_*factor, y_*factor, z_*factor, true);
  }

  Vector unit() const {
    return normalized();
  }

  //------

  Vector &setMagnitude(T magnitude) {
    T factor = 0.0;

    if (normalized_)
      factor = magnitude;
    else {
      T len = length();

      if (len > 0.0)
        factor = magnitude/len;
    }

    x_ *= factor; y_ *= factor; z_ *= factor;

    normalized_ = false;

    return *this;
  }

  //------

  T getDistance(const Vector &vector) const {
    Vector diff = *this - vector;

    return diff.length();
  }

  T getDistanceSqr(const Vector &vector) const {
    Vector diff = *this - vector;

    return diff.lengthSqr();
  }

  //------

  void incX(T x = 1.0) { x_ += x; normalized_ = false; }
  void incY(T y = 1.0) { y_ += y; normalized_ = false; }
  void incZ(T z = 1.0) { z_ += z; normalized_ = false; }

  void decX(T x = 1.0) { x_ -= x; normalized_ = false; }
  void decY(T y = 1.0) { y_ -= y; normalized_ = false; }
  void decZ(T z = 1.0) { z_ -= z; normalized_ = false; }

  void scaleX(T x = 1.0) { x_ *= x; normalized_ = false; }
  void scaleY(T y = 1.0) { y_ *= y; normalized_ = false; }
  void scaleZ(T z = 1.0) { z_ *= z; normalized_ = false; }

  //------

  T minComponent() { return std::min(std::min(x_, y_), z_); }
  T maxComponent() { return std::max(std::max(x_, y_), z_); }

  T minAbsComponent() { return std::min(std::min(::fabs(x_), ::fabs(y_)), ::fabs(z_)); }
  T maxAbsComponent() { return std::max(std::max(::fabs(x_), ::fabs(y_)), ::fabs(z_)); }

  //------

  static Vector min(const Vector &lhs, const Vector &rhs) {
    return Vector(min(lhs.x_, rhs.x_), min(lhs.y_, rhs.y_), min(lhs.z_, rhs.z_));
  }

  static Vector max(const Vector &lhs, const Vector &rhs) {
    return Vector(max(lhs.x_, rhs.x_), max(lhs.y_, rhs.y_), max(lhs.z_, rhs.z_));
  }

  //------

  Vector &operator=(const Point &point) {
    x_ = point.x; y_ = point.y; z_ = point.z;

    normalized_ = false;

    return *this;
  }

  //------

  // operators

  // unary +/-
  Vector operator+() const {
    return Vector(x_, y_, z_, normalized_);
  }

  Vector operator-() const {
    return Vector(-x_, -y_, -z_, normalized_);
  }

  // addition
  Vector &operator+=(const Vector &rhs) {
    x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_;

    normalized_ = false;

    return *this;
  }

  Vector operator+(const Vector &rhs) const {
    return Vector(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
  }

  // subtraction
  Vector &operator-=(const Vector &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_;

    normalized_ = false;

    return *this;
  }

  Vector operator-(const Vector &rhs) const {
    return Vector(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
  }

  // scalar multiplication/division
  Vector &operator*=(T rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs;

    normalized_ = false;

    return *this;
  }

  Vector operator*(const Vector &rhs) {
    Vector t(*this);

    t *= rhs;

    return t;
  }

  friend Vector operator*(const Vector &lhs, T rhs) {
    return Vector(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs);
  }

  friend Vector operator*(T lhs, const Vector &rhs) {
    return Vector(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_);
  }

  Vector &operator/=(T rhs) {
    T irhs = 1.0/rhs;

    x_ *= irhs; y_ *= irhs; z_ *= irhs;

    normalized_ = false;

    return *this;
  }

  Vector operator/(T rhs) const {
    T irhs = 1.0/rhs;

    return Vector(x_*irhs, y_*irhs, z_*irhs);
  }

  //------

  // dot product
  T dotProduct(const Vector &v) const {
    return (x_*v.x_ + y_*v.y_ + z_*v.z_);
  }

//T dotProduct(const Normal &normal) const {
//  return (x_*normal.getX() + y_*normal.getY() + z_*normal.getZ());
//}

  static T dotProduct(const Vector &v1, const Vector &v2) {
    return v1.dotProduct(v2);
  }

  static T dotProduct(const Point &v1, const Vector &v2) {
    return Vector(v1.x, v1.y, v1.z).dotProduct(Vector(v2));
  }

  T dotProduct(const Point &point) const {
    return dotProduct(Vector(point.x, point.y, point.z));
  }

  T dotProduct(T x, T y, T z) const {
    return dotProduct(Vector(x, y, z));
  }

  static T dotProduct(const Vector &v1, T x2, T y2, T z2) {
    return v1.dotProduct(Vector(x2, y2, z2));
  }

  T dotProductSelf() const {
    return dotProduct(*this);
  }

  static T absDotProduct(const Vector &v1, const Vector &v2) {
    return ::fabs(dotProduct(v1, v2));
  }

  //------

  // cross product
  Vector crossProduct(const Vector &v) const {
    return Vector(y_*v.z_ - z_*v.y_, z_*v.x_ - x_*v.z_, x_*v.y_ - y_*v.x_);
  }

  static Vector crossProduct(const Vector &v1, const Vector &v2) {
    return v1.crossProduct(v2);
  }

  Vector crossProduct(const Point &point) const {
    return crossProduct(Vector(point.x, point.y, point.z));
  }

  Vector crossProduct(T x, T y, T z) const {
    return crossProduct(Vector(x, y, z));
  }

  Vector unitCrossProduct(const Vector &v) const {
    return crossProduct(v).normalize();
  }

  // Note: area of parallelogram is
  // CVector v1(x2 - x1, y2 - y1), v2(x4 - x1, y4 - y1);
  // area = v1.crossProduct(v2).length();

  //------

  friend Point operator+(const Point &lhs, const Vector &rhs) {
    return Point(lhs.x + rhs.x_, lhs.y + rhs.y_, lhs.z + rhs.z_);
  }

  friend Point operator+=(Point &lhs, const Vector &rhs) {
    lhs.x += rhs.x_; lhs.y += rhs.y_; lhs.z += rhs.z_;

    return lhs;
  }

  friend Point operator-(const Point &lhs, const Vector &rhs) {
    return Point(lhs.x - rhs.x_, lhs.y - rhs.y_, lhs.z - rhs.z_);
  }

  friend Point operator-=(Point &lhs, const Vector &rhs) {
    lhs.x -= rhs.x_; lhs.y -= rhs.y_; lhs.z -= rhs.z_;

    return lhs;
  }

#if 0
  friend Vector operator-(const Point &lhs, const Point &rhs) {
    return Vector(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }
#endif

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
    vector1.cosIncluded(vector2);
  }

  //------

  void directionCosines(T *x, T *y, T *z) {
    Vector vector1 = normalized();

    *x = vector1.x;
    *y = vector1.y;
    *z = vector1.z;
   }

  //------

  static void coordinateSystem(const Vector &vector1, Vector *vector2, Vector *vector3) {
    if (! vector1.normalized_) {
      coordinateSystem(vector1.normalized(), vector2, vector3);
      return;
    }

    if (::fabs(vector1.x_) > ::fabs(vector1.y_)) {
      T ilen = 1.0/::sqrt(vector1.x_*vector1.x_ + vector1.z_*vector1.z_);

      *vector2 = Vector(-vector1.z_*ilen, 0.0, vector1.x_*ilen);
    }
    else {
      T ilen = 1.0/::sqrt(vector1.y_*vector1.y_ + vector1.z_*vector1.z_);

      *vector2 = Vector(0.0, vector1.z_*ilen, -vector1.y_*ilen);
    }

    *vector3 = vector1.crossProduct(*vector2);
  }

  //------

  Vector perp() const {
    return perp(*this);
  }

  static Vector perp(const Vector &u) {
    Vector u1 = u.unit();
    Vector v1 = Vector(0,1,0);

    T vu = v1.dotProduct(u1);

    v1 = (v1 - vu*u1).unit();

    if (v1.isZero()) {
      v1 = Vector(0,0,1);

      vu = v1.dotProduct(u1);

      v1 = (v1 - vu*u1).unit();
    }

    Vector w1 = u1.crossProduct(v1);

    v1 = w1.crossProduct(u1);

    return v1;
  }

  //------

#if 0
  void getBarycentrics(const Vector &vector1, const Vector &vector2,
                       const Vector &vector3, const Vector &vector4,
                       T barycentric[4]) const {
    // compute the vectors relative to V3 of the tetrahedron
    Vector diff[4] = { vector1 - vector4, vector2 - vector4, vector3 - vector4, *this - vector4 };

    // If the vertices have large magnitude, the linear system of equations
    // for computing barycentric coordinates can be ill-conditioned.  To avoid
    // this, uniformly scale the tetrahedron edges to be of order 1.  The
    // scaling of all differences does not change the barycentric coordinates.
    T maxval = 0.0;

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        maxval = max(maxval, ::fabs(diff[i][j]));

    // scale down only large data
    if (maxval > 1.0) {
      T imaxval = 1.0/maxval;

      for (int i = 0; i < 4; i++)
        diff[i] *= imaxval;
    }

    T det = diff[0].dotProduct(diff[1].crossProduct(diff[2]));

    Vector c12 = diff[1].crossProduct(diff[2]);
    Vector c20 = diff[2].crossProduct(diff[0]);
    Vector c01 = diff[0].crossProduct(diff[1]);

    if (::fabs(det) > 1E-6 ) {
      T idet = 1.0/det;

      barycentric[0] = diff[3].dotProduct(c12)*idet;
      barycentric[1] = diff[3].dotProduct(c20)*idet;
      barycentric[2] = diff[3].dotProduct(c01)*idet;
      barycentric[3] = 1.0 -
                       barycentric[0] - barycentric[1] - barycentric[2];
    }
    else {
      // The tetrahedron is potentially flat.  Determine the face of
      // maximum area and compute barycentric coordinates with respect
      // to that face
      Vector d13     = vector1 - vector3;
      Vector d23     = vector2 - vector3;
      Vector d13cd23 = d13.crossProduct(d23);

      T max_length = d13cd23.lengthSqr();

      int max_ind = 3;

      T length2 = c01.lengthSqr();

      if (length2 > max_length) {
        max_ind    = 0;
        max_length = length2;
      }

      length2 = c12.lengthSqr();

      if (length2 > max_length) {
        max_ind    = 1;
        max_length = length2;
      }

      length2 = c20.lengthSqr();

      if (length2 > max_length) {
        max_ind    = 2;
        max_length = length2;
      }

      if (max_length > 1E-6) {
        T imax_length = 1.0/max_length;

        Vector tmp;

        if      (max_ind == 0) {
          tmp = diff[3].crossProduct(diff[1]);

          barycentric[0] = c01.dotProduct(tmp)*imax_length;

          tmp = diff[0].crossProduct(diff[3]);

          barycentric[1] = c01.dotProduct(tmp)*imax_length;
          barycentric[2] = 0.0;
          barycentric[3] = 1.0 - barycentric[0] - barycentric[1];
        }
        else if (max_ind == 1) {
          barycentric[0] = 0.0;

          tmp = diff[3].crossProduct(diff[2]);

          barycentric[1] = c12.dotProduct(tmp)*imax_length;

          tmp = diff[1].crossProduct(diff[3]);

          barycentric[2] = c12.dotProduct(tmp)*imax_length;
          barycentric[3] = 1.0 - barycentric[1] - barycentric[2];
        }
        else if (max_ind == 2) {
          tmp = diff[2].crossProduct(diff[3]);

          barycentric[0] = c20.dotProduct(tmp)*imax_length;
          barycentric[1] = 0.0;

          tmp = diff[3].crossProduct(diff[0]);

          barycentric[2] = c20.dotProduct(tmp)*imax_length;
          barycentric[3] = 1.0 - barycentric[0] - barycentric[2];
        }
        else {
          diff[3] = *this - rkV2;

          tmp = diff[3].crossProduct(d23);

          barycentric[0] = d13cd23.dotProduct(tmp)*imax_length;

          tmp = d13.crossProduct(diff[3]);

          barycentric[1] = d13cd23.dotProduct(tmp)*imax_length;
          barycentric[2] = 1.0 - barycentric[0] - barycentric[1];
          barycentric[3] = 0.0;
        }
      }
      else {
        // The tetrahedron is potentially a sliver.  Determine the edge of
        // maximum length and compute barycentric coordinates with respect
        // to that edge.
        T max_length2 = diff[0].lengthSqr();

        max_ind = 0;  // <V0,V3>

        T length2 = diff[1].lengthSqr();

        if (length2 > max_length2) {
          max_ind     = 1;  // <V1,V3>
          max_length2 = length2;
        }

        length2 = diff[2].lengthSqr();

        if (length2 > max_length2) {
          max_ind     = 2;  // <V2,V3>
          max_length2 = length2;
        }

        length2 = d13.lengthSqr();

        if (length2 > max_length2) {
          max_ind     = 3;  // <V0,V2>
          max_length2 = length2;
        }

        length2 = d23.lengthSqr();

        if (length2 > max_length2) {
          max_ind     = 4;  // <V1,V2>
          max_length2 = length2;
        }

        Vector v12 = vector1 - vector2;

        length2 = v12.lengthSqr();

        if (length2 > max_length2) {
          max_ind     = 5;  // <V0,V1>
          max_length2 = length2;
        }

        if ( max_length2 > 1E-6 ) {
          T imax_length2 = 1.0/max_length2;

          if ( max_ind == 0 ) {
            // P-V3 = t*(V0-V3)
            barycentric[0] = diff[3].dotProduct(diff[0])*imax_length2;
            barycentric[1] = 0.0;
            barycentric[2] = 0.0;
            barycentric[3] = 1.0 - barycentric[0];
          }
          else if ( max_ind == 1 ) {
            // P-V3 = t*(V1-V3)
            barycentric[0] = 0.0;
            barycentric[1] = diff[3].dotProduct(diff[1])*imax_length2;
            barycentric[2] = 0.0;
            barycentric[3] = 1.0 - barycentric[1];
          }
          else if ( max_ind == 2 ) {
            // P-V3 = t*(V2-V3)
            barycentric[0] = 0.0;
            barycentric[1] = 0.0;
            barycentric[2] = diff[3].dotProduct(diff[2])*imax_length2;
            barycentric[3] = 1.0 - barycentric[2];
          }
          else if ( max_ind == 3 ) {
            // P-V2 = t*(V0-V2)
            diff[3] = *this - vector2;

            barycentric[0] = diff[3].dotProduct(d13)*imax_length2;
            barycentric[1] = 0.0;
            barycentric[2] = 1.0 - barycentric[0];
            barycentric[3] = 0.0;
          }
          else if ( max_ind == 4 ) {
            // P-V2 = t*(V1-V2)
            diff[3] = *this - vector2;

            barycentric[0] = 0.0;
            barycentric[1] = diff[3].dotProduct(d23)*imax_length2;
            barycentric[2] = 1.0 - barycentric[1];
            barycentric[3] = 0.0;
          }
          else {
            // P-V1 = t*(V0-V1)
            diff[3] = *this - vector2;

            barycentric[0] = diff[3].dotProduct(v12)*imax_length2;
            barycentric[1] = 1.0 - barycentric[0];
            barycentric[2] = 0.0;
            barycentric[3] = 0.0;
          }
        }
        else {
            // tetrahedron is a nearly a point, just return equal weights
            barycentric[0] = 0.25;
            barycentric[1] = barycentric[0];
            barycentric[2] = barycentric[0];
            barycentric[3] = barycentric[0];
        }
      }
    }
  }
#endif

  //------

 private:
  CVector3DT(T x, T y, T z, bool normalized) :
   x_(x), y_(y), z_(z), normalized_(normalized) {
  }
};

typedef CVector3DT<double> CVector3D;
typedef CVector3DT<float>  CVector3DF;

#endif
