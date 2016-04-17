#ifndef CVECTOR_3DH_H
#define CVECTOR_3DH_H

// TODO: use CMath::sqrt() and CMath::isqrt() to allow replacement

#include <CMathGen.h>
#include <CPoint3D.h>

template<typename T>
class CPoint3DT;

template<typename T>
class CVector3DHT {
 private:
  typedef CPoint3DT<T>   Point;
  typedef CVector3DHT<T> Vector;

 private:
  T x_, y_, z_, w_;

  bool normalized_;

 public:
  enum Type {
    ZERO,
    UNIT,
    UNIT_X,
    UNIT_Y,
    UNIT_Z,
    UNIT_W
  };

 public:
  // constructor/destructor
  CVector3DHT() { }

  CVector3DHT(T x, T y, T z, T w=1.0) :
   x_(x), y_(y), z_(z), w_(w), normalized_(false) {
  }

 ~CVector3DHT() { }

  //------

  // copy operations
  CVector3DHT(const Vector &vector) :
    x_(vector.x_), y_(vector.y_), z_(vector.z_), w_(vector.w_),
    normalized_(vector.normalized_) {
  }

  Vector &operator=(const Vector &vector) {
    x_ = vector.x_; y_ = vector.y_; z_ = vector.z_; w_ = vector.w_;

    normalized_ = vector.normalized_;

    return *this;
  }

  //------

  // output
  void print(ostream &os) const {
    os << "(" << x_ << "," << y_ << "," << z_ << "," << w_ << ")";
  }

  friend ostream &operator<<(ostream &os, const Vector &vector) {
    vector.print(os);

    return os;
  }

  //------

  // accessors

  // get
  T x() const { return x_; }
  T y() const { return y_; }
  T z() const { return z_; }
  T w() const { return w_; }

  T getX() const { return x_; }
  T getY() const { return y_; }
  T getZ() const { return z_; }
  T getW() const { return w_; }

  void getXYZ(T *x, T *y, T *z) const {
    T w1 = 1.0/w_;

    *x = x_*w1; *y = y_*w1; *z = z_*w1;
  }

  void getXYZW(T *x, T *y, T *z, T *w) const {
    *x = x_; *y = y_; *z = z_; *w = w_;
  }

  // Reference routine would break encapsulation
  // T &getX() const { return x_; }
  // T &getY() const { return y_; }
  // T &getZ() const { return z_; }
  // T &getW() const { return w_; }

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

  void setW(T w) {
    w_ = w;

    normalized_ = false;
  }

  void setXYZ(T x, T y, T z) {
    x_ = x; y_ = y; z_ = z; w_ = 1.0;

    normalized_ = false;
  }

  void setXYZW(T x, T y, T z, T w) {
    x_ = x; y_ = y; z_ = z; w_ = w;

    normalized_ = false;
  }

  // Reference routine would break encapsulation

//T &operator[](uint i) {
//  assert(i < 3);
//
//  return (&x_)[i];
//}

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

    return ::sqrt(x_*x_ + y_*y_ + z_*z_ + w_*w_);
  }

  T modulus() const {
    return length();
  }

  T lengthSqr() const {
    if (normalized_)
      return 1.0;

    return (x_*x_ + y_*y_ + z_*z_ + w_*w_);
  }

  // bool isZero();
  // bool isUnit();

  //------

  // comparison
  int cmp(const Vector &v) const {
    if      (x_ < v.x_) return -1;
    else if (x_ > v.x_) return  1;
    else if (y_ < v.y_) return -1;
    else if (y_ > v.y_) return  1;
    else if (z_ < v.z_) return -1;
    else if (z_ > v.z_) return  1;
    else if (w_ < v.w_) return -1;
    else if (w_ > v.w_) return  1;
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

  //------

  explicit CVector3DHT(const Point &point) :
    x_(point.x), y_(point.y), z_(point.z), w_(1.0), normalized_(false) {
  }

  CVector3DHT(const Point &point1, const Point &point2) :
   x_(point2.x - point1.x), y_(point2.y - point1.y),
   z_(point2.z - point1.z), w_(1.0), normalized_(false) {
  }

  CVector3DHT(const Vector &vector1, const Vector &vector2) :
   x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_),
   x_(vector2.z_ - vector1.z_), w_(1.0), normalized_(false) {
  }

  //------

  CVector3DHT(Type type) {
    if      (type == ZERO  ) {
      x_ = 0.0; y_ = 0.0; z_ = 0.0; w_ = 0.0;
    }
    else if (type == UNIT  ) {
      x_ = 1.0; y_ = 1.0; z_ = 1.0; w_ = 1.0;
    }
    else if (type == UNIT_X) {
      x_ = 1.0; y_ = 0.0; z_ = 0.0; w_ = 0.0;
    }
    else if (type == UNIT_Y) {
      x_ = 0.0; y_ = 1.0; z_ = 0.0; w_ = 0.0;
    }
    else if (type == UNIT_Z) {
      x_ = 0.0; y_ = 0.0; z_ = 1.0; w_ = 0.0;
    }
    else if (type == UNIT_W) {
      x_ = 0.0; y_ = 0.0; z_ = 0.0; w_ = 1.0;
    }
  }

  //------

  Point point() const {
    T w1 = 1.0/w_;

    return Point(x_*w1, y_*w1, z_*w1);
  }

  //------

  Vector &zero() {
    x_ = 0.0; y_ = 0.0; z_ = 0.0; w_ = 0.0;

    normalized_ = false;

    return *this;
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
    z_ *= factor;
    w_ *= factor;

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

    return Vector(x_*factor, y_*factor, z_*factor, w_*factor, true);
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

    x_ *= factor; y_ *= factor; z_ *= factor; w_ *= factor;

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

  void incZ(T z = 1.0) {
    z_ += z;

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

  void decZ(T z = 1.0) {
    z_ -= z;

    normalized_ = false;
  }

  void scaleX(T x = 1.0) {
    x_ *= x;

    normalized_ = false;
  }

  void scaleY(T y = 1.0) {
    y_ *= y;

    normalized_ = false;
  }

  void scaleZ(T z = 1.0) {
    z_ *= z;

    normalized_ = false;
  }

  //------

  T minComponent() {
    return min(min(min(x_, y_), z_), w_);
  }

  T maxComponent() {
    return max(max(max(x_, y_), z_), w_);
  }

  T minAbsComponent() {
    return min(min(min(::fabs(x_), ::fabs(y_)), ::fabs(z_)), ::fabs(w_));
  }

  T maxAbsComponent() {
    return max(max(max(::fabs(x_), ::fabs(y_)), ::fabs(z_)), ::fabs(w_));
  }

  //------

  static Vector min(const Vector &lhs, const Vector &rhs) {
    return Vector(min(lhs.x_, rhs.x_),
                  min(lhs.y_, rhs.y_),
                  min(lhs.z_, rhs.z_),
                  min(lhs.w_, rhs.w_));
  }

  static Vector max(const Vector &lhs, const Vector &rhs) {
    return Vector(max(lhs.x_, rhs.x_),
                  max(lhs.y_, rhs.y_),
                  max(lhs.z_, rhs.z_),
                  max(lhs.w_, rhs.w_));
  }

  //------

  Vector &operator=(const Point &point) {
    x_ = point.x; y_ = point.y; z_ = point.z; w_ = 1.0;

    normalized_ = false;

    return *this;
  }

  //------

  // operators

  // unary +/-
  Vector operator+() const {
    return Vector(x_, y_, z_, w_);
  }

  Vector operator-() const {
    return Vector(-x_, -y_, -z_, -w_);
  }

  // addition
  Vector &operator+=(const Vector &rhs) {
    x_ += rhs.x_;
    y_ += rhs.y_;
    z_ += rhs.z_;
    w_ += rhs.w_;

    normalized_ = false;

    return *this;
  }

  Vector operator+(const Vector &rhs) const {
    return Vector(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_, w_ + rhs.w_);
  }

  // subtraction
  Vector &operator-=(const Vector &rhs) {
    x_ -= rhs.x_;
    y_ -= rhs.y_;
    z_ -= rhs.z_;
    w_ -= rhs.w_;

    normalized_ = false;

    return *this;
  }

  Vector operator-(const Vector &rhs) const {
    return Vector(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_, w_ - rhs.w_);
  }

  //------

  // scalar multiplication/division
  Vector &operator*=(T rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs; w_ *= rhs;

    normalized_ = false;

    return *this;
  }

  friend Vector operator*(const Vector &lhs, T rhs) {
    return Vector(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs, lhs.w_*rhs);
  }

  friend Vector operator*(T lhs, const Vector &rhs) {
    return Vector(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_, lhs*rhs.w_);
  }

  Vector &operator/=(T rhs) {
    T irhs = 1.0/rhs;

    x_ *= irhs; y_ *= irhs; z_ *= irhs; w_ *= irhs;

    normalized_ = false;

    return *this;
  }

  Vector operator/(T rhs) const {
    T irhs = 1.0/rhs;

    return Vector(x_*irhs, y_*irhs, z_*irhs, w_*irhs);
  }

  //------

  // dot product
  T dotProduct(const Vector &v) const {
    return (x_*v.x_ + y_*v.y_ + z_*v.z_ + w_*v.w_);
  }

//T dotProduct(const Normal &normal) const {
//  return (x_*normal.getX() + y_*normal.getY() + z_*normal.getZ());
//}

  static T dotProduct(const Vector &v1, const Vector &v2) {
    return (v1.x_*v2.x_ + v1.y_*v2.y_ + v1.z_*v2.z_ + v1.w_*v2.w_);
  }

  T dotProduct(T x, T y, T z, T w) const {
    return (x_*x + y_*y + z_*z + w_*w);
  }

  static T dotProduct(const Vector &v1, T x2, T y2, T z2, T w2) {
    return (v1.x_*x2 + v1.y_*y2 + v1.z_*z2, v1.w_*w2);
  }

  T dotProductSelf() const {
    return (x_*x_ + y_*y_ + z_*z_ + w_*w_);
  }

  static T absDotProduct(const Vector &v1, const Vector &v2) {
    return ::fabs(dotProduct(v1, v2));
  }

  //------

  // cross product
  Vector crossProduct(const Vector &v) const {
    return Vector(y_*v.z_ - z_*v.y_,
                  z_*v.w_ - w_*v.z_,
                  w_*v.x_ - x_*v.w_,
                  x_*v.y_ - y_*v.x_,
                  normalized_ && v.normalized_);
  }

  static Vector crossProduct(const Vector &v1, const Vector &v2) {
    return Vector(v1.y_*v2.z_ - v1.z_*v2.y_,
                  v1.z_*v2.w_ - v1.w_*v2.z_,
                  v1.w_*v2.x_ - v1.x_*v2.w_,
                  v1.x_*v2.y_ - v1.y_*v2.x_,
                  v1.normalized_ && v2.normalized_);
  }

  Vector crossProduct(const Point &point) const {
    return Vector(y_*point.z - z_*point.y,
                  z_*point.w - w_*point.z,
                  w_*point.x - x_*point.w,
                  x_*point.y - y_*point.x);
  }

  Vector crossProduct(T x, T y, T z, T w) const {
    return Vector(y_*z - z_*y,
                  z_*w - w_*z,
                  w_*x - x_*w,
                  x_*y - y_*x);
  }

  Vector unitCrossProduct(const Vector &v) const {
    return crossProduct(v).normalize();
  }

  // Note: area of parallelogram is
  // CVector v1(x2 - x1, y2 - y1), v2(x4 - x1, y4 - y1);
  // area = v1.crossProduct(v2).length();

  //------

  friend Point operator+(const Point &lhs, const Vector &rhs) {
    T x, y, z;

    rhs.getXYZ(&x, &y, &z);

    return Point(lhs.x + x, lhs.y + y, lhs.z + z);
  }

  friend Point operator+=(Point &lhs, const Vector &rhs) {
    T x, y, z;

    rhs.getXYZ(&x, &y, &z);

    lhs.x += x; lhs.y += y; lhs.z += z;

    return lhs;
  }

  friend Point operator-(const Point &lhs, const Vector &rhs) {
    T x, y, z;

    rhs.getXYZ(&x, &y, &z);

    return Point(lhs.x - x, lhs.y - y, lhs.z - z);
  }

  friend Point operator-=(Point &lhs, const Vector &rhs) {
    T x, y, z;

    rhs.getXYZ(&x, &y, &z);

    lhs.x -= x; lhs.y -= y; lhs.z -= z;

    return lhs;
  }

  friend Vector operator-(const Point &lhs, const Point &rhs) {
    return Vector(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z, 1.0);
  }

  //------

  Vector normal(const Vector &vector2) {
    return crossProduct(vector2).normalize();
  }

  static Vector normal(const Vector &vector1, const Vector &vector2) {
    return crossProduct(vector1, vector2).normalize();
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
    T dot = vector1.dotProduct(vector2);

    if (! vector1.normalized_)
      dot /= vector1.length();

    if (! vector2.normalized_)
      dot /= vector2.length();

    return dot;
  }

  //------

  void directionCosines(T *x, T *y, T *z) {
    Vector vector1 = normalized();

    *x = vector1.x;
    *y = vector1.y;
    *z = vector1.z;
   }

  //------

  static void coordinateSystem(const Vector &vector1,
                               Vector *vector2, Vector *vector3) {
    if (! vector1.normalized_) {
      coordinateSystem(vector1.normalized(), vector2, vector3);
      return;
    }

    if      (::fabs(vector1.x_) > ::fabs(vector1.y_)) {
      T ilen = 1.0/::sqrt(vector1.x_*vector1.x_ + vector1.z_*vector1.z_);

      *vector2 = Vector(-vector1.z_*ilen, 0.0, vector1.x_*ilen, 0.0);
    }
    else if (::fabs(vector1.y_) > ::fabs(vector1.z_)) {
      T ilen = 1.0/::sqrt(vector1.y_*vector1.y_ + vector1.z_*vector1.z_);

      *vector2 = Vector(0.0, vector1.z_*ilen, -vector1.y_*ilen, 0.0);
    }
    else {
      T ilen = 1.0/::sqrt(vector1.z_*vector1.z_ + vector1.w_*vector1.w_);

      *vector2 = Vector(0.0, 0.0, vector1.w_*ilen, -vector1.z_*ilen);
    }

    *vector3 = vector1.crossProduct(*vector2);
  }

  //------

#if 0
  void getBarycentrics(const Vector &vector1, const Vector &vector2,
                       const Vector &vector3, const Vector &vector4,
                       T barycentric[4]) const {
    // compute the vectors relative to V3 of the tetrahedron
    Vector diff[4] = {
      vector1 - vector4,
      vector2 - vector4,
      vector3 - vector4,
      *this   - vector4
    };

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
  CVector3DHT(T x, T y, T z, T w, bool normalized) :
   x_(x), y_(y), z_(z), w_(w), normalized_(normalized) {
  }
};

typedef CVector3DHT<double> CVector3D;
typedef CVector3DHT<float>  CVector3DF;

#endif
