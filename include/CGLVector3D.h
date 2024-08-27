#ifndef CGLVector3D_H
#define CGLVector3D_H

// TODO: use CMath::sqrt() and CMath::isqrt() to allow replacement

#include <CMathGen.h>
#include <CMathUtil.h>
#include <CPoint3D.h>
#include <CMathMacros.h>

class CGLVector3D {
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
  CGLVector3D() { }

  CGLVector3D(float x, float y, float z) :
   x_(x), y_(y), z_(z) {
  }

 ~CGLVector3D() { }

  //------

  // copy operations
  CGLVector3D(const CGLVector3D &vector) :
   x_(vector.x_), y_(vector.y_), z_(vector.z_) {
  }

  CGLVector3D &operator=(const CGLVector3D &v) {
    x_ = v.x_; y_ = v.y_; z_ = v.z_;

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << "," << z_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CGLVector3D &vector) {
    vector.print(os);

    return os;
  }

  //------

  // accessors

  // get
  float x() const { return x_; }
  float y() const { return y_; }
  float z() const { return z_; }

  float getX() const { return x_; }
  float getY() const { return y_; }
  float getZ() const { return z_; }

  void getXYZ(float *x, float *y, float *z) const {
    *x = x_; *y = y_; *z = z_;
  }

  const float *getValues() const { return &x_; }

  // Reference routine would break encapsulation

  float operator[](uint i) const { assert(i < 3); return (&x_)[i]; }

  // set
  CGLVector3D &setX(float x) {
    x_ = x;

    return *this;
  }

  CGLVector3D &setY(float y) {
    y_ = y;

    return *this;
  }

  CGLVector3D &setZ(float z) {
    z_ = z;

    return *this;
  }

  CGLVector3D &setXYZ(float x, float y, float z) {
    x_ = x; y_ = y; z_ = z;

    return *this;
  }

  void iset(uint i, float v) {
    assert(i < 3);

    (&x_)[i] = v;
  }

  // more get accessors
  bool getNormalized() const {
    return fabs(length() - 1.0) < 1E-6;
  }

  float length() const {
    return float(std::sqrt(lengthSqr()));
  }

  float fastLength() const {
    return float(CMathGen::fastDistance(double(x_), double(y_), double(z_)));
  }

  float modulus() const {
    return length();
  }

  float lengthSqr() const {
    return (x_*x_ + y_*y_ + z_*z_);
  }

  //------

  // comparison
  // TODO: tolerance ? use eq()
  int cmp(const CGLVector3D &v) const {
    if      (x_ < v.x_) return -1;
    else if (x_ > v.x_) return  1;
    else if (y_ < v.y_) return -1;
    else if (y_ > v.y_) return  1;
    else if (z_ < v.z_) return -1;
    else if (z_ > v.z_) return  1;
    else                return  0;
  }

  friend bool operator==(const CGLVector3D &lhs, const CGLVector3D &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const CGLVector3D &lhs, const CGLVector3D &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const CGLVector3D &lhs, const CGLVector3D &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const CGLVector3D &lhs, const CGLVector3D &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const CGLVector3D &lhs, const CGLVector3D &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const CGLVector3D &lhs, const CGLVector3D &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  bool eq(const CGLVector3D &rhs) const {
    return fabs(x_ - rhs.x_) < 1E-6 &&
           fabs(y_ - rhs.y_) < 1E-6 &&
           fabs(z_ - rhs.z_) < 1E-6;
  }

  //------

  explicit CGLVector3D(const CPoint3D &point) :
    x_(float(point.x)), y_(float(point.y)), z_(float(point.z)) {
  }

  // v = point2 - point1
  CGLVector3D(const CPoint3D &point1, const CPoint3D &point2) :
   x_(float(point2.x - point1.x)),
   y_(float(point2.y - point1.y)),
   z_(float(point2.z - point1.z)) {
  }

  // v = vector2 - vector1
  CGLVector3D(const CGLVector3D &vector1, const CGLVector3D &vector2) :
   x_(float(vector2.x_ - vector1.x_)),
   y_(float(vector2.y_ - vector1.y_)),
   z_(float(vector2.z_ - vector1.z_)) {
  }

  //------

  // parametric creation
  CGLVector3D(Type type, float s=1.0) {
    if      (type == ZERO  ) { x_ = 0.0; y_ = 0.0; z_ = 0.0; }
    else if (type == UNIT  ) { x_ =   s; y_ =   s; z_ =   s; }
    else if (type == UNIT_X) { x_ =   s; y_ = 0.0; z_ = 0.0; }
    else if (type == UNIT_Y) { x_ = 0.0; y_ =   s; z_ = 0.0; }
    else if (type == UNIT_Z) { x_ = 0.0; y_ = 0.0; z_ =   s; }
  }

  //------

  // to point
  CPoint3D point() const {
    return CPoint3D(x_, y_, z_);
  }

  //------

  CGLVector3D &zero() {
    x_ = 0.0; y_ = 0.0; z_ = 0.0;

    return *this;
  }

  bool isZero() const {
    return CMathUtil::realEq(lengthSqr(), 0);
  }

  bool isUnit() const {
    return CMathUtil::realEq(lengthSqr(), 1);
  }

  //------

  CGLVector3D &normalize() {
    float len = length();

    // assert on len == 0 ?

    float factor = 0.0;

    if (len > 0.0)
      factor = 1.0f/len;

    x_ *= factor;
    y_ *= factor;
    z_ *= factor;

    return *this;
  }

  CGLVector3D normalized() const {
    float len = length();

    // assert on len == 0 ?

    float factor = 0.0;

    if (len > 0.0)
      factor = 1.0f/len;

    return CGLVector3D(x_*factor, y_*factor, z_*factor);
  }

  CGLVector3D unit() const {
    return normalized();
  }

  //------

  CGLVector3D &setMagnitude(float magnitude) {
    float factor = 0.0;

    float len = length();

    if (len > 0.0)
      factor = magnitude/len;

    x_ *= factor; y_ *= factor; z_ *= factor;

    return *this;
  }

  //------

  float getDistance(const CGLVector3D &vector) const {
    CGLVector3D diff = *this - vector;

    return diff.length();
  }

  float getDistanceSqr(const CGLVector3D &vector) const {
    CGLVector3D diff = *this - vector;

    return diff.lengthSqr();
  }

  //------

  void incX(float x = 1.0) { x_ += x; }
  void incY(float y = 1.0) { y_ += y; }
  void incZ(float z = 1.0) { z_ += z; }

  void decX(float x = 1.0) { x_ -= x; }
  void decY(float y = 1.0) { y_ -= y; }
  void decZ(float z = 1.0) { z_ -= z; }

  void scaleX(float x = 1.0) { x_ *= x; }
  void scaleY(float y = 1.0) { y_ *= y; }
  void scaleZ(float z = 1.0) { z_ *= z; }

  //------

  float minComponent() { return std::min(std::min(x_, y_), z_); }
  float maxComponent() { return std::max(std::max(x_, y_), z_); }

  float minAbsComponent() { return float(std::min(std::min(::fabs(x_), ::fabs(y_)), ::fabs(z_))); }
  float maxAbsComponent() { return float(std::max(std::max(::fabs(x_), ::fabs(y_)), ::fabs(z_))); }

  //------

  static CGLVector3D min(const CGLVector3D &lhs, const CGLVector3D &rhs) {
    return CGLVector3D(std::min(lhs.x_, rhs.x_),
                       std::min(lhs.y_, rhs.y_),
                       std::min(lhs.z_, rhs.z_));
  }

  static CGLVector3D max(const CGLVector3D &lhs, const CGLVector3D &rhs) {
    return CGLVector3D(std::max(lhs.x_, rhs.x_),
                       std::max(lhs.y_, rhs.y_),
                       std::max(lhs.z_, rhs.z_));
  }

  //------

  CGLVector3D &operator=(const CPoint3D &point) {
    x_ = float(point.x); y_ = float(point.y); z_ = float(point.z);

    return *this;
  }

  //------

  // operators

  // unary +/-
  CGLVector3D operator+() const {
    return CGLVector3D(x_, y_, z_);
  }

  CGLVector3D operator-() const {
    return CGLVector3D(-x_, -y_, -z_);
  }

  // addition
  CGLVector3D &operator+=(const CGLVector3D &rhs) {
    x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_;

    return *this;
  }

  CGLVector3D operator+(const CGLVector3D &rhs) const {
    return CGLVector3D(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
  }

  // subtraction
  CGLVector3D &operator-=(const CGLVector3D &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_;

    return *this;
  }

  CGLVector3D operator-(const CGLVector3D &rhs) const {
    return CGLVector3D(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
  }

  // scalar multiplication/division
  CGLVector3D &operator*=(const CGLVector3D &rhs) {
    x_ *= rhs.x_; y_ *= rhs.y_; z_ *= rhs.z_;

    return *this;
  }

  CGLVector3D &operator*=(float rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs;

    return *this;
  }

  CGLVector3D operator*(const CGLVector3D &rhs) {
    CGLVector3D t(*this);

    t *= rhs;

    return t;
  }

  friend CGLVector3D operator*(const CGLVector3D &lhs, float rhs) {
    return CGLVector3D(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs);
  }

  friend CGLVector3D operator*(float lhs, const CGLVector3D &rhs) {
    return CGLVector3D(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_);
  }

  CGLVector3D &operator/=(float rhs) {
    float irhs = 1.0f/rhs;

    x_ *= irhs; y_ *= irhs; z_ *= irhs;

    return *this;
  }

  CGLVector3D operator/(float rhs) const {
    float irhs = 1.0f/rhs;

    return CGLVector3D(x_*irhs, y_*irhs, z_*irhs);
  }

  //------

  // dot product
  float dotProduct(const CGLVector3D &v) const {
    return (x_*v.x_ + y_*v.y_ + z_*v.z_);
  }

//float dotProduct(const Normal &normal) const {
//  return (x_*normal.getX() + y_*normal.getY() + z_*normal.getZ());
//}

  static float dotProduct(const CGLVector3D &v1, const CGLVector3D &v2) {
    return v1.dotProduct(v2);
  }

  static float dotProduct(const CPoint3D &v1, const CGLVector3D &v2) {
    return CGLVector3D(float(v1.x), float(v1.y), float(v1.z)).dotProduct(CGLVector3D(v2));
  }

  float dotProduct(const CPoint3D &point) const {
    return dotProduct(CGLVector3D(float(point.x), float(point.y), float(point.z)));
  }

  float dotProduct(float x, float y, float z) const {
    return dotProduct(CGLVector3D(x, y, z));
  }

  static float dotProduct(const CGLVector3D &v1, float x2, float y2, float z2) {
    return v1.dotProduct(CGLVector3D(x2, y2, z2));
  }

  float dotProductSelf() const {
    return dotProduct(*this);
  }

  static float absDotProduct(const CGLVector3D &v1, const CGLVector3D &v2) {
    return float(std::fabs(dotProduct(v1, v2)));
  }

  //------

  // cross product
  CGLVector3D crossProduct(const CGLVector3D &v) const {
    return CGLVector3D(y_*v.z_ - z_*v.y_, z_*v.x_ - x_*v.z_, x_*v.y_ - y_*v.x_);
  }

  static CGLVector3D crossProduct(const CGLVector3D &v1, const CGLVector3D &v2) {
    return v1.crossProduct(v2);
  }

  CGLVector3D crossProduct(const CPoint3D &point) const {
    return crossProduct(CGLVector3D(float(point.x), float(point.y), float(point.z)));
  }

  CGLVector3D crossProduct(float x, float y, float z) const {
    return crossProduct(CGLVector3D(x, y, z));
  }

  CGLVector3D unitCrossProduct(const CGLVector3D &v) const {
    return crossProduct(v).normalize();
  }

  // Note: area of parallelogram is
  // CVector v1(x2 - x1, y2 - y1), v2(x4 - x1, y4 - y1);
  // area = v1.crossProduct(v2).length();

  //------

  friend CPoint3D operator+(const CPoint3D &lhs, const CGLVector3D &rhs) {
    return CPoint3D(float(lhs.x + rhs.x_), float(lhs.y + rhs.y_), float(lhs.z + rhs.z_));
  }

  friend CPoint3D operator+=(CPoint3D &lhs, const CGLVector3D &rhs) {
    lhs.x += rhs.x_; lhs.y += rhs.y_; lhs.z += rhs.z_;

    return lhs;
  }

  friend CPoint3D operator-(const CPoint3D &lhs, const CGLVector3D &rhs) {
    return CPoint3D(lhs.x - rhs.x_, lhs.y - rhs.y_, lhs.z - rhs.z_);
  }

  friend CPoint3D operator-=(CPoint3D &lhs, const CGLVector3D &rhs) {
    lhs.x -= rhs.x_; lhs.y -= rhs.y_; lhs.z -= rhs.z_;

    return lhs;
  }

#if 0
  friend CGLVector3D operator-(const CPoint3D &lhs, const CPoint3D &rhs) {
    return CGLVector3D(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }
#endif

  //------

  CGLVector3D normal(const CGLVector3D &vector2) const {
    return unitCrossProduct(vector2);
  }

  static CGLVector3D normal(const CGLVector3D &vector1, const CGLVector3D &vector2) {
    return vector1.normal(vector2);
  }

  //------

  float cosIncluded(const CGLVector3D &vector1) const {
    float dot = dotProduct(vector1);

    dot /= length();
    dot /= vector1.length();

    return dot;
  }

  static float cosIncluded(const CGLVector3D &vector1, const CGLVector3D &vector2) {
    return vector1.cosIncluded(vector2);
  }

  //------

  void directionCosines(float *x, float *y, float *z) {
    CGLVector3D vector1 = normalized();

    *x = vector1.x_;
    *y = vector1.y_;
    *z = vector1.z_;
   }

  //------

  static void coordinateSystem(const CGLVector3D &vector1, CGLVector3D *vector2,
                               CGLVector3D *vector3) {
    auto v1 = vector1.normalized();

    if (::fabs(v1.x_) > ::fabs(v1.y_)) {
      float ilen = 1.0f/float(std::sqrt(v1.x_*v1.x_ + v1.z_*v1.z_));

      *vector2 = CGLVector3D(-v1.z_*ilen, 0.0, v1.x_*ilen);
    }
    else {
      float ilen = 1.0f/float(std::sqrt(v1.y_*v1.y_ + v1.z_*v1.z_));

      *vector2 = CGLVector3D(0.0, v1.z_*ilen, -v1.y_*ilen);
    }

    *vector3 = v1.crossProduct(*vector2);
  }

  //------

  CGLVector3D perp() const {
    return perp(*this);
  }

  static CGLVector3D perp(const CGLVector3D &u) {
    CGLVector3D u1 = u.unit();
    CGLVector3D v1 = CGLVector3D(0,1,0);

    float vu = v1.dotProduct(u1);

    v1 = (v1 - vu*u1).unit();

    if (v1.isZero()) {
      v1 = CGLVector3D(0,0,1);

      vu = v1.dotProduct(u1);

      v1 = (v1 - vu*u1).unit();
    }

    CGLVector3D w1 = u1.crossProduct(v1);

    v1 = w1.crossProduct(u1);

    return v1;
  }

  //------

#if 0
  void getBarycentrics(const CGLVector3D &vector1, const CGLVector3D &vector2,
                       const CGLVector3D &vector3, const CGLVector3D &vector4,
                       float barycentric[4]) const {
    // compute the vectors relative to V3 of the tetrahedron
    CGLVector3D diff[4] = { vector1 - vector4, vector2 - vector4,
                          vector3 - vector4, *this - vector4 };

    // If the vertices have large magnitude, the linear system of equations
    // for computing barycentric coordinates can be ill-conditioned.  To avoid
    // this, uniformly scale the tetrahedron edges to be of order 1.  The
    // scaling of all differences does not change the barycentric coordinates.
    float maxval = 0.0;

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        maxval = max(maxval, ::fabs(diff[i][j]));

    // scale down only large data
    if (maxval > 1.0) {
      float imaxval = 1.0/maxval;

      for (int i = 0; i < 4; i++)
        diff[i] *= imaxval;
    }

    float det = diff[0].dotProduct(diff[1].crossProduct(diff[2]));

    CGLVector3D c12 = diff[1].crossProduct(diff[2]);
    CGLVector3D c20 = diff[2].crossProduct(diff[0]);
    CGLVector3D c01 = diff[0].crossProduct(diff[1]);

    if (::fabs(det) > 1E-6 ) {
      float idet = 1.0/det;

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
      CGLVector3D d13     = vector1 - vector3;
      CGLVector3D d23     = vector2 - vector3;
      CGLVector3D d13cd23 = d13.crossProduct(d23);

      float max_length = d13cd23.lengthSqr();

      int max_ind = 3;

      float length2 = c01.lengthSqr();

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
        float imax_length = 1.0/max_length;

        CGLVector3D tmp;

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
        float max_length2 = diff[0].lengthSqr();

        max_ind = 0;  // <V0,V3>

        float length2 = diff[1].lengthSqr();

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

        CGLVector3D v12 = vector1 - vector2;

        length2 = v12.lengthSqr();

        if (length2 > max_length2) {
          max_ind     = 5;  // <V0,V1>
          max_length2 = length2;
        }

        if ( max_length2 > 1E-6 ) {
          float imax_length2 = 1.0/max_length2;

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
  float x_ { 0.0f };
  float y_ { 0.0f };
  float z_ { 0.0f };
};

#endif
