#ifndef CVECTOR_3D_H
#define CVECTOR_3D_H

// TODO: use CMath::sqrt() and CMath::isqrt() to allow replacement

#include <CMathGen.h>
#include <CPoint3D.h>
#include <CMathMacros.h>

class CVector3D {
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
  CVector3D() { }

  CVector3D(double x, double y, double z) :
   x_(x), y_(y), z_(z) {
  }

 ~CVector3D() { }

  //------

  // copy operations
  CVector3D(const CVector3D &vector) :
   x_(vector.x_), y_(vector.y_), z_(vector.z_), normalized_(vector.normalized_) {
  }

  CVector3D &operator=(const CVector3D &v) {
    x_ = v.x_; y_ = v.y_; z_ = v.z_;

    normalized_ = v.normalized_;

    return *this;
  }

  //------

  // output
  void print(std::ostream &os) const {
    os << "(" << x_ << "," << y_ << "," << z_ << ")";
  }

  friend std::ostream &operator<<(std::ostream &os, const CVector3D &vector) {
    vector.print(os);

    return os;
  }

  //------

  // accessors

  // get
  double x() const { return x_; }
  double y() const { return y_; }
  double z() const { return z_; }

  double getX() const { return x_; }
  double getY() const { return y_; }
  double getZ() const { return z_; }

  void getXYZ(double *x, double *y, double *z) const {
    *x = x_; *y = y_; *z = z_;
  }

  const double *getValues() const { return &x_; }

  // Reference routine would break encapsulation

  double operator[](uint i) const { assert(i < 3); return (&x_)[i]; }

  // set
  CVector3D &setX(double x) {
    x_ = x;

    normalized_ = false;

    return *this;
  }

  CVector3D &setY(double y) {
    y_ = y;

    normalized_ = false;

    return *this;
  }

  CVector3D &setZ(double z) {
    z_ = z;

    normalized_ = false;

    return *this;
  }

  CVector3D &setXYZ(double x, double y, double z) {
    x_ = x; y_ = y; z_ = z;

    normalized_ = false;

    return *this;
  }

  void iset(uint i, double v) {
    assert(i < 3);

    (&x_)[i] = v;

    normalized_ = false;
  }

  // more get accessors
  bool getNormalized() const {
    return normalized_;
  }

  double length() const {
    if (normalized_)
      return 1.0;

    return ::sqrt(lengthSqr());
  }

  double fastLength() const {
    if (normalized_)
      return 1.0;

    return CMathGen::fastDistance(x_, y_, z_);
  }

  double modulus() const {
    return length();
  }

  double lengthSqr() const {
    if (normalized_)
      return 1.0;

    return (x_*x_ + y_*y_ + z_*z_);
  }

  //------

  // comparison
  // TODO: tolerance ? use eq()
  int cmp(const CVector3D &v) const {
    if      (x_ < v.x_) return -1;
    else if (x_ > v.x_) return  1;
    else if (y_ < v.y_) return -1;
    else if (y_ > v.y_) return  1;
    else if (z_ < v.z_) return -1;
    else if (z_ > v.z_) return  1;
    else                return  0;
  }

  friend bool operator==(const CVector3D &lhs, const CVector3D &rhs) {
    return lhs.cmp(rhs) == 0;
  }

  friend bool operator!=(const CVector3D &lhs, const CVector3D &rhs) {
    return lhs.cmp(rhs) != 0;
  }

  friend bool operator< (const CVector3D &lhs, const CVector3D &rhs) {
    return lhs.cmp(rhs) <  0;
  }

  friend bool operator<=(const CVector3D &lhs, const CVector3D &rhs) {
    return lhs.cmp(rhs) <= 0;
  }

  friend bool operator> (const CVector3D &lhs, const CVector3D &rhs) {
    return lhs.cmp(rhs) >  0;
  }

  friend bool operator>=(const CVector3D &lhs, const CVector3D &rhs) {
    return lhs.cmp(rhs) >= 0;
  }

  bool eq(const CVector3D &rhs) const {
    return fabs(x_ - rhs.x_) < 1E-6 &&
           fabs(y_ - rhs.y_) < 1E-6 &&
           fabs(z_ - rhs.z_) < 1E-6;
  }

  //------

  explicit CVector3D(const CPoint3D &point) :
    x_(point.x), y_(point.y), z_(point.z) {
  }

  // v = point2 - point1
  CVector3D(const CPoint3D &point1, const CPoint3D &point2) :
   x_(point2.x - point1.x), y_(point2.y - point1.y), z_(point2.z - point1.z) {
  }

  // v = vector2 - vector1
  CVector3D(const CVector3D &vector1, const CVector3D &vector2) :
   x_(vector2.x_ - vector1.x_), y_(vector2.y_ - vector1.y_), z_(vector2.z_ - vector1.z_) {
  }

  //------

  // parametric creation
  CVector3D(Type type, double s=1.0) {
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

  CVector3D &zero() {
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

  CVector3D &normalize() {
    if (normalized_)
      return *this;

    double len = length();

    // assert on len == 0 ?

    double factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    x_ *= factor;
    y_ *= factor;
    z_ *= factor;

    normalized_ = true;

    return *this;
  }

  CVector3D normalized() const {
    if (normalized_)
      return *this;

    double len = length();

    // assert on len == 0 ?

    double factor = 0.0;

    if (len > 0.0)
      factor = 1.0/len;

    return CVector3D(x_*factor, y_*factor, z_*factor, true);
  }

  CVector3D unit() const {
    return normalized();
  }

  //------

  CVector3D &setMagnitude(double magnitude) {
    double factor = 0.0;

    if (normalized_)
      factor = magnitude;
    else {
      double len = length();

      if (len > 0.0)
        factor = magnitude/len;
    }

    x_ *= factor; y_ *= factor; z_ *= factor;

    normalized_ = false;

    return *this;
  }

  //------

  double getDistance(const CVector3D &vector) const {
    CVector3D diff = *this - vector;

    return diff.length();
  }

  double getDistanceSqr(const CVector3D &vector) const {
    CVector3D diff = *this - vector;

    return diff.lengthSqr();
  }

  //------

  void incX(double x = 1.0) { x_ += x; normalized_ = false; }
  void incY(double y = 1.0) { y_ += y; normalized_ = false; }
  void incZ(double z = 1.0) { z_ += z; normalized_ = false; }

  void decX(double x = 1.0) { x_ -= x; normalized_ = false; }
  void decY(double y = 1.0) { y_ -= y; normalized_ = false; }
  void decZ(double z = 1.0) { z_ -= z; normalized_ = false; }

  void scaleX(double x = 1.0) { x_ *= x; normalized_ = false; }
  void scaleY(double y = 1.0) { y_ *= y; normalized_ = false; }
  void scaleZ(double z = 1.0) { z_ *= z; normalized_ = false; }

  //------

  double minComponent() { return std::min(std::min(x_, y_), z_); }
  double maxComponent() { return std::max(std::max(x_, y_), z_); }

  double minAbsComponent() { return std::min(std::min(::fabs(x_), ::fabs(y_)), ::fabs(z_)); }
  double maxAbsComponent() { return std::max(std::max(::fabs(x_), ::fabs(y_)), ::fabs(z_)); }

  //------

  static CVector3D min(const CVector3D &lhs, const CVector3D &rhs) {
    return CVector3D(std::min(lhs.x_, rhs.x_), std::min(lhs.y_, rhs.y_), std::min(lhs.z_, rhs.z_));
  }

  static CVector3D max(const CVector3D &lhs, const CVector3D &rhs) {
    return CVector3D(std::max(lhs.x_, rhs.x_), std::max(lhs.y_, rhs.y_), std::max(lhs.z_, rhs.z_));
  }

  //------

  CVector3D &operator=(const CPoint3D &point) {
    x_ = point.x; y_ = point.y; z_ = point.z;

    normalized_ = false;

    return *this;
  }

  //------

  // operators

  // unary +/-
  CVector3D operator+() const {
    return CVector3D(x_, y_, z_, normalized_);
  }

  CVector3D operator-() const {
    return CVector3D(-x_, -y_, -z_, normalized_);
  }

  // addition
  CVector3D &operator+=(const CVector3D &rhs) {
    x_ += rhs.x_; y_ += rhs.y_; z_ += rhs.z_;

    normalized_ = false;

    return *this;
  }

  CVector3D operator+(const CVector3D &rhs) const {
    return CVector3D(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
  }

  // subtraction
  CVector3D &operator-=(const CVector3D &rhs) {
    x_ -= rhs.x_; y_ -= rhs.y_; z_ -= rhs.z_;

    normalized_ = false;

    return *this;
  }

  CVector3D operator-(const CVector3D &rhs) const {
    return CVector3D(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
  }

  // scalar multiplication/division
  CVector3D &operator*=(const CVector3D &rhs) {
    x_ *= rhs.x_; y_ *= rhs.y_; z_ *= rhs.z_;

    normalized_ = false;

    return *this;
  }

  CVector3D &operator*=(double rhs) {
    x_ *= rhs; y_ *= rhs; z_ *= rhs;

    normalized_ = false;

    return *this;
  }

  CVector3D operator*(const CVector3D &rhs) {
    CVector3D t(*this);

    t *= rhs;

    return t;
  }

  friend CVector3D operator*(const CVector3D &lhs, double rhs) {
    return CVector3D(lhs.x_*rhs, lhs.y_*rhs, lhs.z_*rhs);
  }

  friend CVector3D operator*(double lhs, const CVector3D &rhs) {
    return CVector3D(lhs*rhs.x_, lhs*rhs.y_, lhs*rhs.z_);
  }

  CVector3D &operator/=(double rhs) {
    double irhs = 1.0/rhs;

    x_ *= irhs; y_ *= irhs; z_ *= irhs;

    normalized_ = false;

    return *this;
  }

  CVector3D operator/(double rhs) const {
    double irhs = 1.0/rhs;

    return CVector3D(x_*irhs, y_*irhs, z_*irhs);
  }

  //------

  // dot product
  double dotProduct(const CVector3D &v) const {
    return (x_*v.x_ + y_*v.y_ + z_*v.z_);
  }

//double dotProduct(const Normal &normal) const {
//  return (x_*normal.getX() + y_*normal.getY() + z_*normal.getZ());
//}

  static double dotProduct(const CVector3D &v1, const CVector3D &v2) {
    return v1.dotProduct(v2);
  }

  static double dotProduct(const CPoint3D &v1, const CVector3D &v2) {
    return CVector3D(v1.x, v1.y, v1.z).dotProduct(CVector3D(v2));
  }

  double dotProduct(const CPoint3D &point) const {
    return dotProduct(CVector3D(point.x, point.y, point.z));
  }

  double dotProduct(double x, double y, double z) const {
    return dotProduct(CVector3D(x, y, z));
  }

  static double dotProduct(const CVector3D &v1, double x2, double y2, double z2) {
    return v1.dotProduct(CVector3D(x2, y2, z2));
  }

  double dotProductSelf() const {
    return dotProduct(*this);
  }

  static double absDotProduct(const CVector3D &v1, const CVector3D &v2) {
    return ::fabs(dotProduct(v1, v2));
  }

  //------

  // cross product
  CVector3D crossProduct(const CVector3D &v) const {
    return CVector3D(y_*v.z_ - z_*v.y_, z_*v.x_ - x_*v.z_, x_*v.y_ - y_*v.x_);
  }

  static CVector3D crossProduct(const CVector3D &v1, const CVector3D &v2) {
    return v1.crossProduct(v2);
  }

  CVector3D crossProduct(const CPoint3D &point) const {
    return crossProduct(CVector3D(point.x, point.y, point.z));
  }

  CVector3D crossProduct(double x, double y, double z) const {
    return crossProduct(CVector3D(x, y, z));
  }

  CVector3D unitCrossProduct(const CVector3D &v) const {
    return crossProduct(v).normalize();
  }

  // Note: area of parallelogram is
  // CVector v1(x2 - x1, y2 - y1), v2(x4 - x1, y4 - y1);
  // area = v1.crossProduct(v2).length();

  //------

  friend CPoint3D operator+(const CPoint3D &lhs, const CVector3D &rhs) {
    return CPoint3D(lhs.x + rhs.x_, lhs.y + rhs.y_, lhs.z + rhs.z_);
  }

  friend CPoint3D operator+=(CPoint3D &lhs, const CVector3D &rhs) {
    lhs.x += rhs.x_; lhs.y += rhs.y_; lhs.z += rhs.z_;

    return lhs;
  }

  friend CPoint3D operator-(const CPoint3D &lhs, const CVector3D &rhs) {
    return CPoint3D(lhs.x - rhs.x_, lhs.y - rhs.y_, lhs.z - rhs.z_);
  }

  friend CPoint3D operator-=(CPoint3D &lhs, const CVector3D &rhs) {
    lhs.x -= rhs.x_; lhs.y -= rhs.y_; lhs.z -= rhs.z_;

    return lhs;
  }

#if 0
  friend CVector3D operator-(const CPoint3D &lhs, const CPoint3D &rhs) {
    return CVector3D(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }
#endif

  //------

  CVector3D normal(const CVector3D &vector2) const {
    return unitCrossProduct(vector2);
  }

  static CVector3D normal(const CVector3D &vector1, const CVector3D &vector2) {
    return vector1.normal(vector2);
  }

  //------

  double cosIncluded(const CVector3D &vector1) const {
    double dot = dotProduct(vector1);

    if (! normalized_)
      dot /= length();

    if (! vector1.normalized_)
      dot /= vector1.length();

    return dot;
  }

  static double cosIncluded(const CVector3D &vector1, const CVector3D &vector2) {
    return vector1.cosIncluded(vector2);
  }

  //------

  void directionCosines(double *x, double *y, double *z) {
    CVector3D vector1 = normalized();

    *x = vector1.x_;
    *y = vector1.y_;
    *z = vector1.z_;
   }

  //------

  static void coordinateSystem(const CVector3D &vector1, CVector3D *vector2, CVector3D *vector3) {
    if (! vector1.normalized_) {
      coordinateSystem(vector1.normalized(), vector2, vector3);
      return;
    }

    if (::fabs(vector1.x_) > ::fabs(vector1.y_)) {
      double ilen = 1.0/::sqrt(vector1.x_*vector1.x_ + vector1.z_*vector1.z_);

      *vector2 = CVector3D(-vector1.z_*ilen, 0.0, vector1.x_*ilen);
    }
    else {
      double ilen = 1.0/::sqrt(vector1.y_*vector1.y_ + vector1.z_*vector1.z_);

      *vector2 = CVector3D(0.0, vector1.z_*ilen, -vector1.y_*ilen);
    }

    *vector3 = vector1.crossProduct(*vector2);
  }

  //------

  CVector3D perp() const {
    return perp(*this);
  }

  static CVector3D perp(const CVector3D &u) {
    CVector3D u1 = u.unit();
    CVector3D v1 = CVector3D(0,1,0);

    double vu = v1.dotProduct(u1);

    v1 = (v1 - vu*u1).unit();

    if (v1.isZero()) {
      v1 = CVector3D(0,0,1);

      vu = v1.dotProduct(u1);

      v1 = (v1 - vu*u1).unit();
    }

    CVector3D w1 = u1.crossProduct(v1);

    v1 = w1.crossProduct(u1);

    return v1;
  }

  //------

#if 0
  void getBarycentrics(const CVector3D &vector1, const CVector3D &vector2,
                       const CVector3D &vector3, const CVector3D &vector4,
                       double barycentric[4]) const {
    // compute the vectors relative to V3 of the tetrahedron
    CVector3D diff[4] = { vector1 - vector4, vector2 - vector4, vector3 - vector4, *this - vector4 };

    // If the vertices have large magnitude, the linear system of equations
    // for computing barycentric coordinates can be ill-conditioned.  To avoid
    // this, uniformly scale the tetrahedron edges to be of order 1.  The
    // scaling of all differences does not change the barycentric coordinates.
    double maxval = 0.0;

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        maxval = max(maxval, ::fabs(diff[i][j]));

    // scale down only large data
    if (maxval > 1.0) {
      double imaxval = 1.0/maxval;

      for (int i = 0; i < 4; i++)
        diff[i] *= imaxval;
    }

    double det = diff[0].dotProduct(diff[1].crossProduct(diff[2]));

    CVector3D c12 = diff[1].crossProduct(diff[2]);
    CVector3D c20 = diff[2].crossProduct(diff[0]);
    CVector3D c01 = diff[0].crossProduct(diff[1]);

    if (::fabs(det) > 1E-6 ) {
      double idet = 1.0/det;

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
      CVector3D d13     = vector1 - vector3;
      CVector3D d23     = vector2 - vector3;
      CVector3D d13cd23 = d13.crossProduct(d23);

      double max_length = d13cd23.lengthSqr();

      int max_ind = 3;

      double length2 = c01.lengthSqr();

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
        double imax_length = 1.0/max_length;

        CVector3D tmp;

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
        double max_length2 = diff[0].lengthSqr();

        max_ind = 0;  // <V0,V3>

        double length2 = diff[1].lengthSqr();

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

        CVector3D v12 = vector1 - vector2;

        length2 = v12.lengthSqr();

        if (length2 > max_length2) {
          max_ind     = 5;  // <V0,V1>
          max_length2 = length2;
        }

        if ( max_length2 > 1E-6 ) {
          double imax_length2 = 1.0/max_length2;

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
  CVector3D(double x, double y, double z, bool normalized) :
   x_(x), y_(y), z_(z), normalized_(normalized) {
  }

 private:
  double x_ { 0 }, y_ { 0 }, z_ { 0 };
  bool   normalized_ { false };
};

#endif
