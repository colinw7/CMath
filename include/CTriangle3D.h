#ifndef CTRIANGLE_3D_H
#define CTRIANGLE_3D_H

#include <CShape3D.h>
#include <CLine3D.h>
#include <CVector2D.h>
#include <COptVal.h>
#include <CPolygonOrientation.h>
#include <CMathGeom2D.h>

template<typename T>
class CTriangle3DT : public CShape3DT<T> {
 private:
  typedef typename CShape3DT<T>::BBox BBox;
  typedef CVector3DT<T>               Vector;
  typedef CLine3DT<T>                 Line;
  typedef CPoint3DT<T>                Point;
  typedef CTriangle3DT<T>             Triangle;

 private:
  Point            point1_, point2_, point3_;
  COptValT<Vector> normal_;
  COptValT<Vector> normal1_, normal2_, normal3_;

 public:
  // constructor/destructor
  CTriangle3DT(const Point &point1=Point(0,0,0),
               const Point &point2=Point(1,0,0),
               const Point &point3=Point(0.5,1,0)) :
   point1_(point1), point2_(point2), point3_(point3),
   normal_(), normal1_(), normal2_(), normal3_() {
  }

 ~CTriangle3DT() { }

  //------

  // copy operations
  CTriangle3DT(const Triangle &rhs) :
   CShape3DT<T>(rhs), point1_(rhs.point1_), point2_(rhs.point2_),
   point3_(rhs.point3_), normal_(rhs.normal_), normal1_(rhs.normal1_),
   normal2_(rhs.normal2_), normal3_(rhs.normal3_) {
  }

  CTriangle3DT &operator=(const Triangle &rhs) {
    point1_  = rhs.point1_;
    point2_  = rhs.point2_;
    point3_  = rhs.point3_;
    normal_  = rhs.normal_;
    normal1_ = rhs.normal1_;
    normal2_ = rhs.normal2_;
    normal3_ = rhs.normal3_;

    return *this;
  }

  //------

  // accessors
  const Point &getPoint1() const { return point1_; }
  const Point &getPoint2() const { return point2_; }
  const Point &getPoint3() const { return point3_; }

  void setNormals(const Vector &n1, const Vector &n2, const Vector &n3) {
    normal1_ = n1;
    normal2_ = n2;
    normal3_ = n3;
  }

  //------

  Point centroid() { return 0.333333333*(point1_ + point2_ + point3_); }

  T area() { return 0.5*area2(); }

  T area2() {
    Vector p21 = Vector(point1_, point2_);
    Vector p31 = Vector(point1_, point3_);

    return p21.crossProduct(p31).length();
  }

  BBox getBBox() const {
    Point p1(CPoint3D::min(CPoint3D::min(point1_, point2_), point3_));
    Point p2(CPoint3D::max(CPoint3D::max(point1_, point2_), point3_));

    return BBox(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  CPolygonOrientation orientationXY() {
    return CMathGeom2D::PolygonOrientation(point1_.x, point1_.y,
                                           point2_.x, point2_.y,
                                           point3_.x, point3_.y);
  }

  bool inside(const Point &point) const {
    Vector v0(point1_, point2_);
    Vector v1(point2_, point3_);

    Vector n = v0.crossProduct(v1);

    Vector wTest = v0.crossProduct(Vector(point1_, point));

    if (wTest.dotProduct(n) < 0.0)
      return false;

    wTest = v1.crossProduct(Vector(point2_, point));

    if (wTest.dotProduct(n) < 0.0)
      return false;

    Vector v2(point3_ - point1_);

    wTest = v2.crossProduct(point3_, point);

    if (wTest.dotProduct(n) < 0.0)
      return false;

    return true;
  }

  bool intersect(const Line &line, T *t) const {
    Point p1 = CShape3D::transformTo(line.start());
    Point p2 = CShape3D::transformTo(line.end  ());

    Vector ld(p1, p2); // p2 - p1

    Vector e1(point1_, point2_); // point2_ - point1_
    Vector e2(point1_, point3_); // point3_ - point1_

    // compute first barycentric coordinate
    Vector s1 = ld.crossProduct(e2);

    T denom = s1.dotProduct(e1);

    if (fabs(denom) < 1E-6) return false;

    T idenom = 1.0/denom;

    Vector d(point1_, p1); // p1 - point1_

    T beta = d.dotProduct(s1)*idenom;

    if (beta < 0.0 || beta > 1.0)
      return false;

    // Compute second barycentric coordinate
    Vector s2 = d.crossProduct(e1);

    T gamma = ld.dotProduct(s2)*idenom;

    if (gamma < 0.0 || beta + gamma > 1.0)
      return false;

    // Compute intersection param
    *t = e2.dotProduct(s2)*idenom;

    return true;
  }

  void getBarycentrics(const Point &point, T *u, T *v) const {
    Vector v1 = Vector(point  , point3_);
    Vector v2 = Vector(point3_, point1_);
    Vector v3 = Vector(point3_, point2_);

    T a, b, c, d, e, f, g, h, i;

    v2.getXYZ(&a, &d, &g);
    v3.getXYZ(&b, &e, &h);
    v1.getXYZ(&c, &f, &i);

    if (a == 0 && b == 0) {
      std::swap(a, d);
      std::swap(b, e);
      std::swap(c, f);
   }

    T a1 = b*(f + i) - c*(e + h);
    T a2 = a*(e + h) - b*(d + g);
    T b1 = a*(f + i) - c*(d + g);
    T b2 = b*(d + g) - a*(e + h);

    *u = (a2 != 0 ? a1/a2 : 0.0);
    *v = (b2 != 0 ? b1/b2 : 0.0);
  }

  Vector pointNormal(const Point &point) const {
    Vector n;

    if (normal1_.isValid() && normal2_.isValid() && normal3_.isValid()) {
      T u, v;

      getBarycentrics(point, &u, &v);

      T w = 1 - u - v;

      n = u*normal1_.getValue() + v*normal2_.getValue() + w*normal3_.getValue();
    }
    else
      n = getNormal();

    return CShape3D::transformFrom(n);
  }

  const Vector &getNormal() const {
    if (! normal_.isValid()) {
      CTriangle3DT *th = const_cast<CTriangle3DT *>(this);

      Vector p01 = Vector(point1_, point2_);
      Vector p02 = Vector(point1_, point3_);

      th->normal_.setValue(Vector::crossProduct(p01, p02).unit());
    }

    return normal_.getValue();
  }

  CVector2D pointToSurfaceVector(const Point &p) const {
    T u, v;

    getBarycentrics(p, &u, &v);

    return CVector2D(u, v);
  }
};

typedef CTriangle3DT<double> CTriangle3D;

#endif
