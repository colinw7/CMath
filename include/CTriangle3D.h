#ifndef CTRIANGLE_3D_H
#define CTRIANGLE_3D_H

#include <CShape3D.h>
#include <CLine3D.h>
#include <CVector2D.h>
#include <COptVal.h>
#include <CPolygonOrientation.h>
#include <CMathGeom2D.h>

class CTriangle3D : public CShape3D {
 private:
  CPoint3D            point1_, point2_, point3_;
  COptValT<CVector3D> normal_;
  COptValT<CVector3D> normal1_, normal2_, normal3_;

 public:
  // constructor/destructor
  CTriangle3D(const CPoint3D &point1=CPoint3D(0,0,0),
              const CPoint3D &point2=CPoint3D(1,0,0),
              const CPoint3D &point3=CPoint3D(0.5,1,0)) :
   point1_(point1), point2_(point2), point3_(point3) {
  }

 ~CTriangle3D() { }

  //------

  // copy operations
  CTriangle3D(const CTriangle3D &rhs) :
   CShape3D(rhs), point1_(rhs.point1_), point2_(rhs.point2_),
   point3_(rhs.point3_), normal_(rhs.normal_), normal1_(rhs.normal1_),
   normal2_(rhs.normal2_), normal3_(rhs.normal3_) {
  }

  CTriangle3D &operator=(const CTriangle3D &rhs) {
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
  const CPoint3D &getPoint1() const { return point1_; }
  const CPoint3D &getPoint2() const { return point2_; }
  const CPoint3D &getPoint3() const { return point3_; }

  void setNormals(const CVector3D &n1, const CVector3D &n2, const CVector3D &n3) {
    normal1_ = n1;
    normal2_ = n2;
    normal3_ = n3;
  }

  //------

  CPoint3D centroid() { return 0.333333333*(point1_ + point2_ + point3_); }

  double area() { return 0.5*area2(); }

  double area2() {
    CVector3D p21 = CVector3D(point1_, point2_);
    CVector3D p31 = CVector3D(point1_, point3_);

    return p21.crossProduct(p31).length();
  }

  CBBox3D getBBox() const {
    CPoint3D p1(CPoint3D::min(CPoint3D::min(point1_, point2_), point3_));
    CPoint3D p2(CPoint3D::max(CPoint3D::max(point1_, point2_), point3_));

    return CBBox3D(CShape3D::transformFrom(p1), CShape3D::transformFrom(p2));
  }

  CPolygonOrientation orientationXY() {
    return CMathGeom2D::PolygonOrientation(point1_.x, point1_.y,
                                           point2_.x, point2_.y,
                                           point3_.x, point3_.y);
  }

  bool inside(const CPoint3D &point) const {
    CVector3D v0(point1_, point2_);
    CVector3D v1(point2_, point3_);

    CVector3D n = v0.crossProduct(v1);

    CVector3D wTest = v0.crossProduct(CVector3D(point1_, point));

    if (wTest.dotProduct(n) < 0.0)
      return false;

    wTest = v1.crossProduct(CVector3D(point2_, point));

    if (wTest.dotProduct(n) < 0.0)
      return false;

    CVector3D v2(point3_ - point1_);
    CVector3D v3(point   - point3_);

    wTest = v2.crossProduct(v3);

    if (wTest.dotProduct(n) < 0.0)
      return false;

    return true;
  }

  bool intersect(const CLine3D &line, double *t) const {
    CPoint3D p1 = CShape3D::transformTo(line.start());
    CPoint3D p2 = CShape3D::transformTo(line.end  ());

    CVector3D ld(p1, p2); // p2 - p1

    CVector3D e1(point1_, point2_); // point2_ - point1_
    CVector3D e2(point1_, point3_); // point3_ - point1_

    // compute first barycentric coordinate
    CVector3D s1 = ld.crossProduct(e2);

    double denom = s1.dotProduct(e1);

    if (fabs(denom) < 1E-6) return false;

    double idenom = 1.0/denom;

    CVector3D d(point1_, p1); // p1 - point1_

    double beta = d.dotProduct(s1)*idenom;

    if (beta < 0.0 || beta > 1.0)
      return false;

    // Compute second barycentric coordinate
    CVector3D s2 = d.crossProduct(e1);

    double gamma = ld.dotProduct(s2)*idenom;

    if (gamma < 0.0 || beta + gamma > 1.0)
      return false;

    // Compute intersection param
    *t = e2.dotProduct(s2)*idenom;

    return true;
  }

  void getBarycentrics(const CPoint3D &point, double *u, double *v) const {
    CVector3D v1 = CVector3D(point  , point3_);
    CVector3D v2 = CVector3D(point3_, point1_);
    CVector3D v3 = CVector3D(point3_, point2_);

    double a, b, c, d, e, f, g, h, i;

    v2.getXYZ(&a, &d, &g);
    v3.getXYZ(&b, &e, &h);
    v1.getXYZ(&c, &f, &i);

    if (a == 0 && b == 0) {
      std::swap(a, d);
      std::swap(b, e);
      std::swap(c, f);
   }

    double a1 = b*(f + i) - c*(e + h);
    double a2 = a*(e + h) - b*(d + g);
    double b1 = a*(f + i) - c*(d + g);
    double b2 = b*(d + g) - a*(e + h);

    *u = (a2 != 0 ? a1/a2 : 0.0);
    *v = (b2 != 0 ? b1/b2 : 0.0);
  }

  CVector3D pointNormal(const CPoint3D &point) const {
    CVector3D n;

    if (normal1_.isValid() && normal2_.isValid() && normal3_.isValid()) {
      double u, v;

      getBarycentrics(point, &u, &v);

      double w = 1 - u - v;

      n = u*normal1_.getValue() + v*normal2_.getValue() + w*normal3_.getValue();
    }
    else
      n = getNormal();

    return CShape3D::transformFrom(n);
  }

  const CVector3D &getNormal() const {
    if (! normal_.isValid()) {
      CTriangle3D *th = const_cast<CTriangle3D *>(this);

      CVector3D p01 = CVector3D(point1_, point2_);
      CVector3D p02 = CVector3D(point1_, point3_);

      th->normal_.setValue(CVector3D::crossProduct(p01, p02).unit());
    }

    return normal_.getValue();
  }

  CVector2D pointToSurfaceVector(const CPoint3D &p) const {
    double u, v;

    getBarycentrics(p, &u, &v);

    return CVector2D(u, v);
  }
};

#endif
