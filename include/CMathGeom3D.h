#ifndef CMATH_GEOM_3D_H
#define CMATH_GEOM_3D_H

#include <CPoint3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CTriangle3D.h>
#include <CPolygonOrientation.h>

namespace CMathGeom3D {
  CPolygonOrientation
    PolygonOrientation(double x1, double y1, double z1, double x2, double y2, double z2,
                       double x3, double y3, double z3, double eye_x, double eye_y, double eye_z);

  CPolygonOrientation
    PolygonOrientation(const CPoint3D &point1, const CPoint3D &point2,
                       const CPoint3D &point3, const CPoint3D &eye);

  bool LinePlaneIntersect(double line_x1, double line_y1, double line_z1,
                          double line_x2, double line_y2, double line_z2,
                          double plane_x, double plane_y, double plane_z, double plane_c,
                          double *intersect_x, double *intersect_y,
                          double *intersect_z, double *intersect_param);
  bool LinePlaneIntersect(const CVector3D &vector1, const CVector3D &vector2,
                          const CVector3D &normal, double normalc,
                          CVector3D &intersect_point, double &intersect_param);

  bool LineBetweenLines(const CLine3D &line1, const CLine3D &line2, CLine3D &linei);

  void FaceNormal(double x1, double y1, double z1, double x2, double y2, double z2,
                  double x3, double y3, double z3, double *nx, double *ny, double *nz, double *nc);

  bool Collinear(double x1, double y1, double z1, double x2, double y2, double z2,
                 double x3, double y3, double z3);

  double DotProduct(double x1, double y1, double z1, double x2, double y2, double z2);
  double DotProduct(double x1, double y1, double z1, const CPoint3D &point2);
  double DotProduct(const CPoint3D &point1, double x2, double y2, double z2);
  double DotProduct(const CPoint3D &point1, const CPoint3D &point2);

  void CrossProduct(double x1, double y1, double z1, double x2, double y2, double z2,
                    double *x3, double *y3, double *z3);
  CPoint3D CrossProduct(const CPoint3D &point1, const CPoint3D &point2);

  bool PointLineDistance(const CPoint3D &point, const CLine3D &line, double *dist);

  double PointPointDistance(const CPoint3D &point1, const CPoint3D &point2);

  bool TrianglesCentroid(CTriangle3D *triangles, int num_triangles, CPoint3D &cpoint);

  bool PointInside(double x, double y, double z, double *px, double *py, double *pz, int np);

  double PointAngleSum(double x, double y, double z, double *px, double *py, double *pz, int np);

  bool SphereLineIntersect(double xc, double yc, double zc, double r,
                           double x1, double y1, double z1, double x2, double y2, double z2,
                           double *xi1, double *yi1, double *zi1,
                           double *xi2, double *yi2, double *zi2, int *num_i);

  bool FourPointSphere(double x1, double y1, double z1, double x2, double y2, double z2,
                       double x3, double y3, double z3, double x4, double y4, double z4,
                       double *xc, double *yc, double *zc, double *r);

  double TetrahedronVolume(double x1, double y1, double z1, double x2, double y2, double z2,
                           double x3, double y3, double z3, double x4, double y4, double z4);

  void orthogonalize(const std::vector<CVector3D> &ivectors, std::vector<CVector3D> &ovectors);
}

class CPlane3D;

namespace CMathGeom3D {
  bool LinePlaneIntersect(const CLine3D &line, const CPlane3D &plane,
                          CPoint3D &ipoint, double &iparam);

  bool PlaneIntersect(const CPlane3D &plane1, const CPlane3D &plane2, CLine3D &iline);
  bool PlaneIntersect(const CPlane3D &plane1, const CPlane3D &plane2,
                      const CPlane3D &plane3, CPoint3D &ipoint);

  bool PointPlaneDistance(const CPoint3D &point, const CPlane3D &plane, double *dist);
}

class CNPlane3D;

namespace CMathGeom3D {
  bool LinePlaneIntersect(const CPoint3D &point1, const CPoint3D &point2,
                          const CNPlane3D &normal, CPoint3D &intersect_point,
                          double &intersect_param);
  bool LinePlaneIntersect(const CPoint3D &point, const CVector3D &direction,
                          const CNPlane3D &normal, CPoint3D &intersect_point,
                          double &intersect_param);

  void FaceNormal(const CPoint3D &point1, const CPoint3D &point2,
                  const CPoint3D &point3, CNPlane3D &normal);
  void FaceNormal(double x1, double y1, double z1, double x2, double y2, double z2,
                  double x3, double y3, double z3, CNPlane3D &normal);
}

#endif
