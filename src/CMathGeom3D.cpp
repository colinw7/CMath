#include <CMathGeom3D.h>
#include <CPlane3D.h>
#include <CNPlane3D.h>
#include <CMathGen.h>
#include <CPoint3D.h>
#include <CVector3D.h>
#include <CLine3D.h>
#include <CTriangle3D.h>
#include <CMatrix3D.h>
#include <CMatrix3DH.h>

//! polygon orientation - clockwise or anti-clockwise
CPolygonOrientation
CMathGeom3D::
PolygonOrientation(double x1, double y1, double z1,
                   double x2, double y2, double z2,
                   double x3, double y3, double z3,
                   double eye_x, double eye_y, double eye_z)
{
  double off_x = eye_x - x1;
  double off_y = eye_y - y1;
  double off_z = eye_z - z1;

  double xd1d2, yd1d2, zd1d2;

  CrossProduct(x2 - x1, y2 - y1, z2 - z1, x3 - x2, y3 - y2, z3 - z2, &xd1d2, &yd1d2, &zd1d2);

  double dotprod = DotProduct(xd1d2, yd1d2, zd1d2, off_x, off_y, off_z);

  return (CPolygonOrientation) CMathGen::sign(dotprod);
}

//! polygon orientation - clockwise or anti-clockwise
CPolygonOrientation
CMathGeom3D::
PolygonOrientation(const CPoint3D &point1, const CPoint3D &point2,
                   const CPoint3D &point3, const CPoint3D &eye)
{
  CVector3D off(point1, eye);

  CVector3D d1d2 =
    CVector3D(point1, point2).crossProduct(CVector3D(point2, point3));

  double dotprod = d1d2.dotProduct(off);

  return (CPolygonOrientation) CMathGen::sign(dotprod);
}

//------

//! intersection of line and plane
bool
CMathGeom3D::
LinePlaneIntersect(double line_x1, double line_y1, double line_z1,
                   double line_x2, double line_y2, double line_z2,
                   double plane_x, double plane_y, double plane_z,
                   double plane_c, double *intersect_x,
                   double *intersect_y, double *intersect_z, double *iparam)
{
  double line_dx = line_x2 - line_x1;
  double line_dy = line_y2 - line_y1;
  double line_dz = line_z2 - line_z1;

  double d1 = DotProduct(line_dx, line_dy, line_dz, plane_x, plane_y, plane_z);

  if (fabs(d1) < CMathGen::EPSILON_E6)
    return false;

  double d2 = DotProduct(line_x1, line_y1, line_z1, plane_x, plane_y, plane_z);

  *iparam = (plane_c - d2)/d1;

  *intersect_x = line_x1 + (*iparam)*line_dx;
  *intersect_y = line_y1 + (*iparam)*line_dy;
  *intersect_z = line_z1 + (*iparam)*line_dz;

  return true;
}

//! intersection of line and plane
bool
CMathGeom3D::
LinePlaneIntersect(const CPoint3D &point1, const CPoint3D &point2,
                   const CNPlane3D &normal, CPoint3D &ipoint, double &iparam)
{
  CVector3D dvector(point1, point2);

  return LinePlaneIntersect(point1, dvector, normal, ipoint, iparam);
}

//! intersection of line and plane
bool
CMathGeom3D::
LinePlaneIntersect(const CPoint3D &point, const CVector3D &direction,
                   const CNPlane3D &normal, CPoint3D &ipoint, double &iparam)
{
  double d1 = direction.dotProduct(normal.getDirection());

  if (fabs(d1) < CMathGen::EPSILON_E6)
    return false;

  CVector3D vector1(point);

  double d2 = vector1.dotProduct(normal.getDirection());

  iparam = (normal.getScalar() - d2)/d1;
  ipoint = point + direction*iparam;

  return true;
}

//! intersection of line and plane
bool
CMathGeom3D::
LinePlaneIntersect(const CVector3D &vector1, const CVector3D &vector2,
                   const CVector3D &normal, double normalc, CVector3D &ipoint, double &iparam)
{
  CVector3D direction = vector2 - vector1;

  double d1 = direction.dotProduct(normal);

  if (fabs(d1) < CMathGen::EPSILON_E6)
    return false;

  double d2 = vector1.dotProduct(normal);

  iparam = (normalc - d2)/d1;
  ipoint = vector1 + direction*iparam;

  return true;
}

//! intersection of line and plane
bool
CMathGeom3D::
LinePlaneIntersect(const CLine3D &line, const CPlane3D &plane, CPoint3D &ipoint, double &iparam)
{
  const CVector3D &normal = plane.getNormal();

  double d1 = normal.dotProduct(line.vector());

  if (fabs(d1) < CMathGen::EPSILON_E6)
    return false;

  double d2 = normal.dotProduct(CVector3D(line.start(), plane.getPoint()));

  iparam = d2/d1;
  ipoint = line.point(iparam);

  return true;
}

//--------

//! perpendicular connecting line between two lines
bool
CMathGeom3D::
LineBetweenLines(const CLine3D &line1, const CLine3D &line2, CLine3D &linei)
{
  const CVector3D &v1 = line1.vector();
  const CVector3D &v2 = line2.vector();

  const CVector3D vd(line2.start(), line1.start());

  double v22 = v2.dotProduct(v2);

  if (fabs(v22) < CMathGen::EPSILON_E6)
    return false;

  double v11 = v1.dotProduct(v1);
  double v21 = v2.dotProduct(v1);

  double d = v11*v22 - v21*v21;

  if (fabs(d) < CMathGen::EPSILON_E6)
    return false;

  double vd1 = vd.dotProduct(v1);
  double vd2 = vd.dotProduct(v2);

  double mu1 = (vd2*v21 - vd1*v22)/d;

  double mu2 = (vd2 + mu1*v21)/v22;

  linei = CLine3D(line1.point(mu1), line2.point(mu2));

  return true;
}

//--------

//! line intersection between two planes
bool
CMathGeom3D::
PlaneIntersect(const CPlane3D &plane1, const CPlane3D &plane2, CLine3D &iline)
{
  const CVector3D &normal1 = plane1.getNormal();
  const CVector3D &normal2 = plane2.getNormal();

  double n11 = normal1.dotProduct(normal1);
  double n12 = normal1.dotProduct(normal2);
  double n22 = normal2.dotProduct(normal2);

  double det = n11*n22 - n12*n12;

  if (fabs(det) < CMathGen::EPSILON_E6)
    return false;

  double idet = 1.0/det;

  double c1 = (plane1.getConstant()*n22 - plane2.getConstant())*n12*idet;
  double c2 = (plane2.getConstant()*n11 - plane1.getConstant())*n12*idet;

  CVector3D p = normal1.crossProduct(normal2);

  CPoint3D p1 = (c1*normal1 + c2*normal2).point();
  CPoint3D p2 = p1 + p;

  iline = CLine3D(p1, p2);

  return true;
}

bool
CMathGeom3D::
PlaneIntersect(const CPlane3D &plane1, const CPlane3D &plane2,
               const CPlane3D &plane3, CPoint3D &ipoint)
{
  const CVector3D &normal1 = plane1.getNormal();
  const CVector3D &normal2 = plane2.getNormal();
  const CVector3D &normal3 = plane3.getNormal();

  double d = normal1.dotProduct(normal2.crossProduct(normal3));

  if (fabs(d) < CMathGen::EPSILON_E6)
    return false;

  double id = 1.0/d;

  ipoint = (plane1.getConstant()*normal2.crossProduct(normal3) +
            plane2.getConstant()*normal3.crossProduct(normal1) +
            plane3.getConstant()*normal1.crossProduct(normal2)).point();

  ipoint *= id;

  return true;
}

void
CMathGeom3D::
FaceNormal(double x1, double y1, double z1, double x2, double y2, double z2,
           double x3, double y3, double z3, double *nx, double *ny, double *nz, double *nc)
{
  double xdiff1 = x2 - x1, ydiff1 = y2 - y1, zdiff1 = z2 - z1;
  double xdiff2 = x3 - x2, ydiff2 = y3 - y2, zdiff2 = z3 - z2;

  CrossProduct(xdiff1, ydiff1, zdiff1, xdiff2, ydiff2, zdiff2, nx, ny, nz);

  *nc = DotProduct(*nx, *ny, *nz, x1, y1, z1);
}

void
CMathGeom3D::
FaceNormal(double x1, double y1, double z1, double x2, double y2, double z2,
           double x3, double y3, double z3, CNPlane3D &normal)
{
  CVector3D diff1(x2 - x1, y2 - y1, z2 - z1);
  CVector3D diff2(x3 - x2, y3 - y2, z3 - z2);

  normal.setDirection(CVector3D::crossProduct(diff1, diff2));

  normal.setScalar(CVector3D::dotProduct(normal.getDirection(), x1, y1, z1));
}

void
CMathGeom3D::
FaceNormal(const CPoint3D &point1, const CPoint3D &point2,
           const CPoint3D &point3, CNPlane3D &normal)
{
  CVector3D diff1(point1, point2);
  CVector3D diff2(point2, point3);

  normal.setDirection(CVector3D::crossProduct(diff1, diff2));

  normal.setScalar(CVector3D::dotProduct(normal.getDirection(), point1.x, point1.y, point1.z));
}

bool
CMathGeom3D::
Collinear(double x1, double y1, double z1, double x2, double y2, double z2,
          double x3, double y3, double z3)
{
  double xdiff1 = x2 - x1, ydiff1 = y2 - y1, zdiff1 = z2 - z1;
  double xdiff2 = x3 - x2, ydiff2 = y3 - y2, zdiff2 = z3 - z2;

  double nx, ny, nz;

  CrossProduct(xdiff1, ydiff1, zdiff1, xdiff2, ydiff2, zdiff2, &nx, &ny, &nz);

  return (fabs(nx) < CMathGen::EPSILON_E6 ||
          fabs(ny) < CMathGen::EPSILON_E6 ||
          fabs(nz) < CMathGen::EPSILON_E6);
}

double
CMathGeom3D::
DotProduct(const CPoint3D &point1, const CPoint3D &point2)
{
  return (point1.x*point2.x + point1.y*point2.y + point1.z*point2.z);
}

double
CMathGeom3D::
DotProduct(double x1, double y1, double z1, const CPoint3D &point2)
{
  return (x1*point2.x + y1*point2.y + z1*point2.z);
}

double
CMathGeom3D::
DotProduct(const CPoint3D &point1, double x2, double y2, double z2)
{
  return (point1.x*x2 + point1.y*y2 + point1.z*z2);
}

double
CMathGeom3D::
DotProduct(double x1, double y1, double z1, double x2, double y2, double z2)
{
  return (x1*x2 + y1*y2 + z1*z2);
}

void
CMathGeom3D::
CrossProduct(double x1, double y1, double z1, double x2, double y2, double z2,
             double *x3, double *y3, double *z3)
{
  *x3 = y1*z2 - z1*y2;
  *y3 = z1*x2 - x1*z2;
  *z3 = x1*y2 - y1*x2;
}

CPoint3D
CMathGeom3D::
CrossProduct(const CPoint3D &point1, const CPoint3D &point2)
{
  CPoint3D point;

  point.x = point1.y*point2.z - point1.z*point2.y;
  point.y = point1.z*point2.x - point1.x*point2.z;
  point.z = point1.x*point2.y - point1.y*point2.x;

  return point;
}

bool
CMathGeom3D::
PointPlaneDistance(const CPoint3D &point, const CPlane3D &plane, double *dist)
{
  const CVector3D &normal = plane.getNormal();
  const CPoint3D  &ppoint = plane.getPoint();

  *dist = CVector3D(ppoint, point).dotProduct(normal)/normal.length();

  return true;
}

bool
CMathGeom3D::
PointLineDistance(const CPoint3D &point, const CLine3D &line, double *dist)
{
  CVector3D pl(point, line.start());

  CVector3D l = line.vector();

  double u1 = pl.dotProduct(l);
  double u2 = l .lengthSqr();

  if (u2 <= 0.0)
    return false;

  double u = u1/u2;

  if (u < 0.0 || u > 1.0)
    return false;

  CPoint3D intersection(line.start() + u*l);

  *dist = PointPointDistance(point, intersection);

  return true;
}

double
CMathGeom3D::
PointPointDistance(const CPoint3D &point1, const CPoint3D &point2)
{
  return sqrt((point1.x - point2.x)*(point1.x - point2.x) +
              (point1.y - point2.y)*(point1.y - point2.y) +
              (point1.z - point2.z)*(point1.z - point2.z));
}

bool
CMathGeom3D::
TrianglesCentroid(CTriangle3D *triangles, int num_triangles, CPoint3D &cpoint)
{
  CPoint3D c;
  double   a, d = 0;

  cpoint.zero();

  for (int i = 0; i < num_triangles; ++i) {
    c = triangles[i].centroid();
    a = triangles[i].area2();

    cpoint += a*c;

    d += a;
  }

  if (fabs(d) < CMathGen::EPSILON_E6)
    return false;

  double id = 1.0/d;

  cpoint *= id;

  return true;
}

bool
CMathGeom3D::
PointInside(double x, double y, double z, double *px, double *py, double *pz, int np)
{
  double anglesum = PointAngleSum(x, y, z, px, py, pz, np);

  if (fabs(anglesum - 2.0*M_PI) < CMathGen::EPSILON_E6)
    return true;

  return false;
}

double
CMathGeom3D::
PointAngleSum(double x, double y, double z, double *px, double *py, double *pz, int np)
{
  double m1, m2;
  double x1, y1, z1, x2, y2, z2;
  double anglesum = 0, costheta;

  int i = np - 1;

  for (int j = 0; j < np; i = j++) {
    x1 = px[i] - x;
    y1 = py[i] - y;
    z1 = pz[i] - z;

    x2 = px[j] - x;
    y2 = py[j] - y;
    z2 = pz[j] - z;

    m1 = x1*x1 + y1*y1 + z1*z1;
    m2 = x2*x2 + y2*y2 + z2*z2;

    if (m1*m2 < CMathGen::EPSILON_E6)
      return 2.0*M_PI;

    costheta = (x1*x2 + y1*y2 + z1*z2)/(m1*m2);

    anglesum += acos(costheta);
  }

  return anglesum;
}

bool
CMathGeom3D::
SphereLineIntersect(double xc, double yc, double zc, double r,
                    double x1, double y1, double z1, double x2, double y2, double z2,
                    double *xi1, double *yi1, double *zi1, double *xi2, double *yi2, double *zi2,
                    int *num_i)
{
  double dx = x2 - x1;
  double dy = y2 - y1;
  double dz = z2 - z1;

  double x1c = x1 - xc;
  double y1c = y1 - yc;
  double z1c = z1 - zc;

  double a = dx*dx + dy*dy + dz*dz;

  if (fabs(a) < CMathGen::EPSILON_E6) {
    *num_i = 0;
    return false;
  }

  double i2a = 0.5/a;

  double b = 2.0*(dx*x1c + dy*y1c + dz*z1c);

  double mb_i2a = -b*i2a;

  if (mb_i2a < 0.0 || mb_i2a > 1.0) // Not on line
    return false;

  double c = xc*xc + yc*yc + zc*zc + x1*x1 + y1*y1 + z1*z1 -
             2.0*(xc*x1 + yc*y1 + zc*z1) - r*r;

  double b2_4ac = b*b - 4.0*a*c;

  if (b2_4ac < 0.0) {
    *num_i = 0;
    return false;
  }

  if (fabs(b2_4ac) < CMathGen::EPSILON_E6) {
    *num_i = 1;

    double mu = mb_i2a;

    *xi1 = x1 + mu*dx;
    *yi1 = y1 + mu*dy;
    *zi1 = z1 + mu*dz;

    return true;
  }

  double r_b2_4ac_i2a = sqrt(b2_4ac)*i2a;

  *num_i = 2;

  double mu1 = mb_i2a + r_b2_4ac_i2a;

  *xi1 = x1 + mu1*dx;
  *yi1 = y1 + mu1*dy;
  *xi1 = z1 + mu1*dz;

  double mu2 = mb_i2a - r_b2_4ac_i2a;

  *xi2 = x1 + mu2*dx;
  *yi2 = y1 + mu2*dy;
  *zi2 = z1 + mu2*dz;

  return true;
}

bool
CMathGeom3D::
FourPointSphere(double x1, double y1, double z1, double x2, double y2, double z2,
                double x3, double y3, double z3, double x4, double y4, double z4,
                double *xc, double *yc, double *zc, double *r)
{
  double d1 = x1*x1 + y1*y1 + z1*z1;
  double d2 = x2*x2 + y2*y2 + z2*z2;
  double d3 = x3*x3 + y3*y3 + z3*z3;
  double d4 = x4*x4 + y4*y4 + z4*z4;

  CMatrix3DH m11(x1, y1, z1, 1.0, x2, y2, z2, 1.0, x3, y3, z3, 1.0, x4, y4, z4, 1.0);

  double dm11 = m11.determinant();

  if (fabs(dm11) < CMathGen::EPSILON_E6)
    return false;

  double idm11 = 1.0/dm11;

  CMatrix3DH m12(d1, y1, z1, 1.0, d2, y2, z2, 1.0, d3, y3, z3, 1.0, d4, y4, z4, 1.0);

  double dm12 = m12.determinant();

  CMatrix3DH m13(x1, d1, z1, 1.0, x2, d2, z2, 1.0, x3, d3, z3, 1.0, x4, d4, z4, 1.0);

  double dm13 = m13.determinant();

  CMatrix3DH m14(x1, y1, d1, 1.0, x2, y2, d2, 1.0, x3, y3, d3, 1.0, x4, y4, d4, 1.0);

  double dm14 = m14.determinant();

  CMatrix3DH m15(d1, x1, y1, z1, d2, x2, y2, z2, d3, x3, y3, z3, d4, x4, y4, z4);

  double dm15 = m15.determinant();

  *xc = 0.5*dm12*idm11;
  *yc = 0.5*dm13*idm11;
  *zc = 0.5*dm14*idm11;
  *r  = sqrt((*xc)*(*xc) + (*yc)*(*yc) + (*zc)*(*zc) - dm15*idm11);

  return true;
}

//--------

double
CMathGeom3D::
TetrahedronVolume(double x1, double y1, double z1, double x2, double y2, double z2,
                  double x3, double y3, double z3, double x4, double y4, double z4)
{
  double ax = x1 - x4, ay = y1 - y4, az = z1 - z4;
  double bx = x2 - x4, by = y2 - y4, bz = z2 - z4;
  double cx = x3 - x4, cy = y3 - y4, cz = z3 - z4;

  return ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx);
}

//--------

void
CMathGeom3D::
orthogonalize(const std::vector<CVector3D> &ivectors, std::vector<CVector3D> &ovectors)
{
  ovectors = ivectors;

  uint num = ivectors.size();

  for (uint j = 1; j < num; ++j) {
    for (uint i = 0; i < j - 1; ++i) {
      CVector3D p = ovectors[i]*((ovectors[i].dotProduct(ovectors[j]))/
                                 (ovectors[i].dotProduct(ovectors[i])));

      ovectors[j] = ovectors[j] - p;
    }
  }
}
