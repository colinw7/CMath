#ifndef CMathGeom2D_H
#define CMathGeom2D_H

#include <list>
#include <map>
#include <vector>

#include <CPoint2D.h>
#include <CBBox2D.h>
#include <CPolygonOrientation.h>
#include <CLineJoinType.h>

struct CPoint2DParam {
  CPoint2D p;
  double   t;

  CPoint2DParam(const CPoint2D &tp, double tt) :
   p(tp), t(tt) {
  }
};

struct CPoint2DParam2 {
  CPoint2D p;
  double   t1;
  double   t2;

  CPoint2DParam2(const CPoint2D &tp, double tt1, double tt2) :
   p(tp), t1(tt1), t2(tt2) {
  }
};

namespace CMathGeom2D {
  enum ClipZone {
    CLIP_ZONE_LEFT  = (1<<0),
    CLIP_ZONE_RIGHT = (1<<1),
    CLIP_ZONE_TOP   = (1<<2),
    CLIP_ZONE_BOT   = (1<<3),

    CLIP_ZONE_0  = 0,
    CLIP_ZONE_1  = CLIP_ZONE_LEFT,
    CLIP_ZONE_2  = CLIP_ZONE_RIGHT,
    CLIP_ZONE_4  = CLIP_ZONE_TOP,
    CLIP_ZONE_8  = CLIP_ZONE_BOT,
    CLIP_ZONE_5  = (CLIP_ZONE_LEFT  | CLIP_ZONE_TOP),
    CLIP_ZONE_6  = (CLIP_ZONE_RIGHT | CLIP_ZONE_TOP),
    CLIP_ZONE_9  = (CLIP_ZONE_LEFT  | CLIP_ZONE_BOT),
    CLIP_ZONE_10 = (CLIP_ZONE_RIGHT | CLIP_ZONE_BOT)
  };

  enum PointPosition {
    POINT_POSITION_NONE  = -999,
    POINT_POSITION_LEFT  = -1,
    POINT_POSITION_RIGHT = 1,
    POINT_POSITION_ON    = 0
  };

  //------

  //! get clip zone for point
  inline ClipZone getClipZone(double x, double y, double xmin, double ymin,
                              double xmax, double ymax) {
    ClipZone zone;

    if      (x < xmin) {
      if      (y < ymin)
        zone = CLIP_ZONE_9;
      else if (y <= ymax)
        zone = CLIP_ZONE_1;
      else
        zone = CLIP_ZONE_5;
    }
    else if (x <= xmax) {
      if      (y < ymin)
        zone = CLIP_ZONE_8;
      else if (y <= ymax)
        zone = CLIP_ZONE_0;
      else
        zone = CLIP_ZONE_4;
    }
    else {
      if      (y < ymin)
        zone = CLIP_ZONE_10;
      else if (y <= ymax)
        zone = CLIP_ZONE_2;
      else
        zone = CLIP_ZONE_6;
    }

    return zone;
  }

  //! get clip sector for start and end points of line
  inline bool clipBySector(double xmin, double ymin, double xmax, double ymax,
                           double *x1, double *y1, double *x2, double *y2,
                           ClipZone *zone1, ClipZone *zone2, bool *intersect) {
    static const double EPSILON_E6 = 1E-6;

    if (xmin >= xmax || ymin >= ymax) {
      *intersect = false;

      return true;
    }

    if (*x1 == *x2) {
      if (*x1 < xmin || *x1 > xmax) {
        *intersect = false;

        return true;
      }

      if      (*y1 < ymin)
        *y1 = ymin;
      else if (*y1 > ymax)
        *y1 = ymax;

      if      (*y2 < ymin)
        *y2 = ymin;
      else if (*y2 > ymax)
        *y2 = ymax;

      if (fabs(*x2 - *x1) < EPSILON_E6 && fabs(*y2 - *y1) < EPSILON_E6) {
        *intersect = false;

        return true;
      }

      *intersect = true;

      return true;
    }

    if (*y1 == *y2) {
      if (*y1 < ymin || *y1 > ymax) {
        *intersect = false;

        return true;
      }

      if      (*x1 < xmin)
        *x1 = xmin;
      else if (*x1 > xmax)
        *x1 = xmax;

      if      (*x2 < xmin)
        *x2 = xmin;
      else if (*x2 > xmax)
        *x2 = xmax;

      if (fabs(*x2 - *x1) < EPSILON_E6 && fabs(*y2 - *y1) < EPSILON_E6) {
        *intersect = false;

        return true;
      }

      *intersect = true;

      return true;
    }

    *zone1 = getClipZone(*x1, *y1, xmin, ymin, xmax, ymax);
    *zone2 = getClipZone(*x2, *y2, xmin, ymin, xmax, ymax);

    if (*zone1 & *zone2) {
      *intersect = false;

      return true;
    }

    return false;
  }

  //! Clip line (x1,y1) -> (x2,y2) in region (xmin,xmax,ymin,ymax)
  inline bool clipLine(double xmin, double ymin, double xmax, double ymax,
                       double *x1, double *y1, double *x2, double *y2) {
    static const double EPSILON_E6 = 1E-6;

    bool     intersect;
    ClipZone zone1, zone2;

    if (clipBySector(xmin, ymin, xmax, ymax, x1, y1, x2, y2, &zone1, &zone2, &intersect))
      return intersect;

    double x, y;

    if      (zone1 & CLIP_ZONE_BOT) {
      y = ymin;
      x = *x1 + (*x2 - *x1)*(y - *y1)/(*y2 - *y1);

      *x1 = x;
      *y1 = y;

      zone1 = getClipZone(*x1, *y1, xmin, ymin, xmax, ymax);
    }
    else if (zone1 & CLIP_ZONE_TOP) {
      y = ymax;
      x = *x1 + (*x2 - *x1)*(y - *y1)/(*y2 - *y1);

      *x1 = x;
      *y1 = y;

      zone1 = getClipZone(*x1, *y1, xmin, ymin, xmax, ymax);
    }

    if      (zone1 & CLIP_ZONE_LEFT) {
      x = xmin;
      y = *y1 + (*y2 - *y1)*(x - *x1)/(*x2 - *x1);

      *y1 = y;
      *x1 = x;
    }
    else if (zone1 & CLIP_ZONE_RIGHT) {
      x = xmax;
      y = *y1 + (*y2 - *y1)*(x - *x1)/(*x2 - *x1);

      *y1 = y;
      *x1 = x;
    }

    if      (zone2 & CLIP_ZONE_BOT) {
      y = ymin;
      x = *x2 + (*x1 - *x2)*(y - *y2)/(*y1 - *y2);

      *x2 = x;
      *y2 = y;

      zone2 = getClipZone(*x2, *y2, xmin, ymin, xmax, ymax);
    }
    else if (zone2 & CLIP_ZONE_TOP) {
      y = ymax;
      x = *x2 + (*x1 - *x2)*(y - *y2)/(*y1 - *y2);

      *x2 = x;
      *y2 = y;

      zone2 = getClipZone(*x2, *y2, xmin, ymin, xmax, ymax);
    }

    if      (zone2 & CLIP_ZONE_LEFT) {
      x = xmin;
      y = *y2 + (*y1 - *y2)*(x - *x2)/(*x1 - *x2);

      *y2 = y;
      *x2 = x;
    }
    else if (zone2 & CLIP_ZONE_RIGHT) {
      x = xmax;
      y = *y2 + (*y1 - *y2)*(x - *x2)/(*x1 - *x2);

      *y2 = y;
      *x2 = x;
    }

    if (fabs(*x2 - *x1) < EPSILON_E6 && fabs(*y2 - *y1) < EPSILON_E6)
      return false;

    return true;
  }

  inline bool CircleLineIntersect(double xc, double yc, double r,
                                  double x1, double y1, double x2, double y2,
                                  double *xi1, double *yi1,
                                  double *xi2, double *yi2, uint *num_i) {
    // Transform to origin
    x1 -= xc; y1 -= yc;
    x2 -= xc; y2 -= yc;

    // Calculate discriminant
    double dx = x2 - x1;
    double dy = y2 - y1;

    double dr2 = dx*dx + dy*dy;

    if (dr2 <= 0.0)
      return false;

    double idr2 = 1.0/dr2;

    double dd = x1*y2 - x2*y1;

    double dd2 = dd*dd;

    double r2 = r*r;

    double dis = r2*dr2 - dd2;

    if (dis < 0.0) {
      *num_i = 0;

      return false;
    }

      // Calculate intersection
    if (dis != 0) {
      double sdy = (dy < 0.0 ? -1.0 : 1.0);

      dis = sqrt(dis);

      double dd_idr2 = dd*idr2;

      double dddx_idr2 = dd_idr2*dx;
      double dddy_idr2 = dd_idr2*dy;

      double sdydis_idr2 = sdy*dis*idr2;

      double sdydxdis_idr2 = sdydis_idr2*dx;
      double sdydydis_idr2 = sdydis_idr2*dy;

      if (xi1 != NULL)
        *xi1 =  dddy_idr2 - sdydxdis_idr2 + xc;

      if (yi1 != NULL)
        *yi1 = -dddx_idr2 - sdydydis_idr2 + yc;

      if (xi2 != NULL)
        *xi2 =  dddy_idr2 + sdydxdis_idr2 + xc;

      if (yi2 != NULL)
        *yi2 = -dddx_idr2 + sdydydis_idr2 + yc;

      *num_i = 2;
    }
    else {
      *num_i = 1;

      double dd_idr2 = dd*idr2;

      if (xi1 != NULL)
        *xi1 =  dd_idr2*dy + xc;

      if (yi1 != NULL)
        *yi1 = -dd_idr2*dx + yc;

      if (xi2 != NULL)
        *xi2 = *xi1;

      if (yi2 != NULL)
        *yi2 = *yi1;
    }

    return true;
  }
}

#include <CPolygonOrientation.h>
#include <CVector2D.h>

namespace CMathGeom2D {
  inline CPolygonOrientation PolygonOrientation(double x1, double y1, double x2, double y2,
                                                double x3, double y3) {
    double dx1 = x2 - x1;
    double dy1 = y2 - y1;

    double dx2 = x3 - x2;
    double dy2 = y3 - y2;

    return (CPolygonOrientation) CMathGen::sign(dx1*dy2 - dy1*dx2);
  }

  inline CPolygonOrientation PolygonOrientation(const CPoint2D &point1, const CPoint2D &point2,
                                                const CPoint2D &point3) {
    CVector2D d1(point1, point2);
    CVector2D d2(point2, point3);

    return (CPolygonOrientation) CMathGen::sign(d1.getX()*d2.getY() - d1.getY()*d2.getX());
  }

  inline CPolygonOrientation PolygonOrientation(const double *x, const double *y, uint num_xy) {
    int i = 2;

    while (i < (int) num_xy) {
      CPolygonOrientation orient =
        PolygonOrientation(x[i - 2], y[i - 2], x[i - 1], y[i - 1], x[i], y[i]);

      if (orient != CPOLYGON_ORIENTATION_UNKNOWN)
        return orient;

      ++i;
    }

    return CPOLYGON_ORIENTATION_UNKNOWN;
  }

  inline CPolygonOrientation PolygonOrientation(const std::vector<CPoint2D> &points) {
    int i = 2;

    uint num_points = points.size();

    while (i < (int) num_points) {
      CPolygonOrientation orient = PolygonOrientation(points[i - 2], points[i - 1], points[i]);

      if (orient != CPOLYGON_ORIENTATION_UNKNOWN)
        return orient;

      ++i;
    }

    return CPOLYGON_ORIENTATION_UNKNOWN;
  }

  //------

  inline bool IntersectPolygons(const std::vector<CPoint2D> &points1,
                                const std::vector<CPoint2D> &points2,
                                std::vector<CPoint2D> &ipoints) {
    static const double EPSILON_E6 = 1E-6;

    static CPoint2D *f[2];
    static uint      num_f;

    ipoints.clear();

    uint num_points1 = points1.size();
    uint num_points2 = points2.size();

    // fail if polygons are degenerate
    if (num_points1 < 3 || num_points2 < 3)
      return false;

    CPolygonOrientation orient1 = PolygonOrientation(points1);
    CPolygonOrientation orient2 = PolygonOrientation(points2);

    // max number of intersection
    uint ni = num_points1*num_points2;

    // make sure intersection buffer is large enough
    if (num_f < ni) {
      num_f = ni;

      delete [] f[0];
      delete [] f[1];

      f[0] = new CPoint2D [num_f];
      f[1] = new CPoint2D [num_f];
    }

    // store polygon one in start point array
    // Note: if orients don't match we invert the first polygon's point order
    int l1 = 0;

    ni = num_points1;

    if (orient1 == orient2) {
      for (uint i = 0; i < ni; ++i)
        f[l1][i] = points1[i];
    }
    else {
      for (uint i = 0, j = ni - 1; i < ni; ++i, --j)
        f[l1][i] = points1[j];
    }

    // intersect current set of points with each line (end1, end2)
    // of the second polygon (points2)
    CPoint2D end1 = points2[num_points2 - 1];

    for (uint i = 0; i < num_points2; ++i) {
      CPoint2D end2 = points2[i];

      // l2 is destination point index (inverse of current l1)
      int l2 = 1 - l1;

      // calc line coefficients
      double ca = end2.x - end1.x; // (x2 - x1), (y2 - y1)
      double cb = end1.y - end2.y; // (x2 - x1), (y2 - y1)
      double cc = -end1.x*cb - end1.y*ca; // -x1*(y2 - y1) - y1*(x2 - x1)

      // calc side of line for first point
      CPoint2D v1     = f[l1][ni - 1];
      double   fv1    = ca*v1.y + cb*v1.x + cc;
      double   absfv1 = fabs(fv1);

      int index1 = 0;

      if (absfv1 >= EPSILON_E6)
        index1 = CMathGen::sign(fv1)*orient2;

      int ni1 = 0;

      for (uint j = 0; j < ni; ++j) {
        // calc side of line for second point
        CPoint2D v2     = f[l1][j];
        double   fv2    = ca*v2.y + cb*v2.x + cc;
        double   absfv2 = fabs(fv2);

        int index2 = 0;

        if (absfv2 >= EPSILON_E6)
          index2 = CMathGen::sign(fv2)*orient2;

        // add start point
        if (index1 >= 0)
          f[l2][ni1++] = v1;

        // add intersection point (if changed sides)
        if (index1 != 0 && index1 != index2 && index2 != 0) {
          double delta = absfv1 + absfv2;

          double xi = (absfv2*v1.x + absfv1*v2.x)/delta;
          double yi = (absfv2*v1.y + absfv1*v2.y)/delta;

          f[l2][ni1++] = CPoint2D(xi, yi);
        }

        // move to next line
        v1     = v2;
        absfv1 = absfv2;
        index1 = index2;
      }

      // degenerate result so fail
      if (ni1 < 3)
        return false;

      l1   = l2;
      end1 = end2;
      ni   = ni1;
    }

    ipoints.resize(ni);

    for (uint i = 0; i < ni; ++i)
      ipoints[i] = f[l1][i];

    return true;
  }

  inline bool IntersectPolygons(const double *x1, const double *y1, uint n1,
                                const double *x2, const double *y2, uint n2,
                                double **xi, double **yi, uint *ni) {
    *xi = NULL;
    *yi = NULL;
    *ni = 0;

    // convert input polygons to vector of points
    std::vector<CPoint2D> vpoints1, vpoints2;

    vpoints1.resize(n1);
    vpoints2.resize(n2);

    for (uint i = 0; i < n1; ++i) vpoints1[i] = CPoint2D(x1[i], y1[i]);
    for (uint i = 0; i < n2; ++i) vpoints2[i] = CPoint2D(x2[i], y2[i]);

    // call actual implementation
    std::vector<CPoint2D> vipoints;

    if (! IntersectPolygons(vpoints1, vpoints2, vipoints))
      return false;

    // convert result point vector to return array
    *ni = vipoints.size();
    *xi = new double [*ni + 1];
    *yi = new double [*ni + 1];

    for (uint i = 0; i < *ni; ++i) {
      (*xi)[i] = vipoints[i].x;
      (*yi)[i] = vipoints[i].y;
    }

    return true;
  }

  bool PolygonLineIntersect(const double *x, const double *y, uint nxy,
                            double x1, double y1, double x2, double y2,
                            double *xi, double *yi, uint *num_i);

  inline double TriangleArea2(const CPoint2D &point1, const CPoint2D &point2,
                              const CPoint2D &point3) {
    return (point2.x - point1.x)*(point3.y - point1.y) -
           (point3.x - point1.x)*(point2.y - point1.y);
  }

  inline bool PointLineLeft(const CPoint2D &lpoint1, const CPoint2D &lpoint2,
                            const CPoint2D &point) {
    return TriangleArea2(lpoint1, lpoint2, point) > 0.0;
  }

  inline bool PointLineRight(const CPoint2D &lpoint1, const CPoint2D &lpoint2,
                             const CPoint2D &point) {
    return TriangleArea2(lpoint1, lpoint2, point) < 0.0;
  }

  inline bool PointLineOn(const CPoint2D &lpoint1, const CPoint2D &lpoint2,
                          const CPoint2D &point) {
    return TriangleArea2(lpoint1, lpoint2, point) == 0.0;
  }

  inline bool PointLineLeftOn(const CPoint2D &lpoint1, const CPoint2D &lpoint2,
                              const CPoint2D &point) {
    return TriangleArea2(lpoint1, lpoint2, point) >= 0.0;
  }

  inline bool PointLineRightOn(const CPoint2D &lpoint1, const CPoint2D &lpoint2,
                               const CPoint2D &point) {
    return TriangleArea2(lpoint1, lpoint2, point) >= 0.0;
  }

  void PolygonCentroid(const double *x, const double *y, int num_xy, double *xc, double *yc);
  void PolygonCentroid(const std::vector<CPoint2D> &points, CPoint2D &p);

  inline double IncludedAngle(const CPoint2D &point1, const CPoint2D &point2,
                              const CPoint2D &point3) {
    double s1 = (point2.x - point1.x)*(point2.x - point1.x) +
                (point2.y - point1.y)*(point2.y - point1.y);
    double s2 = (point3.x - point2.x)*(point3.x - point2.x) +
                (point3.y - point2.y)*(point3.y - point2.y);

    if (s1 == 0.0 || s2 == 0.0)
      return 0.0;

    double d1 = sqrt(s1);
    double d2 = sqrt(s2);

    double s3 = (point3.x - point1.x)*(point3.x - point1.x) +
                (point3.y - point1.y)*(point3.y - point1.y);

    double cs = fabs((s1 + s2 - s3)/(2.0*d1*d2));

    double theta = acos(cs);

    return theta;
  }
}

namespace CMathGeom2D {
  inline bool ArcThrough(double x1, double y1, double x2, double y2, double x3, double y3,
                         double xr, double yr, double *xc, double *yc, double *xt1, double *yt1,
                         double *xt2, double *yt2) {
    if (xr <= 0.0 || yr <= 0.0)
      return false;

    double s1 = (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1);
    double s2 = (x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2);
    double s3 = (x3 - x1)*(x3 - x1) + (y3 - y1)*(y3 - y1);

    if (s1 == 0.0 || s2 == 0.0)
      return false;

    double d1 = sqrt(s1);
    double d2 = sqrt(s2);

    double cs = fabs((s1 + s2 - s3)/(2.0*d1*d2));

    double theta = acos(cs);

    if (s1 + s2 - s3 < 0)
      theta = M_PI - theta;

    if (theta == 0.0)
      return false;

    double xd = 0.0;
    double yd = 0.0;

    if (theta != M_PI)  {
      xd = xr/tan(theta/2.0);
      yd = yr/tan(theta/2.0);
    }

    *xt1 = x2 + xd*(x1 - x2)/d1;
    *yt1 = y2 + yd*(y1 - y2)/d1;

    *xt2 = x2 + xd*(x3 - x2)/d2;
    *yt2 = y2 + yd*(y3 - y2)/d2;

    double dc = CMathGen::Hypot(xr, xd);

    double xtm = (*xt1 + *xt2)/2.0;
    double ytm = (*yt1 + *yt2)/2.0;

    double dm = CMathGen::Hypot(xtm - x2, ytm - y2);

    *xc = x2 + dc*(xtm - x2)/dm;
    *yc = y2 + dc*(ytm - y2)/dm;

    return true;
  }

  inline bool ConvertFromSVGArc(double x1, double y1, double x2, double y2,
                                double phi, double rx, double ry, int fa, int fs,
                                double *cx, double *cy, double *theta, double *delta) {
    if (x1 == x2 && y1 == y2) return false;

    if (rx == 0.0 || ry == 0.0) return false;

    // Step 1

    double dx2 = (x1 - x2)/2.0;
    double dy2 = (y1 - y2)/2.0;

    double phi1 = CMathGen::DegToRad(phi);

    double c = cos(phi1);
    double s = sin(phi1);

    double xx =  c*dx2 + s*dy2;
    double yy = -s*dx2 + c*dy2;

    rx = fabs(rx);
    ry = fabs(ry);

    double rx2 = rx*rx;
    double ry2 = ry*ry;
    double xx2 = xx*xx;
    double yy2 = yy*yy;

    double rc = xx2/rx2 + yy2/ry2;

    if (rc > 1) {
      double rc2 = sqrt(rc);

      rx *= rc2;
      ry *= rc2;

      rx2 = rx*rx;
      ry2 = ry*ry;
    }

    // Step 2

    int sign = (fa == fs) ? -1 : 1;

    double sq = (rx2*ry2 - rx2*yy2 - ry2*xx2)/(rx2*yy2 + ry2*xx2);

    sq = (sq < 0) ? 0 : sq;

    double coef = sign*sqrt(sq);

    double cx1 =  coef*((rx*yy)/ry);
    double cy1 = -coef*((ry*xx)/rx);

    // Step 3

    double sx2 = (x1 + x2)/2.0;
    double sy2 = (y1 + y2)/2.0;

    *cx = sx2 + c*cx1 - s*cy1;
    *cy = sy2 + s*cx1 + c*cy1;

    // Step 4

    double ux = ( xx - cx1)/rx;
    double uy = ( yy - cy1)/ry;

    double vx = (-xx - cx1)/rx;
    double vy = (-yy - cy1)/ry;

    double mod_u = CMathGen::Hypot(ux, uy);
    double mod_v = ux;

    sign = (uy < 0) ? -1 : 1;

    *theta = sign*acos(mod_v/mod_u);

    *theta = CMathGen::RadToDeg(*theta);

    while (*theta >  360)
      *theta -= 360;

    while (*theta < -360)
      *theta += 360;

    mod_u = sqrt((ux*ux + uy*uy) * (vx*vx + vy*vy));
    mod_v = ux*vx + uy*vy;

    sign = ((ux*vy - uy*vx) < 0) ? -1 : 1;

    *delta = sign*acos(mod_v/mod_u);

    *delta = CMathGen::RadToDeg(*delta);

    if      (fs == 0 && *delta > 0)
      *delta -= 360;
    else if (fs == 1 && *delta < 0)
      *delta += 360;

    while (*delta > 360)
      *delta -= 360;

    while (*delta < -360)
      *delta += 360;

    return true;
  }

  inline bool ConvertToSVGArc(double cx, double cy, double rx, double ry, double theta,
                              double delta, double phi, double *x0, double *y0,
                              double *x1, double *y1, int *fa, int *fs) {
    // Convert angles to radians

    double theta1 = CMathGen::DegToRad(theta);
    double theta2 = CMathGen::DegToRad(theta + delta);
    double phi_r  = CMathGen::DegToRad(phi);

    // Figure out the coordinates of the beginning and ending points

    *x0 = cx + cos(phi_r) * rx * cos(theta1) + sin(-phi_r) * ry * sin(theta1);
    *y0 = cy + sin(phi_r) * rx * cos(theta1) + cos( phi_r) * ry * sin(theta1);

    *x1 = cx + cos(phi_r) * rx * cos(theta2) + sin(-phi_r) * ry * sin(theta2);
    *y1 = cy + sin(phi_r) * rx * cos(theta2) + cos( phi_r) * ry * sin(theta2);

    *fa = (delta > 180) ? 1 : 0;
    *fs = (delta > 0  ) ? 1 : 0;

    return true;
  }
}

#include <CBBox2D.h>

namespace CMathGeom2D {
  inline bool PointInsideConvex(double x, double y, const double *px, const double *py, uint np) {
    CPolygonOrientation orient = PolygonOrientation(px, py, np);

    if (orient == CPOLYGON_ORIENTATION_UNKNOWN) return false;

    if (orient == CPOLYGON_ORIENTATION_ANTICLOCKWISE) {
      int i1 = np - 1;

      // iterate forwards through polygon lines (i1 -> i2)
      for (int i2 = 0; i2 < (int) np; i1 = i2++) {
        // test orientation of triangle made from point and
        // line matches polygon orientation
        double f = (px[i2] - px[i1])*(y - py[i1]) - (py[i2] - py[i1])*(x - px[i1]);

        if (f < 0)
          return false;
      }
    }
    else {
      int i1 = 0;

      // iterate backwards through polygon lines (i1 -> i2)
      for (int i2 = np - 1; i2 >= 0; i1 = i2--) {
        // test orientation of triangle made from point and
        // line matches polygon orientation
        double f = (px[i2] - px[i1])*(y - py[i1]) - (py[i2] - py[i1])*(x - px[i1]);

        if (f < 0)
          return false;
      }
    }

    return true;
  }

  void EllipsePointAtAngle(double cx, double cy, double xr, double yr, double a,
                           double *x, double *y);

  inline bool PointInsideEvenOdd(const CPoint2D &point, const CPoint2D *points, uint num_points) {
    double xinters;

    int counter = 0;

    int             i2     = num_points - 1;
    const CPoint2D *point2 = &points[i2];

    // iterate through all lines of the polygon
    for (int i1 = 0; i1 < (int) num_points; ++i1) {
      const CPoint2D *point1 = &points[i1];

      // intersect current line with horizontal line at inside point
      if (point.y > std::min(point1->y, point2->y)) {
        if (point.y <= std::max(point1->y, point2->y)) {
          if (point.x <= std::max(point1->x, point2->x)) {
            if (point1->y != point2->y) {
              // if we have an intersection, increase count
              xinters = (point.y   - point1->y)*(point2->x - point1->x)/
                        (point2->y - point1->y) + point1->x;

              if (point1->x == point2->x || point.x <= xinters)
                ++counter;
            }
          }
        }
      }

      // next line
      i2     = i1;
      point2 = point1;
    }

    // if odd then success
    return ((counter % 2) != 0);
  }

  inline bool PointInsideEvenOdd(const CPoint2D &point, const std::vector<CPoint2D> &points) {
    return PointInsideEvenOdd(point, &points[0], points.size());
  }

  inline double PolygonArea2(const double *x, const double *y, uint num_xy) {
    double area = 0.0;

    int i1 = num_xy - 1;
    int i2 = 0;

    for ( ; i2 < (int) num_xy; i1 = i2++)
      area += x[i1]*y[i2] - y[i1]*x[i2];

    return area;
  }

  inline double PolygonArea(const double *x, const double *y, uint num_xy) {
    return fabs(0.5*PolygonArea2(x, y, num_xy));
  }

  inline bool PointBetween(const CPoint2D &lpoint1, const CPoint2D &lpoint2,
                           const CPoint2D &point) {
    if (lpoint1.x != lpoint2.x)
      return (lpoint1.x <= point.x && point.x <= lpoint2.x) ||
             (lpoint1.x >= point.x && point.x >= lpoint2.x);
    else
      return (lpoint1.x <= point.y && point.x <= lpoint2.y) ||
             (lpoint1.x >= point.y && point.x >= lpoint2.y);
  }

  inline bool Intersects(const CPoint2D &l1point1, const CPoint2D &l1point2,
                         const CPoint2D &l2point1, const CPoint2D &l2point2) {
     double area11 = TriangleArea2(l1point1, l1point2, l2point1);

    if (area11 == 0.0 && PointBetween(l1point1, l1point2, l2point1))
      return true;

    double area12 = TriangleArea2(l1point1, l1point2, l2point2);

    if (area12 == 0.0 && PointBetween(l1point1, l1point2, l2point2))
      return true;

    double area21 = TriangleArea2(l2point1, l2point2, l1point1);

    if (area21 == 0.0 && PointBetween(l2point1, l2point2, l1point1))
      return true;

    double area22 = TriangleArea2(l2point1, l2point2, l1point2);

    if (area22 == 0.0 && PointBetween(l2point1, l2point2, l1point1))
      return true;

    return ((area11*area12 < 0.0) && (area21*area22 < 0.0));
  }
}

#include <CLine2D.h>

namespace CMathGeom2D {
  bool PolygonIsConvex(const double *x, const double *y, int num_xy);

  bool PointLineDistance(const CPoint2D &point, const CLine2D &line, double *dist);

  inline bool IntersectLine(double x11, double y11, double x21, double y21,
                            double x12, double y12, double x22, double y22,
                            double *xi, double *yi, double *mu1=0, double *mu2=0) {
    double offset_x1 = x21 - x11;
    double offset_y1 = y21 - y11;
    double offset_x2 = x22 - x12;
    double offset_y2 = y22 - y12;

    double delta = offset_x1*offset_y2 - offset_y1*offset_x2;

    if (fabs(delta) < 1E-6) { // parallel
      if (xi != NULL) *xi = 0.0;
      if (yi != NULL) *yi = 0.0;

      if (mu1 != NULL) *mu1 = 0.0;
      if (mu2 != NULL) *mu2 = 0.0;

      return false;
    }

    double idelta = 1.0/delta;

    double dx = x12 - x11;
    double dy = y12 - y11;

    double m1 = (dx*offset_y1 - dy*offset_x1)*idelta;
    double m2 = (dx*offset_y2 - dy*offset_x2)*idelta;

    if (xi != NULL) *xi = x11 + m2*offset_x1;
    if (yi != NULL) *yi = y12 + m1*offset_y2;

    if (mu1 != NULL) *mu1 = m2;
    if (mu2 != NULL) *mu2 = m1;

    return true;
  }

  inline bool IntersectLine(const CPoint2D &line1s, const CPoint2D &line1e,
                            const CPoint2D &line2s, const CPoint2D &line2e,
                            CPoint2D *point, double *mu1, double *mu2) {
    if (point != NULL)
      return IntersectLine(line1s.x, line1s.y, line1e.x, line1e.y,
                           line2s.x, line2s.y, line2e.x, line2e.y,
                           &point->x, &point->y, mu1, mu2);
    else
      return IntersectLine(line1s.x, line1s.y, line1e.x, line1e.y,
                           line2s.x, line2s.y, line2e.x, line2e.y,
                           NULL, NULL, mu1, mu2);
  }

  inline bool IntersectLine(const CLine2D &line1, const CLine2D &line2,
                            CPoint2D *point, double *mu1, double *mu2) {
    return IntersectLine(line1.start(), line1.end(), line2.start(), line2.end(), point, mu1, mu2);
  }

  bool SlicePolygonByLines(const std::vector<CPoint2D> &poly, const CLine2D &line,
                           std::vector< std::vector<CPoint2D> > &opolys);

  //---

  inline void PointsRange(const std::vector<CPoint2D> &points, CPoint2D &min_point,
                          CPoint2D &max_point) {
    uint num_points = points.size();

    assert(num_points > 0);

    min_point = points[0];
    max_point = min_point;

    for (uint i = 1; i < num_points; ++i) {
      min_point.x = std::min(min_point.x, points[i].x);
      min_point.y = std::min(min_point.y, points[i].y);
      max_point.x = std::max(max_point.x, points[i].x);
      max_point.y = std::max(max_point.y, points[i].y);
    }
  }

  inline void PointsRange(const std::vector<CIPoint2D> &points, CIPoint2D &min_point,
                          CIPoint2D &max_point) {
    uint num_points = points.size();

    assert(num_points > 0);

    min_point = points[0];
    max_point = min_point;

    for (uint i = 1; i < num_points; ++i) {
      min_point.x = std::min(min_point.x, points[i].x);
      min_point.y = std::min(min_point.y, points[i].y);
      max_point.x = std::max(max_point.x, points[i].x);
      max_point.y = std::max(max_point.y, points[i].y);
    }
  }
}

namespace CMathGeom2D {
  inline bool PointInsideRect(double x, double y, double xmin, double ymin,
                              double xmax, double ymax) {
    // ensure rectangle is properly ordered (assert ?)
    if (xmin > xmax) std::swap(xmin, xmax);
    if (ymin > ymax) std::swap(ymin, ymax);

    return (x >= xmin && x <= xmax && y >= ymin && y <= ymax);
  }

  inline bool PointInsideRect(const CPoint2D &point, const CBBox2D &rect) {
    // TODO: needed - CBBox2D has this code already
    return PointInsideRect(point.x, point.y, rect.getXMin(), rect.getYMin(),
                           rect.getXMax(), rect.getYMax());
  }
}

namespace CMathGeom2D {
  bool CutPolygons(const std::vector<CPoint2D> &points1, const std::vector<CPoint2D> &points2,
                   std::vector<CPoint2D> &opoints);
  bool CutPolygons(const std::vector<CPoint2D> &points1, const std::vector<CPoint2D> &points2,
                   std::vector< std::vector<CPoint2D> > &opoints);

  bool AddPolygons(const std::vector<CPoint2D> &points1, const std::vector<CPoint2D> &points2,
                   std::vector<CPoint2D> &opoints);

  //------

  bool PointInsideConvex(const CPoint2D &point, const std::vector<CPoint2D> &points);
  bool PointInsideConvex(const CPoint2D &point, const CPoint2D *points, uint num_points);
  bool PointInsideConvex(double x, double y, const CPoint2D *points, uint num_points);

  bool PointInsideEvenOdd(double x, double y, const CPoint2D *points, uint num_points);
  bool PointInsideEvenOdd(double x, double y, const double *xp, const double *yp, uint np);

  bool PointInsideByAngle(const CPoint2D &point, const std::vector<CPoint2D> &points);
  bool PointInsideByAngle(const CPoint2D &point, const CPoint2D *points, uint num_points);
  bool PointInsideByAngle(double x, double y, const CPoint2D *points, uint num_points);
  bool PointInsideByAngle(double x, double y, const double *xp, const double *yp, uint np);

  bool PointInsideConvex1(double x, double y, const double *xp, const double *yp, uint np);

  //------

  bool ArcThrough(const CPoint2D &point1, const CPoint2D &point2, const CPoint2D &point3,
                  double r, CPoint2D &center, CPoint2D &t1, CPoint2D &t2);
  bool ArcThrough(double x1, double y1, double x2, double y2, double x3, double y3,
                  double xr, double yr, double *xc, double *yc, double *xt1,
                  double *yt1, double *xt2, double *yt2);

  bool LinesAreCoincident(double x11, double y11, double x21, double y21,
                          double x12, double y12, double x22, double y22);
  bool LinesAreCoincident(const CPoint2D &p11, const CPoint2D &p12,
                          const CPoint2D &p21, const CPoint2D &p22);

  double PolygonArea(double x1, double y1, double x2, double y2, double x3, double y3);
  double PolygonArea(double x1, double y1, double x2, double y2,
                     double x3, double y3, double x4, double y4);
  double PolygonArea(const CPoint2D *points, uint num_points);
  double PolygonArea(const std::vector<CPoint2D> &points);

  double PolygonArea2(double x1, double y1, double x2, double y2,
                      double x3, double y3, double x4, double y4);
  double PolygonArea2(const CPoint2D *points, uint num_points);
  double PolygonArea2(const std::vector<CPoint2D> &points);
  double PolygonArea2(const std::list<CPoint2D> &points);

  double TriangleArea(const CPoint2D &point1, const CPoint2D &point2, const CPoint2D &point3);
  double TriangleArea(double x1, double y1, double x2, double y2, double x3, double y3);
  double TriangleArea2(double x1, double y1, double x2, double y2, double x3, double y3);

  double IncludedAngle(double x1, double y1, double x2, double y2, double x3, double y3);
  double IncludedAngle(double x1, double y1, double x2, double y2);

  bool clipLine(double xmin, double ymin, double xmax, double ymax,
                double *x1, double *y1, double *x2, double *y2);
  bool clipLine1(double xmin, double ymin, double xmax, double ymax,
                 double *x1, double *y1, double *x2, double *y2);

  double PointPointDistance(const CPoint2D &point1, const CPoint2D &point2);

  double Hypot(double a, double b);

  bool ThreePointCircle(const CPoint2D &point1, const CPoint2D &point2,
                        const CPoint2D &point3, CPoint2D &center, double &radius);
  bool ThreePointCircle(double x1, double y1, double x2, double y2,
                        double x3, double y3, double *xc, double *yc, double *r);

  PointPosition PointLinePosition(const CPoint2D &lpoint1, const CPoint2D &lpoint2,
                                  const CPoint2D &point);

  bool Collinear(const CPoint2D &point1, const CPoint2D &point2, const CPoint2D &point3);

  bool IntersectsProperly(const CPoint2D &l1point1, const CPoint2D &l1point2,
                          const CPoint2D &l2point1, const CPoint2D &l2point2);

  CPoint2D RotatePoint(const CPoint2D &point, double angle, const CPoint2D &o);

  bool IsPerpendicular(double x21, double y21, double x32, double y32);

  bool ThreePointCircle1(double x1, double y1, double x2, double y2,
                         double x3, double y3, double *xc, double *yc, double *r);

  bool CircleCircleIntersect(double x1, double y1, double r1, double x2, double y2, double r2,
                             double *xi1, double *yi1, double *xi2, double *yi2);
}

#include <CLine2D.h>

namespace CMathGeom2D {
  bool PolygonIsConvex(const std::vector<CLine2D> &lines);
  bool PolygonIsConvex(const std::list<CLine2D> &lines);

  bool clipLine(const CLine2D &line, const CBBox2D &bbox, CLine2D &cline);

  bool PointLineDistance(const CPoint2D &point, const CLine2D &line, double *dist);

  bool PolygonLineIntersect(const std::vector<CPoint2D> &points, const CLine2D &line,
                            std::vector<CPoint2D> &ipoints);

  bool PolygonLineIntersect(const std::vector<CPoint2D> &points, const CLine2D &line,
                            std::map<uint,CPoint2D> &ipoints);

  bool LinesAreCoincident(const CLine2D &line1, const CLine2D &line2);

  PointPosition PointLinePosition(const CLine2D &line, const CPoint2D &point);

  bool PointLineLeft   (const CLine2D &line, const CPoint2D &point);
  bool PointLineRight  (const CLine2D &line, const CPoint2D &point);
  bool PointLineOn     (const CLine2D &line, const CPoint2D &point);
  bool PointLineLeftOn (const CLine2D &line, const CPoint2D &point);
  bool PointLineRightOn(const CLine2D &line, const CPoint2D &point);
}

#include <C3Bezier2D.h>
#include <C2Bezier2D.h>

namespace CMathGeom2D {
  void ArcToBeziers(double x, double y, double xr, double yr, double angle1, double angle2,
                    std::vector<C3Bezier2D> &beziers);
  void ArcNToBeziers(double x, double y, double xr, double yr, double angle1, double angle2,
                     std::vector<C3Bezier2D> &beziers);

  void BezierToLines(const C3Bezier2D &bezier, std::vector<CPoint2D> &points, double tol=-1.0);
  void BezierToLines(const C2Bezier2D &bezier, std::vector<CPoint2D> &points, double tol=-1.0);

  bool BezierLineIntersect(const C2Bezier2D &bezier, const CLine2D &line,
                           std::vector<CPoint2D> &ipoints, double tol);
  bool BezierLineIntersect(const C2Bezier2D &bezier, const CLine2D &line,
                           std::vector<CPoint2DParam2> &ibpoints, double tol);
  bool BezierLineIntersect(const C2Bezier2D &bezier, const CLine2D &line, double t1, double t2,
                           std::vector<CPoint2DParam2> &ibpoints, double tol);

  bool BezierBezierIntersect(const C2Bezier2D &bezier1, const C2Bezier2D &bezier2,
                             std::vector<CPoint2D> &ipoints, double tol);
  bool BezierBezierIntersect(const C2Bezier2D &bezier1, const C2Bezier2D &bezier2,
                             std::vector<CPoint2DParam2> &ipoints, double tol);
  bool BezierBezierIntersect(const C2Bezier2D &bezier1, double t11, double t12,
                             const C2Bezier2D &bezier2, double t21, double t22,
                             std::vector<CPoint2DParam2> &ipoints, double tol);

  bool BezierLineIntersect(const C3Bezier2D &bezier, const CLine2D &line,
                           std::vector<CPoint2D> &ipoints, double tol);
  bool BezierLineIntersect(const C3Bezier2D &bezier, const CLine2D &line,
                           std::vector<CPoint2DParam2> &ibpoints, double tol);
  bool BezierLineIntersect(const C3Bezier2D &bezier, const CLine2D &line, double t1, double t2,
                           std::vector<CPoint2DParam2> &ibpoints, double tol);

  bool BezierBezierIntersect(const C3Bezier2D &bezier1, const C3Bezier2D &bezier2,
                             std::vector<CPoint2D> &ipoints, double tol);
  bool BezierBezierIntersect(const C3Bezier2D &bezier1, const C3Bezier2D &bezier2,
                             std::vector<CPoint2DParam2> &ipoints, double tol);
  bool BezierBezierIntersect(const C3Bezier2D &bezier1, double t11, double t12,
                             const C3Bezier2D &bezier2, double t21, double t22,
                             std::vector<CPoint2DParam2> &ipoints, double tol);

  void AddUniquePoint(std::vector<CPoint2DParam2> &points, const CPoint2DParam2 &p);

  bool LineToPolygon(const CPoint2D &p1, const CPoint2D &p2,
                     std::vector<CPoint2D> &points, double line_width=1.0);

  bool LineJoinToPolygon(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3,
                         CLineJoinType lineJoin, std::vector<CPoint2D> &points,
                         double line_width=1.0);

  bool MitreLineJoinToPolygon(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3,
                              std::vector<CPoint2D> &points, double line_width=1.0,
                              double mitre_limit=1.0);
  bool RoundLineJoinToPolygon(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3,
                              std::vector<CPoint2D> &points);
  bool BevelLineJoinToPolygon(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3,
                              std::vector<CPoint2D> &points, double line_width=1.0);

  void AddWidthToPoint(const CPoint2D &p, double g, double line_width, CPoint2D &p1, CPoint2D &p2);
}

struct COffsetBezierData {
  std::vector<C3Bezier2D> beziers;
  std::vector<CPoint2D>   points;
};

namespace CMathGeom2D {
  C3Bezier2D offsetBezier(const C3Bezier2D &b, double offset, uint num_samples=32);

  void offsetBezier(const C3Bezier2D &b, double offset, COffsetBezierData &data,
                    uint num_samples=32, int num_curves=1);
}

#endif
