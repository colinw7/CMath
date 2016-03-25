#ifndef CMathGeom2D_H
#define CMathGeom2D_H

#include <CMathGen.h>
#include <vector>

#define EPSILON_E6 1E-6

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

  bool CircleLineIntersect(double xc, double yc, double r,
                           double x1, double y1, double x2, double y2,
                           double *xi1, double *yi1,
                           double *xi2, double *yi2, uint *num_i);
}

namespace CMathGeom2D {
  bool IntersectPolygons(const double *x1, const double *y1, uint n1,
                         const double *x2, const double *y2, uint n2,
                         double **xi, double **yi, uint *ni);

  bool PolygonLineIntersect(const double *x, const double *y, uint nxy,
                            double x1, double y1, double x2, double y2,
                            double *xi, double *yi, uint *num_i);

  bool PointLineLeft   (const CPoint2D &lpoint1, const CPoint2D &lpoint2, const CPoint2D &point);
  bool PointLineRight  (const CPoint2D &lpoint1, const CPoint2D &lpoint2, const CPoint2D &point);
  bool PointLineOn     (const CPoint2D &lpoint1, const CPoint2D &lpoint2, const CPoint2D &point);
  bool PointLineLeftOn (const CPoint2D &lpoint1, const CPoint2D &lpoint2, const CPoint2D &point);
  bool PointLineRightOn(const CPoint2D &lpoint1, const CPoint2D &lpoint2, const CPoint2D &point);

  double TriangleArea2(const CPoint2D &point1, const CPoint2D &point2, const CPoint2D &point3);

  void PolygonCentroid(const double *x, const double *y, int num_xy, double *xc, double *yc);
  void PolygonCentroid(const std::vector<CPoint2D> &points, CPoint2D &c);

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
#include <CPolygonOrientation.h>

namespace CMathGeom2D {
  inline CPolygonOrientation PolygonOrientation(double x1, double y1, double x2, double y2,
                                                double x3, double y3) {
    double dx1 = x2 - x1;
    double dy1 = y2 - y1;

    double dx2 = x3 - x2;
    double dy2 = y3 - y2;

    return (CPolygonOrientation) CMathGen::sign(dx1*dy2 - dy1*dx2);
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

  bool PointInsideEvenOdd(const CPoint2D &point, const std::vector<CPoint2D> &points);
  bool PointInsideEvenOdd(const CPoint2D &point, const CPoint2D *points, uint num_points);

  double PolygonArea(const double *x, const double *y, uint num_xy);

  bool Intersects(const CPoint2D &l1point1, const CPoint2D &l1point2,
                  const CPoint2D &l2point1, const CPoint2D &l2point2);
}

#include <CLine2D.h>

namespace CMathGeom2D {
  bool PolygonIsConvex(const std::vector<CLine2D> &lines);

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

#endif
