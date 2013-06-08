#ifndef CMathGeom2D_H
#define CMathGeom2D_H

#define EPSILON_E6 1E-6

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

  CPolygonOrientation PolygonOrientation(double x1, double y1, double x2, double y2,
                                         double x3, double y3);
}

#endif
