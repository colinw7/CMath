#ifndef CBresenhamShapes_H
#define CBresenhamShapes_H

#include <cmath>

/********************************************************************
*                                                                   *
*                    Curve Rasterizing Algorithm                    *
*                                                                   *
********************************************************************/

/**
 * @author Zingl Alois
 * @date 16.5.2012
 * @version 1.0
*/

class CBresenhamShapes {
 public:
  virtual ~CBresenhamShapes() { }

 public:
  void plotLine(int x0, int y0, int x1, int y1) {
    int dx =  abs(x1 - x0), sx = (x0 < x1 ? 1 : -1);
    int dy = -abs(y1 - y0), sy = (y0 < y1 ? 1 : -1);

    int err = dx + dy, e2; /* error value e_xy */

    for (;;) { /* loop */
      drawPoint(x0, y0);

      e2 = 2*err;

      if (e2 >= dy) { /* e_xy + e_x > 0 */
        if (x0 == x1) break;
        err += dy; x0 += sx;
      }
      if (e2 <= dx) { /* e_xy + e_y < 0 */
        if (y0 == y1) break;
        err += dx; y0 += sy;
      }
    }
  }

 public:
  void plotEllipse(int xm, int ym, int a, int b) {
    int x = -a, y = 0;   /* II. quadrant from bottom left to top right */

    long e2 = (long)b*b, err = (long)x*(2*e2 + x) + e2; /* error of 1.step */

    do {
      drawPoint(xm - x, ym + y); /*   I. Quadrant */
      drawPoint(xm + x, ym + y); /*  II. Quadrant */
      drawPoint(xm + x, ym - y); /* III. Quadrant */
      drawPoint(xm - x, ym - y); /*  IV. Quadrant */

      e2 = 2*err;

      if (e2 >= (x*2 + 1)*(long)b*b) /* e_xy + e_x > 0 */
        err += (++x*2 + 1)*(long)b*b;

      if (e2 <= (y*2 + 1)*(long)a*a) /* e_xy + e_y < 0 */
        err += (++y*2 + 1)*(long)a*a;
    } while (x <= 0);

    while (y++ < b) { /* too early stop of flat ellipses a=1, */
      drawPoint(xm, ym + y); /* -> finish tip of ellipse */
      drawPoint(xm, ym - y);
    }
  }

 public:
  void plotOptimizedEllipse(int xm, int ym, int a, int b) {
    long x = -a, y = 0; /* II. quadrant from bottom left to top right */
    long e2 = b, dx = (1 + 2*x)*e2*e2; /* error increment  */
    long dy = x*x, err = dx + dy;      /* error of 1.step */

    do {
      drawPoint(xm - x, ym + y); /*   I. Quadrant */
      drawPoint(xm + x, ym + y); /*  II. Quadrant */
      drawPoint(xm + x, ym - y); /* III. Quadrant */
      drawPoint(xm - x, ym - y); /*  IV. Quadrant */

      e2 = 2*err;
      if (e2 >= dx) { x++; err += dx += 2*(long)b*b; } /* x step */
      if (e2 <= dy) { y++; err += dy += 2*(long)a*a; } /* y step */
    } while (x <= 0);

    while (y++ < b) { /* too early stop for flat ellipses with a=1, */
      drawPoint(xm, ym + y); /* -> finish tip of ellipse */
      drawPoint(xm, ym - y);
    }
  }

 public:
  void plotCircle(int xm, int ym, int r) {
    int x = -r, y = 0, err = 2 - 2*r; /* bottom left to top right */
    do {
      drawPoint(xm - x, ym + y); /*   I. Quadrant +x +y */
      drawPoint(xm - y, ym - x); /*  II. Quadrant -x +y */
      drawPoint(xm + x, ym - y); /* III. Quadrant -x -y */
      drawPoint(xm + y, ym + x); /*  IV. Quadrant +x -y */

      r = err;
      if (r <= y) err += ++y*2 + 1; /* e_xy + e_y < 0 */
      if (r > x || err > y)         /* e_xy + e_x > 0 or no 2nd y - step */
        err += ++x*2 + 1;           /* -> x - step now */
    } while (x < 0);
  }

 private:
  /* rectangular parameter enclosing the ellipse */
  void plotEllipseRect(int x0, int y0, int x1, int y1) {
    long a = abs(x1 - x0), b = abs(y1 - y0), b1 = b&1; /* diameter */
    double dx = 4*(1.0 - a)*b*b, dy = 4*(b1 + 1)*a*a;  /* error increment */
    double err = dx + dy + b1*a*a, e2;                 /* error of 1.step */

    if (x0 > x1) { x0 = x1; x1 += a; } /* if called with swapped points */
    if (y0 > y1) y0 = y1;              /* .. exchange them */
    y0 += (b + 1)/2; y1 = y0 - b1;     /* starting pixel */
    a = 8*a*a; b1 = 8*b*b;

    do {
      drawPoint(x1, y0); /*   I. Quadrant */
      drawPoint(x0, y0); /*  II. Quadrant */
      drawPoint(x0, y1); /* III. Quadrant */
      drawPoint(x1, y1); /*  IV. Quadrant */

      e2 = 2*err;
      if (e2 <= dy) { y0++; y1--; err += dy += a; }                /* y step */
      if (e2 >= dx || 2*err > dy) { x0++; x1--; err += dx += b1; } /* x */
    } while (x0 <= x1);

    while (y0 - y1 <= b) {     /* too early stop of flat ellipses a=1 */
      drawPoint(x0 - 1, y0);   /* -> finish tip of ellipse */
      drawPoint(x1 + 1, y0++);
      drawPoint(x0 - 1, y1);
      drawPoint(x1 + 1, y1--);
    }
  }

 public:
  void plotBasicQuadBezier(int x0, int y0, int x1, int y1, int x2, int y2) {
    int sx = x0<x2 ? 1 : -1, sy = y0<y2 ? 1 : -1; /* step direction */

    double x = x0 - 2*x1 + x2, y = y0 - 2*y1 + y2, xy = 2*x*y*sx*sy;

    double cur = sx*sy*(x*(y2 - y0) - y*(x2 - x0))/2; /* curvature */

    /* compute error increments of P0 */
    double dx = (1 - 2*abs(x0 - x1))*y*y + abs(y0 - y1)*xy - 2*cur*abs(y0 - y2);
    double dy = (1 - 2*abs(y0 - y1))*x*x + abs(x0 - x1)*xy + 2*cur*abs(x0 - x2);

    /* compute error increments of P2 */
    double ex = (1 - 2*abs(x2 - x1))*y*y + abs(y2 - y1)*xy + 2*cur*abs(y0 - y2);
    double ey = (1 - 2*abs(y2 - y1))*x*x + abs(x2 - x1)*xy - 2*cur*abs(x0 - x2);

    /* sign of gradient must not change */
    assert((x0 - x1)*(x2 - x1) <= 0 && (y0 - y1)*(y2 - y1) <= 0);

    if (cur==0) { plotLine(x0, y0, x2, y2); return; } /* straight line */

    x *= 2*x; y *= 2*y;

    if (cur < 0) { /* negated curvature */
      x = -x; dx = -dx; ex = -ex; xy = -xy;
      y = -y; dy = -dy; ey = -ey;
    }

    /* algorithm fails for almost straight line, check error values */
    if (dx >= -y || dy <= -x || ex <= -y || ey >= -x) {
      x1 = (x0 + 4*x1 + x2)/6; y1 = (y0 + 4*y1 + y2)/6; /* approximation */

      plotLine(x0, y0, x1, y1);
      plotLine(x1, y1, x2, y2);

      return;
    }

    dx -= xy; ex = dx + dy; dy -= xy; /* error of 1.step */

    for (;;) { /* plot curve */
      drawPoint(x0, y0);

      ey = 2*ex - dy; /* save value for test of y step */

      if (2*ex >= dx) { /* x step */
        if (x0 == x2) break;
        x0 += sx; dy -= xy; ex += dx += y;
      }
      if (ey <= 0) {    /* y step */
        if (y0 == y2) break;
        y0 += sy; dx -= xy; ex += dy += x;
      }
    }
  }

 public:
  void plotFineQuadBezier(int x0, int y0, int x1, int y1, int x2, int y2) {
    int sx = x0<x2 ? 1 : -1, sy = y0<y2 ? 1 : -1; /* step direction */

    long f = 1, fx = x0 - 2*x1 + x2, fy = y0 - 2*y1 + y2;

    long long x = 2*fx*fx, y = 2*fy*fy, xy = 2*fx*fy*sx*sy;

    long long cur = sx*sy*(fx*(y2 - y0) - fy*(x2 - x0)); /* curvature */

    /* compute error increments of P0 */
    long long dx = abs(y0 - y1)*xy - abs(x0 - x1)*y - cur*abs(y0 - y2);
    long long dy = abs(x0 - x1)*xy - abs(y0 - y1)*x + cur*abs(x0 - x2);

    /* compute error increments of P2 */
    long long ex = abs(y2 - y1)*xy - abs(x2 - x1)*y + cur*abs(y0 - y2);
    long long ey = abs(x2 - x1)*xy - abs(y2 - y1)*x - cur*abs(x0 - x2);

    /* sign of gradient must not change */
    assert((x0 - x1)*(x2 - x1) <= 0 && (y0 - y1)*(y2 - y1) <= 0);

    if (cur == 0) { plotLine(x0, y0, x2, y2); return; } /* straight line */

    /* compute required minimum resolution factor */
    if (dx == 0 || dy == 0 || ex == 0 || ey == 0)
      f = abs(xy/cur)/2 + 1; /* division by zero: use curvature */
    else {
      fx = 2*y/dx; if (fx > f) f = fx; /* increase resolution */
      fx = 2*x/dy; if (fx > f) f = fx;
      fx = 2*y/ex; if (fx > f) f = fx;
      fx = 2*x/ey; if (fx > f) f = fx;
    } /* negated curvature? */

    if (cur < 0) { x = -x; y = -y; dx = -dx; dy = -dy; xy = -xy; }

    dx = f*dx + y/2 - xy; dy = f*dy + x/2 - xy; ex = dx + dy + xy; /* error 1.step */

    for (fx = fy = f; ; ) { /* plot curve */
      drawPoint(x0, y0);

      if (x0 == x2 && y0 == y2) break;

      do { /* move f sub-pixel */
        ey = 2*ex - dy; /* save value for test of y step */

        if (2*ex >= dx) { fx--; dy -= xy; ex += dx += y; } /* x step */
        if (  ey <= 0 ) { fy--; dx -= xy; ex += dy += x; } /* y step */
      } while (fx > 0 && fy > 0);                          /* pixel complete? */

      if (2*fx <= f) { x0 += sx; fx += f; } /* sufficient sub-steps.. */
      if (2*fy <= f) { y0 += sy; fy += f; } /* .. for a pixel? */
    }
  }

 private:
  /* plot a limited quadratic Bezier segment */
  void plotQuadBezierSeg(int x0, int y0, int x1, int y1, int x2, int y2) {
    int sx = x2 - x1, sy = y2 - y1;
    long xx = x0 - x1, yy = y0 - y1, xy;     /* relative values for checks */
    double dx, dy, err, cur = xx*sy - yy*sx; /* curvature */

    assert(xx*sx <= 0 && yy*sy <= 0); /* sign of gradient must not change */

    if (sx*(long)sx + sy*(long)sy > xx*xx + yy*yy) { /* begin with longer part */
      x2 = x0; x0 = sx + x1; y2 = y0; y0 = sy + y1; cur = -cur; /* swap P0 P2 */
    }
    if (cur != 0) { /* no straight line */
      xx += sx; xx *= sx = x0 < x2 ? 1 : -1; /* x step direction */
      yy += sy; yy *= sy = y0 < y2 ? 1 : -1; /* y step direction */
      xy = 2*xx*yy; xx *= xx; yy *= yy;      /* differences 2nd degree */

      if (cur*sx*sy < 0) { /* negated curvature? */
        xx = -xx; yy = -yy; xy = -xy; cur = -cur;
      }

      dx = 4.0*sy*cur*(x1 - x0) + xx - xy; /* differences 1st degree */
      dy = 4.0*sx*cur*(y0 - y1) + yy - xy;

      xx += xx; yy += yy; err = dx + dy + xy; /* error 1st step */

      do {
        drawPoint(x0, y0); /* plot curve */

        if (x0 == x2 && y0 == y2) return; /* last pixel -> curve finished */

        y1 = 2*err < dx; /* save value for test of y step */

        if (2*err > dy) { x0 += sx; dx -= xy; err += dy += yy; } /* x step */
        if (    y1    ) { y0 += sy; dy -= xy; err += dx += xx; } /* y step */
      } while (dy < 0 && dx > 0); /* gradient negates -> algorithm fails */
    }
    plotLine(x0, y0, x2, y2); /* plot remaining part to end */
  }

 public:
  /* plot any quadratic Bezier curve */
  void plotQuadBezier(int x0, int y0, int x1, int y1, int x2, int y2) {
    int    x = x0 - x1, y = y0 - y1;
    double t = x0 - 2*x1 + x2, r;

    if ((long)x*(x2 - x1) > 0) { /* horizontal cut at P4? */
      if ((long)y*(y2 - y1) > 0) /* vertical cut at P6 too? */
        if (fabs((y0 - 2*y1 + y2)/t*x) > abs(y)) { /* which first? */
          x0 = x2; x2 = x + x1; y0 = y2; y2 = y + y1; /* swap points */
        } /* now horizontal cut at P4 comes first */

      t = (x0 - x1)/t;
      r = (1 - t)*((1 - t)*y0 + 2.0*t*y1) + t*t*y2; /* By(t=P4) */
      t = (x0*x2 - x1*x1)*t/(x0 - x1);              /* gradient dP4/dx=0 */
      x = std::floor(t + 0.5); y = std::floor(r + 0.5);
      r = (y1 - y0)*(t - x0)/(x1 - x0) + y0;        /* intersect P3 | P0 P1 */

      plotQuadBezierSeg(x0, y0, x, std::floor(r + 0.5), x, y);

      r = (y1 - y2)*(t - x2)/(x1 - x2) + y2; /* intersect P4 | P1 P2 */

      x0 = x1 = x; y0 = y; y1 = std::floor(r + 0.5); /* P0 = P4, P1 = P8 */
    }
    if ((long)(y0 - y1)*(y2 - y1) > 0) {            /* vertical cut at P6? */
      t = y0 - 2*y1 + y2; t = (y0 - y1)/t;
      r = (1 - t)*((1 - t)*x0 + 2.0*t*x1) + t*t*x2; /* Bx(t=P6) */
      t = (y0*y2 - y1*y1)*t/(y0 - y1);              /* gradient dP6/dy=0 */
      x = std::floor(r + 0.5); y = std::floor(t + 0.5);
      r = (x1 - x0)*(t - y0)/(y1 - y0) + x0;        /* intersect P6 | P0 P1 */

      plotQuadBezierSeg(x0, y0, std::floor(r + 0.5), y, x, y);

      r = (x1 - x2)*(t - y2)/(y1 - y2) + x2;         /* intersect P7 | P1 P2 */

      x0 = x; x1 = std::floor(r + 0.5); y0 = y1 = y; /* P0 = P6, P1 = P7 */
    }
    plotQuadBezierSeg(x0, y0, x1, y1, x2, y2);       /* remaining part */
  }

 private:
  /* plot a limited rational Bezier segment, squared weight */
  void plotQuadRationalBezierSeg(int x0, int y0, int x1, int y1, int x2, int y2, float w) {
    int sx = x2 - x1, sy = y2 - y1; /* relative values for checks */

    double dx = x0 - x2, dy = y0 - y2, xx = x0 - x1, yy = y0 - y1;
    double xy = xx*sy + yy*sx, cur = xx*sy - yy*sx, err; /* curvature */

    assert(xx*sx <= 0.0 && yy*sy <= 0.0); /* sign of gradient must not change */

    if (cur != 0.0 && w > 0.0) { /* no straight line */
      if (sx*(long)sx + sy*(long)sy > xx*xx + yy*yy) {    /* begin with longer part */
        x2 = x0; x0 -= dx; y2 = y0; y0 -= dy; cur = -cur; /* swap P0 P2 */
      }

      xx = 2.0*(4.0*w*sx*xx + dx*dx); /* differences 2nd degree */
      yy = 2.0*(4.0*w*sy*yy + dy*dy);

      sx = x0 < x2 ? 1 : -1; /* x step direction */
      sy = y0 < y2 ? 1 : -1; /* y step direction */

      xy = -2.0*sx*sy*(2.0*w*xy + dx*dy);

      if (cur*sx*sy < 0.0) { /* negated curvature? */
        xx = -xx; yy = -yy; xy = -xy; cur = -cur;
      }

      dx = 4.0*w*(x1 - x0)*sy*cur + xx/2.0 + xy; /* differences 1st degree */
      dy = 4.0*w*(y0 - y1)*sx*cur + yy/2.0 + xy;

      if (w < 0.5 && (dy > xy || dx < xy)) { /* flat ellipse, algorithm fails */
        cur = (w + 1.0)/2.0; w = sqrt(w); xy = 1.0/(w + 1.0);

        sx = std::floor((x0 + 2.0*w*x1 + x2)*xy/2.0 + 0.5); /* subdivide curve in half */
        sy = std::floor((y0 + 2.0*w*y1 + y2)*xy/2.0 + 0.5);

        dx = std::floor((w*x1 + x0)*xy + 0.5); dy = std::floor((y1*w + y0)*xy + 0.5);

        plotQuadRationalBezierSeg(x0, y0, dx, dy, sx, sy, cur); /* plot separately */

        dx = std::floor((w*x1 + x2)*xy + 0.5); dy = std::floor((y1*w + y2)*xy + 0.5);

        plotQuadRationalBezierSeg(sx, sy, dx, dy, x2, y2, cur);

        return;
      }

      err = dx + dy - xy; /* error 1.step */

      do {
        drawPoint(x0, y0); /* plot curve */

        if (x0 == x2 && y0 == y2) return; /* last pixel -> curve finished */

        x1 = 2*err > dy; y1 = 2*(err + yy) < -dy; /* save value for test of x step */

        if (2*err < dx || y1) { y0 += sy; dy += xy; err += dx += xx; } /* y step */
        if (2*err > dx || x1) { x0 += sx; dx += xy; err += dy += yy; } /* x step */

      } while (dy <= xy && dx >= xy); /* gradient negates -> algorithm fails */
    }

    plotLine(x0, y0, x2, y2); /* plot remaining needle to end */
  }

 private:
  /* plot any quadratic rational Bezier curve */
  void plotQuadRationalBezier(int x0, int y0, int x1, int y1, int x2, int y2, float w) {
    int x = x0 - 2*x1 + x2, y = y0 - 2*y1 + y2;
    double xx = x0 - x1, yy = y0 - y1, ww, t, q;

    assert(w >= 0.0);

    if (xx*(x2 - x1) > 0) { /* horizontal cut at P4? */
      if (yy*(y2 - y1) > 0) /* vertical cut at P6 too? */
        if (fabs(xx*y) > fabs(yy*x)) { /* which first? */
          x0 = x2; x2 = xx + x1; y0 = y2; y2 = yy + y1; /* swap points */
        } /* now horizontal cut at P4 comes first */
      if (x0 == x2 || w == 1.0) t = (x0 - x1)/(double)x;
      else { /* non-rational or rational case */
        q = sqrt(4.0*w*w*(x0 - x1)*(x2 - x1)+(x2 - x0)*(long)(x2 - x0));
        if (x1 < x0) q = -q;
        t = (2.0*w*(x0 - x1) - x0 + x2 + q)/(2.0*(1.0 - w)*(x2 - x0)); /* t at P4 */
      }

      q = 1.0/(2.0*t*(1.0 - t)*(w - 1.0) + 1.0); /* sub-divide at t */

      xx = (t*t*(x0 - 2.0*w*x1 + x2) + 2.0*t*(w*x1 - x0) + x0)*q; /* = P4 */
      yy = (t*t*(y0 - 2.0*w*y1 + y2) + 2.0*t*(w*y1 - y0) + y0)*q;
      ww = t*(w - 1.0) + 1.0; ww *= ww*q; /* squared weight P3 */

      w = ((1.0 - t)*(w - 1.0) + 1.0)*sqrt(q);            /* weight P8 */
      x = std::floor(xx + 0.5); y = std::floor(yy + 0.5); /* P4 */
      yy = (xx - x0)*(y1 - y0)/(x1 - x0) + y0;            /* intersect P3 | P0 P1 */

      plotQuadRationalBezierSeg(x0, y0, x, std::floor(yy + 0.5), x, y, ww);

      yy = (xx - x2)*(y1 - y2)/(x1 - x2) + y2;        /* intersect P4 | P1 P2 */
      y1 = std::floor(yy + 0.5); x0 = x1 = x; y0 = y; /* P0 = P4, P1 = P8 */
    }
    if ((y0 - y1)*(long)(y2 - y1) > 0) { /* vertical cut at P6? */
      if (y0 == y2 || w == 1.0) t = (y0 - y1)/(y0 - 2.0*y1 + y2);
      else { /* non-rational or rational case */
        q = sqrt(4.0*w*w*(y0 - y1)*(y2 - y1)+(y2 - y0)*(long)(y2 - y0));
        if (y1 < y0) q = -q;
        t = (2.0*w*(y0 - y1) - y0 + y2 + q)/(2.0*(1.0 - w)*(y2 - y0)); /* t at P6 */
      }

      q = 1.0/(2.0*t*(1.0 - t)*(w - 1.0) + 1.0); /* sub-divide at t */

      xx = (t*t*(x0 - 2.0*w*x1 + x2) + 2.0*t*(w*x1 - x0) + x0)*q; /* = P6 */
      yy = (t*t*(y0 - 2.0*w*y1 + y2) + 2.0*t*(w*y1 - y0) + y0)*q;
      ww = t*(w - 1.0) + 1.0; ww *= ww*q; /* squared weight P5 */

      w = ((1.0 - t)*(w - 1.0) + 1.0)*sqrt(q);            /* weight P7 */
      x = std::floor(xx + 0.5); y = std::floor(yy + 0.5); /* P6 */
      xx = (x1 - x0)*(yy - y0)/(y1 - y0) + x0;            /* intersect P6 | P0 P1 */

      plotQuadRationalBezierSeg(x0, y0, std::floor(xx + 0.5), y, x, y, ww);

      xx = (x1 - x2)*(yy - y2)/(y1 - y2) + x2;        /* intersect P7 | P1 P2 */
      x1 = std::floor(xx + 0.5); x0 = x; y0 = y1 = y; /* P0 = P6, P1 = P7 */
    }

    plotQuadRationalBezierSeg(x0, y0, x1, y1, x2, y2, w*w); /* remaining */
  }

 public:
  /* plot ellipse rotated by angle (radian) */
  void plotRotatedEllipse(int x, int y, int a, int b, float angle) {
    float xd = (long)a*a, yd = (long)b*b;
    float s  = sin(angle), zd = (xd - yd)*s; /* ellipse rotation */

    xd = sqrt(xd - zd*s), yd = sqrt(yd + zd*s); /* surrounding rectangle */

    a = xd + 0.5; b = yd + 0.5; zd = zd*a*b/(xd*yd); /* scale to integer */

    plotRotatedEllipseRect(x - a, y - b, x + a, y + b, (long)(4*zd*cos(angle)));
  }

 private:
  /* rectangle enclosing the ellipse, integer rotation angle */
  void plotRotatedEllipseRect(int x0, int y0, int x1, int y1, long zd) {
    int xd = x1 - x0, yd = y1 - y0;

    float w = xd*(long)yd;

    if (zd == 0) return plotEllipseRect(x0, y0, x1, y1); /* looks nicer */

    if (w != 0.0) w = (w - zd)/(w + w); /* squared weight of P1 */
    assert(w <= 1.0 && w >= 0.0);       /* limit angle to |zd|<=xd*yd */

    xd = std::floor(xd*w + 0.5); yd = std::floor(yd*w + 0.5); /* snap xe, ye to int */

    plotQuadRationalBezierSeg(x0, y0 + yd, x0, y0, x0 + xd, y0, 1.0 - w);
    plotQuadRationalBezierSeg(x0, y0 + yd, x0, y1, x1 - xd, y1, w);
    plotQuadRationalBezierSeg(x1, y1 - yd, x1, y1, x1 - xd, y1, 1.0 - w);
    plotQuadRationalBezierSeg(x1, y1 - yd, x1, y0, x0 + xd, y0, w);
  }

 private:
  /* plot limited cubic Bezier segment */
  void plotCubicBezierSeg(int x0, int y0, float x1, float y1, float x2, float y2, int x3, int y3) {
    int f, fx, fy, leg = 1;
    int sx = x0 < x3 ? 1 : -1, sy = y0 < y3 ? 1 : -1; /* step direction */
    float xc = -fabs(x0 + x1 - x2 - x3), xa = xc - 4*sx*(x1 - x2), xb = sx*(x0 - x1 - x2 + x3);
    float yc = -fabs(y0 + y1 - y2 - y3), ya = yc - 4*sy*(y1 - y2), yb = sy*(y0 - y1 - y2 + y3);
    double ab, ac, bc, cb, xx, xy, yy, dx, dy, ex, *pxy, EP = 0.01;
    /* check for curve restrains */
    /* slope P0-P1 == P2-P3 and (P0-P3 == P1-P2 or no slope change) */
    assert((x1 - x0)*(x2 - x3) < EP && ((x3 - x0)*(x1 - x2) < EP || xb*xb < xa*xc + EP));
    assert((y1 - y0)*(y2 - y3) < EP && ((y3 - y0)*(y1 - y2) < EP || yb*yb < ya*yc + EP));

    if (xa == 0 && ya == 0) { /* quadratic Bezier */
      sx = std::floor((3*x1 - x0 + 1)/2); sy = std::floor((3*y1 - y0 + 1)/2); /* new midpoint */
      return plotQuadBezierSeg(x0, y0, sx, sy, x3, y3);
    }
    x1 = (x1 - x0)*(x1 - x0)+(y1 - y0)*(y1 - y0); /* line lengths */
    x2 = (x2 - x3)*(x2 - x3)+(y2 - y3)*(y2 - y3);
    do { /* loop over both ends */
      ab = xa*yb - xb*ya; ac = xa*yc - xc*ya; bc = xb*yc - xc*yb;
      ex = ab*(ab + ac - 3*bc) + ac*ac; /* P0 part of self-intersection loop? */

      f = ex > 0 ? 1 : sqrt(1 + 1024/x1); /* calculate resolution */

      ab *= f; ac *= f; bc *= f; ex *= f*f; /* increase resolution */

      xy = 9*(ab + ac + bc)/8; cb = 8*(xa - ya); /* init differences of 1st degree */

      dx = 27*(8*ab*(yb*yb - ya*yc) + ex*(ya + 2*yb + yc))/64 - ya*ya*(xy - ya);
      dy = 27*(8*ab*(xb*xb - xa*xc) - ex*(xa + 2*xb + xc))/64 - xa*xa*(xy + xa);
      /* init differences of 2nd degree */
      xx = 3*(3*ab*(3*yb*yb - ya*ya - 2*ya*yc) - ya*(3*ac*(ya + yb) + ya*cb))/4;
      yy = 3*(3*ab*(3*xb*xb - xa*xa - 2*xa*xc) - xa*(3*ac*(xa + xb) + xa*cb))/4;

      xy = xa*ya*(6*ab + 6*ac - 3*bc + cb); ac = ya*ya; cb = xa*xa;
      xy = 3*(xy + 9*f*(cb*yb*yc - xb*xc*ac) - 18*xb*yb*ab)/8;

      if (ex < 0) { /* negate values if inside self-intersection loop */
         dx = -dx; dy = -dy; xx = -xx; yy = -yy; xy = -xy; ac = -ac; cb = -cb;
      } /* init differences of 3rd degree */
      ab = 6*ya*ac; ac = -6*xa*ac; bc = 6*ya*cb; cb = -6*xa*cb;
      dx += xy; ex = dx + dy; dy += xy; /* error of 1st step */

      for (pxy = &xy, fx = fy = f; x0 != x3 && y0 != y3; ) {
        drawPoint(x0, y0); /* plot curve */
        do { /* move sub-steps of one pixel */
          if (dx > *pxy || dy < *pxy) goto exit; /* confusing values */
          y1 = 2*ex - dy;   /* save value for test of y step */
          if (2*ex >= dx) { /* x sub-step */
            fx--; ex += dx += xx; dy += xy += ac; yy += bc; xx += ab;
          }
          if (y1 <= 0) { /* y sub-step */
            fy--; ex += dy += yy; dx += xy += bc; xx += ac; yy += cb;
          }
        } while (fx > 0 && fy > 0); /* pixel complete? */
        if (2*fx <= f) { x0 += sx; fx += f; } /* x step */
        if (2*fy <= f) { y0 += sy; fy += f; } /* y step */
        if (pxy == &xy && dx < 0 && dy > 0) pxy = &EP; /* pixel ahead valid */
      }
  exit:
      xx = x0; x0 = x3; x3 = xx; sx = -sx; xb = -xb;          /* swap legs */
      yy = y0; y0 = y3; y3 = yy; sy = -sy; yb = -yb; x1 = x2;
    } while (leg--);                                          /* try other end */
    plotLine(x0, y0, x3, y3); /* remaining part in case of cusp or crunode */
  }

 public:
  /* plot any cubic Bezier curve */
  void plotCubicBezier(int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3) {
    int n = 0, i = 0;

    long xc = x0 + x1 - x2 - x3, xa = xc - 4*(x1 - x2);
    long xb = x0 - x1 - x2 + x3, xd = xb + 4*(x1 + x2);
    long yc = y0 + y1 - y2 - y3, ya = yc - 4*(y1 - y2);
    long yb = y0 - y1 - y2 + y3, yd = yb + 4*(y1 + y2);

    float fx0 = x0, fx1, fx2, fx3, fy0 = y0, fy1, fy2, fy3;

    double t1 = xb*xb - xa*xc, t2, t[5];

    /* sub-divide curve at gradient sign changes */
    if (xa == 0) { /* horizontal */
      if (abs(xc) < 2*abs(xb)) t[n++] = xc/(2.0*xb); /* one change */
    } else if (t1 > 0.0) {                           /* two changes */
      t2 = sqrt(t1);
      t1 = (xb - t2)/xa; if (fabs(t1) < 1.0) t[n++] = t1;
      t1 = (xb + t2)/xa; if (fabs(t1) < 1.0) t[n++] = t1;
    }

    t1 = yb*yb - ya*yc;

    if (ya == 0) { /* vertical */
      if (abs(yc) < 2*abs(yb)) t[n++] = yc/(2.0*yb); /* one change */
    } else if (t1 > 0.0) {                           /* two changes */
      t2 = sqrt(t1);
      t1 = (yb - t2)/ya; if (fabs(t1) < 1.0) t[n++] = t1;
      t1 = (yb + t2)/ya; if (fabs(t1) < 1.0) t[n++] = t1;
    }

    for (i = 1; i < n; i++) /* bubble sort of 4 points */
      if ((t1 = t[i - 1]) > t[i]) { t[i - 1] = t[i]; t[i] = t1; i = 0; }

    t1 = -1.0; t[n] = 1.0; /* begin / end point */

    for (i = 0; i <= n; i++) { /* plot each segment separately */
      t2 = t[i]; /* sub-divide at t[i - 1], t[i] */

      fx1 = (t1*(t1*xb - 2*xc) - t2*(t1*(t1*xa - 2*xb) + xc) + xd)/8 - fx0;
      fy1 = (t1*(t1*yb - 2*yc) - t2*(t1*(t1*ya - 2*yb) + yc) + yd)/8 - fy0;
      fx2 = (t2*(t2*xb - 2*xc) - t1*(t2*(t2*xa - 2*xb) + xc) + xd)/8 - fx0;
      fy2 = (t2*(t2*yb - 2*yc) - t1*(t2*(t2*ya - 2*yb) + yc) + yd)/8 - fy0;

      fx0 -= fx3 = (t2*(t2*(3*xb - t2*xa) - 3*xc) + xd)/8;
      fy0 -= fy3 = (t2*(t2*(3*yb - t2*ya) - 3*yc) + yd)/8;

      x3 = std::floor(fx3 + 0.5); y3 = std::floor(fy3 + 0.5); /* scale bounds to int */

      if (fx0 != 0.0) { fx1 *= fx0 = (x0 - x3)/fx0; fx2 *= fx0; }
      if (fy0 != 0.0) { fy1 *= fy0 = (y0 - y3)/fy0; fy2 *= fy0; }

      if (x0 != x3 || y0 != y3) /* segment t1 - t2 */
        plotCubicBezierSeg(x0, y0, x0 + fx1, y0 + fy1, x0 + fx2, y0 + fy2, x3, y3);

      x0 = x3; y0 = y3; fx0 = fx3; fy0 = fy3; t1 = t2;
    }
  }

 public:
  /* plot quadratic spline, destroys input arrays x, y */
  void plotQuadSpline(int n, int x[], int y[]) {
    #define M_MAX 6

    float mi = 1, m[M_MAX]; /* diagonal constants of matrix */

    int i, x0, y0, x1, y1, x2 = x[n], y2 = y[n];

    assert(n > 1); /* need at least 3 points P[0]..P[n] */

    x[1] = x0 = 8*x[1] - 2*x[0]; /* first row of matrix */
    y[1] = y0 = 8*y[1] - 2*y[0];

    for (i = 2; i < n; i++) { /* forward sweep */
      if (i - 2 < M_MAX) m[i - 2] = mi = 1.0/(6.0 - mi);

      x[i] = x0 = std::floor(8*x[i] - x0*mi + 0.5); /* store yi */
      y[i] = y0 = std::floor(8*y[i] - y0*mi + 0.5);
    }

    x1 = std::floor((x0 - 2*x2)/(5.0 - mi) + 0.5); /* correction last row */
    y1 = std::floor((y0 - 2*y2)/(5.0 - mi) + 0.5);

    for (i = n - 2; i > 0; i--) { /* back substitution */
      if (i <= M_MAX) mi = m[i - 1];

      x0 = std::floor((x[i] - x1)*mi + 0.5); /* next corner */
      y0 = std::floor((y[i] - y1)*mi + 0.5);

      plotQuadBezier((x0 + x1)/2, (y0 + y1)/2, x1, y1, x2, y2);

      x2 = (x0 + x1)/2; x1 = x0;
      y2 = (y0 + y1)/2; y1 = y0;
    }

    plotQuadBezier(x[0], y[0], x1, y1, x2, y2);
  }

 public:
  /* plot cubic spline, destroys input arrays x, y */
  void plotCubicSpline(int n, int x[], int y[]) {
    #define M_MAX 6

    float mi = 0.25, m[M_MAX]; /* diagonal constants of matrix */

    int x3 = x[n - 1], y3 = y[n - 1], x4 = x[n], y4 = y[n];

    int i, x0, y0, x1, y1, x2, y2;

    assert(n > 2); /* need at least 4 points P[0]..P[n] */

    x[1] = x0 = 12*x[1] - 3*x[0]; /* first row of matrix */
    y[1] = y0 = 12*y[1] - 3*y[0];

    for (i = 2; i < n; i++) { /* foreward sweep */
       if (i - 2 < M_MAX) m[i - 2] = mi = 0.25/(2.0 - mi);

       x[i] = x0 = std::floor(12*x[i] - 2*x0*mi + 0.5);
       y[i] = y0 = std::floor(12*y[i] - 2*y0*mi + 0.5);
    }

    x2 = std::floor((x0 - 3*x4)/(7 - 4*mi) + 0.5); /* correct last row */
    y2 = std::floor((y0 - 3*y4)/(7 - 4*mi) + 0.5);

    plotCubicBezier(x3, y3, (x2 + x4)/2, (y2 + y4)/2, x4, y4, x4, y4);

    if (n - 3 < M_MAX) mi = m[n - 3];

    x1 = std::floor((x[n - 2] - 2*x2)*mi + 0.5);
    y1 = std::floor((y[n - 2] - 2*y2)*mi + 0.5);

    for (i = n - 3; i > 0; i--) { /* back substitution */
       if (i <= M_MAX) mi = m[i - 1];

       x0 = std::floor((x[i] - 2*x1)*mi + 0.5);
       y0 = std::floor((y[i] - 2*y1)*mi + 0.5);
       x4 = std::floor((x0 + 4*x1 + x2 + 3)/6.0); /* reconstruct P[i] */
       y4 = std::floor((y0 + 4*y1 + y2 + 3)/6.0);

       plotCubicBezier(x4, y4,
                       std::floor((2*x1 + x2)/3 + 0.5), std::floor((2*y1 + y2)/3 + 0.5),
                       std::floor((x1 + 2*x2)/3 + 0.5), std::floor((y1 + 2*y2)/3 + 0.5),
                       x3, y3);

       x3 = x4; y3 = y4; x2 = x1; y2 = y1; x1 = x0; y1 = y0;
    }

    x0 = x[0]; x4 = std::floor((3*x0 + 7*x1 + 2*x2 + 6)/12.0); /* reconstruct P[1] */
    y0 = y[0]; y4 = std::floor((3*y0 + 7*y1 + 2*y2 + 6)/12.0);

    plotCubicBezier(x4, y4, std::floor((2*x1 + x2)/3 + 0.5), std::floor((2*y1 + y2)/3 + 0.5),
                    std::floor((x1 + 2*x2)/3 + 0.5), std::floor((y1 + 2*y2)/3 + 0.5), x3, y3);

    plotCubicBezier(x0, y0, x0, y0, (x0 + x1)/2, (y0 + y1)/2, x4, y4);
  }

 protected:
  virtual void drawPoint(int x, int y) = 0;
};

#endif
