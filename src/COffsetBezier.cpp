#include <COffsetBezier.h>

struct ParamPoint {
  CPoint2D p;
  double   t;

  ParamPoint(const CPoint2D &p1, double t1) :
   p(p1), t(t1) {
  }
};

struct PointFunc {
  PointFunc(const ParamPoint *points, uint n) :
   points_(points), n_(n) {
  }

  CPoint2D operator()(double t) {
    CPoint2D p;

    for (uint i = 0; i < n_ - 1; ++i) {
      double t1 = points_[i    ].t;
      double t2 = points_[i + 1].t;

      assert(t1 <= t2);

      if (t >= t1 && t <= t2) {
        double x1 = points_[i    ].p.x, y1 = points_[i    ].p.y;
        double x2 = points_[i + 1].p.x, y2 = points_[i + 1].p.y;

        if (t1 != t2) {
          p.x = (t - t1)*(x2 - x1)/(t2 - t1) + x1;
          p.y = (t - t1)*(y2 - y1)/(t2 - t1) + y1;
        }
        else {
          p.x = x1;
          p.y = y1;
        }

        return p;
      }
    }

    //----

    assert(points_[0].t < points_[n_ - 1].t);

    if      (t <= points_[0].t) {
      double t1 = points_[0].t;
      double t2 = points_[1].t;

      double x1 = points_[0].p.x, y1 = points_[0].p.y;
      double x2 = points_[1].p.x, y2 = points_[1].p.y;

      if (t1 != t2) {
        p.x = (t - t1)*(x2 - x1)/(t2 - t1) + x1;
        p.y = (t - t1)*(y2 - y1)/(t2 - t1) + y1;
      }
      else {
        p.x = x1;
        p.y = y1;
      }
    }
    else if (t >= points_[n_ - 1].t) {
      double t1 = points_[n_ - 2].t;
      double t2 = points_[n_ - 1].t;

      double x1 = points_[n_ - 2].p.x, y1 = points_[n_ - 2].p.y;
      double x2 = points_[n_ - 1].p.x, y2 = points_[n_ - 1].p.y;

      if (t1 != t2) {
        p.x = (t - t1)*(x2 - x1)/(t2 - t1) + x1;
        p.y = (t - t1)*(y2 - y1)/(t2 - t1) + y1;
      }
      else {
        p.x = x2;
        p.y = y2;
      }
    }

    return p;
  }

  const ParamPoint *points_;
  uint              n_;
};

C3Bezier2D
CMathGeom2D::
offsetBezier(const C3Bezier2D &b, double offset, uint num_samples)
{
  COffsetBezierData data;

  offsetBezier(b, offset, data, num_samples, 1);

  return data.beziers.front();
}

void
CMathGeom2D::
offsetBezier(const C3Bezier2D &b, double offset, COffsetBezierData &data,
             uint num_samples, int num_curves)
{
  if (num_curves > 1) {
    double dt = 1.0/num_curves;

    C3Bezier2D bb = b;

    C3Bezier2D b1, b2;

    for (int i = 1; i < num_curves; ++i) {
      bb.split(dt, b1, b2);

      offsetBezier(b1, offset, data, num_samples, 1);

      bb = b2;
    }

    offsetBezier(bb, offset, data, num_samples, 1);

    return;
  }

  //----

  assert(num_samples >= 4);

  std::vector<ParamPoint> points;

  double dt = 1.0/(num_samples - 1);

  for (uint i = 0; i < num_samples; ++i) {
    double t = i*dt;

    CPoint2D p;

    b.calc(t, p);

    points.push_back(ParamPoint(p, t));
  }

  std::vector<ParamPoint> opoints;

  for (uint i = 0; i < num_samples; ++i) {
    double a;

    const CPoint2D &p2 = points[i].p;

    if      (i > 0 && i < num_samples - 1) {
      const CPoint2D &p1 = points[i - 1].p;
      //const CPoint2D &p3 = points[i + 1].p;

      double a1 = atan2(p2.y - p1.y, p2.x - p1.x);
      //double a2 = atan2(p3.y - p2.y, p3.x - p2.x);

      //a = (a1 + a2)/2;
      a = a1;
    }
    else if (i == 0) {
      const CPoint2D &p3 = points[i + 1].p;

      a = atan2(p3.y - p2.y, p3.x - p2.x);
    }
    else {
      const CPoint2D &p1 = points[i - 1].p;

      a = atan2(p2.y - p1.y, p2.x - p1.x);
    }

    CPoint2D p4, p5;

    CMathGeom2D::AddWidthToPoint(p2, a, offset, p4, p5);

    opoints.push_back(ParamPoint(p4, points[i].t));

    data.points.push_back(p4);
  }

  PointFunc f(&opoints[0], num_samples);

  C3Bezier2D b1 = C3Bezier2D::bestParamFit(f);

  data.beziers.push_back(b1);
}
