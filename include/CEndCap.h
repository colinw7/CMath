#ifndef C_END_CAP_H
#define C_END_CAP_H

#include <CIPoint2D.h>
#include <cmath>
#include <cassert>

template<typename T>
int sign(const T &t) {
  return (t >= T(0) ? (t > T(0) ? 1 : 0) : -1);
}

class CEndCap {
 public:
  enum Direction {
    DIR_NONE,
    DIR_HORIZONTAL,
    DIR_VERTICAL,
    DIR_45_UP,
    DIR_45_DOWN
  };

  enum EndCap {
    END_CAP_NONE,
    END_CAP_FLAT,
    END_CAP_EXTEND,
    END_CAP_OCTAGON
  };

  enum PointFlags {
    POINT_NONE       = 0,
    POINT_FLIP_X     = (1<<0),
    POINT_FLIP_Y     = (1<<1),
    POINT_ROTATE_45  = (1<<3),
    POINT_ROTATE_90  = (1<<4),
    POINT_ROTATE_135 = (1<<5),
    POINT_OPPOSITE   = (POINT_FLIP_X|POINT_FLIP_Y)
  };

  typedef std::vector<CIPoint2D> Points;

 public:
  CEndCap(int width) {
    hw_ = width/2.0;
    s2_ = sqrt(2.0);

    a1_ = hw_/s2_;
    a2_ = hw_*s2_;
    a3_ = hw_*(s2_ - 1.0);

    initPoints();
  }

  void calcPolygon(const CIPoint2D &p1, const CIPoint2D &p2, EndCap startCap, EndCap stopCap,
                   Points &polygon) {
    Points points;

    points.push_back(p1);
    points.push_back(p2);

    calcPolygon(points, startCap, stopCap, polygon);
  }

  void calcPolygon(const Points &points, EndCap startCap, EndCap stopCap, Points &polygon) {
    uint np = points.size();

    assert(np >= 2);

    EndCap midCap = END_CAP_EXTEND;

    if (startCap == END_CAP_OCTAGON || stopCap == END_CAP_OCTAGON)
      midCap = END_CAP_OCTAGON;

    polygon.clear();

    polygon.reserve(4*np + 4);

    hpoints_.clear(); hpoints_.reserve(2*np + 2);
    lpoints_.clear(); lpoints_.reserve(2*np + 2);

    addStartPoints(points[0], points[1], startCap);

    for (uint i = 1; i < np - 1; ++i) {
      const CIPoint2D &p0 = points[i - 1];
      const CIPoint2D &p1 = points[i    ];
      const CIPoint2D &p2 = points[i + 1];

      addMidPoints(p0, p1, p2, midCap);
    }

    addEndPoints(points[np - 1], points[np - 2], stopCap);

    std::copy(hpoints_. begin(), hpoints_. end(), std::back_inserter(polygon));
    std::copy(lpoints_.rbegin(), lpoints_.rend(), std::back_inserter(polygon));
  }

 private:
  void addStartPoints(const CIPoint2D &p1, const CIPoint2D &p2, EndCap startCap) {
    uint pf = calcPointFlags(p1, p2);

    if      (startCap == END_CAP_OCTAGON) {
      addHPoint(getOctagonPoint(p1, 0, pf));
      addHPoint(getOctagonPoint(p1, 1, pf));
      addLPoint(getOctagonPoint(p1, 7, pf));
      addLPoint(getOctagonPoint(p1, 6, pf));
    }
    else if (startCap == END_CAP_EXTEND) {
      addHPoint(getExtendPoint(p1, 1, pf));
      addLPoint(getExtendPoint(p1, 7, pf));
    }
    else {
      addHPoint(getFlatPoint(p1, 2, pf));
      addLPoint(getFlatPoint(p1, 6, pf));
    }
  }

  void addEndPoints(const CIPoint2D &p2, const CIPoint2D &p1, EndCap stopCap) {
    uint pf = calcPointFlags(p1, p2);

    if      (stopCap == END_CAP_OCTAGON) {
      addHPoint(getOctagonPoint(p2, 2, pf));
      addHPoint(getOctagonPoint(p2, 3, pf));
      addLPoint(getOctagonPoint(p2, 5, pf));
      addLPoint(getOctagonPoint(p2, 4, pf));
    }
    else if (stopCap == END_CAP_EXTEND) {
      addHPoint(getExtendPoint(p2, 3, pf));
      addLPoint(getExtendPoint(p2, 5, pf));
    }
    else {
      addHPoint(getFlatPoint(p2, 2, pf));
      addLPoint(getFlatPoint(p2, 6, pf));
    }
  }

  void addMidPoints(const CIPoint2D &p1, const CIPoint2D &p2, const CIPoint2D &p3, EndCap midCap) {
    uint pf = POINT_NONE;

    int angle = calcLineAngle(p1, p2, p3, pf);

    if      (angle == 0 || angle == 180) {
      addHPoint(getFlatPoint(p2, 2, pf));
      addLPoint(getFlatPoint(p2, 6, pf));
    }
    else if (angle == 90) {
      addHPoint(getExtendPoint(p2, 1, pf));

      if (midCap != END_CAP_OCTAGON)
        addLPoint(getExtendPoint(p2, 5, pf));
      else {
        if (! isFlipped(pf)) {
          addLPoint(getOctagonPoint(p2, 5, pf));
          addLPoint(getOctagonPoint(p2, 4, pf));
        }
        else {
          addLPoint(getOctagonPoint(p2, 4, pf));
          addLPoint(getOctagonPoint(p2, 5, pf));
        }
      }
    }
    else if (angle == 270) {
      if (midCap != END_CAP_OCTAGON)
        addHPoint(getExtendPoint(p2, 3, pf));
      else {
        if (! isFlipped(pf)) {
          addHPoint(getOctagonPoint(p2, 2, pf));
          addHPoint(getOctagonPoint(p2, 3, pf));
        }
        else {
          addHPoint(getOctagonPoint(p2, 3, pf));
          addHPoint(getOctagonPoint(p2, 2, pf));
        }
      }

      addLPoint(getExtendPoint(p2, 7, pf));
    }
    else if (angle == 45) {
      addHPoint(getOctagonPoint(p2, 1, pf));
      addLPoint(getOctagonPoint(p2, 5, pf));
    }
    else if (angle == 315) {
      addHPoint(getOctagonPoint(p2, 2, pf));
      addLPoint(getOctagonPoint(p2, 6, pf));
    }
    else if (angle == 135) {
      addHPoint(getIntersectPoint(p2, 0, pf));
      addLPoint(getOctagonPoint(p2, 5, pf));
      addLPoint(getOctagonPoint(p2, 4, pf));
      addLPoint(getOctagonPoint(p2, 3, pf));
    }
    else if (angle == 225) {
      addHPoint(getOctagonPoint(p2, 2, pf));
      addHPoint(getOctagonPoint(p2, 3, pf));
      addHPoint(getOctagonPoint(p2, 4, pf));
      addLPoint(getIntersectPoint(p2, 7, pf));
    }
    else
      assert(false);
  }

 public:
  void getEndBoundary(const CIPoint2D &p1, const CIPoint2D &p2, const Direction &direction,
                      const EndCap &endCap, Points &points) {
    if      (endCap == END_CAP_OCTAGON) {
      for (uint i = 0; i < 8; ++i)
        points.push_back(getOctagonPoint(p1, i));
    }
    else {
      uint pf = POINT_NONE;

      if (direction == DIR_45_UP || direction == DIR_45_DOWN)
        pf = POINT_ROTATE_45;

      CIPoint2D p;

      if (endCap == END_CAP_FLAT) {
        if (direction == DIR_HORIZONTAL || direction == DIR_45_UP) {
          int sx = sign(p2.getX() - p1.getX());

          p = p1 + CIPoint2D(sx*hw_, 0);
        }
        else {
          int sy = sign(p2.getY() - p1.getY());

          p = p1 + CIPoint2D(0, sy*hw_);
        }
      }
      else
        p = p1;

      points.push_back(getExtendPoint(p, 7, pf));
      points.push_back(getExtendPoint(p, 5, pf));
      points.push_back(getExtendPoint(p, 3, pf));
      points.push_back(getExtendPoint(p, 1, pf));
    }
  }

  //----

 private:
  CIPoint2D getExtendPoint(const CIPoint2D &o, uint i, uint flags=POINT_NONE) const {
    if (flags & POINT_ROTATE_45 ) i = (i + 7) % 8;
    if (flags & POINT_ROTATE_90 ) i = (i + 6) % 8;
    if (flags & POINT_ROTATE_135) i = (i + 5) % 8;
    if (flags & POINT_FLIP_X    ) i = (i > 4 ? 12 - i : 4 - i);
    if (flags & POINT_FLIP_Y    ) i = (8 - i) % 8;

    return CIPoint2D(round(o.getX() + ex_[i]), round(o.getY() + ey_[i]));
  }

  CIPoint2D getFlatPoint(const CIPoint2D &o, uint i, uint flags=POINT_NONE) const {
    if (flags & POINT_ROTATE_45 ) i = (i + 7) % 8;
    if (flags & POINT_ROTATE_90 ) i = (i + 6) % 8;
    if (flags & POINT_ROTATE_135) i = (i + 5) % 8;
    if (flags & POINT_FLIP_X    ) i = (i > 4 ? 12 - i : 4 - i);
    if (flags & POINT_FLIP_Y    ) i = (8 - i) % 8;

    return CIPoint2D(round(o.getX() + fx_[i]), round(o.getY() + fy_[i]));
  }

  CIPoint2D getOctagonPoint(const CIPoint2D &o, uint i, uint flags=POINT_NONE) const {
    if (flags & POINT_ROTATE_45 ) i = (i + 7) % 8;
    if (flags & POINT_ROTATE_90 ) i = (i + 6) % 8;
    if (flags & POINT_ROTATE_135) i = (i + 5) % 8;
    if (flags & POINT_FLIP_X    ) i = (i > 3 ? 11 - i : 3 - i);
    if (flags & POINT_FLIP_Y    ) i = 7 - i;

    return CIPoint2D(round(o.getX() + hx_[i]), round(o.getY() + hy_[i]));
  }

  CIPoint2D getIntersectPoint(const CIPoint2D &o, uint i, uint flags=POINT_NONE) const {
    if (flags & POINT_ROTATE_45 ) i = (i + 7) % 8;
    if (flags & POINT_ROTATE_90 ) i = (i + 6) % 8;
    if (flags & POINT_ROTATE_135) i = (i + 5) % 8;
    if (flags & POINT_FLIP_X    ) i = (i > 3 ? 11 - i : 3 - i);
    if (flags & POINT_FLIP_Y    ) i = 7 - i;

    return CIPoint2D(round(o.getX() + ix_[i]), round(o.getY() + iy_[i]));
  }

  //----

  void initPoints() {
    // extend rect points:
    //     2
    //   1   3
    // 0       4
    //   7   5
    //     6
    ex_[0] = -a2_;
    ex_[1] = -hw_; ex_[7] = ex_[1];
    ex_[2] =    0; ex_[6] = ex_[2];
    ex_[3] =  hw_; ex_[5] = ex_[3];
    ex_[4] =  a2_;

    ey_[6] = -a2_;
    ey_[7] = -hw_; ey_[5] = ey_[7];
    ey_[0] =    0; ey_[4] = ey_[0];
    ey_[1] =  hw_; ey_[3] = ey_[1];
    ey_[2] =  a2_;

    // flat rect points:
    //     2
    //   1   3
    // 0       4
    //   7   5
    //     6
    fx_[0] = -hw_;
    fx_[1] = -a1_; fx_[7] = fx_[1];
    fx_[2] =    0; fx_[6] = fx_[2];
    fx_[3] =  a1_; fx_[5] = fx_[3];
    fx_[4] =  hw_;

    fy_[6] = -hw_;
    fy_[7] = -a1_; fy_[5] = fy_[7];
    fy_[0] =    0; fy_[4] = fy_[0];
    fy_[1] =  a1_; fy_[3] = fy_[1];
    fy_[2] =  hw_;

    // Octagon Points are stored as:
    //     1 2
    //    0   3
    //    7   4
    //     6 5
    hx_[0] = -hw_; hx_[7] = hx_[0];
    hx_[1] = -a3_; hx_[6] = hx_[1];
    hx_[2] =  a3_; hx_[5] = hx_[2];
    hx_[3] =  hw_; hx_[4] = hx_[3];

    hy_[0] =  a3_; hy_[3] = hy_[0];
    hy_[1] =  hw_; hy_[2] = hy_[1];
    hy_[6] = -hw_; hy_[5] = hy_[6];
    hy_[7] = -a3_; hy_[4] = hy_[7];

    // Intersect Octagon Points are stored as:
    //     1 2
    //    0   3
    //    7   4
    //     6 5
    double s = 2*hw_ + a3_;

    ix_[0] =  -s ; ix_[7] = ix_[0];
    ix_[1] = -hw_; ix_[6] = ix_[1];
    ix_[2] =  hw_; ix_[5] = ix_[2];
    ix_[3] =   s ; ix_[4] = ix_[3];

    iy_[0] =  hw_; iy_[3] = iy_[0];
    iy_[1] =   s ; iy_[2] = iy_[1];
    iy_[6] =  -s ; iy_[5] = iy_[6];
    iy_[7] = -hw_; iy_[4] = iy_[7];
  }

  //----

 public:
  static uint calcPointFlags(const CIPoint2D &p1, const CIPoint2D &p2) {
    Direction direction;
    int       xs, ys;

    calcDirection(p1, p2, direction, xs, ys);

    uint pf = POINT_NONE;

    if (direction == DIR_HORIZONTAL || direction == DIR_45_UP) {
      if (xs < 0) pf |= POINT_OPPOSITE;
    }
    else {
      if (ys < 0) pf |= POINT_OPPOSITE;
    }

    if      (direction == DIR_45_UP   ) pf |= POINT_ROTATE_45;
    else if (direction == DIR_VERTICAL) pf |= POINT_ROTATE_90;
    else if (direction == DIR_45_DOWN ) pf |= POINT_ROTATE_135;

    return pf;
  }

  static void calcDirection(const CIPoint2D &p1, const CIPoint2D &p2, Direction &direction,
                            int &xs, int &ys) {
    int dx = p2.getX() - p1.getX();
    int dy = p2.getY() - p1.getY();

    xs = sign(dx);
    ys = sign(dy);

    if      (dx && ! dy)
      direction = DIR_HORIZONTAL;
    else if (dy && ! dx)
      direction = DIR_VERTICAL;
    else if (abs(dx) == abs(dy)) {
      if (sign(dx) == sign(dy))
        direction = DIR_45_UP;
      else
        direction = DIR_45_DOWN;
    }
    else
      direction = DIR_NONE;
  }

  static int calcLineAngle(const CIPoint2D &p1, const CIPoint2D &p2,
                           const CIPoint2D &p3, uint &pf) {
    Direction direction1, direction2;
    int       xs1, ys1, xs2, ys2;

    calcDirection(p1, p2, direction1, xs1, ys1);
    calcDirection(p2, p3, direction2, xs2, ys2);

    if (direction1 == direction2) {
      pf = calcPointFlags(p1, p2);

      return 0;
    }

    if      (direction1 == DIR_HORIZONTAL) {
      if (xs1 < 0) pf |= POINT_OPPOSITE;

      if      (direction2 == DIR_VERTICAL  ) {
        return (xs1*ys2 > 0 ? 90 : 270);
      }
      else if (direction2 == DIR_45_UP     ) {
        return (xs1*xs2 > 0 ? 45 : 225);
      }
      else if (direction2 == DIR_45_DOWN   ) {
        return (xs1*ys2 > 0 ? 135 : 315);
      }
    }
    else if (direction1 == DIR_VERTICAL) {
      if      (direction2 == DIR_HORIZONTAL) {
        if (xs2 > 0) pf |= POINT_FLIP_X;
        else         pf |= POINT_FLIP_Y;

        return (ys1*xs2 > 0 ? 270 : 90);
      }
      else if (direction2 == DIR_45_UP     ) {
        pf |= POINT_ROTATE_90;

        if (ys1 < 0) pf |= POINT_OPPOSITE;

        return (ys1*xs2 > 0 ? 315 : 135);
      }
      else if (direction2 == DIR_45_DOWN   ) {
        pf |= POINT_ROTATE_90;

        if (ys1 < 0) pf |= POINT_OPPOSITE;

        return (ys1*ys2 > 0 ? 45 : 225);
      }
    }
    else if (direction1 == DIR_45_UP) {
      if      (direction2 == DIR_HORIZONTAL) {
        if (xs1     < 0) pf |= POINT_FLIP_Y;
        if (xs2     > 0) pf |= POINT_FLIP_X;
        if (xs1*xs2 < 0) pf |= POINT_ROTATE_45;

        return (xs1*xs2 > 0 ? 315 : 135);
      }
      else if (direction2 == DIR_VERTICAL  ) {
        pf |= POINT_ROTATE_45;

        if (xs1 < 0) pf |= POINT_OPPOSITE;

        return (xs1*ys2 > 0 ? 45 : 225);
      }
      else if (direction2 == DIR_45_DOWN   ) {
        pf |= POINT_ROTATE_45;

        if (xs1 < 0) pf |= POINT_OPPOSITE;

        return (xs1*ys2 > 0 ? 90 : 270);
      }
    }
    else if (direction1 == DIR_45_DOWN) {
      if      (direction2 == DIR_HORIZONTAL) {
        if (xs2     < 0) pf |= POINT_FLIP_Y;
        if (ys1     < 0) pf |= POINT_FLIP_X;
        if (ys1*xs2 > 0) pf |= POINT_ROTATE_135;

        return (ys1*xs2 > 0 ? 225 : 45);
      }
      else if (direction2 == DIR_VERTICAL  ) {
        pf |= POINT_ROTATE_135;

        if (ys1 < 0) pf |= POINT_OPPOSITE;

        return (ys1*ys2 > 0 ? 315 : 135);
      }
      else if (direction2 == DIR_45_UP     ) {
        pf |= POINT_ROTATE_135;

        if (ys1 < 0) pf |= POINT_OPPOSITE;

        return (ys1*xs2 > 0 ? 270 : 90);
      }
    }

    assert(false);
  }

  static bool isFlipped(uint pf) {
    bool flip_x = (pf & POINT_FLIP_X);
    bool flip_y = (pf & POINT_FLIP_Y);

    if (flip_x && flip_y) return false;

    return (flip_x || flip_y);
  }

 private:
  void addHPoint(const CIPoint2D &p) { hpoints_.push_back(p); }
  void addLPoint(const CIPoint2D &p) { lpoints_.push_back(p); }

 private:
  double hw_;
  double s2_;
  double a1_, a2_, a3_;
  double hx_[8], hy_[8];
  double ix_[8], iy_[8];
  double ex_[8], ey_[8];
  double fx_[8], fy_[8];
  Points hpoints_, lpoints_;
};

#endif
