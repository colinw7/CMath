#ifndef C_END_CAP_H
#define C_END_CAP_H

#include <CIPoint2D.h>
#include <cmath>
#include <cassert>

template<typename T>
int sign(const T &t) {
  return (t >= T(0) ? (t > T(0) ? 1 : 0) : -1);
}

class CPathOutline {
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
    POINT_NONE     = 0,
    POINT_FLIP_X   = (1<<0),
    POINT_FLIP_Y   = (1<<1),
    POINT_OPPOSITE = (1<<2)
  };

 public:
  CPathOutline(int width) {
    hw_ = width/2.0;
    s2_ = sqrt(2.0);

    a1_ = hw_/s2_;
    a2_ = hw_*s2_;
    a3_ = hw_*(s2_ - 1.0);

    initRect     ();
    initOctagon  ();
    initIntersect();
  }

  void calcPolygon(const std::vector<CIPoint2D> &points, EndCap startCap, EndCap stopCap,
                   std::vector<CIPoint2D> &polygon) {
    uint np = points.size();

    assert(np >= 2);

    EndCap midCap = END_CAP_EXTEND;

    if (startCap == END_CAP_OCTAGON || stopCap == END_CAP_OCTAGON)
      midCap = END_CAP_OCTAGON;

    polygon.clear();

    polygon.reserve(4*np + 4);

    addStartPoints(points[0], points[1], startCap);

    for (uint i = 1; i < np - 1; ++i) {
      const CIPoint2D &p0 = points[i - 1];
      const CIPoint2D &p1 = points[i    ];
      const CIPoint2D &p2 = points[i + 1];

      addMidPoints(p0, p1, p2, midCap);
    }

    addEndPoints(points[np - 1], points[np - 2], stopCap);
  }

 private:
  void addStartPoints(const CIPoint2D &p1, const CIPoint2D &p2, EndCap startCap) {
    int angle = calcAngle(p1, p2);

    uint pf = POINT_NONE;

    if      (angle == 0 || angle == 180) {
      if (angle == 180) pf |= POINT_OPPOSITE;

      if (startCap != END_CAP_OCTAGON) {
        addHPoint(getHRectPoint(p1, 1, startCap, pf));
        addLPoint(getHRectPoint(p1, 0, startCap, pf));
      }
      else {
        addHPoint(getOctagonPoint(p1, 0, pf));
        addHPoint(getOctagonPoint(p1, 1, pf));
        addLPoint(getOctagonPoint(p1, 7, pf));
        addLPoint(getOctagonPoint(p1, 6, pf));
      }
    }
    else if (angle == 90 || angle == -90) {
      if (angle == -90) pf |= POINT_OPPOSITE;

      if (startCap != END_CAP_OCTAGON) {
        addHPoint(getVRectPoint(p1, 0, startCap, pf));
        addLPoint(getVRectPoint(p1, 3, startCap, pf));
      }
      else {
        addHPoint(getOctagonPoint(p1, 6, pf));
        addHPoint(getOctagonPoint(p1, 7, pf));
        addLPoint(getOctagonPoint(p1, 5, pf));
        addLPoint(getOctagonPoint(p1, 4, pf));
      }
    }
    else if (angle == 45 || angle == 135) {
      if (angle == 135) pf |= POINT_OPPOSITE;

      if (startCap != END_CAP_OCTAGON) {
        addHPoint(get45EndPoint(p1, 0, 3, startCap, pf));
        addLPoint(get45EndPoint(p1, 3, 0, startCap, pf));
      }
      else {
        addHPoint(getOctagonPoint(p1, 7, pf));
        addHPoint(getOctagonPoint(p1, 0, pf));
        addLPoint(getOctagonPoint(p1, 6, pf));
        addLPoint(getOctagonPoint(p1, 5, pf));
      }
    }
    else if (angle == -45 || angle == -135) {
      if (xs == -135) pf |= POINT_OPPOSITE;

      if (startCap != END_CAP_OCTAGON) {
        addHPoint(get45EndPoint(p1, 1, 0, startCap, pf));
        addLPoint(get45EndPoint(p1, 0, 1, startCap, pf));
      }
      else {
        addHPoint(getOctagonPoint(p1, 1, pf));
        addHPoint(getOctagonPoint(p1, 2, pf));
        addLPoint(getOctagonPoint(p1, 0, pf));
        addLPoint(getOctagonPoint(p1, 7, pf));
      }
    }
    else
      assert(false);
  }

  void addEndPoints(const CIPoint2D &p2, const CIPoint2D &p1, EndCap stopCap) {
    int angle = calcAngle(p1, p2);

    uint pf = POINT_NONE;

    if      (direction == DIR_HORIZONTAL) {
      if (xs < 0) pf |= POINT_OPPOSITE;

      if (stopCap != END_CAP_OCTAGON) {
        addHPoint(getHRectPoint(p2, 2, stopCap, pf));
        addLPoint(getHRectPoint(p2, 3, stopCap, pf));
      }
      else {
        addHPoint(getOctagonPoint(p2, 2, pf));
        addHPoint(getOctagonPoint(p2, 3, pf));
        addLPoint(getOctagonPoint(p2, 5, pf));
        addLPoint(getOctagonPoint(p2, 4, pf));
      }
    }
    else if (direction == DIR_VERTICAL) {
      if (ys < 0) pf |= POINT_OPPOSITE;

      if (stopCap != END_CAP_OCTAGON) {
        addHPoint(getVRectPoint(p2, 1, stopCap, pf));
        addLPoint(getVRectPoint(p2, 2, stopCap, pf));
      }
      else {
        addHPoint(getOctagonPoint(p2, 0, pf));
        addHPoint(getOctagonPoint(p2, 1, pf));
        addLPoint(getOctagonPoint(p2, 3, pf));
        addLPoint(getOctagonPoint(p2, 2, pf));
      }
    }
    else if (direction == DIR_45_UP) {
      if (xs < 0) pf |= POINT_OPPOSITE;

      if (stopCap != END_CAP_OCTAGON) {
        addHPoint(get45EndPoint(p2, 1, 2, stopCap, pf));
        addLPoint(get45EndPoint(p2, 2, 1, stopCap, pf));
      }
      else {
        addHPoint(getOctagonPoint(p2, 1, pf));
        addHPoint(getOctagonPoint(p2, 2, pf));
        addLPoint(getOctagonPoint(p2, 4, pf));
        addLPoint(getOctagonPoint(p2, 3, pf));
      }
    }
    else if (direction == DIR_45_DOWN) {
      if (xs < 0) pf |= POINT_OPPOSITE;

      if (stopCap != END_CAP_OCTAGON) {
        addHPoint(get45EndPoint(p2, 2, 3, stopCap, pf));
        addLPoint(get45EndPoint(p2, 3, 2, stopCap, pf));
      }
      else {
        addHPoint(getOctagonPoint(p2, 3, pf));
        addHPoint(getOctagonPoint(p2, 4, pf));
        addLPoint(getOctagonPoint(p2, 6, pf));
        addLPoint(getOctagonPoint(p2, 5, pf));
      }
    }
    else
      assert(false);
  }

  void addMidPoints(const CIPoint2D &p1, const CIPoint2D &p2, const CIPoint2D &p3, EndCap midCap) {
    int angle1 = calcAngle(p1, p2);
    int angle2 = calcAngle(p2, p3);

    uint pf = POINT_NONE;

    if      (direction1 == direction2) {
      if      (direction1 == DIR_HORIZONTAL) {
        if (xs1 < 0) pf |= POINT_OPPOSITE;

        addHPoint(getHRectPoint(p2, 2, END_CAP_FLAT, pf));
        addLPoint(getHRectPoint(p2, 0, END_CAP_FLAT, pf));
      }
      else if (direction1 == DIR_VERTICAL) {
        if (ys1 < 0) pf |= POINT_OPPOSITE;

        addHPoint(getVRectPoint(p2, 0, END_CAP_FLAT, pf));
        addLPoint(getVRectPoint(p2, 2, END_CAP_FLAT, pf));
      }
      else if (direction1 == DIR_45_UP) {
        if (xs1 < 0) pf |= POINT_OPPOSITE;

        addHPoint(get45EndPoint(p2, 0, 3, END_CAP_FLAT, pf));
        addLPoint(get45EndPoint(p2, 3, 0, END_CAP_FLAT, pf));
      }
      else if (direction1 == DIR_45_DOWN) {
        if (xs1 < 0) pf |= POINT_OPPOSITE;

        addHPoint(get45EndPoint(p2, 1, 0, END_CAP_FLAT, pf));
        addLPoint(get45EndPoint(p2, 0, 1, END_CAP_FLAT, pf));
      }
      else
        assert(false);
    }
    else if (direction1 == DIR_HORIZONTAL) {
      if      (direction2 == DIR_VERTICAL) {
        if (midCap != END_CAP_OCTAGON) {
          if (ys2 < 0) pf |= POINT_FLIP_X;
          if (xs1 < 0) pf |= POINT_FLIP_Y;

          addHPoint(getHRectPoint(p2, 1, END_CAP_EXTEND, pf));
          addLPoint(getHRectPoint(p2, 3, END_CAP_EXTEND, pf));
        }
        else {
          if (xs1 < 0) pf |= POINT_OPPOSITE;

          if (xs1*ys2 > 0) {
            addHPoint(getHRectPoint(p2, 1, END_CAP_EXTEND, pf));
            addLPoint(getOctagonPoint(p2, 5, pf));
            addLPoint(getOctagonPoint(p2, 4, pf));
          }
          else {
            addLPoint(getHRectPoint(p2, 0, END_CAP_EXTEND, pf));
            addHPoint(getOctagonPoint(p2, 2, pf));
            addHPoint(getOctagonPoint(p2, 3, pf));
          }
        }
      }
      else if (direction2 == DIR_45_UP) {
        if (xs1 < 0) pf |= POINT_OPPOSITE;

        if (xs1*xs2 > 0) {
          addHPoint(getOctagonPoint(p2, 1, pf));
          addLPoint(getOctagonPoint(p2, 5, pf));
        }
        else {
          addHPoint(getOctagonPoint(p2, 2, pf));
          addHPoint(getOctagonPoint(p2, 3, pf));
          addHPoint(getOctagonPoint(p2, 4, pf));
          addLPoint(getIntersectPoint(p2, 7, pf));
        }
      }
      else if (direction2 == DIR_45_DOWN) {
        if (xs1 < 0) pf |= POINT_OPPOSITE;

        if (xs1*xs2 > 0) {
          addHPoint(getOctagonPoint(p2, 2, pf));
          addLPoint(getOctagonPoint(p2, 6, pf));
        }
        else {
          addHPoint(getIntersectPoint(p2, 0, pf));
          addLPoint(getOctagonPoint(p2, 5, pf));
          addLPoint(getOctagonPoint(p2, 4, pf));
          addLPoint(getOctagonPoint(p2, 3, pf));
        }
      }
      else
        assert(false);
    }
    else if (direction1 == DIR_VERTICAL) {
      if      (direction2 == DIR_HORIZONTAL) {
        if (midCap != END_CAP_OCTAGON) {
          if (xs2 < 0) pf |= POINT_FLIP_Y;
          if (ys1 < 0) pf |= POINT_FLIP_X;

          addLPoint(getHRectPoint(p2, 3, END_CAP_EXTEND, pf));
          addHPoint(getHRectPoint(p2, 1, END_CAP_EXTEND, pf));
        }
        else {
          if (ys1 < 0) pf |= POINT_OPPOSITE;

          if (ys1*xs2 > 0) {
            addLPoint(getHRectPoint(p2, 3, END_CAP_EXTEND, pf));
            addHPoint(getOctagonPoint(p2, 0, pf));
            addHPoint(getOctagonPoint(p2, 1, pf));
          }
          else {
            addHPoint(getHRectPoint(p2, 0, END_CAP_EXTEND, pf));
            addLPoint(getOctagonPoint(p2, 3, pf));
            addLPoint(getOctagonPoint(p2, 2, pf));
          }
        }
      }
      else if (direction2 == DIR_45_UP) {
        if (ys1 < 0) pf |= POINT_OPPOSITE;

        if (ys1*xs2 > 0) {
          addLPoint(getOctagonPoint(p2, 4, pf));
          addHPoint(getOctagonPoint(p2, 0, pf));
        }
        else {
          addLPoint(getOctagonPoint(p2, 3, pf));
          addLPoint(getOctagonPoint(p2, 2, pf));
          addLPoint(getOctagonPoint(p2, 1, pf));
          addHPoint(getIntersectPoint(p2, 6, pf));
        }
      }
      else if (direction2 == DIR_45_DOWN) {
        if (ys1 < 0) pf |= POINT_OPPOSITE;

        if (ys1*xs2 > 0) {
          addHPoint(getOctagonPoint(p2, 0, pf));
          addHPoint(getOctagonPoint(p2, 1, pf));
          addHPoint(getOctagonPoint(p2, 2, pf));
          addLPoint(getIntersectPoint(p2, 5, pf));
        }
        else {
          addHPoint(getOctagonPoint(p2, 7, pf));
          addLPoint(getOctagonPoint(p2, 3, pf));
        }
      }
      else
        assert(false);
    }
    else if (direction1 == DIR_45_UP) {
      if (xs1 < 0) pf |= POINT_OPPOSITE;

      if      (direction2 == DIR_HORIZONTAL) {
        if (xs1*xs2 > 0) {
          addHPoint(getOctagonPoint(p2, 1, pf));
          addLPoint(getOctagonPoint(p2, 5, pf));
        }
        else {
          addHPoint(getIntersectPoint(p2, 7, pf));
          addLPoint(getOctagonPoint(p2, 4, pf));
          addLPoint(getOctagonPoint(p2, 3, pf));
          addLPoint(getOctagonPoint(p2, 2, pf));
        }
      }
      else if (direction2 == DIR_VERTICAL) {
        if (xs1*ys2 > 0) {
          addHPoint(getOctagonPoint(p2, 0, pf));
          addLPoint(getOctagonPoint(p2, 4, pf));
        }
        else {
          addHPoint(getOctagonPoint(p2, 1, pf));
          addHPoint(getOctagonPoint(p2, 2, pf));
          addHPoint(getOctagonPoint(p2, 3, pf));
          addLPoint(getIntersectPoint(p2, 6, pf));
        }
      }
      else if (direction2 == DIR_45_DOWN) {
        if (xs1*xs2 > 0) {
          if (midCap != END_CAP_OCTAGON) {
            addHPoint(get45MidPoint(p2, 1, pf));
            addLPoint(get45MidPoint(p2, 3, pf));
          }
          else {
            addHPoint(getOctagonPoint(p2, 1, pf));
            addHPoint(getOctagonPoint(p2, 2, pf));
            addLPoint(get45MidPoint(p2, 3, pf));
          }
        }
        else {
          if (midCap != END_CAP_OCTAGON) {
            addHPoint(get45MidPoint(p2, 0, pf));
            addLPoint(get45MidPoint(p2, 2, pf));
          }
          else {
            addHPoint(get45MidPoint(p2, 0, pf));
            addLPoint(getOctagonPoint(p2, 4, pf));
            addLPoint(getOctagonPoint(p2, 3, pf));
          }
        }
      }
    }
    else if (direction1 == DIR_45_DOWN) {
      if (xs1 < 0) pf |= POINT_OPPOSITE;

      if      (direction2 == DIR_HORIZONTAL) {
        if (xs1*xs2 > 0) {
          addLPoint(getOctagonPoint(p2, 6, pf));
          addHPoint(getOctagonPoint(p2, 2, pf));
        }
        else {
          addLPoint(getIntersectPoint(p2, 0, pf));
          addHPoint(getOctagonPoint(p2, 3, pf));
          addHPoint(getOctagonPoint(p2, 4, pf));
          addHPoint(getOctagonPoint(p2, 5, pf));
        }
      }
      else if (direction2 == DIR_VERTICAL) {
        if (xs1*ys2 > 0) {
          addLPoint(getOctagonPoint(p2, 6, pf));
          addLPoint(getOctagonPoint(p2, 5, pf));
          addLPoint(getOctagonPoint(p2, 4, pf));
          addHPoint(getIntersectPoint(p2, 1, pf));
        }
        else {
          addLPoint(getOctagonPoint(p2, 7, pf));
          addHPoint(getOctagonPoint(p2, 3, pf));
        }
      }
      else if (direction2 == DIR_45_UP) {
        if (xs1*xs2 > 0) {
          if (midCap != END_CAP_OCTAGON) {
            addLPoint(get45MidPoint(p2, 3, pf));
            addHPoint(get45MidPoint(p2, 1, pf));
          }
          else {
            addLPoint(getOctagonPoint(p2, 6, pf));
            addLPoint(getOctagonPoint(p2, 5, pf));
            addHPoint(get45MidPoint(p2, 1, pf));
          }
        }
        else {
          if (midCap != END_CAP_OCTAGON) {
            addLPoint(get45MidPoint(p2, 0, pf));
            addHPoint(get45MidPoint(p2, 2, pf));
          }
          else {
            addLPoint(get45MidPoint(p2, 0, pf));
            addHPoint(getOctagonPoint(p2, 3, pf));
            addHPoint(getOctagonPoint(p2, 4, pf));
          }
        }
      }
    }
    else
      assert(false);
  }

 public:
  void getEndBoundary(const CIPoint2D &p1, const CIPoint2D &p2, const Direction &direction,
                      const EndCap &endCap, std::vector<CIPoint2D> &points) {
    if      (endCap == END_CAP_OCTAGON) {
      for (uint i = 0; i < 8; ++i)
        points.push_back(getOctagonPoint(p1, i));
    }
    else {
      double x[4], y[4];

      if (endCap == END_CAP_FLAT) {
        if (direction == DIR_HORIZONTAL || direction == DIR_45_UP) {
          int sx = sign(p2.getX() - p1.getX());

          x[0] = p1.getX()           ; x[3] = x[0];
          x[2] = p1.getX() + sx*2*hw_; x[1] = x[2];

          y[0] = p1.getY() - hw_; y[1] = y[0];
          y[2] = p1.getY() + hw_; y[3] = y[2];
        }
        else {
          int sy = sign(p2.getY() - p1.getY());

          y[0] = p1.getY()           ; y[1] = y[0];
          y[2] = p1.getY() + sy*2*hw_; y[3] = y[2];

          x[0] = p1.getX() - hw_; x[3] = x[0];
          x[2] = p1.getX() + hw_; x[1] = x[2];
        }
      }
      else {
        x[0] = p1.getX() - hw_; x[3] = x[0];
        x[2] = p1.getX() + hw_; x[1] = x[2];
        y[0] = p1.getY() - hw_; y[1] = y[0];
        y[2] = p1.getY() + hw_; y[3] = y[2];
      }

      if (direction == DIR_45_UP || direction == DIR_45_DOWN)
        rotatePolygon45(x, y, 4, p1.getX(), p1.getY());

      for (uint i = 0; i < 4; ++i)
        points.push_back(CIPoint2D(round(x[i]), round(y[i])));
    }
  }

  //----

 private:
  CIPoint2D getHRectPoint(const CIPoint2D &o, uint i, EndCap endCap,
                          uint flags=POINT_NONE) const {
    int es = (endCap == END_CAP_FLAT ? 0 : 1);

    if (flags & POINT_FLIP_X  ) i = 3 - i;
    if (flags & POINT_FLIP_Y  ) i = (i > 1 ? 5 - i : 1 - i);
    if (flags & POINT_OPPOSITE) i = (i + 2) % 4;

    return CIPoint2D(round(o.getX() + es*rx_[i]), round(o.getY() + ry_[i]));
  }

  CIPoint2D getVRectPoint(const CIPoint2D &o, uint i, EndCap endCap,
                          uint flags=POINT_NONE) const {
    int es = (endCap == END_CAP_FLAT ? 0 : 1);

    if (flags & POINT_FLIP_X  ) i = 3 - i;
    if (flags & POINT_FLIP_Y  ) i = (i > 1 ? 5 - i : 1 - i);
    if (flags & POINT_OPPOSITE) i = (i + 2) % 4;

    return CIPoint2D(round(o.getX() + rx_[i]), round(o.getY() + es*ry_[i]));
  }

  // Rectangle Points are:
  //  1  2
  //  0  3
  void initRect() {
    rx_[0] = -hw_; rx_[1] = rx_[0];
    rx_[2] =  hw_; rx_[3] = rx_[2];

    ry_[0] = -hw_; ry_[3] = ry_[0];
    ry_[1] =  hw_; ry_[2] = ry_[1];
  }

  //----

  // End 45 Rectangle Points are:
  //    1
  //  0   2
  //    3
  CIPoint2D get45EndPoint(const CIPoint2D &o, uint i1, uint i2, EndCap endCap,
                          uint flags=POINT_NONE) const {
    if (flags & POINT_OPPOSITE) { i1 = (i1 + 2) % 4; i2 = (i2 + 2) % 4; }

    if (endCap == END_CAP_FLAT) {
      if      (i1 == 0 && i2 == 1) return CIPoint2D(round(o.getX() - a1_), round(o.getY() - a1_));
      else if (i1 == 3 && i2 == 2) return CIPoint2D(round(o.getX() - a1_), round(o.getY() - a1_));
      else if (i1 == 1 && i2 == 0) return CIPoint2D(round(o.getX() + a1_), round(o.getY() + a1_));
      else if (i1 == 2 && i2 == 3) return CIPoint2D(round(o.getX() + a1_), round(o.getY() + a1_));

      else if (i1 == 0 && i2 == 3) return CIPoint2D(round(o.getX() - a1_), round(o.getY() + a1_));
      else if (i1 == 1 && i2 == 2) return CIPoint2D(round(o.getX() - a1_), round(o.getY() + a1_));
      else if (i1 == 2 && i2 == 1) return CIPoint2D(round(o.getX() + a1_), round(o.getY() - a1_));
      else if (i1 == 3 && i2 == 0) return CIPoint2D(round(o.getX() + a1_), round(o.getY() - a1_));
      else                         assert(false);
    }
    else {
      if      (i1 == 0) return CIPoint2D(round(o.getX() - a2_), round(o.getY()      ));
      else if (i1 == 1) return CIPoint2D(round(o.getX()      ), round(o.getY() + a2_));
      else if (i1 == 2) return CIPoint2D(round(o.getX() + a2_), round(o.getY()      ));
      else if (i1 == 3) return CIPoint2D(round(o.getX()      ), round(o.getY() - a2_));
      else              assert(false);
    }
  }

  // Rotated Rectangle Points are:
  //    1
  //  0   2
  //    3
  CIPoint2D get45MidPoint(const CIPoint2D &o, uint i, uint flags=POINT_NONE) const {
    if (flags & POINT_FLIP_X  ) i = 3 - i;
    if (flags & POINT_FLIP_Y  ) i = (i > 1 ? 5 - i : 1 - i);
    if (flags & POINT_OPPOSITE) i = (i + 2) % 4;

    if      (i == 0) return CIPoint2D(round(o.getX() - hw_ - a3_), round(o.getY()            ));
    else if (i == 1) return CIPoint2D(round(o.getX()            ), round(o.getY() + hw_ + a3_));
    else if (i == 2) return CIPoint2D(round(o.getX() + hw_ + a3_), round(o.getY()            ));
    else if (i == 3) return CIPoint2D(round(o.getX()            ), round(o.getY() - hw_ - a3_));
    else             assert(false);
  }

  //----

  CIPoint2D getOctagonPoint(const CIPoint2D &o, uint i, uint flags=POINT_NONE) const {
    if (flags & POINT_OPPOSITE) i = (i + 4) % 8;

    return CIPoint2D(round(o.getX() + hx_[i]), round(o.getY() + hy_[i]));
  }

  // Octagon Points are stored as:
  //     1 2
  //    0   3
  //    7   4
  //     6 5
  void initOctagon() {
    hx_[0] = -hw_; hx_[7] = hx_[0];
    hx_[1] = -a3_; hx_[6] = hx_[1];
    hx_[2] =  a3_; hx_[5] = hx_[2];
    hx_[3] =  hw_; hx_[4] = hx_[3];

    hy_[0] =  a3_; hy_[3] = hy_[0];
    hy_[1] =  hw_; hy_[2] = hy_[1];
    hy_[6] = -hw_; hy_[5] = hy_[6];
    hy_[7] = -a3_; hy_[4] = hy_[7];
  }

  //----

  // if 'i' is negative then use opposite point
  CIPoint2D getIntersectPoint(const CIPoint2D &o, uint i, uint flags=POINT_NONE) const {
    if (flags & POINT_OPPOSITE) i = (i + 4) % 8;

    return CIPoint2D(round(o.getX() + ix_[i]), round(o.getY() + iy_[i]));
  }

  // Intersect Octagon Points are stored as:
  //     1 2
  //    0   3
  //    7   4
  //     6 5
  void initIntersect() {
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

#if 0
  double getHVExtent(EndCap endCap) const {
    if      (endCap == END_CAP_FLAT   ) return 0;
    else if (endCap == END_CAP_EXTEND ) return hw_;
    else if (endCap == END_CAP_OCTAGON) return hw_;
    else    assert(false);
  }
#endif

 public:
  static int calcAngle(const CIPoint2D &p1, const CIPoint2D &p2) {
    int dx = p2.getX() - p1.getX();
    int dy = p2.getY() - p1.getY();

    if      (dx > 0) {
      if      (dy > 0) return  45;
      else if (dy < 0) return -45;
      else             return   0;
    }
    else if (dx < 0) {
      if      (dy > 0) return  135;
      else if (dy < 0) return -135;
      else             return  180;
    }
    else {
      if      (dy > 0) return  90;
      else if (dy < 0) return -90;
      else             return  0;
    }
  }

  static void rotatePolygon45(double *x, double *y, uint num_xy, double xc, double yc) {
    for (uint i = 0; i < num_xy; ++i)
      rotatePoint45(x[i], y[i], xc, yc, &x[i], &y[i]);
  }

 private:
  void addHPoint(const CIPoint2D &p) { hpoints_.push_back(p); }
  void addLPoint(const CIPoint2D &p) { lpoints_.push_back(p); }

 private:
  static void rotatePoint45(double x, double y, double xc, double yc, double *rx, double *ry) {
    static double c = 1.0/sqrt(2.0);
    static double s = c;

    double x1 = x - xc;
    double y1 = y - yc;

    double x2 = x1*c - y1*s;
    double y2 = x1*s + y1*c;

    x2 += xc;
    y2 += yc;

    *rx = x2;
    *ry = y2;
  }

 private:
  double                 hw_;
  double                 s2_;
  double                 a1_, a2_, a3_;
  double                 rx_[4], ry_[4];
  double                 hx_[8], hy_[8];
  double                 ix_[8], iy_[8];
  std::vector<CIPoint2D> hpoints_, lpoints_;
};

#endif
