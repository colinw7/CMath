#ifndef CMathGen_H
#define CMathGen_H

namespace CMathGen {
  enum AxisType2D {
    X_AXIS_2D = (1<<0),
    Y_AXIS_2D = (1<<1),

    XY_AXIS_2D = (X_AXIS_2D | Y_AXIS_2D)
  };

  enum IntersectType {
    INTERSECT_NONE    = 0,
    INTERSECT_ALL     = (1<<0),
    INTERSECT_INSIDE  = (1<<1),
    INTERSECT_OUTSIDE = (1<<2),
    INTERSECT_VALID   = (INTERSECT_INSIDE | INTERSECT_OUTSIDE)
  };
}

#endif
