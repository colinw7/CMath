#ifndef CPolygonOrientation_H
#define CPolygonOrientation_H

enum class CPolygonOrientation {
  UNKNOWN          =  0,
  CLOCKWISE        = -1,
  ANTICLOCKWISE    =  1,
  COUNTERCLOCKWISE = ANTICLOCKWISE
};

#endif
