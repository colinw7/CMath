#ifndef CIRange2D_H
#define CIRange2D_H

#include <iostream>

class CIRange2D {
 public:
  explicit CIRange2D(int x11 = 0, int y11 = 0, int x21 = 0, int y21 = 0) {
    set(x11, y11, x21, y21);
  }

  void set(int x11, int y11, int x21, int y21) {
    x1 = x11, y1 = y11, x2 = x21, y2 = y21;
  }

  void get(int *x11, int *y11, int *x21, int *y21) const {
    *x11 = x1, *y11 = y1, *x21 = x2, *y21 = y2;
  }

  int dx() const { return x2 - x1; }
  int dy() const { return y2 - y1; }

  int xmid() const { return (x2 + x1) >> 1; }
  int ymid() const { return (y2 + y1) >> 1; }

  int xmin() const { return std::min(x1, x2); }
  int ymin() const { return std::min(y1, y2); }
  int xmax() const { return std::max(x1, x2); }
  int ymax() const { return std::max(y1, y2); }

  int xsize() const { return abs(x2 - x1) + 1; }
  int ysize() const { return abs(y2 - y1) + 1; }

  void inc(int dx, int dy) { x1 += dx; y1 += dy; x2 += dx; y2 += dy; }

  void incX(int dx) { x1 += dx; x2 += dx; }
  void incY(int dy) { y1 += dy; y2 += dy; }

  CIRange2D &operator=(const CIRange2D &range) {
    x1 = range.x1; y1 = range.y1;
    x2 = range.x2; y2 = range.y2;

    return *this;
  }

  void print(std::ostream &os) const {
    os << x1 << "," << y1 << "," << x2 << "," << y2;
  }

  friend std::ostream &operator<<(std::ostream &os, const CIRange2D &range) {
    range.print(os);

    return os;
  }

 public:
  int x1 { 0 }, y1 { 0 }, x2 { 0 }, y2 { 0 };
};

#endif
