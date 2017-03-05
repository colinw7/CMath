#ifndef CISIZE_2D_H
#define CISIZE_2D_H

#include <iostream>
#include <sys/types.h>

class CISize2D {
 public:
  int width, height;

  CISize2D() { }

  CISize2D(int w, int h) :
   width(w), height(h) {
  }

  CISize2D(const CISize2D &size) :
   width(size.width), height(size.height) {
  }

  void set(int w, int h) {
    width  = w;
    height = h;
  }

  void get(int *w, int *h) const {
    *w = width;
    *h = height;
  }

  int getWidth () const { return width ; }
  int getHeight() const { return height; }

  void setWidth (int w) { width  = w; }
  void setHeight(int h) { height = h; }

  int area() const { return width*height; }

  friend CISize2D operator*(int m, const CISize2D &size) {
    return CISize2D(m*size.width, m*size.height);
  }

  friend CISize2D operator*(const CISize2D &size, int m) {
    return CISize2D(m*size.width, m*size.height);
  }

  friend CISize2D operator/(const CISize2D &size, int m) {
    return CISize2D(size.width/m, size.height/m);
  }

  friend CISize2D operator*(double m, const CISize2D &size) {
    return CISize2D(int(m*size.width), int(m*size.height));
  }

  friend CISize2D operator*(const CISize2D &size, double m) {
    return CISize2D(int(m*size.width), int(m*size.height));
  }

  friend CISize2D operator/(const CISize2D &size, double m) {
    return CISize2D(int(size.width/m), int(size.height/m));
  }

  bool operator==(const CISize2D &size) const {
    return (width == size.width && height == size.height);
  }

  bool operator!=(const CISize2D &size) const {
    return (width != size.width || height != size.height);
  }

  friend std::ostream &operator<<(std::ostream &os, const CISize2D &size) {
    return os << "(" << size.width << "," << size.height << ")";
  }
};

#endif
