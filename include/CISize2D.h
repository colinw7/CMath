#ifndef CISIZE_2D_H
#define CISIZE_2D_H

#include <iostream>

template<typename T>
class CISize2DT {
 public:
  T width, height;

  CISize2DT() { }

  CISize2DT(T w, T h) :
   width(w), height(h) {
  }

  CISize2DT(const CISize2DT &size) :
   width(size.width), height(size.height) {
  }

  void set(T w, T h) {
    width  = w;
    height = h;
  }

  void get(T *w, T *h) const {
    *w = width;
    *h = height;
  }

  T getWidth () const { return width ; }
  T getHeight() const { return height; }

  void setWidth (T w) { width  = w; }
  void setHeight(T h) { height = h; }

  T area() const { return width*height; }

  friend CISize2DT operator*(T m, const CISize2DT &size) {
    return CISize2DT(m*size.width, m*size.height);
  }

  friend CISize2DT operator*(const CISize2DT &size, T m) {
    return CISize2DT(m*size.width, m*size.height);
  }

  friend CISize2DT operator/(const CISize2DT &size, T m) {
    return CISize2DT(size.width/m, size.height/m);
  }

  friend CISize2DT operator*(double m, const CISize2DT &size) {
    return CISize2DT(int(m*size.width), int(m*size.height));
  }

  friend CISize2DT operator*(const CISize2DT &size, double m) {
    return CISize2DT(int(m*size.width), int(m*size.height));
  }

  friend CISize2DT operator/(const CISize2DT &size, double m) {
    return CISize2DT(int(size.width/m), int(size.height/m));
  }

  bool operator==(const CISize2DT &size) const {
    return (width == size.width && height == size.height);
  }

  bool operator!=(const CISize2DT &size) const {
    return (width != size.width || height != size.height);
  }

  friend std::ostream &operator<<(std::ostream &os, const CISize2DT &size) {
    return os << "(" << size.width << "," << size.height << ")";
  }
};

typedef CISize2DT<int>  CISize2D;
typedef CISize2DT<uint> CUISize2D;

#endif
