#ifndef CSIZE_2D_H
#define CSIZE_2D_H

#include <CISize2D.h>

template<typename T>
class CSize2DT {
 private:
  typedef CSize2DT<T>  Size;
  typedef CPoint2DT<T> Point;

 public:
  T width, height;

  CSize2DT() :
   width(0), height() {
  }

  CSize2DT(T w, T h) :
   width(w), height(h) {
  }

  CSize2DT(const CISize2D &size) :
   width(size.getWidth()), height(size.getHeight()) {
  }

  CSize2DT(const CSize2DT &size) :
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

  friend Size operator*(T m, const Size &size) {
    return Size(m*size.width, m*size.height);
  }

  friend Size operator*(const Size &size, T m) {
    return Size(m*size.width, m*size.height);
  }

  friend Size operator/(const Size &size, T m) {
    return Size(size.width/m, size.height/m);
  }

  friend std::ostream &operator<<(std::ostream &os, const Size &size) {
    return os << "(" << size.width << "," << size.height << ")";
  }

  friend Point operator+(const Size &s, const Point &p) {
    return Point(p.x + s.width, p.y + s.height);
  }

  friend Point operator+(const Point &p, const Size &s) {
    return (s + p);
  }
};

typedef CSize2DT<double> CSize2D;

#endif
