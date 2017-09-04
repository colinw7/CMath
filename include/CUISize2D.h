#ifndef CUISIZE_2D_H
#define CUISIZE_2D_H

#include <iostream>
#include <sys/types.h>

class CUISize2D {
 public:
  unsigned int width, height;

  CUISize2D() { }

  CUISize2D(unsigned int w, unsigned int h) :
   width(w), height(h) {
  }

  CUISize2D(const CUISize2D &size) :
   width(size.width), height(size.height) {
  }

  void set(unsigned int w, unsigned int h) {
    width  = w;
    height = h;
  }

  void get(unsigned int *w, unsigned int *h) const {
    *w = width;
    *h = height;
  }

  unsigned int getWidth () const { return width ; }
  unsigned int getHeight() const { return height; }

  void setWidth (unsigned int w) { width  = w; }
  void setHeight(unsigned int h) { height = h; }

  unsigned int area() const { return width*height; }

  friend CUISize2D operator*(unsigned int m, const CUISize2D &size) {
    return CUISize2D(m*size.width, m*size.height);
  }

  friend CUISize2D operator*(const CUISize2D &size, unsigned int m) {
    return CUISize2D(m*size.width, m*size.height);
  }

  friend CUISize2D operator/(const CUISize2D &size, unsigned int m) {
    return CUISize2D(size.width/m, size.height/m);
  }

  friend CUISize2D operator*(double m, const CUISize2D &size) {
    return CUISize2D((unsigned int)(m*size.width), (unsigned int)(m*size.height));
  }

  friend CUISize2D operator*(const CUISize2D &size, double m) {
    return CUISize2D((unsigned int)(m*size.width), (unsigned int)(m*size.height));
  }

  friend CUISize2D operator/(const CUISize2D &size, double m) {
    return CUISize2D((unsigned int)(size.width/m), (unsigned int)(size.height/m));
  }

  bool operator==(const CUISize2D &size) const {
    return (width == size.width && height == size.height);
  }

  bool operator!=(const CUISize2D &size) const {
    return (width != size.width || height != size.height);
  }

  friend std::ostream &operator<<(std::ostream &os, const CUISize2D &size) {
    return os << "(" << size.width << "," << size.height << ")";
  }
};

#endif
