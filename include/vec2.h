#ifndef vec2_H
#define vec2_H

#include <cmath>
#include <iostream>

struct vec2 {
  vec2() { }

  explicit vec2(double r) :
   x(r), y(r) {
  }

  vec2(double x1, double y1) :
   x(x1), y(y1) {
  }

  vec2 operator-() const {
    return vec2(-x, -y);
  }

  //---

  friend vec2 operator-(const vec2 &lhs, const vec2 &rhs) {
    return vec2(lhs.x - rhs.x, lhs.y - rhs.y);
  }

  friend vec2 operator-(const vec2 &lhs, double rhs) {
    return vec2(lhs.x - rhs, lhs.y - rhs);
  }

  //---

  friend vec2 operator+(const vec2 &lhs, const vec2 &rhs) {
    return vec2(lhs.x + rhs.x, lhs.y + rhs.y);
  }

  friend vec2 operator+(const vec2 &lhs, double rhs) {
    return vec2(lhs.x + rhs, lhs.y + rhs);
  }

  friend vec2 operator*(double f, const vec2 &v) {
    return vec2(v.x*f, v.y*f);
  }

  friend vec2 operator*(const vec2 &v, double f) {
    return vec2(v.x*f, v.y*f);
  }

  friend vec2 operator*(const vec2 &lhs, const vec2 &rhs) {
    return vec2(lhs.x*rhs.x, lhs.y*rhs.y);
  }

  friend vec2 operator/(const vec2 &v, double f) {
    return vec2(v.x/f, v.y/f);
  }

  //---

  double &operator[](int i) { return (&x)[i]; }

  vec2 xy() const { return vec2(x, y); }
  vec2 yx() const { return vec2(y, x); }

  //---

  vec2 min(const vec2 &v) const {
    return vec2(std::min(x, v.x), std::min(y, v.y));
  }

  vec2 min(double f) const { return min(vec2(f, f)); }

  vec2 max(const vec2 &v) const {
    return vec2(std::max(x, v.x), std::max(y, v.y));
  }

  vec2 max(double f) const { return max(vec2(f, f)); }

  //---

  vec2 abs() const { return vec2(std::abs(x), std::abs(y)); }

  vec2 floor() const { return vec2(std::floor(x), std::floor(y)); }

  //---

  double length() const { return std::hypot(x, y); }

  double dot(const vec2 &v) const { return (x*v.x + y*v.y); }

  //---

  friend std::ostream &operator<<(std::ostream &os, const vec2 &v) {
    os << v.x << " " << v.y;
    return os;
  }

  //---

  double x { 0.0 };
  double y { 0.0 };
};

#endif
