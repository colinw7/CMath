#ifndef vec3_H
#define vec3_H

#include <vec2.h>

struct vec3 {
  vec3() { }

  explicit vec3(double r) :
   x(r), y(r), z(r) {
  }

  vec3(double x1, double y1, double z1) :
   x(x1), y(y1), z(z1) {
  }

  //---

  vec3 &operator+=(const vec3 &v) {
    x += v.x; y += v.y; z += v.z; return *this;
  }

  vec3 &operator-=(const vec3 &v) {
    x -= v.x; y -= v.y; z -= v.z; return *this;
  }

  //---

  friend vec3 operator+(const vec3 &lhs, double rhs) {
    return vec3(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs);
  }

  friend vec3 operator+(double lhs, const vec3 &rhs) {
    return vec3(lhs + rhs.x, lhs + rhs.y, lhs + rhs.z);
  }

  friend vec3 operator+(const vec3 &lhs, const vec3 &rhs) {
    return vec3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
  }

  //---

  friend vec3 operator-(const vec3 &lhs, double rhs) {
    return vec3(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs);
  }

  friend vec3 operator-(double lhs, const vec3 &rhs) {
    return vec3(lhs - rhs.x, lhs - rhs.y, lhs - rhs.z);
  }

  friend vec3 operator-(const vec3 &lhs, const vec3 &rhs) {
    return vec3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }

  //---

  vec3 &operator*=(double f) {
    x *= f; y *= f; z *= f; return *this;
  }

  friend vec3 operator*(const vec3 &v, double f) {
    return vec3(v.x*f, v.y*f, v.z*f);
  }

  friend vec3 operator*(double f, const vec3 &v) {
    return vec3(f*v.x, f*v.y, f*v.z);
  }

  //---

  vec3 &operator/=(double f) {
    x /= f; y /= f; z /= f; return *this;
  }

  friend vec3 operator/(const vec3 &v, double f) {
    return vec3(v.x/f, v.y/f, v.z/f);
  }

  //---

  vec2 xy() const { return vec2(x, y); }
  vec2 xz() const { return vec2(x, z); }
  vec2 yz() const { return vec2(y, z); }

  vec3 xyz() const { return vec3(x, y, z); }
  vec3 yzx() const { return vec3(y, z, x); }
  vec3 zxy() const { return vec3(z, x, y); }

  void setXY(const vec2 &v) { x = v.x; y = v.y; }

  //---

  vec3 abs() const { return vec3(std::abs(x), std::abs(y), std::abs(z)); }

  //---

  vec3 min(const vec3 &v) const {
    return vec3(std::min(x, v.x), std::min(y, v.y), std::min(z, v.z));
  }

  vec3 min(double f) const { return min(vec3(f, f, f)); }

  vec3 max(const vec3 &v) const {
    return vec3(std::max(x, v.x), std::max(y, v.y), std::max(z, v.z));
  }

  vec3 max(double f) const { return max(vec3(f, f, f)); }

  //---

  double length() const { return std::sqrt(x*x + y*y + z*z); }

  //---

  // dot product
  double dot(const vec3 &v) const { return (x*v.x + y*v.y + z*v.z); }

  //---

  friend vec3 cos(const vec3 &v) { return v.cos(); }
  friend vec3 sin(const vec3 &v) { return v.sin(); }

  vec3 cos() const { return vec3(std::cos(x), std::cos(y), std::cos(z)); }
  vec3 sin() const { return vec3(std::sin(x), std::sin(y), std::sin(z)); }

  //---

  vec3 mix(const vec3 &v, double f) const {
    auto rmix = [](double x, double y, double f) {
      return (x*(1.0 - f) + y*f);
    };

    return vec3(rmix(x, v.x, f), rmix(y, v.y, f), rmix(z, v.z, f));
  }

  //---

  vec3 normalize() const {
    auto len = length();
    if (len <= 0.0) return *this;

    auto factor = 1.0/len;

    return vec3(x*factor, y*factor, z*factor);
  }

  //---

  double x { 0.0 };
  double y { 0.0 };
  double z { 0.0 };
};

#endif
