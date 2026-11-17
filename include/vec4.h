#ifndef vec4_H
#define vec4_H

#include <vec2.h>
#include <vec3.h>

struct vec4 {
  vec4() { }

  explicit vec4(double r) :
   x(r), y(r), z(r), w(1) {
  }

  vec4(double x1, double y1, double z1) :
   x(x1), y(y1), z(z1), w(1) {
  }

  vec4(double x1, double y1, double z1, double w1) :
   x(x1), y(y1), z(z1), w(w1) {
  }

  //---

  vec4 &operator+=(const vec4 &v) {
    x += v.x; y += v.y; z += v.z; return *this;
  }

  vec4 &operator-=(const vec4 &v) {
    x -= v.x; y -= v.y; z -= v.z; return *this;
  }

  //---

  friend vec4 operator+(const vec4 &lhs, double rhs) {
    return vec4(lhs.x + rhs, lhs.y + rhs, lhs.z + rhs);
  }

  friend vec4 operator+(double lhs, const vec4 &rhs) {
    return vec4(lhs + rhs.x, lhs + rhs.y, lhs + rhs.z);
  }

  friend vec4 operator+(const vec4 &lhs, const vec4 &rhs) {
    return vec4(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
  }

  //---

  vec4 operator-() {
    return vec4(-x, -y, -z);
  }

  friend vec4 operator-(const vec4 &lhs, double rhs) {
    return vec4(lhs.x - rhs, lhs.y - rhs, lhs.z - rhs);
  }

  friend vec4 operator-(double lhs, const vec4 &rhs) {
    return vec4(lhs - rhs.x, lhs - rhs.y, lhs - rhs.z);
  }

  friend vec4 operator-(const vec4 &lhs, const vec4 &rhs) {
    return vec4(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }

  //---

  vec4 &operator*=(double f) {
    x *= f; y *= f; z *= f; return *this;
  }

  friend vec4 operator*(const vec4 &v, double f) {
    return vec4(v.x*f, v.y*f, v.z*f);
  }

  friend vec4 operator*(double f, const vec4 &v) {
    return vec4(f*v.x, f*v.y, f*v.z);
  }

  friend vec4 operator*(const vec4 &lhs, const vec4 &rhs) {
    return vec4(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
  }

  //---

  vec4 &operator/=(double f) {
    x /= f; y /= f; z /= f; return *this;
  }

  friend vec4 operator/(const vec4 &v, double f) {
    return vec4(v.x/f, v.y/f, v.z/f);
  }

  friend vec4 operator/(const vec4 &lhs, const vec4 &rhs) {
    return vec4(lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z);
  }

  //---

  vec2 xy() const { return vec2(x, y); }
  vec2 xz() const { return vec2(x, z); }
  vec2 yz() const { return vec2(y, z); }

  vec2 xx() const { return vec2(x, x); }
  vec2 yy() const { return vec2(y, y); }
  vec2 zz() const { return vec2(z, z); }

  vec3 xyz() const { return vec3(x, y, z); }
  vec3 yzx() const { return vec3(y, z, x); }
  vec3 zxy() const { return vec3(z, x, y); }

  vec3 www() const { return vec3(w, w, w); }

  void setXYZ(const vec3 &v) { x = v.x; y = v.y; z = v.y; }

  //---

  vec4 abs() const { return vec4(std::abs(x), std::abs(y), std::abs(z), 1); }

  //---

  vec4 min(const vec4 &v) const {
    return vec4(std::min(x, v.x), std::min(y, v.y), std::min(z, v.z));
  }

  vec4 min(double f) const { return min(vec4(f, f, f)); }

  vec4 max(const vec4 &v) const {
    return vec4(std::max(x, v.x), std::max(y, v.y), std::max(z, v.z));
  }

  vec4 max(double f) const { return max(vec4(f, f, f)); }

  //---

  double length() const { return std::sqrt(x*x + y*y + z*z); }

  //---

  // dot product
  double dot(const vec4 &v) const { return (x*v.x + y*v.y + z*v.z); }

  //---

  vec4 sqrt() { return vec4(std::sqrt(x), std::sqrt(y), std::sqrt(z)); }

  //---

  friend vec4 cos(const vec4 &v) { return v.cos(); }
  friend vec4 sin(const vec4 &v) { return v.sin(); }

  vec4 cos() const { return vec4(std::cos(x), std::cos(y), std::cos(z)); }
  vec4 sin() const { return vec4(std::sin(x), std::sin(y), std::sin(z)); }

  //---

  vec4 mix(const vec4 &v, double f) const {
    auto rmix = [](double x1, double y1, double f1) {
      return (x1*(1.0 - f1) + y1*f1);
    };

    return vec4(rmix(x, v.x, f), rmix(y, v.y, f), rmix(z, v.z, f));
  }

  //---

  vec4 normalize() const {
    auto len = length();
    if (len <= 0.0) return *this;

    auto factor = 1.0/len;

    return vec4(x*factor, y*factor, z*factor);
  }

  //---

  friend std::ostream &operator<<(std::ostream &os, const vec4 &v) {
    os << v.x << " " << v.y << " " << v.z << " " << v.w;
    return os;
  }

  //---

  double x { 0.0 };
  double y { 0.0 };
  double z { 0.0 };
  double w { 1.0 };
};

#endif
