#ifndef CAngle_H
#define CAngle_H

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>

// Angle type
class CAngle {
 public:
  enum class Type {
    DEGREES,
    RADIANS
  };

 public:
  static CAngle makeDegrees(double r) { return CAngle(Type::DEGREES, r); }
  static CAngle makeRadians(double r) { return CAngle(Type::RADIANS, r); }

  static double Deg2Rad(double d) { return M_PI*d/180.0; }
  static double Rad2Deg(double d) { return 180.0*d/M_PI; }

  CAngle(Type type=Type::DEGREES, double angle=0.0) :
   type_(type), angle_(angle) {
  }

  double degrees() const {
    if      (type_ == Type::DEGREES) return angle_;
    else if (type_ == Type::RADIANS) return Rad2Deg(angle_);
    else    assert(false);
  }

  double radians() const {
    if      (type_ == Type::RADIANS) return angle_;
    else if (type_ == Type::DEGREES) return Deg2Rad(angle_);
    else    assert(false);
  }

  void setDegrees(double a) {
    type_ = Type::DEGREES; angle_ = a;
  }

  void setRadians(double a) {
    type_ = Type::RADIANS; angle_ = a;
  }

  std::string toString() const {
    std::ostringstream ostr;

    ostr << degrees();

    return ostr.str();
  }

  void fromString(const std::string &str) {
    std::istringstream istr(str);

    double r;

    istr >> r;

    type_  = Type::DEGREES;
    angle_ = r;
  }

  void print(std::ostream &os) const {
    os << degrees();
  }

  friend std::ostream &operator<<(std::ostream &os, const CAngle &a) {
    a.print(os);

    return os;
  }

 private:
  Type   type_  { Type::DEGREES };
  double angle_ { 0.0 };
};

#endif
