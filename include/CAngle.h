#ifndef CAngle_H
#define CAngle_H

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>

// Angle type
//  . Value in degress, radians or gradians
//
//  degress->radians  = PI/180.0
//  degress->gradians = 100.0/90.0
class CAngle {
 public:
  enum class Type {
    DEGREES,
    RADIANS,
    GRADS
  };

 public:
  static CAngle makeDegrees(double r) { return CAngle(Type::DEGREES, r); }
  static CAngle makeRadians(double r) { return CAngle(Type::RADIANS, r); }
  static CAngle makeGrads  (double r) { return CAngle(Type::GRADS  , r); }

  //---

  // degrees to radians or gradians
  static double Deg2Rad (double d) { return M_PI*d/180.0; }
  static double Deg2Grad(double d) { return 100.0*d/90.0; }
  // radians to degrees or gradians
  static double Rad2Deg (double r) { return 180.0*r/M_PI; }
  static double Rad2Grad(double r) { return Deg2Grad(Rad2Deg(r)); }
  // gradians to degrees or radians
  static double Grad2Deg(double g) { return 90.0*g/100.0; }
  static double Grad2Rad(double g) { return Deg2Rad(Grad2Deg(g)); }

  //---

  explicit CAngle(double angle=0.0) :
   type_(Type::DEGREES), angle_(angle) {
  }

  explicit CAngle(Type type, double angle=0.0) :
   type_(type), angle_(angle) {
  }

  //---

  const Type &type() const { return type_; }

  double value() const { return angle_; }

  //---

  double degrees() const {
    if      (type_ == Type::DEGREES) return angle_;
    else if (type_ == Type::RADIANS) return Rad2Deg (angle_);
    else if (type_ == Type::GRADS  ) return Grad2Deg(angle_);
    else    assert(false);
  }

  double radians() const {
    if      (type_ == Type::RADIANS) return angle_;
    else if (type_ == Type::DEGREES) return Deg2Rad (angle_);
    else if (type_ == Type::GRADS  ) return Grad2Rad(angle_);
    else    assert(false);
  }

  double grads() const {
    if      (type_ == Type::GRADS  ) return angle_;
    else if (type_ == Type::DEGREES) return Deg2Grad(angle_);
    else if (type_ == Type::RADIANS) return Rad2Grad(angle_);
    else    assert(false);
  }

  //---

  void setDegrees(double a) {
    type_ = Type::DEGREES; angle_ = a;
  }

  void setRadians(double a) {
    type_ = Type::RADIANS; angle_ = a;
  }

  void setGrads(double a) {
    type_ = Type::GRADS; angle_ = a;
  }

  //---

  std::string toString() const {
    std::ostringstream ostr;

    // always degrees ?
    ostr << degrees();

    return ostr.str();
  }

  void fromString(const std::string &str) {
    std::istringstream istr(str);

    // handle postfix type ?
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

  //---

  friend bool operator==(const CAngle &lhs, const CAngle &rhs) {
    if (lhs.type_ != rhs.type_)
      return lhs.degrees() == rhs.degrees();
    else
      return lhs.angle_ == rhs.angle_;
  }

  friend bool operator<(const CAngle &lhs, const CAngle &rhs) {
    if (lhs.type_ != rhs.type_)
      return lhs.degrees() < rhs.degrees();
    else
      return lhs.angle_ < rhs.angle_;
  }

  //---

 private:
  Type   type_  { Type::DEGREES };
  double angle_ { 0.0 };
};

#endif
