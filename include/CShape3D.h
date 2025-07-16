#ifndef CSHAPE_3D_H
#define CSHAPE_3D_H

#include <CBBox3D.h>
#include <CMatrix3D.h>
#include <CPoint3D.h>
#include <CVector3D.h>

#include <optional>

class CLine3D;

class CShape3D {
 public:
  CShape3D() { }

  virtual ~CShape3D() { }

  virtual CBBox3D getBBox() const = 0;

  void translate(const CPoint3D &d) {
    translate_ += d;

    transform_  = OptMatrix();
    itransform_ = OptMatrix();
  }

  void scale(const CPoint3D &s) {
    scale_ *= s;

    transform_  = OptMatrix();
    itransform_ = OptMatrix();
  }

  void rotate(const CPoint3D &r) {
    rotate_ += r;

    transform_  = OptMatrix();
    itransform_ = OptMatrix();
  }

  const CPoint3D &getTranslate() const { return translate_; }
  const CPoint3D &getScale    () const { return scale_    ; }
  const CPoint3D &getRotate   () const { return rotate_   ; }

  CPoint3D transformFrom(const CPoint3D &point) const {
    return getTransform().multiplyPoint(point);
  }

  CVector3D transformFrom(const CVector3D &vector) const {
    return getTransform().multiplyVector(vector);
  }

  CPoint3D transformTo(const CPoint3D &point) const {
    return getITransform().multiplyPoint(point);
  }

  CVector3D transformTo(const CVector3D &vector) const {
    return getITransform().multiplyVector(vector);
  }

  const CMatrix3D &getTransform() const {
    if (! transform_) {
      auto *th = const_cast<CShape3D *>(this);

      th->transform_ = CMatrix3D::translation(translate_.x, translate_.y, translate_.z)*
                       CMatrix3D::rotation   (rotate_   .x, rotate_   .y, rotate_   .z)*
                       CMatrix3D::scale      (scale_    .x, scale_    .y, scale_    .z);
    }

    return transform_.value();
  }

  const CMatrix3D &getITransform() const {
    if (! itransform_) {
      const auto &transform = getTransform();

      auto *th = const_cast<CShape3D *>(this);

      th->itransform_ = transform.inverse();
    }

    return itransform_.value();
  }

  virtual bool intersect(const CLine3D &, double *, double *) const { return false; }

 protected:
  struct TRange {
    bool   set  { false };
    double tmin { 0.0 };
    double tmax { 0.0 };

    TRange() { }

    bool isOutside(double t) {
      return (! set || (t < tmin || t > tmax));
    }

    void update(double t) {
      if (! set) {
        tmin = t;
        tmax = t;

        set = true;
      }
      else {
        tmin = std::min(t, tmin);
        tmax = std::max(t, tmax);
      }
    }
  };

  using OptMatrix = std::optional<CMatrix3D>;

  CPoint3D  translate_ { 0.0, 0.0, 0.0 }; //! Translate
  CPoint3D  scale_     { 1.0, 1.0, 1.0 }; //! Scale
  CPoint3D  rotate_    { 0.0, 0.0, 0.0 }; //! Rotate
  OptMatrix transform_;                   //! Transform
  OptMatrix itransform_;                  //! Inverse Transform
};

#endif
