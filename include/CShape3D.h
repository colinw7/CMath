#ifndef CSHAPE_3D_H
#define CSHAPE_3D_H

#include <CBBox3D.h>
#include <CMatrix3D.h>
#include <CPoint3D.h>
#include <CVector3D.h>
#include <COptVal.h>

class CShape3D {
 protected:
  struct TRange {
    bool   set;
    double tmin;
    double tmax;

    TRange() { set = false; }

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

  CPoint3D            translate_;  //! Translate
  CPoint3D            scale_    ;  //! Scale
  CPoint3D            rotate_   ;  //! Rotate
  COptValT<CMatrix3D> transform_;  //! Transform
  COptValT<CMatrix3D> itransform_; //! Inverse Transform

 public:
  CShape3D() :
   translate_(0,0,0), scale_(1,1,1), rotate_(0,0,0), transform_(), itransform_() {
    transform_ .setInvalid();
    itransform_.setInvalid();
  }

  virtual ~CShape3D() { }

  virtual CBBox3D getBBox() const = 0;

  void translate(const CPoint3D &d) {
    translate_ += d;

    transform_ .setInvalid();
    itransform_.setInvalid();
  }

  void scale(const CPoint3D &s) {
    scale_ *= s;

    transform_ .setInvalid();
    itransform_.setInvalid();
  }

  void rotate(const CPoint3D &r) {
    rotate_ += r;

    transform_ .setInvalid();
    itransform_.setInvalid();
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
    if (! transform_.isValid()) {
      CShape3D *th = const_cast<CShape3D *>(this);

      th->transform_.setValue(CMatrix3D::translation(translate_.x, translate_.y, translate_.z)*
                              CMatrix3D::rotation   (rotate_   .x, rotate_   .y, rotate_   .z)*
                              CMatrix3D::scale      (scale_    .x, scale_    .y, scale_    .z));
    }

    return transform_.getValue();
  }

  const CMatrix3D &getITransform() const {
    if (! itransform_.isValid()) {
      const CMatrix3D &transform = getTransform();

      CShape3D *th = const_cast<CShape3D *>(this);

      th->itransform_.setValue(transform.inverse());
    }

    return itransform_.getValue();
  }
};

#endif
