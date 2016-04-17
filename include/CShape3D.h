#ifndef CSHAPE_3D_H
#define CSHAPE_3D_H

#include <CBBox3D.h>
#include <CMatrix3D.h>
#include <CPoint3D.h>
#include <CVector3D.h>
#include <COptVal.h>

template<typename T>
class CShape3DT {
 public:
  typedef CBBox3DT<T> BBox;

 protected:
  struct TRange {
    bool set;
    T    tmin;
    T    tmax;

    TRange() { set = false; }

    bool isOutside(T t) {
      return (! set || (t < tmin || t > tmax));
    }

    void update(T t) {
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

  typedef CMatrix3DT<T> Matrix;
  typedef CPoint3DT<T>  Point;
  typedef CVector3DT<T> Vector;

  Point            translate_;  //! Translate
  Point            scale_    ;  //! Scale
  Point            rotate_   ;  //! Rotate
  COptValT<Matrix> transform_;  //! Transform
  COptValT<Matrix> itransform_; //! Inverse Transform

 public:
  CShape3DT() :
   translate_(0,0,0), scale_(1,1,1), rotate_(0,0,0), transform_(), itransform_() {
    transform_ .setInvalid();
    itransform_.setInvalid();
  }

  virtual ~CShape3DT() { }

  virtual BBox getBBox() const = 0;

  void translate(const Point &d) {
    translate_ += d;

    transform_ .setInvalid();
    itransform_.setInvalid();
  }

  void scale(const Point &s) {
    scale_ *= s;

    transform_ .setInvalid();
    itransform_.setInvalid();
  }

  void rotate(const Point &r) {
    rotate_ += r;

    transform_ .setInvalid();
    itransform_.setInvalid();
  }

  const Point &getTranslate() const { return translate_; }
  const Point &getScale    () const { return scale_    ; }
  const Point &getRotate   () const { return rotate_   ; }

  Point transformFrom(const Point &point) const {
    return getTransform().multiplyPoint(point);
  }

  Vector transformFrom(const Vector &vector) const {
    return getTransform().multiplyVector(vector);
  }

  Point transformTo(const Point &point) const {
    return getITransform().multiplyPoint(point);
  }

  Vector transformTo(const Vector &vector) const {
    return getITransform().multiplyVector(vector);
  }

  const Matrix &getTransform() const {
    if (! transform_.isValid()) {
      CShape3DT *th = const_cast<CShape3DT *>(this);

      th->transform_.setValue(Matrix::translation(translate_.x, translate_.y, translate_.z)*
                              Matrix::rotation   (rotate_   .x, rotate_   .y, rotate_   .z)*
                              Matrix::scale      (scale_    .x, scale_    .y, scale_    .z));
    }

    return transform_.getValue();
  }

  const Matrix &getITransform() const {
    if (! itransform_.isValid()) {
      const Matrix &transform = getTransform();

      CShape3DT *th = const_cast<CShape3DT *>(this);

      th->itransform_.setValue(transform.inverse());
    }

    return itransform_.getValue();
  }
};

typedef CShape3DT<double> CShape3D;

#endif
