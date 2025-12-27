#ifndef CMATRIX_STACK_3D_H
#define CMATRIX_STACK_3D_H

#include <CMatrix3D.h>
#include <CPoint3D.h>
#include <vector>
#include <sstream>
#include <cassert>

enum class CMatrix3DTransformType {
  NONE,
  TRANSLATE,
  SCALE,
  ROTATE,
  MATRIX
};

class CMatrixStack3D {
 public:
  class Transform {
   public:
    enum { NUM_VALUES = 12 };

   public:
    Transform(CMatrix3DTransformType type, double v) :
     type_(type) {
      memset(&v_[0], 0, NUM_VALUES*sizeof(double));

      v_[0] = v;
    }

    Transform(CMatrix3DTransformType type, double v1, double v2, double v3) :
     type_(type) {
      memset(&v_[0], 0, NUM_VALUES*sizeof(double));

      v_[0] = v1;
      v_[1] = v2;
      v_[2] = v3;
    }

    Transform(const CMatrix3D &m) :
     type_(CMatrix3DTransformType::MATRIX) {
      m.getValues(v_, 6);
    }

    CMatrix3DTransformType type() const { return type_; }

    //---

    double xangle() const {
      assert(type_ == CMatrix3DTransformType::ROTATE);
      return value(0);
    }

    double yangle() const {
      assert(type_ == CMatrix3DTransformType::ROTATE);
      return value(1);
    }

    double zangle() const {
      assert(type_ == CMatrix3DTransformType::ROTATE);
      return value(2);
    }

    //---

    double xscale() const {
      assert(type_ == CMatrix3DTransformType::SCALE);
      return value(0);
    }

    double yscale() const {
      assert(type_ == CMatrix3DTransformType::SCALE);
      return value(1);
    }

    double zscale() const {
      assert(type_ == CMatrix3DTransformType::SCALE);
      return value(2);
    }

    //---

    double dx() const {
      assert(type_ == CMatrix3DTransformType::TRANSLATE);
      return value(0);
    }

    double dy() const {
      assert(type_ == CMatrix3DTransformType::TRANSLATE);
      return value(1);
    }

    double dz() const {
      assert(type_ == CMatrix3DTransformType::TRANSLATE);
      return value(2);
    }

    //---

    double value(int i) const { return v_[i]; }
    void setValue(int i, double v) { v_[i] = v; }

    const double *values() const { return &v_[0]; }

    //---

    void setDx(double v) {
      assert(type_ == CMatrix3DTransformType::TRANSLATE);
      setValue(0, v);
    }

    void setDy(double v) {
      assert(type_ == CMatrix3DTransformType::TRANSLATE);
      setValue(1, v);
    }

    void setXScale(double v) {
      assert(type_ == CMatrix3DTransformType::SCALE);
      setValue(0, v);
    }

    void setYScale(double v) {
      assert(type_ == CMatrix3DTransformType::SCALE);
      setValue(1, v);
    }

    //---

    void setXAngle(double v) {
      assert(type_ == CMatrix3DTransformType::ROTATE);

      setValue(0, v);
    }

    void setYAngle(double v) {
      assert(type_ == CMatrix3DTransformType::ROTATE);

      setValue(1, v);
    }

    void setZAngle(double v) {
      assert(type_ == CMatrix3DTransformType::ROTATE);

      setValue(2, v);
    }

    //---

    CMatrix3D valueMatrix() const {
      CMatrix3D m;

      m.setValues(v_[0], v_[ 1], v_[ 2],
                  v_[3], v_[ 4], v_[ 5],
                  v_[6], v_[ 7], v_[ 8],
                  v_[9], v_[10], v_[11]);

      return m;
    }

    CMatrix3D calcMatrix() const {
      switch (type_) {
        case CMatrix3DTransformType::TRANSLATE: return CMatrix3D::translation(v_[0], v_[1], v_[2]);
        case CMatrix3DTransformType::SCALE    : return CMatrix3D::scale      (v_[0], v_[1], v_[2]);
        case CMatrix3DTransformType::ROTATE   : return CMatrix3D::rotation   (v_[0], v_[1], v_[2]);
        case CMatrix3DTransformType::MATRIX   : return valueMatrix           ();
        default                               : assert(false); return CMatrix3D();
      }
    }

    //---

    std::string name() const {
      switch (type_) {
        case CMatrix3DTransformType::TRANSLATE: return "translate";
        case CMatrix3DTransformType::SCALE    : return "scale";
        case CMatrix3DTransformType::ROTATE   : return "rotate";
        case CMatrix3DTransformType::MATRIX   : return "matrix";
        default                               : assert(false); return "";
      }
    }

    //---

    void printValues(std::ostream &os, int n) const {
      for (int i = 0; i < n; ++i) {
        if (i > 0) os << ",";

        os << v_[i];
      }
    }

    void printParts(std::ostream &os) const {
      switch (type_) {
        case CMatrix3DTransformType::TRANSLATE: printValues(os, 3); return;
        case CMatrix3DTransformType::SCALE    : printValues(os, 3); return;
        case CMatrix3DTransformType::ROTATE   : printValues(os, 4); return;
        case CMatrix3DTransformType::MATRIX   : printValues(os, 6); return;
        default                               : assert(false)     ; return;
      }
    }

    void printMatrix(std::ostream &os) const {
      os << calcMatrix();
    }

    friend std::ostream &operator<<(std::ostream &os, const Transform &rhs) {
      os << rhs.name() << "("; rhs.printParts(os); os << ")";

      return os;
    }

   private:
    CMatrix3DTransformType type_;
    double               v_[NUM_VALUES];
  };

  using TransformStack = std::vector<Transform>;

  //---

 public:
  static Transform translateTransform(const CPoint3D &d) {
    return Transform(CMatrix3DTransformType::TRANSLATE, d.x, d.y, d.z);
  }

  static Transform translateTransform(double dx, double dy, double dz) {
    return Transform(CMatrix3DTransformType::TRANSLATE, dx, dy, dz);
  }

  static Transform scaleTransform(double s) {
    return Transform(CMatrix3DTransformType::SCALE, s, s, s);
  }

  static Transform scaleTransform(double sx, double sy, double sz) {
    return Transform(CMatrix3DTransformType::SCALE, sx, sy, sz);
  }

  static Transform rotateTransform(double ax, double ay, double az) {
    return Transform(CMatrix3DTransformType::ROTATE, ax, ay, az);
  }

  //---

  static CMatrixStack3D translation(double dx, double dy, double dz) {
    CMatrixStack3D m;
    m.addTranslation(dx, dy, dz);
    return m;
  }

  static CMatrixStack3D scale(double sx, double sy, double sz) {
    CMatrixStack3D m;
    m.addScale(sx, sy, sz);
    return m;
  }

  static CMatrixStack3D rotation(double ax, double ay, double az) {
    CMatrixStack3D m;
    m.addRotation(ax, ay, az);
    return m;
  }

  //---

  // default is identity
  CMatrixStack3D() { }

  CMatrixStack3D(const CMatrix3D &m) {
    transformStack_.push_back(Transform(m));
  }

  //---

  const TransformStack &transformStack() const { return transformStack_; }

  CMatrixStack3D &addTranslation(const CPoint3D &d) {
    transformStack_.push_back(translateTransform(d));

    mValid_ = false;

    return *this;
  }

  CMatrixStack3D &addTranslation(double dx, double dy, double dz) {
    transformStack_.push_back(translateTransform(dx, dy, dz));

    mValid_ = false;

    return *this;
  }

  //---

  CMatrixStack3D &addScale(double s) {
    transformStack_.push_back(scaleTransform(s, s, s));

    mValid_ = false;

    return *this;
  }

  CMatrixStack3D &addScale(double sx, double sy, double sz) {
    transformStack_.push_back(scaleTransform(sx, sy, sz));

    mValid_ = false;

    return *this;
  }

  //---

  CMatrixStack3D &addRotation(double ax, double ay, double az) {
    transformStack_.push_back(rotateTransform(ax, ay, az));

    mValid_ = false;

    return *this;
  }

  //---

  CMatrixStack3D &addMatrix(double m00, double m01, double m02,
                            double m10, double m11, double m12,
                            double m20, double m21, double m22,
                            double tx, double ty, double tz) {
    CMatrix3D m;

    m.setValues(m00, m01, m02, m10, m11, m12, m20, m21, m22, tx, ty, tz);

    transformStack_.push_back(Transform(m));

    mValid_ = false;

    return *this;
  }

  CMatrixStack3D &addMatrix(const CMatrix3D &m) {
    transformStack_.push_back(Transform(m));

    mValid_ = false;

    return *this;
  }

  //---

  CMatrixStack3D &reset() {
    transformStack_.clear();

    mValid_ = false;

    return *this;
  }

  bool isEmpty() const {
    return transformStack_.empty();
  }

  size_t length() const {
    return transformStack_.size();
  }

  //---

  const Transform &transform(uint i) const {
    return transformStack_[i];
  }

  void setTransform(uint i, const Transform &t) {
    assert(size_t(i) < length());

    transformStack_[i] = t;

    mValid_ = false;
  }

  //---

  CMatrixStack3D &prepend(const CMatrixStack3D &m) {
    auto transformStack = m.transformStack_;

    for (const auto &t : transformStack_)
      transformStack.push_back(t);

    transformStack_ = transformStack;

    mValid_ = false;

    return *this;
  }

  CMatrixStack3D &append(const CMatrixStack3D &m) {
    for (const auto &t : m.transformStack_)
      transformStack_.push_back(t);

    mValid_ = false;

    return *this;
  }

  //---

  bool moveUp(uint ind) {
    auto len = uint(length());

    assert(ind < len);

    if (ind == 0) return false;

    std::swap(transformStack_[ind], transformStack_[ind - 1]);

    mValid_ = false;

    return true;
  }

  bool moveDown(uint ind) {
    auto len = uint(length());

    assert(ind < len);

    if (ind == len - 1) return false;

    std::swap(transformStack_[ind], transformStack_[ind + 1]);

    mValid_ = false;

    return true;
  }

  void remove(uint ind) {
    auto len = uint(length());

    assert(ind < len);

    for (uint i = ind + 1; i < len; ++i)
      transformStack_[i - 1] = transformStack_[i];

    transformStack_.pop_back();

    mValid_ = false;
  }

  //---

  const CMatrix3D &getMatrix() const {
    if (! mValid_) {
      auto *th = const_cast<CMatrixStack3D *>(this);

      th->m_.setIdentity();

      auto num = uint(length());

      for (uint i = 0; i < num; ++i) {
        const Transform &t = transform(i);

        th->m_ *= t.calcMatrix();
      }

      th->mValid_ = true;
    }

    return m_;
  }

  CMatrix3D getIMatrix() const {
    return getMatrix().inverse();
  }

  //---

  void multiplyPoint(const CPoint3D &point1, CPoint3D &point2) const {
    getMatrix().multiplyPoint(point1, point2);
  }

  void multiplyPoint(double xi, double yi, double zi, double *xo, double *yo, double *zo) const {
    getMatrix().multiplyPoint(xi, yi, zi, xo, yo, zo);
  }

  void preMultiplyPoint(const CPoint3D &point1, CPoint3D &point2) const {
    getMatrix().preMultiplyPoint(point1, point2);
  }

  void preMultiplyPoint(double xi, double yi, double zi, double *xo, double *yo, double *zo) const {
    getMatrix().preMultiplyPoint(xi, yi, zi, xo, yo, zo);
  }

  //---

  friend std::ostream &operator<<(std::ostream &os, const CMatrixStack3D &rhs) {
    rhs.print(os);

    return os;
  }

  std::string toString() const {
    std::stringstream str;

    str << *this;

    return str.str();
  }

  void print(std::ostream &os) const {
    auto num = uint(length());

    for (uint i = 0; i < num; ++i) {
      const auto &t = transform(i);

      if (i > 0) os << " ";

      os << t.name() << "(";

      t.printParts(os);

      os << ")";
    }
  }

 private:
  TransformStack transformStack_;
  bool           mValid_ { false };
  CMatrix3D      m_;
};

#endif
