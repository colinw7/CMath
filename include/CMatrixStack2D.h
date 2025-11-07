#ifndef CMATRIX_STACK_2D_H
#define CMATRIX_STACK_2D_H

#include <CMatrix2D.h>
#include <CPoint2D.h>
#include <vector>
#include <sstream>
#include <cassert>

enum class CMatrixTransformType {
  NONE,
  TRANSLATE,
  SCALE1,
  SCALE2,
  ROTATE,
  ROTATE_ORIGIN,
  SKEWX,
  SKEWY,
  MATRIX
};

class CMatrixStack2D {
 public:
  class Transform {
   public:
    enum { NUM_VALUES = 6 };

   public:
    Transform(CMatrixTransformType type, double v) :
     type_(type) {
      memset(&v_[0], 0, NUM_VALUES*sizeof(double));

      v_[0] = v;
    }

    Transform(CMatrixTransformType type, double v1, double v2) :
     type_(type) {
      memset(&v_[0], 0, NUM_VALUES*sizeof(double));

      v_[0] = v1;
      v_[1] = v2;
    }

    Transform(double a, const CPoint2D &p) :
     type_(CMatrixTransformType::ROTATE_ORIGIN) {
      memset(&v_[0], 0, NUM_VALUES*sizeof(double));

      v_[0] = a;
      v_[1] = p.x;
      v_[2] = p.y;
    }

    Transform(const CMatrix2D &m) :
     type_(CMatrixTransformType::MATRIX) {
      m.getValues(v_, 6);
    }

    CMatrixTransformType type() const { return type_; }

    //---

    double angle() const {
      assert(type_ == CMatrixTransformType::ROTATE ||
             type_ == CMatrixTransformType::ROTATE_ORIGIN ||
             type_ == CMatrixTransformType::SKEWX ||
             type_ == CMatrixTransformType::SKEWY);

      return value(0);
    }

    //---

    double xscale() const {
      assert(type_ == CMatrixTransformType::SCALE1 || type_ == CMatrixTransformType::SCALE2);
      return value(0);
    }

    double yscale() const {
      assert(type_ == CMatrixTransformType::SCALE1 || type_ == CMatrixTransformType::SCALE2);
      return value(1);
    }

    //---

    double dx() const {
      assert(type_ == CMatrixTransformType::TRANSLATE);
      return value(0);
    }

    double dy() const {
      assert(type_ == CMatrixTransformType::TRANSLATE);
      return value(1);
    }

    //---

    double xo() const {
      assert(type_ == CMatrixTransformType::ROTATE_ORIGIN);
      return value(1);
    }

    double yo() const {
      assert(type_ == CMatrixTransformType::ROTATE_ORIGIN);
      return value(2);
    }

    //---

    double value(int i) const { return v_[i]; }
    void setValue(int i, double v) { v_[i] = v; }

    const double *values() const { return &v_[0]; }

    //---

    void setDx(double v) {
      assert(type_ == CMatrixTransformType::TRANSLATE);
      setValue(0, v);
    }

    void setDy(double v) {
      assert(type_ == CMatrixTransformType::TRANSLATE);
      setValue(1, v);
    }

    void setXScale(double v) {
      assert(type_ == CMatrixTransformType::SCALE1 || type_ == CMatrixTransformType::SCALE2);
      setValue(0, v);
    }

    void setYScale(double v) {
      assert(type_ == CMatrixTransformType::SCALE1 || type_ == CMatrixTransformType::SCALE2);
      setValue(1, v);
    }

    void setAngle(double v) {
      assert(type_ == CMatrixTransformType::ROTATE ||
             type_ == CMatrixTransformType::ROTATE_ORIGIN ||
             type_ == CMatrixTransformType::SKEWX ||
             type_ == CMatrixTransformType::SKEWY);

      setValue(0, v);
    }

    //---

    CMatrix2D rotateMatrix() const {
      auto m1 = CMatrix2D::translation(-v_[1], -v_[2]);
      auto m2 = CMatrix2D::rotation   ( v_[0]);
      auto m3 = CMatrix2D::translation( v_[1],  v_[2]);

      return m3*m2*m1;
    }

    CMatrix2D valueMatrix() const {
      CMatrix2D m;

      m.setValues(v_[0], v_[1], v_[2], v_[3], v_[4], v_[5]);

      return m;
    }

    CMatrix2D calcMatrix() const {
      switch (type_) {
        case CMatrixTransformType::TRANSLATE    : return CMatrix2D::translation(v_[0], v_[1]);
        case CMatrixTransformType::SCALE1       : return CMatrix2D::scale      (v_[0], v_[0]);
        case CMatrixTransformType::SCALE2       : return CMatrix2D::scale      (v_[0], v_[1]);
        case CMatrixTransformType::ROTATE       : return CMatrix2D::rotation   (v_[0]    );
        case CMatrixTransformType::ROTATE_ORIGIN: return rotateMatrix          ();
        case CMatrixTransformType::SKEWX        : return CMatrix2D::skewX      (v_[0]);
        case CMatrixTransformType::SKEWY        : return CMatrix2D::skewY      (v_[0]);
        case CMatrixTransformType::MATRIX       : return valueMatrix           ();
        default                                 : assert(false); return CMatrix2D();
      }
    }

    //---

    std::string name() const {
      switch (type_) {
        case CMatrixTransformType::TRANSLATE    : return "translate";
        case CMatrixTransformType::SCALE1       : return "scale";
        case CMatrixTransformType::SCALE2       : return "scale";
        case CMatrixTransformType::ROTATE       : return "rotate";
        case CMatrixTransformType::ROTATE_ORIGIN: return "rotate";
        case CMatrixTransformType::SKEWX        : return "skewX";
        case CMatrixTransformType::SKEWY        : return "skewY";
        case CMatrixTransformType::MATRIX       : return "matrix";
        default                                 : assert(false); return "";
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
        case CMatrixTransformType::TRANSLATE    : printValues(os, 2); return;
        case CMatrixTransformType::SCALE1       : printValues(os, 1); return;
        case CMatrixTransformType::SCALE2       : printValues(os, 2); return;
        case CMatrixTransformType::ROTATE       : printValues(os, 1); return;
        case CMatrixTransformType::ROTATE_ORIGIN: printValues(os, 3); return;
        case CMatrixTransformType::SKEWX        : printValues(os, 1); return;
        case CMatrixTransformType::SKEWY        : printValues(os, 1); return;
        case CMatrixTransformType::MATRIX       : printValues(os, 6); return;
        default                                 : assert(false)     ; return;
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
    CMatrixTransformType type_;
    double               v_[NUM_VALUES];
  };

  using TransformStack = std::vector<Transform>;

  //---

 public:
  static Transform translateTransform(const CPoint2D &d) {
    return Transform(CMatrixTransformType::TRANSLATE, d.x, d.y);
  }

  static Transform translateTransform(double dx, double dy) {
    return Transform(CMatrixTransformType::TRANSLATE, dx, dy);
  }

  static CMatrixStack2D translation(double dx, double dy) {
    CMatrixStack2D m;
    m.addTranslation(dx, dy);
    return m;
  }

  static CMatrixStack2D scale(double sx, double sy) {
    CMatrixStack2D m;
    m.addScale(sx, sy);
    return m;
  }

  static CMatrixStack2D rotation(double a) {
    CMatrixStack2D m;
    m.addRotation(a);
    return m;
  }

  //---

  // default is identity
  CMatrixStack2D() { }

  CMatrixStack2D(const CMatrix2D &m) {
    transformStack_.push_back(Transform(m));
  }

  //---

  const TransformStack &transformStack() const { return transformStack_; }

  CMatrixStack2D &addTranslation(const CPoint2D &d) {
    transformStack_.push_back(translateTransform(d));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &addTranslation(double dx, double dy) {
    transformStack_.push_back(translateTransform(dx, dy));

    mValid_ = false;

    return *this;
  }

  //---

  static Transform scaleTransform(double s) {
    return Transform(CMatrixTransformType::SCALE1, s);
  }

  static Transform scaleTransform(double sx, double sy) {
    return Transform(CMatrixTransformType::SCALE2, sx, sy);
  }

  CMatrixStack2D &addScale(double s) {
    transformStack_.push_back(scaleTransform(s));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &addScale(double sx, double sy) {
    transformStack_.push_back(scaleTransform(sx, sy));

    mValid_ = false;

    return *this;
  }

  //---

  static Transform rotateTransform(double a, const CPoint2D &o) {
    return Transform(a, o);
  }

  static Transform rotateTransform(double a) {
    return Transform(CMatrixTransformType::ROTATE, a);
  }

  CMatrixStack2D &addRotation(double a) {
    transformStack_.push_back(rotateTransform(a));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &addRotation(double a, const CPoint2D &o) {
    transformStack_.push_back(rotateTransform(a, o));

    mValid_ = false;

    return *this;
  }

  //---

  CMatrixStack2D &addSkewX(double a) {
    transformStack_.push_back(Transform(CMatrixTransformType::SKEWX, a));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &addSkewY(double a) {
    transformStack_.push_back(Transform(CMatrixTransformType::SKEWY, a));

    mValid_ = false;

    return *this;
  }

  //---

  CMatrixStack2D &matrix(double m00, double m01, double m10, double m11, double tx, double ty) {
    CMatrix2D m;

    m.setValues(m00, m01, m10, m11, tx, ty);

    transformStack_.push_back(Transform(m));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &matrix(const CMatrix2D &m) {
    transformStack_.push_back(Transform(m));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &reset() {
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

  CMatrixStack2D &prepend(const CMatrixStack2D &m) {
    auto transformStack = m.transformStack_;

    for (const auto &t : transformStack_)
      transformStack.push_back(t);

    transformStack_ = transformStack;

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &append(const CMatrixStack2D &m) {
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

  const CMatrix2D &getMatrix() const {
    if (! mValid_) {
      auto *th = const_cast<CMatrixStack2D *>(this);

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

  CMatrix2D getIMatrix() const {
    return getMatrix().inverse();
  }

  //---

  void multiplyPoint(const CPoint2D &point1, CPoint2D &point2) const {
    getMatrix().multiplyPoint(point1, point2);
  }

  void multiplyPoint(double xi, double yi, double *xo, double *yo) const {
    getMatrix().multiplyPoint(xi, yi, xo, yo);
  }

  void preMultiplyPoint(const CPoint2D &point1, CPoint2D &point2) const {
    getMatrix().preMultiplyPoint(point1, point2);
  }

  void preMultiplyPoint(double xi, double yi, double *xo, double *yo) const {
    getMatrix().preMultiplyPoint(xi, yi, xo, yo);
  }

  //---

  friend std::ostream &operator<<(std::ostream &os, const CMatrixStack2D &rhs) {
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
  CMatrix2D      m_;
};

#endif
