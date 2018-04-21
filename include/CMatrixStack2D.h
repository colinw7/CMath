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

    double angle() const {
      assert(type_ == CMatrixTransformType::ROTATE ||
             type_ == CMatrixTransformType::ROTATE_ORIGIN ||
             type_ == CMatrixTransformType::SKEWX ||
             type_ == CMatrixTransformType::SKEWY);

      return value(0);
    }

    double xscale() const {
      assert(type_ == CMatrixTransformType::SCALE1 ||
             type_ == CMatrixTransformType::SCALE2);

      return value(0);
    }

    double yscale() const {
      assert(type_ == CMatrixTransformType::SCALE1 ||
             type_ == CMatrixTransformType::SCALE2);

      return value(1);
    }

    double dx() const {
      assert(type_ == CMatrixTransformType::TRANSLATE);

      return value(0);
    }

    double dy() const {
      assert(type_ == CMatrixTransformType::TRANSLATE);

      return value(1);
    }

    double xo() const {
      assert(type_ == CMatrixTransformType::ROTATE_ORIGIN);

      return value(1);
    }

    double yo() const {
      assert(type_ == CMatrixTransformType::ROTATE_ORIGIN);

      return value(2);
    }

    double value(int i) const { return v_[i]; }

    const double *values() const { return &v_[0]; }

    CMatrix2D rotateMatrix() const {
      CMatrix2D m1 = CMatrix2D::translation(-v_[1], -v_[2]);
      CMatrix2D m2 = CMatrix2D::rotation   ( v_[0]);
      CMatrix2D m3 = CMatrix2D::translation( v_[1],  v_[2]);

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

  typedef std::vector<Transform> TransformStack;

  //---

 public:
  CMatrixStack2D() :
   transformStack_(), mValid_(false), m_() {
  }

  CMatrixStack2D(const CMatrixStack2D &m) :
   transformStack_(m.transformStack_), mValid_(m.mValid_), m_(m.m_) {
  }

  CMatrixStack2D(const CMatrix2D &m) :
   transformStack_() {
    transformStack_.push_back(Transform(m));

    mValid_ = false;
  }

  const CMatrixStack2D &operator=(const CMatrixStack2D &m) {
    transformStack_ = m.transformStack_;
    mValid_         = m.mValid_;
    m_              = m.m_;

    return *this;
  }

  const TransformStack &transformStack() const { return transformStack_; }

  static Transform translateTransform(const CPoint2D &d) {
    return Transform(CMatrixTransformType::TRANSLATE, d.x, d.y);
  }

  static Transform translateTransform(double dx, double dy) {
    return Transform(CMatrixTransformType::TRANSLATE, dx, dy);
  }

  CMatrixStack2D &translate(const CPoint2D &d) {
    transformStack_.push_back(translateTransform(d));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &translate(double dx, double dy) {
    transformStack_.push_back(translateTransform(dx, dy));

    mValid_ = false;

    return *this;
  }

  static Transform scaleTransform(double s) {
    return Transform(CMatrixTransformType::SCALE1, s);
  }

  static Transform scaleTransform(double sx, double sy) {
    return Transform(CMatrixTransformType::SCALE2, sx, sy);
  }

  CMatrixStack2D &scale(double s) {
    transformStack_.push_back(scaleTransform(s));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &scale(double sx, double sy) {
    transformStack_.push_back(scaleTransform(sx, sy));

    mValid_ = false;

    return *this;
  }

  static Transform rotateTransform(double a, const CPoint2D &o) {
    return Transform(a, o);
  }

  static Transform rotateTransform(double a) {
    return Transform(CMatrixTransformType::ROTATE, a);
  }

  CMatrixStack2D &rotate(double a) {
    transformStack_.push_back(rotateTransform(a));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &rotate(double a, const CPoint2D &o) {
    transformStack_.push_back(rotateTransform(a, o));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &skewX(double a) {
    transformStack_.push_back(Transform(CMatrixTransformType::SKEWX, a));

    mValid_ = false;

    return *this;
  }

  CMatrixStack2D &skewY(double a) {
    transformStack_.push_back(Transform(CMatrixTransformType::SKEWY, a));

    mValid_ = false;

    return *this;
  }

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

  int length() const {
    return transformStack_.size();
  }

  const Transform &transform(int i) const {
    return transformStack_[i];
  }

  void setTransform(int i, const Transform &t) {
    assert(i >= 0 && i < int(transformStack_.size()));

    transformStack_[i] = t;
  }

  CMatrixStack2D &append(const CMatrixStack2D &m) {
    for (const auto &t : m.transformStack_)
      transformStack_.push_back(t);

    mValid_ = false;

    return *this;
  }

  const CMatrix2D &getMatrix() const {
    if (! mValid_) {
      CMatrixStack2D *th = const_cast<CMatrixStack2D *>(this);

      th->m_.setIdentity();

      uint num = transformStack_.size();

      for (uint i = 0; i < num; ++i) {
        const Transform &t = transformStack_[i];

        th->m_ *= t.calcMatrix();
      }

      th->mValid_ = true;
    }

    return m_;
  }

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
    uint num = transformStack_.size();

    for (uint i = 0; i < num; ++i) {
      const Transform &t = transformStack_[i];

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
