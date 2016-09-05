#ifndef CMATRIX_STACK_2D_H
#define CMATRIX_STACK_2D_H

#include <CMatrix2D.h>
#include <CPoint2D.h>
#include <vector>
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

template<typename T>
class CMatrixStack2DT {
 public:
  typedef CMatrix2DT<T> Matrix;
  typedef CPoint2DT<T>  Point;

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

    Transform(double a, const Point &p) :
     type_(CMatrixTransformType::ROTATE_ORIGIN) {
      memset(&v_[0], 0, NUM_VALUES*sizeof(double));

      v_[0] = a;
      v_[1] = p.x;
      v_[2] = p.y;
    }

    Transform(const Matrix &m) :
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

    Matrix rotateMatrix() const {
      Matrix m1 = Matrix::translation(-v_[1], -v_[2]);
      Matrix m2 = Matrix::rotation   ( v_[0]);
      Matrix m3 = Matrix::translation( v_[1],  v_[2]);

      return m3*m2*m1;
    }

    Matrix valueMatrix() const {
      Matrix m;

      m.setValues(v_[0], v_[1], v_[2], v_[3], v_[4], v_[5]);

      return m;
    }

    Matrix calcMatrix() const {
      switch (type_) {
        case CMatrixTransformType::TRANSLATE    : return Matrix::translation(v_[0], v_[1]);
        case CMatrixTransformType::SCALE1       : return Matrix::scale      (v_[0], v_[0]);
        case CMatrixTransformType::SCALE2       : return Matrix::scale      (v_[0], v_[1]);
        case CMatrixTransformType::ROTATE       : return Matrix::rotation   (v_[0]    );
        case CMatrixTransformType::ROTATE_ORIGIN: return rotateMatrix       ();
        case CMatrixTransformType::SKEWX        : return Matrix::skewX      (v_[0]);
        case CMatrixTransformType::SKEWY        : return Matrix::skewY      (v_[0]);
        case CMatrixTransformType::MATRIX       : return valueMatrix        ();
        default                                 : assert(false); return Matrix();
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
  CMatrixStack2DT() :
   transformStack_(), mValid_(false), m_() {
  }

  CMatrixStack2DT(const CMatrixStack2DT &m) :
   transformStack_(m.transformStack_), mValid_(m.mValid_), m_(m.m_) {
  }

  CMatrixStack2DT(const Matrix &m) :
   transformStack_() {
    transformStack_.push_back(Transform(m));

    mValid_ = false;
  }

  const CMatrixStack2DT &operator=(const CMatrixStack2DT &m) {
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

  void translate(const Point &d) {
    transformStack_.push_back(translateTransform(d));

    mValid_ = false;
  }

  void translate(double dx, double dy) {
    transformStack_.push_back(translateTransform(dx, dy));

    mValid_ = false;
  }

  static Transform scaleTransform(double s) {
    return Transform(CMatrixTransformType::SCALE1, s);
  }

  static Transform scaleTransform(double sx, double sy) {
    return Transform(CMatrixTransformType::SCALE2, sx, sy);
  }

  void scale(double s) {
    transformStack_.push_back(scaleTransform(s));

    mValid_ = false;
  }

  void scale(double sx, double sy) {
    transformStack_.push_back(scaleTransform(sx, sy));

    mValid_ = false;
  }

  static Transform rotateTransform(double a, const CPoint2D &o) {
    return Transform(a, o);
  }

  static Transform rotateTransform(double a) {
    return Transform(CMatrixTransformType::ROTATE, a);
  }

  void rotate(double a) {
    transformStack_.push_back(rotateTransform(a));

    mValid_ = false;
  }

  void rotate(double a, const Point &o) {
    transformStack_.push_back(rotateTransform(a, o));

    mValid_ = false;
  }

  void skewX(double a) {
    transformStack_.push_back(Transform(CMatrixTransformType::SKEWX, a));

    mValid_ = false;
  }

  void skewY(double a) {
    transformStack_.push_back(Transform(CMatrixTransformType::SKEWY, a));

    mValid_ = false;
  }

  void matrix(double m00, double m01, double m10, double m11, double tx, double ty) {
    Matrix m;

    m.setValues(m00, m01, m10, m11, tx, ty);

    transformStack_.push_back(Transform(m));

    mValid_ = false;
  }

  void matrix(const Matrix &m) {
    transformStack_.push_back(Transform(m));

    mValid_ = false;
  }

  void reset() {
    transformStack_.clear();

    mValid_ = false;
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

  void append(const CMatrixStack2DT &m) {
    for (const auto &t : m.transformStack_)
      transformStack_.push_back(t);

    mValid_ = false;
  }

  const Matrix &getMatrix() const {
    if (! mValid_) {
      CMatrixStack2DT *th = const_cast<CMatrixStack2DT *>(this);

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

  void multiplyPoint(const Point &point1, Point &point2) const {
    getMatrix().multiplyPoint(point1, point2);
  }

  void multiplyPoint(T xi, T yi, T *xo, T *yo) const {
    getMatrix().multiplyPoint(xi, yi, xo, yo);
  }

  void preMultiplyPoint(const Point &point1, Point &point2) const {
    getMatrix().preMultiplyPoint(point1, point2);
  }

  void preMultiplyPoint(T xi, T yi, T *xo, T *yo) const {
    getMatrix().preMultiplyPoint(xi, yi, xo, yo);
  }

  friend std::ostream &operator<<(std::ostream &os, const CMatrixStack2DT &rhs) {
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
  Matrix         m_;
};

typedef CMatrixStack2DT<double> CMatrixStack2D;

#endif
