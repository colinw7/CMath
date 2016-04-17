#ifndef CMATRIX_STACK_2D_H
#define CMATRIX_STACK_2D_H

#include <cassert>
#include <CMatrix2D.h>

template<typename T>
class CMatrixStack2DT {
 public:
  typedef CMatrix2DT<T> Matrix;

 public:
  enum Type {
    NONE,
    TRANSLATE,
    SCALE,
    ROTATE
  };

  struct Transform {
    Type   type;
    double v1, v2;

    Transform(Type type_, double v1_, double v2_) :
     type(type_), v1(v1_), v2(v2_) {
    }
  };

 public:
  CMatrixStack2DT() :
   transformStack_(), mValid_(false), m_() {
  }

  void translate(double dx, double dy) {
    transformStack_.push_back(Transform(TRANSLATE, dx, dy));

    mValid_ = false;
  }

  void scale(double sx, double sy) {
    transformStack_.push_back(Transform(SCALE, sx, sy));

    mValid_ = false;
  }

  void rotate(double a) {
    transformStack_.push_back(Transform(ROTATE, a, a));

    mValid_ = false;
  }

  const Matrix &getMatrix() {
    if (! mValid_) {
      CMatrixStack2DT *th = const_cast<CMatrixStack2DT *>(this);

      th->m_.setIdentity();

      uint num = transformStack_.size();

      for (uint i = 0; i < num; ++i) {
        const Transform &t = transformStack_[i];

        switch (t.type) {
          case TRANSLATE: th->m_ *= Matrix::translation(t.v1, t.v2); break;
          case SCALE    : th->m_ *= Matrix::scale      (t.v1, t.v2); break;
          case ROTATE   : th->m_ *= Matrix::rotation   (t.v1      ); break;
          default       : assert(false); break;
        }
      }

      mValid_ = true;
    }

    return m_;
  }

 private:
  typedef std::vector<Transform> TransformStack;

  TransformStack transformStack_;
  bool           mValid_;
  CMatrix2D      m_;
};

typedef CMatrixStack2DT<double> CMatrixStack2D;

#endif
