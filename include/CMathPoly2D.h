#ifndef CMathPoly2D_H
#define CMathPoly2D_H

#include <CPolygon2D.h>

enum CBinaryOp {
  CBINARY_OP_OR,
  CBINARY_OP_AND,
  CBINARY_OP_XOR,
  CBINARY_OP_NOT
};

bool
CPolygonBinaryOp(const CPolygon2D &poly1, const CPolygon2D &poly2, CBinaryOp op,
                 std::vector<CPolygon2D> &rpolys);

#endif
