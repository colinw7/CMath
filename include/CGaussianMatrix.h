#ifndef CGAUSSIAN_MATRIX_H
#define CGAUSSIAN_MATRIX_H

template<typename T, uint ROWS, uint COLS>
class CGaussianMatrix {
 private:
  T m[ROWS][COLS];

 public:
  CGaussianMatrix(T a, T b, T c, T d,
                  T e, T f, T g, T h,
                  T i, T j, T k, T l) {
    m[0][0] = a; m[0][1] = b; m[0][2] = c; m[0][3] = d;
    m[1][0] = e; m[1][1] = f; m[1][2] = g; m[1][3] = h;
    m[2][0] = i; m[2][1] = j; m[2][2] = k; m[2][3] = l;
  }

  CGaussianMatrix(T *d) {
    for (int r = 0; r < ROWS; ++r)
      for (int c = 0; c < COLS; ++c)
        m[r][c] = *d++;
  }

  void print(std::ostream &os) const {
    os << m[0][0] <<" "<< m[0][1] <<" "<< m[0][2] <<" "<< m[0][3] << std::endl;
    os << m[1][0] <<" "<< m[1][1] <<" "<< m[1][2] <<" "<< m[1][3] << std::endl;
    os << m[2][0] <<" "<< m[2][1] <<" "<< m[2][2] <<" "<< m[2][3] << std::endl;
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const CGaussianMatrix &matrix) {
    matrix.print(os);

    return os;
  }

  void swapRows(int r1, int r2) {
    for (int c = 0; c < COLS; ++c)
      swap(m[r1][c], m[r2][c]);
  }

  void scaleRow(int r, T f) {
    for (int c = 0; c < COLS; ++c)
      m[r][c] *= f;
  }

  bool solve(T *res) {
    //std::cout << *this << std::endl;

    // Gaussian elimination

    for (int p = 0; p < ROWS; ++p) {
      // Find max row element for column

      int    mx = p;
      T mv = m[p][p];
      T ma = fabs(mv);

      for (int r = p + 1; r < ROWS; ++r) {
        T mv1 = m[r][p];
        T ma1 = fabs(mv1);

        if (ma1 > ma) {
          mx = r;
          mv = mv1;
          ma = ma1;
        }
      }

      // if zero insoluable

      if (ma < 1E-6) {
        std::cerr << "Insoluable" << std::endl;
        return false;
      }

      // Swap max row into current rows (if required)

      if (mx != p)
        swapRows(p, mx);

      //std::cout << *this << std::endl;

      //----------

      // Scale row so current diagonal element is 1

      T f = 1.0/mv;

      scaleRow(p, f);

      //std::cout << *this << std::endl;

      //----------

      // Zero out lower diagonal using current row

      for (int r = p + 1; r < ROWS; ++r) {
        T f = (1.0*m[r][p])/m[p][p];

        for (int c = 0; c < COLS; ++c)
          m[r][c] -= m[p][c]*f;
      }

      //std::cout << *this << std::endl;
    }

    //std::cout << "Back Subst" << std::endl;

    for (int r = ROWS - 2; r >= 0; --r) {
      for (int c = r + 1; c < COLS - 1; ++c) {
        m[r][COLS - 1] -= m[r][c]*m[c][COLS - 1];

        m[r][c] = 0.0;
      }
    }

    //std::cout << *this << std::endl;

    for (int r = 0; r < ROWS; ++r)
      res[r] = m[r][COLS - 1];

    return true;
  }
};

#endif
