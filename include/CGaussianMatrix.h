#ifndef CGAUSSIAN_MATRIX_H
#define CGAUSSIAN_MATRIX_H

template<uint ROWS, uint COLS>
class CGaussianMatrix {
 public:
  CGaussianMatrix(double a, double b, double c, double d,
                  double e, double f, double g, double h,
                  double i, double j, double k, double l) {
    m[0][0] = a; m[0][1] = b; m[0][2] = c; m[0][3] = d;
    m[1][0] = e; m[1][1] = f; m[1][2] = g; m[1][3] = h;
    m[2][0] = i; m[2][1] = j; m[2][2] = k; m[2][3] = l;
  }

  CGaussianMatrix(double *d) {
    for (uint r = 0; r < ROWS; ++r)
      for (uint c = 0; c < COLS; ++c)
        m[r][c] = *d++;
  }

  void print(std::ostream &os) const {
    os << m[0][0] <<" "<< m[0][1] <<" "<< m[0][2] <<" "<< m[0][3] << std::endl;
    os << m[1][0] <<" "<< m[1][1] <<" "<< m[1][2] <<" "<< m[1][3] << std::endl;
    os << m[2][0] <<" "<< m[2][1] <<" "<< m[2][2] <<" "<< m[2][3] << std::endl;
  }

  friend std::ostream &operator<<(std::ostream &os, const CGaussianMatrix &matrix) {
    matrix.print(os);

    return os;
  }

  void swapRows(int r1, int r2) {
    for (uint c = 0; c < COLS; ++c)
      std::swap(m[r1][c], m[r2][c]);
  }

  void scaleRow(int r, double f) {
    for (uint c = 0; c < COLS; ++c)
      m[r][c] *= f;
  }

  bool solve(double *res) {
    //std::cout << *this << std::endl;

    // Gaussian elimination

    for (uint p = 0; p < ROWS; ++p) {
      // Find max row element for column

      uint   mx = p;
      double mv = m[p][p];
      double ma = fabs(mv);

      for (uint r = p + 1; r < ROWS; ++r) {
        double mv1 = m[r][p];
        double ma1 = fabs(mv1);

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

      double f = 1.0/mv;

      scaleRow(p, f);

      //std::cout << *this << std::endl;

      //----------

      // Zero out lower diagonal using current row

      for (uint r = p + 1; r < ROWS; ++r) {
        double f = (1.0*m[r][p])/m[p][p];

        for (uint c = 0; c < COLS; ++c)
          m[r][c] -= m[p][c]*f;
      }

      //std::cout << *this << std::endl;
    }

    //std::cout << "Back Subst" << std::endl;

    for (uint r = ROWS - 2; r >= 0; --r) {
      for (uint c = r + 1; c < COLS - 1; ++c) {
        m[r][COLS - 1] -= m[r][c]*m[c][COLS - 1];

        m[r][c] = 0.0;
      }
    }

    //std::cout << *this << std::endl;

    for (uint r = 0; r < ROWS; ++r)
      res[r] = m[r][COLS - 1];

    return true;
  }

 private:
  double m[ROWS][COLS];
};

#endif
