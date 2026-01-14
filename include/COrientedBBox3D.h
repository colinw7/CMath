#ifndef COrientedBBox3D_H
#define COrientedBBox3D_H

#include <CVector3D.h>

class COrientedBBox3D {
 public:
  COrientedBBox3D() { }

  const CVector3D &center() const { return center_; }
  CVector3D &center() { return center_; }

  const CVector3D *axes() const { return &axis_[0]; }
  CVector3D *axes() { return &axis_[0]; }

  const CVector3D &axis(int i) const { return axis_[i]; }
  CVector3D &axis(int i) { return axis_[i]; }

  const double *extents() const { return &extent_[0]; }
  double *extents() { return &extent_[0]; }

  const double &extent(int i) const { return extent_[i]; }
  double &extent(int i) { return extent_[i]; }

  //---

  static COrientedBBox3D calcOrientBox(const std::vector<CVector3D> &akPoint);

 private:
  CVector3D center_;
  CVector3D axis_[3];
  double    extent_[3];
};

class CEigen {
 public:
  CEigen(int iSize) {
    assert(iSize >= 2);

    iSize_ = iSize;

    aafMat_.resize(iSize_);

    for (int i = 0; i < iSize_; i++)
      aafMat_[i].resize(iSize_);

    afDiag_.resize(iSize_);
    afSubd_.resize(iSize_);
  }

  double &matrix(int iRow, int iCol) {
    return aafMat_[iRow][iCol];
  }

  double getEigenvector(int iRow, int iCol) const { return aafMat_[iRow][iCol]; }

  double getEigenvalue(int i) const { return afDiag_[i]; }

  void calcIncrSortEigenStuff3() {
    calcTridiagonal3(aafMat_, &afDiag_[0], &afSubd_[0]);

    calcQLAlgorithm(iSize_, &afDiag_[0], &afSubd_[0], aafMat_);

    calcIncreasingSort(iSize_, &afDiag_[0], aafMat_);
  }

 private:
  using Reals  = std::vector<double>;
  using Matrix = std::vector<Reals>;

 private:
  void calcTridiagonal3(Matrix &aafMat, double *afDiag, double *afSubd) {
    double fM00 = aafMat[0][0];
    double fM01 = aafMat[0][1];
    double fM02 = aafMat[0][2];
    double fM11 = aafMat[1][1];
    double fM12 = aafMat[1][2];
    double fM22 = aafMat[2][2];

    afDiag[0] = fM00;
    afSubd[2] = 0.0;

    if (fM02 != 0.0) {
      double fLength = std::sqrt(fM01*fM01+fM02*fM02);
      double fInvLength = 1.0/fLength;
      fM01 *= fInvLength;
      fM02 *= fInvLength;

      double fQ = 2.0*fM01*fM12+fM02*(fM22-fM11);
      afDiag[1] = fM11+fM02*fQ;
      afDiag[2] = fM22-fM02*fQ;

      afSubd[0] = fLength;
      afSubd[1] = fM12-fM01*fQ;

      aafMat[0][0] = 1.0; aafMat[0][1] = 0.0;  aafMat[0][2] = 0.0;
      aafMat[1][0] = 0.0; aafMat[1][1] = fM01; aafMat[1][2] = fM02;
      aafMat[2][0] = 0.0; aafMat[2][1] = fM02; aafMat[2][2] = -fM01;
    }
    else {
      afDiag[1] = fM11;
      afDiag[2] = fM22;

      afSubd[0] = fM01;
      afSubd[1] = fM12;

      aafMat[0][0] = 1.0; aafMat[0][1] = 0.0; aafMat[0][2] = 0.0;
      aafMat[1][0] = 0.0; aafMat[1][1] = 1.0; aafMat[1][2] = 0.0;
      aafMat[2][0] = 0.0; aafMat[2][1] = 0.0; aafMat[2][2] = 1.0;
    }
  }

  bool calcQLAlgorithm(int iSize, double *afDiag, double *afSubd, Matrix &aafMat) {
    const int iMaxIter = 32;

    for (int i0 = 0; i0 < iSize; i0++) {
      int i1;

      for (i1 = 0; i1 < iMaxIter; i1++) {
        int i2;
        for (i2 = i0; i2 <= iSize - 2; i2++) {
          double fTmp = std::abs(afDiag[i2]) + std::abs(afDiag[i2 + 1]);
          if (std::abs(afSubd[i2]) + fTmp == fTmp)
            break;
        }
        if (i2 == i0)
          break;

        double fG = (afDiag[i0 + 1] - afDiag[i0])/(2.0*afSubd[i0]);
        double fR = std::sqrt(fG*fG + 1.0);
        if (fG < 0.0)
          fG = afDiag[i2] - afDiag[i0] + afSubd[i0]/(fG - fR);
        else
          fG = afDiag[i2] - afDiag[i0] + afSubd[i0]/(fG + fR);

        double fSin = 1.0, fCos = 1.0, fP = 0.0;
        for (int i3 = i2 - 1; i3 >= i0; i3--) {
          double fF = fSin*afSubd[i3];
          double fB = fCos*afSubd[i3];

          if (std::abs(fF) >= std::abs(fG)) {
            fCos = fG/fF;
            fR = std::sqrt(fCos*fCos + 1.0);

            afSubd[i3 + 1] = fF*fR;

            fSin = 1.0/fR;
            fCos *= fSin;
          }
          else {
            fSin = fF/fG;
            fR = std::sqrt(fSin*fSin + 1.0);

            afSubd[i3 + 1] = fG*fR;

            fCos = 1.0/fR;
            fSin *= fCos;
          }

          fG = afDiag[i3 + 1] - fP;
          fR = (afDiag[i3] - fG)*fSin + 2.0*fB*fCos;
          fP = fSin*fR;

          afDiag[i3 + 1] = fG + fP;

          fG = fCos*fR - fB;

          for (int i4 = 0; i4 < iSize; i4++) {
            fF = aafMat[i4][i3 + 1];

            aafMat[i4][i3 + 1] = fSin*aafMat[i4][i3] + fCos*fF;
            aafMat[i4][i3    ] = fCos*aafMat[i4][i3] - fSin*fF;
          }
        }
        afDiag[i0] -= fP;
        afSubd[i0] = fG;
        afSubd[i2] = 0.0;
      }
      if (i1 == iMaxIter)
        return false;
    }

    return true;
  }

  void calcIncreasingSort(int iSize, double *afEigval, Matrix &aafEigvec) {
    // sort eigenvalues in increasing order, e[0] <= ... <= e[iSize - 1]
    for (int i0 = 0, i1; i0 <= iSize - 2; i0++) {
      // locate minimum eigenvalue
      i1 = i0;
      double fMin = afEigval[i1];
      int i2;
      for (i2 = i0 + 1; i2 < iSize; i2++) {
        if (afEigval[i2] < fMin) {
          i1 = i1;
          fMin = afEigval[i1];
        }
      }

      if (i1 != i0) {
        // swap eigenvalues
        afEigval[i1] = afEigval[i0];
        afEigval[i0] = fMin;

        // swap eigenvectors
        for (i2 = 0; i2 < iSize; i2++) {
          double fTmp = aafEigvec[i2][i0];
          aafEigvec[i2][i0] = aafEigvec[i2][i1];
          aafEigvec[i2][i1] = fTmp;
        }
      }
    }
  }

 private:
  int    iSize_  { 0 };
  Matrix aafMat_;
  Reals  afDiag_;
  Reals  afSubd_;
};

#endif
