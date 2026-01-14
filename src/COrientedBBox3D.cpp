#include <COrientedBBox3D.h>

namespace COrientedBBox3DUtil {

void gaussPointsFit(const std::vector<CVector3D> &akPoint, CVector3D& rkCenter,
                    CVector3D akAxis[3], double afExtent[3]);

}

//----------------------------------------------------------------------------
COrientedBBox3D
COrientedBBox3D::
calcOrientBox(const std::vector<CVector3D> &akPoint)
{
  COrientedBBox3D kBox;
  COrientedBBox3DUtil::gaussPointsFit(akPoint, kBox.center(), kBox.axes(), kBox.extents());

  // Let C be the box center and let U0, U1, and U2 be the box axes.
  // Each input point is of the form X = C + y0*U0 + y1*U1 + y2*U2.
  // The // following code computes min(y0), max(y0), min(y1), max(y1), min(y2), and max(y2).
  // The box center is then adjusted to be
  //   C' = C + 0.5*(min(y0)+max(y0))*U0 + 0.5*(min(y1)+max(y1))*U1 + 0.5*(min(y2)+max(y2))*U2

  CVector3D kDiff = akPoint[0] - kBox.center();

  double fY0Min = kDiff.dotProduct(kBox.axis(0)), fY0Max = fY0Min;
  double fY1Min = kDiff.dotProduct(kBox.axis(1)), fY1Max = fY1Min;
  double fY2Min = kDiff.dotProduct(kBox.axis(2)), fY2Max = fY2Min;

  auto iQuantity = int(akPoint.size());

  for (int i = 1; i < iQuantity; i++) {
    kDiff = akPoint[i] - kBox.center();

    double fY0 = kDiff.dotProduct(kBox.axis(0));
    if      (fY0 < fY0Min) fY0Min = fY0;
    else if (fY0 > fY0Max) fY0Max = fY0;

    double fY1 = kDiff.dotProduct(kBox.axis(1));
    if      (fY1 < fY1Min) fY1Min = fY1;
    else if (fY1 > fY1Max) fY1Max = fY1;

    double fY2 = kDiff.dotProduct(kBox.axis(2));
    if      (fY2 < fY2Min) fY2Min = fY2;
    else if (fY2 > fY2Max) fY2Max = fY2;
  }

  kBox.center() += (0.5*(fY0Min + fY0Max))*kBox.axis(0) +
                   (0.5*(fY1Min + fY1Max))*kBox.axis(1) +
                   (0.5*(fY2Min + fY2Max))*kBox.axis(2);

  kBox.extent(0) = 0.5*(fY0Max - fY0Min);
  kBox.extent(1) = 0.5*(fY1Max - fY1Min);
  kBox.extent(2) = 0.5*(fY2Max - fY2Min);

  return kBox;
}

//---

namespace COrientedBBox3DUtil {

void
gaussPointsFit(const std::vector<CVector3D> &akPoint, CVector3D &rkCenter,
               CVector3D akAxis[3], double afExtent[3])
{
  auto iQuantity = int(akPoint.size());

  // compute mean of points
  rkCenter = akPoint[0];
  for (int i = 1; i < iQuantity; i++)
    rkCenter += akPoint[i];

  double fInvQuantity = 1.0/iQuantity;
  rkCenter *= fInvQuantity;

  // compute covariances of points
  double fSumXX = 0.0, fSumXY = 0.0, fSumXZ = 0.0;
  double fSumYY = 0.0, fSumYZ = 0.0, fSumZZ = 0.0;

  for (int i = 0; i < iQuantity; i++) {
    CVector3D kDiff = akPoint[i] - rkCenter;

    fSumXX += kDiff.x()*kDiff.x();
    fSumXY += kDiff.x()*kDiff.y();
    fSumXZ += kDiff.x()*kDiff.z();
    fSumYY += kDiff.y()*kDiff.y();
    fSumYZ += kDiff.y()*kDiff.z();
    fSumZZ += kDiff.z()*kDiff.z();
  }

  fSumXX *= fInvQuantity;
  fSumXY *= fInvQuantity;
  fSumXZ *= fInvQuantity;
  fSumYY *= fInvQuantity;
  fSumYZ *= fInvQuantity;
  fSumZZ *= fInvQuantity;

  // compute eigenvectors for covariance matrix
  CEigen kES(3);

  kES.matrix(0, 0) = fSumXX;
  kES.matrix(0, 1) = fSumXY;
  kES.matrix(0, 2) = fSumXZ;
  kES.matrix(1, 0) = fSumXY;
  kES.matrix(1, 1) = fSumYY;
  kES.matrix(1, 2) = fSumYZ;
  kES.matrix(2, 0) = fSumXZ;
  kES.matrix(2, 1) = fSumYZ;
  kES.matrix(2, 2) = fSumZZ;

  kES.calcIncrSortEigenStuff3();

  akAxis[0].setX(kES.getEigenvector(0, 0));
  akAxis[0].setY(kES.getEigenvector(1, 0));
  akAxis[0].setZ(kES.getEigenvector(2, 0));
  akAxis[1].setX(kES.getEigenvector(0, 1));
  akAxis[1].setY(kES.getEigenvector(1, 1));
  akAxis[1].setZ(kES.getEigenvector(2, 1));
  akAxis[2].setX(kES.getEigenvector(0, 2));
  akAxis[2].setY(kES.getEigenvector(1, 2));
  akAxis[2].setZ(kES.getEigenvector(2, 2));

  afExtent[0] = kES.getEigenvalue(0);
  afExtent[1] = kES.getEigenvalue(1);
  afExtent[2] = kES.getEigenvalue(2);
}

}
