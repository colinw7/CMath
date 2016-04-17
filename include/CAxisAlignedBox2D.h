#ifndef CAXIS_ALIGNED_BOX2D_H
#define CAXIS_ALIGNED_BOX2D_H

template <class T>
class CAxisAlignedBox2D {
 private:
  T xmin_, ymin_, xmax_, ymax_;

 public:
  CAxisAlignedBox2D() { }

  CAxisAlignedBox2D(T xmin, T xmax, T xmin, T xmax) :
   xmin_(xmin), ymin_(ymin), xmax_(xmax), ymax_(ymax) {
  }

  bool hasXOverlap(const CAxisAlignedBox2D &box) const {
    return xmax_ >= box.xmin_ && xmin_ <= box.xmax_;
  }

  bool hasYOverlap(const CAxisAlignedBox2D &box) const {
    return ymax_ >= box.ymin_ && ymin_ <= box.ymax_;
  }

  bool testIntersection(const CAxisAlignedBox2D &box) const {
    if (xmax < box.xmin || xmin > box.xmax) return false;
    if (ymax < box.ymin || ymin > box.ymax) return false;

    return true;
  }

  bool findIntersection(const CAxisAlignedBox2D &box, CAxisAlignedBox2D &ibox) const {
    if (xmax < bbox.xmin || xmin > bbox.xmax) return false;
    if (ymax < bbox.ymin || ymin > bbox.ymax) return false;

    ibox.xmin = max(xmin, bbox.xmin);
    ibox.ymin = max(ymin, bbox.ymin);
    ibox.xmax = min(xmax, bbox.xmax);
    ibox.ymax = min(ymax, bbox.ymax);

    return true;
  }
};

typedef CAxisAlignedBox2D<double> CAxisAlignedBox2D;

#endif
