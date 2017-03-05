#ifndef CAXIS_ALIGNED_BOX2D_H
#define CAXIS_ALIGNED_BOX2D_H

class CAxisAlignedBox2D {
 public:
  CAxisAlignedBox2D() { }

  CAxisAlignedBox2D(double xmin, double xmax, double xmin, double xmax) :
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

 private:
  double xmin_ { 0 }, ymin_ { 0 }, xmax_ { 0 }, ymax_ { 0 };

};

#endif
