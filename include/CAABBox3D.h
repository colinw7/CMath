#ifndef CAABBox_H
#define CAABBox_H

class CAABBox3D {
 public:
  using Points = std::vector<CPoint3D>;

 public:
  CAABBox3D() { }

  bool isSet() const { return set_; }

  void addPoint(const CPoint3D &p) {
    if (! set_) {
      pmin_ = point;
      pmax_ = point;

      set_ = true;
    }
    else {
      pmin_ = CPoint3D::min(pmin_, point);
      pmax_ = CPoint3D::max(pmax_, point);
    }
  }

  bool intersects(const CAABBox3D &other) const {
    // Overlap happens if there's an overlap on ALL three axes independently
    return (pmin_.x <= other.pmax_.x && pmax_.x >= other.pmin_.x) &&
           (pmin_.y <= other.pmax_.y && pmax_.y >= other.pmin_.y) &&
           (pmin_.z <= other.pmax_.z && pmax_.z >= other.pmin_.z);
  }

 private:
  bool     set_ { false };
  CPoint3D pmin_;
  CPoint3D pmax_;
};

#endif
