#ifndef CBezierPath_H
#define CBezierPath_H

#include <C3Bezier2D.h>
#include <CArcToBezier.h>

#include <algorithm>

class CBezierPath {
 public:
  using Beziers = std::vector<C3Bezier2D>;
  using Points  = std::vector<CPoint2D>;

 public:
  CBezierPath() { }

  CBezierPath(const Beziers &beziers) :
   beziers_(beziers) {
  }

  //---

  const Beziers &beziers() const { return beziers_; }
  void setBeziers(const Beziers &beziers) { beziers_ = beziers; invalidate(); }

  bool isClosed() const { return closed_; }
  void setClosed(bool b) { closed_ = b; invalidate(); }

  //---

  void clear() { beziers_.clear(); points_.clear(); closed_ = false; invalidate(); }

  bool isEmpty() const { return (arcLength() <= 0.0); }

  //---

  void moveTo(const CPoint2D &p) {
    points_.push_back(p);

    invalidate();
  }

  void lineTo(const CPoint2D &p) {
    beziers_.emplace_back(lastPoint(), lastPoint(), p, p);

    points_.push_back(p);

    invalidate();
  }

  void quadTo(const CPoint2D &p1, const CPoint2D &p2) {
    assert(false); // TODO

    points_.push_back(p1);
    points_.push_back(p2);
  }

  void cubicTo(const CPoint2D &p1, const CPoint2D &p2, const CPoint2D &p3) {
    beziers_.emplace_back(lastPoint(), p1, p2, p3);

    points_.push_back(p1);
    points_.push_back(p2);
    points_.push_back(p3);

    invalidate();
  }

  void arc(const CPoint2D &p, double sx, double sy, double a1, double a2) {
    CArcToBezier::BezierList beziers;
    CArcToBezier::ArcToBeziers(p.x, p.y, sx, sy, a1, a2, beziers);

    for (const auto &b : beziers) {
      beziers_.push_back(b);

      points_.push_back(b.getFirstPoint   ());
      points_.push_back(b.getControlPoint1());
      points_.push_back(b.getControlPoint2());
      points_.push_back(b.getLastPoint    ());
    }
  }

  void close() {
    setClosed(true);
  }

  //---

  void addRoundedRect(const CBBox2D &bbox, double xr, double yr) {
    auto p1 = bbox.getMin();
    auto p2 = bbox.getMax();

    xr = std::min(xr, bbox.getWidth ()/2.0);
    yr = std::min(yr, bbox.getHeight()/2.0);

    CBezierPath bezierPath;

    bezierPath.moveTo(CPoint2D(p1.x + xr, p1.y));

    bezierPath.lineTo(CPoint2D(p2.x - xr, p1.y));

    bezierPath.arc(CPoint2D(p2.x - xr, p1.y + yr), xr, yr, -M_PI/2.0, 0);

    bezierPath.lineTo(CPoint2D(p2.x, p2.y - yr));

    bezierPath.arc(CPoint2D(p2.x - xr, p2.y - yr), xr, yr, 0, M_PI/2.0);

    bezierPath.lineTo(CPoint2D(p1.x + xr, p2.y));

    bezierPath.arc(CPoint2D(p1.x + xr, p2.y - yr), xr, yr, M_PI/2.0, M_PI);

    bezierPath.lineTo(CPoint2D(p1.x, p1.y + yr));

    bezierPath.arc(CPoint2D(p1.x + xr, p1.y + yr), xr, yr, M_PI, 3.0*M_PI/2.0);

    bezierPath.close();

    for (const auto &b : bezierPath.beziers_)
      beziers_.push_back(b);
  }

  void addRect(const CBBox2D &bbox) {
    auto p1 = bbox.getMin();
    auto p2 = bbox.getMax();

    CBezierPath bezierPath;

    bezierPath.moveTo(CPoint2D(p1.x, p1.y));
    bezierPath.lineTo(CPoint2D(p2.x, p1.y));
    bezierPath.lineTo(CPoint2D(p2.x, p2.y));
    bezierPath.lineTo(CPoint2D(p1.x, p2.y));
    bezierPath.close();

    for (const auto &b : bezierPath.beziers_)
      beziers_.push_back(b);
  }

  //---

  void addPolyStar(const CPoint2D &center, int numPoints,
                   double innerRadius, double outerRadius,
                   double innerRoundness, double outerRoundness,
                   double rotation) {
    //assert(innerRoundness >= 0.0 && innerRoundness <= 1.0);
    //assert(outerRoundness >= 0.0 && outerRoundness <= 1.0);
    innerRoundness = CMathUtil::clamp(innerRoundness, 0.0, 1.0);
    outerRoundness = CMathUtil::clamp(outerRoundness, 0.0, 1.0);

    CBezierPath bezierPath;

    auto da = (numPoints > 0 ? 2.0*M_PI/(2*numPoints) : 0);

    struct PolyPoint {
      bool     inner { false };
      CPoint2D point;
      double   angle { 0.0 };
      CPoint2D lpoint;
      CPoint2D rpoint;
    };

    std::vector<PolyPoint> points;

    for (int i = 0; i < numPoints; ++i) {
      PolyPoint pi, po;

      po.angle = 2*i*da + M_PI/2.0 + rotation;
      pi.angle = po.angle + da;

      auto ci = std::cos(pi.angle); auto si = std::sin(pi.angle);
      auto co = std::cos(po.angle); auto so = std::sin(po.angle);

      auto xo = center.x + co*outerRadius;
      auto yo = center.y + so*outerRadius;
      auto xi = center.x + ci*innerRadius;
      auto yi = center.y + si*innerRadius;

      pi.point = CPoint2D(xi, yi);
      po.point = CPoint2D(xo, yo);

      pi.inner = true;

      points.push_back(po);
      points.push_back(pi);
    }

    bool rounded = false;

    if (innerRoundness > 0.0 || outerRoundness > 0.0) {
      auto d = points[0].point.distanceTo(points[1].point)/2.0;

      for (auto &p : points) {
        auto g = CVector2D(center, p.point).perpendicular().normalized();

        auto r = (p.inner ? innerRoundness : outerRoundness);

        auto p1 = p.point + g*r*d;
        auto p2 = p.point - g*r*d;

        p.lpoint = p1;
        p.rpoint = p2;
      }

      rounded = true;
    }

    if (rounded) {
      const auto &p = points[points.size() - 1];

      bezierPath.moveTo(p.point);

      for (size_t i1 = points.size() - 1, i2 = 0; i2 < points.size(); i1 = i2++) {
        const auto &p1 = points[i1];
        const auto &p2 = points[i2];

        bezierPath.cubicTo(p1.rpoint, p2.lpoint, p2.point);
      }
    }
    else {
      int i = 0;

      for (const auto &p : points) {
        if (i == 0)
          bezierPath.moveTo(p.point);
        else
          bezierPath.lineTo(p.point);

        ++i;
      }

      bezierPath.lineTo(points[0].point);
    }

    bezierPath.close();

    for (const auto &b : bezierPath.beziers_)
      beziers_.push_back(b);
  }

  void addPolygon(const CPoint2D &center, int numPoints,
                  double outerRadius, double outerRoundness,
                  double rotation) {
    //assert(outerRoundness >= 0.0 && outerRoundness <= 1.0);
    outerRoundness = CMathUtil::clamp(outerRoundness, 0.0, 1.0);

    CBezierPath bezierPath;

    auto da = (numPoints > 0 ? 2.0*M_PI/numPoints : 0);

    struct PolyPoint {
      CPoint2D point;
      double   angle { 0.0 };
      CPoint2D lpoint;
      CPoint2D rpoint;
    };

    std::vector<PolyPoint> points;

    for (int i = 0; i < numPoints; ++i) {
      PolyPoint po;

      po.angle = i*da + M_PI/2.0 + rotation;

      auto co = std::cos(po.angle); auto so = std::sin(po.angle);

      auto xo = center.x + co*outerRadius;
      auto yo = center.y + so*outerRadius;

      po.point = CPoint2D(xo, yo);

      points.push_back(po);
    }

    bool rounded = false;

    if (outerRoundness > 0.0) {
      auto d = points[0].point.distanceTo(points[1].point)/2.0;

      for (auto &p : points) {
        auto g = CVector2D(center, p.point).perpendicular().normalized();

        auto r = outerRoundness;

        auto p1 = p.point + g*r*d;
        auto p2 = p.point - g*r*d;

        p.lpoint = p1;
        p.rpoint = p2;
      }

      rounded = true;
    }

    if (rounded) {
      const auto &p = points[points.size() - 1];

      bezierPath.moveTo(p.point);

      for (size_t i1 = points.size() - 1, i2 = 0; i2 < points.size(); i1 = i2++) {
        const auto &p1 = points[i1];
        const auto &p2 = points[i2];

        bezierPath.cubicTo(p1.rpoint, p2.lpoint, p2.point);
      }
    }
    else {
      int i = 0;

      for (const auto &p : points) {
        if (i == 0)
          bezierPath.moveTo(p.point);
        else
          bezierPath.lineTo(p.point);

        ++i;
      }

      bezierPath.lineTo(points[0].point);
    }

    bezierPath.close();

    for (const auto &b : bezierPath.beziers_)
      beziers_.push_back(b);
  }

  //---

  double arcLength() const {
    if (! lengthValid_) {
      length_ = 0;

      for (const auto &b : beziers_)
        length_ += b.arcLength();

      lengthValid_ = true;
    }

    return length_;
  }

  //---

  CBezierPath split(double s, double e) const {
    bool flip = false;

    if (s > e) {
      flip = true;

      std::swap(s, e);
    }

    // none
    if (s == e)
      return CBezierPath();

    // all
    if (s <= 0.0 && e >= 1.0)
      return *this;

    //---

    std::vector<C3Bezier2D> beziers;

    auto l = arcLength();

    auto s1 = s*l;
    auto e1 = e*l;

    double l1 = 0.0;
    double l2 = 0.0;

    bool inside = false;

    for (const auto &b : beziers_) {
      auto bl = b.arcLength();

      auto isBreak = b.isBreak();

      l2 = l1 + bl;

      if (! inside) {
        if (s1 >= l1 && s1 < l2) {
          inside = true;

          auto ts = (s1 - l1)/bl;

          if (ts > 0.0) {
            C3Bezier2D bezier1, bezier2;
            auto rc = b.split(ts, bezier1, bezier2);
            assert(rc);

            auto bl1 = bezier1.arcLength();

            l1 += bl1;

            if (e1 >= l1 && e1 < l2) {
              bl -= bl1;

              auto te = (e1 - l1)/bl;

              if (te < 1.0) {
                C3Bezier2D bezier3, bezier4;
                auto rc1 = bezier2.split(te, bezier3, bezier4);
                assert(rc1);

                if (isBreak)
                  bezier3.setBreak(true);

                beziers.push_back(bezier3);
              }
              else
                beziers.push_back(bezier2);

              inside = false;

              break;
            }
            else {
              if (isBreak)
                bezier2.setBreak(true);

              beziers.push_back(bezier2);
            }
          }
          else {
            if (e1 >= l1 && e1 < l2) {
              auto te = (e1 - l1)/bl;

              if (te < 1.0) {
                C3Bezier2D bezier3, bezier4;
                auto rc = b.split(te, bezier3, bezier4);
                assert(rc);

                if (isBreak)
                  bezier3.setBreak(true);

                beziers.push_back(bezier3);
              }
              else
                beziers.push_back(b);

              inside = false;

              break;
            }
            else
              beziers.push_back(b);
          }
        }
      }
      else {
        if (e1 > l1 && e1 <= l2) {
          auto te = (e1 - l1)/bl;

          if (te < 1.0) {
            C3Bezier2D bezier1, bezier2;
            auto rc = b.split(te, bezier1, bezier2);
            assert(rc);

            if (isBreak)
              bezier1.setBreak(true);

            beziers.push_back(bezier1);
          }
          else
            beziers.push_back(b);

          inside = false;

          break;
        }
        else {
          beziers.push_back(b);
        }
      }

      l1 = l2;
    }

    if (inside) {
      if (e1 >= l2)
        inside = false;
    }

    assert(! inside);

    if (flip) {
      for (auto &b : beziers)
        b.flip();

      std::reverse(beziers.begin(), beziers.end());
    }

    return CBezierPath(beziers);
  }

  void combine(const CBezierPath &bezierPath) {
    auto n = beziers_.size();

    if (n > 0)
      beziers_[n - 1].setBreak(true);

    for (const auto &b : bezierPath.beziers_)
      beziers_.push_back(b);
  }

  //---

  CBezierPath rounded(double round) const {
    // compute center
    CPoint2D c;

    auto n = points_.size();

    for (const auto &p : points_) {
      c += p;
    }

    c /= double(n);

    //---

    // add round points
    struct PolyPoint {
      CPoint2D point;
      CPoint2D lpoint;
      CPoint2D rpoint;
    };

    std::vector<PolyPoint> points;

    for (const auto &p : points_) {
      auto v = CVector2D(c, p);
      auto r = v.length();

      PolyPoint po;

      po.point = p;

      auto g = v.perpendicular().normalized();

      auto p1 = p - g*r*round;
      auto p2 = p + g*r*round;

      po.lpoint = p1;
      po.rpoint = p2;

      points.push_back(po);
    }

    //---

    // build bezier
    CBezierPath bezierPath;

    const auto &p = points[n - 1];

    bezierPath.moveTo(p.point);

    for (size_t i1 = n - 1, i2 = 0; i2 < n; i1 = i2++) {
      const auto &p1 = points[i1];
      const auto &p2 = points[i2];

      bezierPath.cubicTo(p1.rpoint, p2.lpoint, p2.point);
    }

    return bezierPath;
  }

  //---

  CPoint2D calc(double t) const {
    std::vector<C3Bezier2D> beziers;

    auto l = arcLength();

    auto tl = t*l;

    double l1 = 0.0;
    double l2 = 0.0;

    for (const auto &b : beziers_) {
      auto bl = b.arcLength();

      l2 = l1 + bl;

      if (tl >= l1 && tl < l2) {
        auto ts = (tl - l1)/bl;

        return b.calc(ts);
      }

      l1 = l2;
    }

    return CPoint2D(0, 0);
  }

  double gradient(double t) const {
    std::vector<C3Bezier2D> beziers;

    auto l = arcLength();

    auto tl = t*l;

    double l1 = 0.0;
    double l2 = 0.0;

    for (const auto &b : beziers_) {
      auto bl = b.arcLength();

      l2 = l1 + bl;

      if (tl >= l1 && tl < l2) {
        auto ts = (tl - l1)/bl;

        return b.gradient(ts);
      }

      l1 = l2;
    }

    return 0.0;
  }

  //---

  void invalidate() {
    lengthValid_ = false;
  }

  //---

  CBBox2D bbox() const {
    CBBox2D bbox;

    for (const auto &b : beziers_)
      bbox += b.getBBox();

    return bbox;
  }

  const CPoint2D &lastPoint() const {
    static CPoint2D s_p;
    return (! points_.empty() ? points_.back() : s_p);
  }

  //---

  friend std::ostream &operator<<(std::ostream &os, const CBezierPath &path) {
    os << path.toString();
    return os;
  }

  std::string toString() const {
    std::string str = "(";

    int i = 0;

    for (const auto &b : beziers_) {
      if (i > 0)
        str += ", ";

      str += b.toString();

      ++i;
    }

    str += ")";

    return str;
  }

 private:
  Beziers beziers_;
  bool    closed_ { false };

  Points points_;

  mutable bool   lengthValid_ { false };
  mutable double length_      { 0.0 };
};

#endif
