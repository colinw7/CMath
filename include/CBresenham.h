#ifndef CBRESENHAM_H
#define CBRESENHAM_H

#include <CLineDash.h>

class CBresenham {
 public:
  virtual ~CBresenham() { }

  void drawLine(int x1, int y1, int x2, int y2) {
    CILineDash line_dash = getLineDash();

    int dx = x2 - x1;

    if (dx == 0) {
      if (y2 > y1) {
        for (int y = y1; y <= y2; ++y) {
          if (line_dash.isDraw())
            drawPoint(x1, y);

          line_dash.step();
        }
      }
      else {
        for (int y = y2; y <= y1; ++y) {
          if (line_dash.isDraw())
            drawPoint(x1, y);

          line_dash.step();
        }
      }

      return;
    }

    int dy = y2 - y1;

    if (dy == 0) {
      if (x2 > x1) {
        for (int x = x1; x <= x2; ++x) {
          if (line_dash.isDraw())
            drawPoint(x, y1);

          line_dash.step();
        }
      }
      else {
        for (int x = x2; x <= x1; ++x) {
          if (line_dash.isDraw())
            drawPoint(x, y1);

          line_dash.step();
        }
      }

      return;
    }

    int adx = abs(dx);
    int ady = abs(dy);

    int eps = 0;

    if (adx > ady) {
      int y = y1;

      if (dx > 0) {
        if (dy > 0) {
          for (int x = x1; x <= x2; ++x) {
            if (line_dash.isDraw())
              drawPoint(x, y);

            line_dash.step();

            eps += ady;

            if ((eps << 1) >= (int) adx) {
              ++y;

              eps -= adx;
            }
          }
        }
        else {
          for (int x = x1; x <= x2; ++x) {
            if (line_dash.isDraw())
              drawPoint(x, y);

            line_dash.step();

            eps += ady;

            if ((eps << 1) >= (int) adx) {
              --y;

              eps -= adx;
            }
          }
        }
      }
      else {
        if (dy > 0) {
          for (int x = x1; x >= x2; --x) {
            if (line_dash.isDraw())
              drawPoint(x, y);

            line_dash.step();

            eps += ady;

            if ((eps << 1) >= (int) adx) {
              ++y;

              eps -= adx;
            }
          }
        }
        else {
          for (int x = x1; x >= x2; --x) {
            if (line_dash.isDraw())
              drawPoint(x, y);

            line_dash.step();

            eps += ady;

            if ((eps << 1) >= (int) adx) {
              --y;

              eps -= adx;
            }
          }
        }
      }
    }
    else {
      int x = x1;

      if (dy > 0) {
        if (dx > 0) {
          for (int y = y1; y <= y2; ++y) {
            if (line_dash.isDraw())
              drawPoint(x, y);

            line_dash.step();

            eps += adx;

            if ((eps << 1) >= (int) ady) {
              ++x;

              eps -= ady;
            }
          }
        }
        else {
          for (int y = y1; y <= y2; ++y) {
            if (line_dash.isDraw())
              drawPoint(x, y);

            line_dash.step();

            eps += adx;

            if ((eps << 1) >= (int) ady) {
              --x;

              eps -= ady;
            }
          }
        }
      }
      else {
        if (dx > 0) {
          for (int y = y1; y >= y2; --y) {
            if (line_dash.isDraw())
              drawPoint(x, y);

            line_dash.step();

            eps += adx;

            if ((eps << 1) >= (int) ady) {
              ++x;

              eps -= ady;
            }
          }
        }
        else {
          for (int y = y1; y >= y2; --y) {
            if (line_dash.isDraw())
              drawPoint(x, y);

            line_dash.step();

            eps += adx;

            if ((eps << 1) >= (int) ady) {
              --x;

              eps -= ady;
            }
          }
        }
      }
    }
  }

  void drawCircle(int x, int y, int r) {
    int xx = 0;
    int yy = r;

    int g  = 3 - 2*r;
    int gd = 10 - 4*r;
    int gr = 6;

    while (xx <= yy) {
      if (xx != 0) {
        if (yy != 0) {
          drawPoint(x + xx, y + yy);
          drawPoint(x + xx, y - yy);
          drawPoint(x - xx, y + yy);
          drawPoint(x - xx, y - yy);

          if (xx != yy) {
            drawPoint(x + yy, y + xx);
            drawPoint(x + yy, y - xx);
            drawPoint(x - yy, y + xx);
            drawPoint(x - yy, y - xx);
          }
        }
        else {
          drawPoint(x + xx, y     );
          drawPoint(x - xx, y     );
          drawPoint(x     , y + xx);
          drawPoint(x     , y - xx);
        }
      }
      else {
        if (yy != 0) {
          drawPoint(x     , y + yy);
          drawPoint(x     , y - yy);
          drawPoint(x + yy, y     );
          drawPoint(x - yy, y     );
        }
        else {
          drawPoint(x, y);
        }
      }

      if (g >= 0) {
        g  += gd;
        gd += 8;

        --yy;
      }
      else {
        g  += gr;
        gd += 4;
      }

      gr += 4;

      ++xx;
    }
  }

  // TODO: draw Pixel only once

  void fillCircle(int x, int y, int r) {
    int xx = 0;
    int yy = r;

    int g  = 3 - 2*r;
    int gd = 10 - 4*r;
    int gr = 6;

    while (xx <= yy) {
      if (xx != 0) {
        if (yy != 0) {
          drawHLine(x - xx, x + xx, y + yy);
          drawHLine(x - xx, x + xx, y - yy);

          if (xx != yy) {
            drawHLine(x - yy, x + yy, y + xx);
            drawHLine(x - yy, x + yy, y - xx);
          }
        }
        else {
          drawHLine(x - xx, x + xx, y);

          drawPoint(x, y + xx);
          drawPoint(x, y - xx);
        }
      }
      else {
        if (yy != 0) {
          drawPoint(x, y + yy);
          drawPoint(x, y - yy);

          drawHLine(x - yy, x + yy, y);
        }
        else {
          drawPoint(x, y);
        }
      }

      if (g >= 0) {
        g  += gd;
        gd += 8;

        --yy;
      }
      else {
        g  += gr;
        gd += 4;
      }

      gr += 4;

      ++xx;
    }
  }

  virtual void drawPoint(int x, int y) = 0;

  virtual const CILineDash &getLineDash() const {
    static CILineDash dash;

    return dash;
  }

 private:
  void drawHLine(int x1, int x2, int y) {
    for (int x = x1; x <= x2; ++x)
      drawPoint(x, y);
  }

  void drawVLine(int x, int y1, int y2) {
    for (int y = y1; y <= y2; ++y)
      drawPoint(x, y);
  }
};

#endif
