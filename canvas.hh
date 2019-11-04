#ifndef CANVAS_HH
#define CANVAS_HH

#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

#include <gd.h>

#include "array2D.hh"
#include "scene.hh"

using namespace std;


enum class DIRECTION : unsigned short { NORTH = 0,
                                        EAST,
                                        SOUTH,
                                        WEST };


class canvas {
public:
  unsigned width, height; // [pixels]

private:
  array2D<double> zbuffer;
  string filename;
  gdImagePtr img_ptr = nullptr;

public:
  canvas(string fn, int x, int y): width(x), height(y), zbuffer(x, y, INT_MAX), filename(fn) {
    // allocate mem
    img_ptr = gdImageCreateTrueColor(width, height);
  }

  ~canvas() {
    // actually the file is only opened here
    FILE* png_ptr = fopen(filename.c_str(), "wb");
    // write to disk
    gdImagePng(img_ptr, png_ptr);
    fclose(png_ptr);
    gdImageDestroy(img_ptr);
  }

  // just write the pixel
  void write_pixel(const int x, const int y,
                   int16_t r, int16_t g, int16_t b) {
    const int32_t col = 127 << 24 | r << 16 | g << 8 | b;
    img_ptr->tpixels[y][x] = col; // assuming TrueColor
  }

  // just write the pixel taking into account the zbuffer
  // true if pixel was drawn
  void write_pixel_zb(const int x, const int y, const double z,
                      int16_t r, int16_t g, int16_t b) {
    if (z < zbuffer(x, y)) {
      zbuffer(x, y) = z;
      const int32_t col = 127 << 24 | r << 16 | g << 8 | b;
      img_ptr->tpixels[y][x] = col; // assuming TrueColor
    }
  }

  // just write the pixel taking into account the zbuffer
  // true if pixel was drawn
  bool would_write_pixel_zb(const int x, const int y, const double z) {
    if (z > (zbuffer(x, y)))
      return false;
    else
      return true;
  }

  bool draw_line(const double x1, const double y1,
                 const double x2, const double y2,
                 const double z,
                 int16_t r, int16_t g, int16_t b,
                 bool draw);

  void render_scene(const scene& S);

  void render_test();

  void bucket_fill(const int r, const int g, const int b);
  void draw_tick(int x, int perp_offset, int tick_length, DIRECTION dir, string label);
  void draw_ticks(const scene& S, bool mirror);
  void vertical_separator();
};

#endif
