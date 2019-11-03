
#include <algorithm> // min, max
#include <cmath>     // modf
#include <iostream>
#include <tuple>
#include <vector>

#include <gd.h>

#include "canvas.hh"
#include "scene.hh"

using namespace std;


// draw line if both endpoints are visible,
// return true if drawn
bool canvas::draw_line(const double x1, const double y1,
                       const double x2, const double y2,
                       const double z,
                       int16_t r, int16_t g, int16_t b,
                       bool draw) {
  const int colour = gdImageColorResolve(img_ptr, r, g, b);
  if (z - 30 < zbuffer(x1, y1) && z - 30 < zbuffer(x2, y2)) {
    if (draw)
      gdImageLine(img_ptr, x1, y1, x2, y2, colour);
    return true;
  }
  return false;
}


void canvas::render_scene(const scene& S) {
  ofstream debug("debug-render_scene", ofstream::out | ofstream::app);

  // fill in here

  debug.close();
}

void canvas::render_test() {
  for (size_t y = 0; y < height; y++) {
    for (size_t x = 0; x < width; x++) {
      write_pixel(x, y, x, 0.1 * x, y);
    }
  }
}

void canvas::bucket_fill(const int r, const int g, const int b) {
  for (size_t y = 0; y < height; y++) {
    for (size_t x = 0; x < width; x++) {
      const int32_t col = 127 << 24 | r << 16 | g << 8 | b;
      img_ptr->tpixels[y][x] = col; // assuming TrueColor
    }
  }
}
