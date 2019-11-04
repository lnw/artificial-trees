
#include <algorithm> // min, max
#include <cmath>     // modf, floor
#include <iostream>
#include <tgmath.h> // log10
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


void canvas::vertical_separator() {
  const int black = gdImageColorResolve(img_ptr, 0, 0, 0);

  gdImageLine(img_ptr, width / 2 - 1, 0, width / 2 - 1, height, black);
  gdImageLine(img_ptr, width / 2, 0, width / 2, height, black);
  gdImageLine(img_ptr, width / 2 + 1, 0, width / 2 + 1, height, black);
}


void canvas::draw_ticks(const scene& S, bool mirror) {
  const double x_min = S.get_x_min();
  const double x_max = S.get_x_max();
  const double y_min = S.get_y_min();
  const double y_max = S.get_y_max();
  const double z_min = S.get_z_min();
  const double z_max = S.get_z_max();
  const double delta_x = x_max - x_min;
  const double delta_y = y_max - y_min;
  const double delta_z = z_max - z_min;

  const double scale = 4.0;
  // cout << "delta: " << delta_x << ", " << delta_y << ", " << delta_z << endl;
  // cout << "log delta: " << log(delta_x) / log(scale) << ", " << log(delta_y) / log(scale) << ", " << log(delta_z) / log(scale) << endl;
  // cout << "log delta: " << floor(log(delta_x) / log(scale) - 1) << ", " << floor(log(delta_y) / log(scale) - 1) << ", " << floor(log(delta_z) / log(scale) - 1) << endl;
  // cout << "log delta: " << pow(scale, floor(log(delta_x) / log(scale) - 1)) << ", " << pow(scale, floor(log(delta_y) / log(scale) - 1)) << ", " << pow(scale, floor(log(delta_z) / log(scale) - 1)) << endl;

  const double x_step = pow(scale, floor(log(delta_x) / log(scale) - 1)); // [m]
  const double y_step = pow(scale, floor(log(delta_y) / log(scale) - 1)); // [m]
  const double z_step = pow(scale, floor(log(delta_z) / log(scale) - 1)); // [m]

  const double pixels_per_m_x = (width / 2) / (delta_x /* * 1.2 */); // [px / m], and leave some margin
  const double pixels_per_m_y = (height) / (delta_y /* * 1.2 */);    // [px / m], and leave some margin
  const double pixels_per_m_z = (height) / (delta_z /* * 1.2 */);    // [px / m], and leave some margin

  const int shift_x = x_min * pixels_per_m_x; // [px]
  const int shift_y = y_min * pixels_per_m_y; // [px]
  const int shift_z = z_min * pixels_per_m_z; // [px]

  const int tick_length = 30;
  for (int x = int(ceil(x_min / x_step)); x <= int(x_max / x_step); x++) {
    int x_tick = x * pixels_per_m_x * x_step - shift_x;
    draw_tick(x_tick, height, tick_length, DIRECTION::NORTH, to_string_with_precision(x * x_step, 2));
    draw_tick(x_tick + width / 2, height, tick_length, DIRECTION::NORTH, to_string_with_precision(x * x_step, 2));
    if (mirror) {
      draw_tick(x_tick, 0, tick_length, DIRECTION::SOUTH, to_string_with_precision(x * x_step, 2));
      draw_tick(x_tick + width / 2, 0, tick_length, DIRECTION::SOUTH, to_string_with_precision(x * x_step, 2));
    }
  }

  for (int y = int(ceil(y_min / y_step)); y <= int(y_max / y_step); y++) {
    int y_tick = y * pixels_per_m_y * y_step - shift_y;
    draw_tick(height - y_tick, width / 2, tick_length, DIRECTION::EAST, to_string_with_precision(y * y_step, 2));
    if (mirror) {
      draw_tick(height - y_tick, width, tick_length, DIRECTION::WEST, to_string_with_precision(y * y_step, 2));
    }
  }

  for (int z = int(ceil(z_min / z_step)); z <= int(z_max / z_step); z++) {
    int z_tick = z * pixels_per_m_z * z_step - shift_z;
    draw_tick(z_tick, 0, tick_length, DIRECTION::EAST, to_string_with_precision(z * z_step, 2));
    if (mirror) {
      draw_tick(z_tick, width / 2, tick_length, DIRECTION::WEST, to_string_with_precision(z * z_step, 2));
    }
  }
}


void canvas::draw_tick(int offset, int perp_offset, int tick_length, DIRECTION dir, string label) {
  const int black = gdImageColorResolve(img_ptr, 0, 0, 0);
  const double fontsize = 25.;
  char* font = "./fonts/vera.ttf";
  const double text_orientation = 0; // [rad]

  switch (dir) {
    case DIRECTION::NORTH:
      gdImageLine(img_ptr, offset - 1, perp_offset - tick_length, offset - 1, perp_offset, black);
      gdImageLine(img_ptr, offset, perp_offset - tick_length, offset, perp_offset, black);
      gdImageLine(img_ptr, offset + 1, perp_offset - tick_length, offset + 1, perp_offset, black);
      break;
    case DIRECTION::EAST:
      gdImageLine(img_ptr, perp_offset, offset - 1, perp_offset + tick_length, offset - 1, black);
      gdImageLine(img_ptr, perp_offset, offset, perp_offset + tick_length, offset, black);
      gdImageLine(img_ptr, perp_offset, offset + 1, perp_offset + tick_length, offset + 1, black);
      break;
    case DIRECTION::SOUTH:
      gdImageLine(img_ptr, offset - 1, perp_offset, offset - 1, perp_offset + tick_length, black);
      gdImageLine(img_ptr, offset, perp_offset, offset, perp_offset + tick_length, black);
      gdImageLine(img_ptr, offset + 1, perp_offset, offset + 1, perp_offset + tick_length, black);
      break;
    case DIRECTION::WEST:
      gdImageLine(img_ptr, perp_offset - tick_length, offset - 1, perp_offset, offset - 1, black);
      gdImageLine(img_ptr, perp_offset - tick_length, offset, perp_offset, offset, black);
      gdImageLine(img_ptr, perp_offset - tick_length, offset + 1, perp_offset, offset + 1, black);
      break;
    default:
      assert(false);
  }

  if (!label.empty()) {
    // get bb of string
    int bb[8]; // SW - SE - NE - NW // SW is 0,0
    char* s1 = const_cast<char*>(label.c_str());
    char* err = gdImageStringFT(nullptr, &bb[0], 0, font, fontsize, 0., 0, 0, s1);
    if (err) {
      fprintf(stderr, "%s", err);
      cout << "not good" << endl;
    }
    // cout << bb[0] << " " << bb[1] << " " << bb[2] << " " << bb[3] << " " << bb[4] << " " << bb[5] << " " << bb[6] << " " << bb[7] << endl;

    int xxx, yyy;
    int space = 10;
    switch (dir) {
      case DIRECTION::NORTH:
        xxx = offset - bb[2] / 2;
        yyy = perp_offset - tick_length - space; // + bb[5];
        err = gdImageStringFT(img_ptr, &bb[0],
                              black, font, fontsize, text_orientation,
                              xxx,
                              yyy, s1);
        break;
      case DIRECTION::EAST:
        xxx = perp_offset + tick_length + space;
        yyy = offset - bb[5] / 2;
        err = gdImageStringFT(img_ptr, &bb[0],
                              black, font, fontsize, text_orientation,
                              xxx,
                              yyy, s1);
        break;
      case DIRECTION::SOUTH:
        xxx = offset - bb[2] / 2;
        yyy = perp_offset + tick_length + space - bb[5];
        err = gdImageStringFT(img_ptr, &bb[0],
                              black, font, fontsize, text_orientation,
                              xxx,
                              yyy, s1);
        break;
      case DIRECTION::WEST:
        xxx = perp_offset - tick_length - space - bb[2];
        yyy = offset - bb[5] / 2;
        err = gdImageStringFT(img_ptr, &bb[0],
                              black, font, fontsize, text_orientation,
                              xxx,
                              yyy, s1);
        break;
      default:
        assert(false);
    }

    if (err) {
      fprintf(stderr, "%s", err);
      cout << "not good" << endl;
    }
  }
}
