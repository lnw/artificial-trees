
#include <cassert>

#include "scene.hh"
#include "tree.hh"


void scene::update_bounds() {

  x_min = -2.0;
  x_max = 2.0;
  y_min = -1.0;
  y_max = 1.0;
  z_min = 0.0;
  z_max = 1.0;

  assert(x_max > x_min);
  assert(y_max > y_min);
  assert(z_max > z_min);
}
