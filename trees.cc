
#include "canvas.hh"
#include "geometry3.hh"
#include "scene.hh"

using namespace std;

int main(int ac, char** av) {

  const string filename = "out.png";
  const int view_x(5000), view_y(2500); // pixels
  canvas V(filename, view_x, view_y);

  scene S;
  S.add_tree();

  V.bucket_fill(100, 100, 100);
  V.vertical_separator();
  const bool mirror = true;
  V.draw_ticks(S, mirror);
  V.render_scene(S);

  return 0;
}
