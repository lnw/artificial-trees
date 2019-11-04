#ifndef SCENE_HH
#define SCENE_HH

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream> // for ofstream
#include <unordered_map>
#include <utility> // for pair

#include "tree.hh"


inline bool file_accessable(const std::string& fn) {
  std::ifstream f(fn.c_str());
  return f.good();
}

// everything about the depicted landscape that has nothing to do with pixels yet
class scene {
  std::vector<Tree> trees;
  std::vector<double> h_offsets; // one for each tree

  double x_min, x_max, y_min, y_max, z_min, z_max;

public:
  void add_tree() {
    trees.push_back(Tree());
    update_bounds();
  }

  void update_bounds();
  double get_x_min() const { return x_min; }
  double get_x_max() const { return x_max; }
  double get_y_min() const { return y_min; }
  double get_y_max() const { return y_max; }
  double get_z_min() const { return z_min; }
  double get_z_max() const { return z_max; }
};

#endif
