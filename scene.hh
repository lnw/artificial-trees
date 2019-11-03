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
  std::vector<double> h_offsets;

public:
  void add_tree() {
    trees.push_back(Tree());
  }
};

#endif
