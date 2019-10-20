#ifndef NODE_HH
#define NODE_HH

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>  // for ofstream
#include <memory>
#include <limits.h>
#include <utility>  // for pair

#include "geometry3.hh"

class Node {
  coord3d c3d;
  std::vector<std::shared_ptr<Node>> children;

public:


};

#endif
