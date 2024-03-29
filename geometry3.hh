#ifndef GEOMETRY3_HH
#define GEOMETRY3_HH

#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "auxiliary.hh"

using namespace std;

struct matrix3d;

struct coord3d {
  double x[3];

  coord3d(const double y[3]) {
    x[0] = y[0];
    x[1] = y[1];
    x[2] = y[2];
  }
  coord3d(const double x_, const double y, const double z) {
    x[0] = x_;
    x[1] = y;
    x[2] = z;
  }
  coord3d() {
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;
  }
  coord3d(const double theta, const double phi) {
    x[0] = sin(theta) * cos(phi);
    x[1] = sin(theta) * sin(phi);
    x[2] = cos(theta);
  } // and r=0

  coord3d operator/(const double d) const { return coord3d(*this) /= d; }
  coord3d& operator/=(const double d) {
    x[0] /= d;
    x[1] /= d;
    x[2] /= d;
    return *this;
  }
  coord3d operator*(const double d) const { return coord3d(*this) *= d; }
  coord3d& operator*=(const double d) {
    x[0] *= d;
    x[1] *= d;
    x[2] *= d;
    return *this;
  }
  coord3d operator+(const coord3d& y) const { return coord3d(*this) += y; }
  coord3d& operator+=(const coord3d& y) {
    x[0] += y[0];
    x[1] += y[1];
    x[2] += y[2];
    return *this;
  }
  coord3d operator-(const coord3d& y) const { return coord3d(*this) -= y; }
  coord3d& operator-=(const coord3d& y) {
    x[0] -= y[0];
    x[1] -= y[1];
    x[2] -= y[2];
    return *this;
  }
  coord3d operator-() const {
    coord3d y(-x[0], -x[1], -x[2]);
    return y;
  }
  coord3d operator*(const matrix3d& m) const;

  coord3d cross(const coord3d& y) const {
    return coord3d(x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]);
  }
  matrix3d outer(const coord3d& y) const;
  double dot(const coord3d& y) const { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; }
  double norm() const { return sqrt(dot(*this)); }

  coord3d normalised() const { return *this / norm(); }

  double& operator[](const unsigned int i) { return x[i]; }
  double operator[](const unsigned int i) const { return x[i]; }

  //vector<double>::iterator begin(){return &x[0];}
  //vector<double>::iterator end(){return &x[3];}


  static double dist(const coord3d& x, const coord3d& y) { return (x - y).norm(); }
  // d/dx_i ||x|| = x_i/||x||.
  static coord3d dnorm(const coord3d& x) { return x / x.norm(); }
  // d^2/(dx_i dx_j) ||x|| = -x_i x_j/||x||^3 + [i==j]/||x||
  static void ddnorm(const coord3d& x, vector<double>& H) {
    const double n = 1.0 / x.norm(), n3 = n * n * n;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        H[i * 3 + j] = -x[i] * x[j] * n3 + (i == j ? n : 0);
  }

  void scale(const double f);

  // calculation of the angle beta at b(0,0,0)
  static double angle(const coord3d& a, const coord3d& c);
  // calculation of the derivative of angle beta at b(0,0,0) according to coordinates a and c with fixed b
  static void dangle(const coord3d& a, const coord3d& c, coord3d& da, coord3d& dc);
  // calculation of the dihedral angle theta at a(0,0,0), b, c and d,  the result is an angle between -\pi and +\pi (in radians)
  static double dihedral(const coord3d& b, const coord3d& c, const coord3d& d);
  // calculation of the derivative of dihedral angle theta at a(0,0,0), b, c and d according to coordinates b, c and d with fixed a
  static void ddihedral(const coord3d& b, const coord3d& c, const coord3d& d, coord3d& db, coord3d& dc, coord3d& dd);

  // coordinates of d-a in the dihedral abcd, where c2=b-a and c3=c-a
  // so the call is coord = foo(b-a, c-a, ...) + a
  static coord3d internal2cart(const coord3d& c2, const coord3d& c3, const double r, const double a, const double d);


  friend vector<coord3d>& operator-=(vector<coord3d>& xs, const coord3d& y) {
    for (size_t i = 0; i < xs.size(); i++)
      xs[i] -= y;
    return xs;
  }

  friend vector<coord3d>& operator*=(vector<coord3d>& xs, const double& y) {
    for (size_t i = 0; i < xs.size(); i++)
      xs[i] *= y;
    return xs;
  }

  friend ostream& operator<<(ostream& s, const coord3d& c3d) {
    s << fixed << "{" << c3d[0] << "," << c3d[1] << "," << c3d[2] << "}";
    return s;
  }
};


struct matrix3d {
  double values[9];

  //  matrix3d()                { memset(values,0,9*sizeof(double)); }
  matrix3d(const double* v) { memcpy(values, v, 9 * sizeof(double)); }
  explicit matrix3d(const double r = 0, const double s = 0, const double t = 0, const double u = 0, const double v = 0, const double w = 0, const double x = 0, const double y = 0, const double z = 0) {
    values[0] = r;
    values[1] = s;
    values[2] = t;
    values[3] = u;
    values[4] = v;
    values[5] = w;
    values[6] = x;
    values[7] = y;
    values[8] = z;
  }

  double& operator()(int i, int j) { return values[i * 3 + j]; }
  double operator()(int i, int j) const { return values[i * 3 + j]; }
  matrix3d operator+(const matrix3d& y) const { return matrix3d(*this) += y; }
  matrix3d& operator+=(const matrix3d& y) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        values[3 * i + j] += y(i, j);
      }
    };
    return *this;
  }
  matrix3d operator-(const matrix3d& y) const { return matrix3d(*this) -= y; }
  matrix3d& operator-=(const matrix3d& y) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        values[3 * i + j] -= y(i, j);
      }
    };
    return *this;
  }
  matrix3d operator*(const double s) const { return matrix3d(*this) *= s; }
  matrix3d& operator*=(const double s) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        values[3 * i + j] *= s;
      }
    };
    return *this;
  }
  matrix3d operator-() const {
    matrix3d m(-values[0], -values[1], -values[2], -values[3], -values[4], -values[5], -values[6], -values[7], -values[8]);
    return m;
  }

  matrix3d transpose() const {
    const matrix3d& M(*this);
    matrix3d Mt;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        Mt(i, j) = M(j, i);
    return Mt;
  }

  matrix3d inverse() const;

  double norm() const {
    const matrix3d& M(*this);

    double _norm = 0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        _norm += M(i, j) * M(i, j);

    return sqrt(_norm);
  }

  matrix3d operator*(const matrix3d& B) const {
    const matrix3d& A(*this);
    matrix3d C;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double sum = 0;
        for (int k = 0; k < 3; k++) {
          sum += A(i, k) * B(k, j);
        }
        C(i, j) = sum;
      }
    }
    return C;
  }

  coord3d operator*(const coord3d& x) const {
    coord3d y;
    for (int j = 0; j < 3; j++)
      y += coord3d(values[j] * x[j], values[3 + j] * x[j], values[6 + j] * x[j]);
    return y;
  }

  vector<coord3d> operator*(const vector<coord3d>& xs) const {
    const matrix3d& A(*this);
    vector<coord3d> ys(xs.size());

    for (size_t i = 0; i < xs.size(); i++)
      ys[i] = A * xs[i];
    return ys;
  }

  coord3d eigenvalues() const;
  coord3d eigenvector(const double lambda) const;
  pair<coord3d, matrix3d> eigensystem() const;

  static matrix3d unit_matrix() {
    return matrix3d(1, 0, 0, 0, 1, 0, 0, 0, 1);
  }

  friend ostream& operator<<(ostream& S, const matrix3d& M) {
    S << "{";
    for (int i = 0; i < 3; i++)
      S << vector<double>(&M.values[i * 3], &M.values[(i + 1) * 3]) << (i + 1 < 3 ? "," : "}");
    return S;
  }
};

#endif
