
#include <cmath>
#include <iomanip>
#include <limits>

#include "geometry3.hh"


coord3d coord3d::operator*(const matrix3d& m) const {
  return coord3d(x[0] * m(0, 0) + x[1] * m(1, 0) + x[2] * m(2, 0), x[0] * m(0, 1) + x[1] * m(1, 1) + x[2] * m(2, 1), x[0] * m(0, 2) + x[1] * m(1, 2) + x[2] * m(2, 2));
}

// calculation of the angle beta at b(0,0,0)
double coord3d::angle(const coord3d& a, const coord3d& c) {
  const double L2 = a.dot(a);
  const double R2 = c.dot(c);
  const double M2 = (c - a).dot(c - a);
  //if (abs(M2)<1.e-10) return 0;
  const double den = 2.0 * sqrt(L2 * R2);
  double arg = (L2 + R2 - M2) / den;
  if (arg > 1)
    return 0;
  if (arg < -1)
    return M_PI;
  return acos(arg);
}

// calculation of the derivative of angle beta at b(0,0,0) according to coordinates a and c with fixed b
void coord3d::dangle(const coord3d& a, const coord3d& c, coord3d& da, coord3d& dc) {
  const double L2 = a.dot(a);
  const double R2 = c.dot(c);
  const double M2 = (c - a).dot(c - a);
  const double den = 2.0 * sqrt(L2 * R2);
  const double arg = (L2 + R2 - M2) / den;

  const coord3d dM2__da = (a - c) * 2.0;
  const coord3d dL2__da = a * 2.0;
  const coord3d dden__da = dL2__da * R2 / sqrt(L2 * R2);
  const coord3d darg__da = dL2__da * 1.0 / den - dM2__da * 1.0 / den - dden__da * (L2 + R2 - M2) / (den * den);

  const coord3d dM2__dc = (c - a) * 2.0;
  const coord3d dR2__dc = c * 2.0;
  const coord3d dden__dc = dR2__dc * L2 / sqrt(L2 * R2);
  const coord3d darg__dc = dR2__dc * 1.0 / den - dM2__dc * 1.0 / den - dden__dc * (L2 + R2 - M2) / (den * den);

  da = -darg__da * 1.0 / sqrt(1.0 - arg * arg);
  dc = -darg__dc * 1.0 / sqrt(1.0 - arg * arg);
}

double coord3d::dihedral(const coord3d& b, const coord3d& c, const coord3d& d) {
  const coord3d ab = b; // a=0
  const coord3d bc = c - b;
  const coord3d cd = d - c;

  const coord3d abc = ab.cross(bc);
  const coord3d bcd = bc.cross(cd);

  const coord3d bc1 = bc / bc.norm();
  const coord3d abc1 = abc / abc.norm();
  const coord3d bcd1 = bcd / bcd.norm();
  const coord3d aux = abc1.cross(bc1);

  const double x = abc1.dot(bcd1);
  const double y = aux.dot(bcd1);

  return atan2(y, x);
}

// calculation of the derivative of dihedral angle theta at a(0,0,0), b, c and d  according to coordinates b, c and d with fixed a
void coord3d::ddihedral(const coord3d& b, const coord3d& c, const coord3d& d, coord3d& db, coord3d& dc, coord3d& dd) {
  const coord3d ab = b; // a=0
  const coord3d bc = c - b;
  const coord3d cd = d - c;

  const double bc_length_inv = 1.0 / bc.norm();
  const coord3d bc1 = bc * bc_length_inv;

  const coord3d abc = ab.cross(bc);
  const coord3d bcd = bc.cross(cd);

  const double abc_length_inv = 1.0 / abc.norm();
  const double bcd_length_inv = 1.0 / bcd.norm();
  const coord3d abc1 = abc * abc_length_inv;
  const coord3d bcd1 = bcd * bcd_length_inv;

  const coord3d aux = abc1.cross(bc1);

  const double x = abc1.dot(bcd1);
  const double y = aux.dot(bcd1);

  //  const double dihedral_abcd = atan2(y,x);
  //  cout << "D: "<< dihedral_abcd<<endl;

  const matrix3d dab__db = matrix3d::unit_matrix();
  const matrix3d dbc__db = -matrix3d::unit_matrix();
  const matrix3d dbc__dc = matrix3d::unit_matrix();
  const matrix3d dcd__dc = -matrix3d::unit_matrix();
  const matrix3d dcd__dd = matrix3d::unit_matrix();

  // bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
  const coord3d dbc_length_inv__dbc = -bc * pow(bc_length_inv, 3);

  // bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
  // vec = vec * mtx
  const coord3d dbc_length_inv__db = dbc_length_inv__dbc * dbc__db;
  const coord3d dbc_length_inv__dc = dbc_length_inv__dbc * dbc__dc;

  const matrix3d dbc1__dbc = matrix3d::unit_matrix() * bc_length_inv;

  // bc1_x=bc_x*bc_length_inv
  // bc1_y=bc_y*bc_length_inv
  // bc1_z=bc_z*bc_length_inv
  // mtx = mtx * mtx + vec outer vec
  const matrix3d dbc1__db = dbc1__dbc * dbc__db + bc.outer(dbc_length_inv__db);
  const matrix3d dbc1__dc = dbc1__dbc * dbc__dc + bc.outer(dbc_length_inv__dc);

  // abc_x=ab_y*bc_z - ab_z*bc_y
  // abc_y=ab_z*bc_x - ab_x*bc_z
  // abc_z=ab_x*bc_y - ab_y*bc_x
  //FIXME is there a more elegant way of doing this?
  const matrix3d dabc__dab = matrix3d(0, bc[2], -bc[1], -bc[2], 0, bc[0], bc[1], -bc[0], 0);
  const matrix3d dabc__dbc = matrix3d(0, -ab[2], ab[1], ab[2], 0, -ab[0], -ab[1], ab[0], 0);

  // bcd_x=bc_y*cd_z - bc_z*cd_y
  // bcd_y=bc_z*cd_x - bc_x*cd_z
  // bcd_z=bc_x*cd_y - bc_y*cd_x
  //FIXME is there a more elegant way of doing this?
  const matrix3d dbcd__dbc = matrix3d(0, cd[2], -cd[1], -cd[2], 0, cd[0], cd[1], -cd[0], 0);
  const matrix3d dbcd__dcd = matrix3d(0, -bc[2], bc[1], bc[2], 0, -bc[0], -bc[1], bc[0], 0);

  // abc_x=-ab_y*bc_z + ab_z*bc_y
  // abc_y=-ab_z*bc_x + ab_x*bc_z
  // abc_z=-ab_x*bc_y + ab_y*bc_x
  // mtx = mtx * mtx + mtx * mtx
  const matrix3d dabc__db = dabc__dab * dab__db + dabc__dbc * dbc__db;
  const matrix3d dabc__dc = dabc__dbc * dbc__dc;

  // bcd_x=-bc_y*cd_z + bc_z*cd_y
  // bcd_y=-bc_z*cd_x + bc_x*cd_z
  // bcd_z=-bc_x*cd_y + bc_y*cd_x
  // mtx = mtx * mtx + mtx * mtx
  const matrix3d dbcd__db = dbcd__dbc * dbc__db;
  const matrix3d dbcd__dc = dbcd__dbc * dbc__dc + dbcd__dcd * dcd__dc;
  const matrix3d dbcd__dd = dbcd__dcd * dcd__dd;

  // abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
  // bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
  const coord3d dabc_length_inv__dabc = -abc * pow(abc_length_inv, 3);
  const coord3d dbcd_length_inv__dbcd = -bcd * pow(bcd_length_inv, 3);

  // abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
  // vec = vec * mtx
  const coord3d dabc_length_inv__db = dabc_length_inv__dabc * dabc__db;
  const coord3d dabc_length_inv__dc = dabc_length_inv__dabc * dabc__dc;

  // bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
  // vec = vec * mtx
  const coord3d dbcd_length_inv__db = dbcd_length_inv__dbcd * dbcd__db;
  const coord3d dbcd_length_inv__dc = dbcd_length_inv__dbcd * dbcd__dc;
  const coord3d dbcd_length_inv__dd = dbcd_length_inv__dbcd * dbcd__dd;

  // abc1_x=abc_x*abc_length_inv
  // abc1_y=abc_y*abc_length_inv
  // abc1_z=abc_z*abc_length_inv
  const matrix3d dabc1__dabc = matrix3d::unit_matrix() * abc_length_inv;

  // abc1_x=abc_x*abc_length_inv
  // abc1_y=abc_y*abc_length_inv
  // abc1_z=abc_z*abc_length_inv
  // mtx = mtx * mtx + vec outer vec
  const matrix3d dabc1__db = dabc1__dabc * dabc__db + abc.outer(dabc_length_inv__db);
  const matrix3d dabc1__dc = dabc1__dabc * dabc__dc + abc.outer(dabc_length_inv__dc);

  // bcd1_x=bcd_x*bcd_length_inv
  // bcd1_y=bcd_y*bcd_length_inv
  // bcd1_z=bcd_z*bcd_length_inv
  const matrix3d dbcd1__dbcd = matrix3d::unit_matrix() * bcd_length_inv;

  // bcd1_x=bcd_x*bcd_length_inv
  // bcd1_y=bcd_y*bcd_length_inv
  // bcd1_z=bcd_z*bcd_length_inv
  // mtx = mtx*mtx + vec outer vec
  const matrix3d dbcd1__db = dbcd1__dbcd * dbcd__db + bcd.outer(dbcd_length_inv__db);
  const matrix3d dbcd1__dc = dbcd1__dbcd * dbcd__dc + bcd.outer(dbcd_length_inv__dc);
  const matrix3d dbcd1__dd = dbcd1__dbcd * dbcd__dd + bcd.outer(dbcd_length_inv__dd);

  // aux_x=abc1_y*bc1_z-bc1_y*abc1_z
  // aux_y=abc1_z*bc1_x-bc1_z*abc1_x
  // aux_z=abc1_x*bc1_y-bc1_x*abc1_y
  //FIXME is there a more elegant way of doing this?
  const matrix3d daux__dabc1 = matrix3d(0, bc1[2], -bc1[1], -bc1[2], 0, bc1[0], bc1[1], -bc1[0], 0);
  const matrix3d daux__dbc1 = matrix3d(0, -abc1[2], abc1[1], abc1[2], 0, -abc1[0], -abc1[1], abc1[0], 0);

  // aux_x=abc1_y*bc1_z-bc1_y*abc1_z
  // aux_y=abc1_z*bc1_x-bc1_z*abc1_x
  // aux_z=abc1_x*bc1_y-bc1_x*abc1_y
  // mtx = mtx*mtx + mtx*mtx
  const matrix3d daux__db = daux__dabc1 * dabc1__db + daux__dbc1 * dbc1__db;
  const matrix3d daux__dc = daux__dabc1 * dabc1__dc + daux__dbc1 * dbc1__dc;

  // y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
  // vec = vec * mtx
  const coord3d dy__db = bcd1 * daux__db + aux * dbcd1__db;
  const coord3d dy__dc = bcd1 * daux__dc + aux * dbcd1__dc;
  const coord3d dy__dd = aux * dbcd1__dd;

  // x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
  // vec = vec * mtx
  const coord3d dx__db = bcd1 * dabc1__db + abc1 * dbcd1__db;
  const coord3d dx__dc = bcd1 * dabc1__dc + abc1 * dbcd1__dc;
  const coord3d dx__dd = abc1 * dbcd1__dd;

  // df__dx=-y/(x**2 + y**2)
  // df__dy=x/(x**2 + y**2)
  const double df__dx = -y / (x * x + y * y);
  const double df__dy = x / (x * x + y * y);

  // f=atan2(y,x)
  // vec = vec*sca + vec*sca
  db = dx__db * df__dx + dy__db * df__dy;
  dc = dx__dc * df__dx + dy__dc * df__dy;
  dd = dx__dd * df__dx + dy__dd * df__dy;
}

// coordinates of d-a in the dihedral abcd, where c2=b-a and c3=c-a
// so the call is coord = foo(b-a, c-a, ...) + a
coord3d coord3d::internal2cart(const coord3d& c1, const coord3d& c2, const double r, const double a, const double d) {
  //ofstream debug("debug", ofstream::out | ofstream::app);

  //debug << "i2c: " << r << ", " << a << ", " << d << endl;
  //debug << "i2c: " << r << ", " << a*180/M_PI<< ", " << d*180/M_PI << endl;
  // c0 is at (0,0,0)
  coord3d c0_work(coord3d(0, 0, 0)), c1_work(c1), c2_work(c2);

  // orient c0, c1, c2: move by -c2
  c0_work -= c2;
  c1_work -= c2;
  c2_work -= c2;
  //debug << "geom (after shift) " << c0_work << ", " << c1_work << ", " << c2_work << endl;

  // orient c0, c1, c2: rotate c1 (and c2) such that c1 lies on (-x,0,0)
  double theta_0 = coord3d::angle(c1_work, coord3d(-1, 0, 0));
  coord3d rot_axis_0 = c1_work.cross(coord3d(-1, 0, 0));
  //debug << theta_0 << ", " << abs(theta_0 - M_PI) << endl;
  if (abs(theta_0 - M_PI) < 1.e-7 || abs(theta_0) < 1.e-7) {
    rot_axis_0 = coord3d(0, 0, 1);
    //debug << "overwriting" << endl;
  }
  //debug << "rot axis " << rot_axis_0 << endl;
  rot_axis_0 = (rot_axis_0 / rot_axis_0.norm());
  double u = rot_axis_0[0], v = rot_axis_0[1], w = rot_axis_0[2];
  //debug << "theta_0 " << theta_0 << endl;
  matrix3d rot0 = matrix3d(u * u + (1 - u * u) * cos(theta_0), u * v * (1 - cos(theta_0)) - w * sin(theta_0), u * w * (1 - cos(theta_0)) + v * sin(theta_0),
                           u * v * (1 - cos(theta_0)) + w * sin(theta_0), v * v + (1 - v * v) * cos(theta_0), v * w * (1 - cos(theta_0)) - u * sin(theta_0),
                           u * w * (1 - cos(theta_0)) - v * sin(theta_0), v * w * (1 - cos(theta_0)) + u * sin(theta_0), w * w + (1 - w * w) * cos(theta_0));
  c0_work = rot0 * c0_work;
  c1_work = rot0 * c1_work;
  //debug << "geom (after rot1) " << c0_work << ", " << c1_work << ", " << c2_work << endl;

  // orient c0, c1, c2: rotate c0 to lie on (-x,y,0)
  coord3d rot_axis_1 = coord3d(1, 0, 0);
  //debug << "rot axis " << rot_axis_1 << endl;
  u = rot_axis_1[0], v = rot_axis_1[1], w = rot_axis_1[2];
  double theta_1 = coord3d::angle(coord3d(0, c0_work[1], c0_work[2]), coord3d(0, 1, 0));
  if (c0_work[2] > 0)
    theta_1 *= -1;
  //debug << "theta_1 " << theta_1 << endl;
  matrix3d rot1 = matrix3d(u * u + (1 - u * u) * cos(theta_1), u * v * (1 - cos(theta_1)) - w * sin(theta_1), u * w * (1 - cos(theta_1)) + v * sin(theta_1),
                           u * v * (1 - cos(theta_1)) + w * sin(theta_1), v * v + (1 - v * v) * cos(theta_1), v * w * (1 - cos(theta_1)) - u * sin(theta_1),
                           u * w * (1 - cos(theta_1)) - v * sin(theta_1), v * w * (1 - cos(theta_1)) + u * sin(theta_1), w * w + (1 - w * w) * cos(theta_1));
  c0_work = rot1 * c0_work;
  //debug << "geom (after rot2) " << c0_work << ", " << c1_work << ", " << c2_work << endl;

  // place c3
  coord3d c3_work(-r, 0, 0);

  // rotate c3 by a around 0/0/z
  coord3d rot_axis_2 = coord3d(0, 0, 1);
  //debug << "rot axis " << rot_axis_2 << endl;
  u = rot_axis_2[0], v = rot_axis_2[1], w = rot_axis_2[2];
  //debug << "a " << a << endl;
  matrix3d rot2 = matrix3d(u * u + (1 - u * u) * cos(-a), u * v * (1 - cos(-a)) - w * sin(-a), u * w * (1 - cos(-a)) + v * sin(-a),
                           u * v * (1 - cos(-a)) + w * sin(-a), v * v + (1 - v * v) * cos(-a), v * w * (1 - cos(-a)) - u * sin(-a),
                           u * w * (1 - cos(-a)) - v * sin(-a), v * w * (1 - cos(a)) + u * sin(-a), w * w + (1 - w * w) * cos(-a));
  c3_work = rot2 * c3_work;
  //debug << "geom " << c3_work << endl;

  // rotate c3 by d around 0/0/z
  coord3d rot_axis_3 = coord3d(1, 0, 0);
  //debug << "rot axis " << rot_axis_3 << endl;
  u = rot_axis_3[0], v = rot_axis_3[1], w = rot_axis_3[2];
  //debug << "d " << d << endl;
  matrix3d rot3 = matrix3d(u * u + (1 - u * u) * cos(-d), u * v * (1 - cos(-d)) - w * sin(-d), u * w * (1 - cos(-d)) + v * sin(-d),
                           u * v * (1 - cos(-d)) + w * sin(-d), v * v + (1 - v * v) * cos(-d), v * w * (1 - cos(-d)) - u * sin(-d),
                           u * w * (1 - cos(-d)) - v * sin(-d), v * w * (1 - cos(-d)) + u * sin(-d), w * w + (1 - w * w) * cos(-d));
  c3_work = rot3 * c3_work;
  //debug << "geom " << c3_work << endl;

  // rotate everything back:
  u = rot_axis_1[0], v = rot_axis_1[1], w = rot_axis_1[2];
  matrix3d rot4 = matrix3d(u * u + (1 - u * u) * cos(-theta_1), u * v * (1 - cos(-theta_1)) - w * sin(-theta_1), u * w * (1 - cos(-theta_1)) + v * sin(-theta_1),
                           u * v * (1 - cos(-theta_1)) + w * sin(-theta_1), v * v + (1 - v * v) * cos(-theta_1), v * w * (1 - cos(-theta_1)) - u * sin(-theta_1),
                           u * w * (1 - cos(-theta_1)) - v * sin(-theta_1), v * w * (1 - cos(-theta_1)) + u * sin(-theta_1), w * w + (1 - w * w) * cos(-theta_1));
  c3_work = rot4 * c3_work;

  // rotate everything back:
  u = rot_axis_0[0], v = rot_axis_0[1], w = rot_axis_0[2];
  matrix3d rot5 = matrix3d(u * u + (1 - u * u) * cos(-theta_0), u * v * (1 - cos(-theta_0)) - w * sin(-theta_0), u * w * (1 - cos(-theta_0)) + v * sin(-theta_0),
                           u * v * (1 - cos(-theta_0)) + w * sin(-theta_0), v * v + (1 - v * v) * cos(-theta_0), v * w * (1 - cos(-theta_0)) - u * sin(-theta_0),
                           u * w * (1 - cos(-theta_0)) - v * sin(-theta_0), v * w * (1 - cos(-theta_0)) + u * sin(-theta_0), w * w + (1 - w * w) * cos(-theta_0));
  c3_work = rot5 * c3_work;
  //debug << "geom " << c3_work << endl;

  c3_work += c2;
  //debug << "geom " << c3_work << endl;

  //debug.close();

  return c3_work;
}

void coord3d::scale(const double f) {
  for (int i = 0; i < 3; i++) {
    x[i] *= f;
  }
}


matrix3d coord3d::outer(const coord3d& y) const {
  return matrix3d(x[0] * y[0], x[0] * y[1], x[0] * y[2], x[1] * y[0], x[1] * y[1], x[1] * y[2], x[2] * y[0], x[2] * y[1], x[2] * y[2]);
}


matrix3d matrix3d::inverse() const {

  const double den = (-(values[2] * values[4] * values[6]) + values[1] * values[5] * values[6] + values[2] * values[3] * values[7] - values[0] * values[5] * values[7] - values[1] * values[3] * values[8] + values[0] * values[4] * values[8]);
  const double m1 = (-(values[5] * values[7]) + values[4] * values[8]) / den;
  const double m2 = (values[2] * values[7] - values[1] * values[8]) / den;
  const double m3 = (-(values[2] * values[4]) + values[1] * values[5]) / den;
  const double m4 = (values[5] * values[6] - values[3] * values[8]) / den;
  const double m5 = (-(values[2] * values[6]) + values[0] * values[8]) / den;
  const double m6 = (values[2] * values[3] - values[0] * values[5]) / den;
  const double m7 = (-(values[4] * values[6]) + values[3] * values[7]) / den;
  const double m8 = (values[1] * values[6] - values[0] * values[7]) / den;
  const double m9 = (-(values[1] * values[3]) + values[0] * values[4]) / den;

  return matrix3d(m1, m2, m3, m4, m5, m6, m7, m8, m9);
}


// Eigenvalue solver specialized to symmetric real 3x3 matrices using Viete's
// closed form solution to cubic polynomials with three real roots.
coord3d matrix3d::eigenvalues() const {
  const matrix3d& M(*this);
  // Make sure that matrix is symmetric. TODO: FP comparison, not exact.
  assert(M(0, 1) == M(1, 0) && M(0, 2) == M(2, 0) && M(1, 2) == M(2, 1));

  // Coefficients up to symmetry
  double a(M(0, 0)), b(M(0, 1)), c(M(0, 2)), d(M(1, 1)), e(M(1, 2)), f(M(2, 2));

  // Coefficients of characteristic polynomial, calculated with Mathematica
  long double
      A = -1.L,
      B = a + d + f,
      C = b * b + c * c - a * d + e * e - a * f - d * f,
      D = -c * c * d + 2 * b * c * e - a * e * e - b * b * f + a * d * f;

  if (D == 0) { // Second order equation. TODO: FP comparison
    long double Disc = sqrtl(B * B - 4 * A * C);
    cout << "D = " << Disc << endl;
    return coord3d(0, (-B - Disc) / (2.L * A), (-B + Disc) / (2.L * A));
  }

  // Depress characteristic polynomial - see http://en.wikipedia.org/wiki/Cubic_equation#Reduction_to_a_depressed_cubic
  long double
      p = (3.L * A * C - B * B) / (3.L * A * A),
      q = (2.L * B * B * B - 9.L * A * B * C + 27.L * A * A * D) / (27.L * A * A * A),
      xc = B / (3.L * A);

  // François Viète's solution to cubic polynomials with three real roots.
  coord3d t;
  long double K = 2 * sqrtl(-p / 3.L),
              theta0 = (1.L / 3.L) * acosl((3.L * q) / (2.L * p) * sqrtl(-3.L / p));
  for (int k = 0; k < 3; k++)
    t[k] = K * cosl(theta0 - k * 2.L * M_PI / 3.L);

  // lambda = t - B/(3A)
  return t - coord3d(xc, xc, xc);
}


coord3d matrix3d::eigenvector(const double lambda) const {

  const matrix3d& M(*this);
  coord3d x;
  // cout << "M, l" << M << ", " << lambda << endl;

  // using the first two eqs
  // [ a_12 * a_23 - a_13 * (a_22 - r) ]
  // [ a_12 * a_13 - a_23 * (a_11 - r) ]
  // [ (a_11 - r) * (a_22 - r) - a_12^2 ]
  x = coord3d(M(0, 1) * M(1, 2) - M(0, 2) * (M(1, 1) - lambda),
              M(0, 1) * M(0, 2) - M(1, 2) * (M(0, 0) - lambda),
              (M(0, 0) - lambda) * (M(1, 1) - lambda) - M(0, 1) * M(0, 1));
  if (x.norm() / (M(0, 0) + M(1, 1) + M(2, 2)) > 1.e-12) // not zero-ish
    return x / x.norm();

  // using the first+last eqs
  // [ a_12 * (a_33 - r) - a_13 * a_23 ]
  // [ a_13^2 - (a_11 - r) * (a_33 - r) ]
  // [ a_23 * (a_11 - r) - a_12 * a_13 ]
  x = coord3d(M(0, 1) * (M(2, 2) - lambda) - M(0, 2) * M(1, 2),
              M(0, 2) * M(0, 2) - (M(0, 0) - lambda) * (M(2, 2) - lambda),
              M(1, 2) * (M(0, 0) - lambda) - M(0, 1) * M(0, 2));
  if (x.norm() / (M(0, 0) + M(1, 1) + M(2, 2)) > 1.e-12) // not zero-ish
    return x / x.norm();

  // using the last two eqs
  // [ a_23^2 - (a_22 - r) * (a_33 - r) ]
  // [ a_12 * (a_33 - r) - a_13 * a_23 ]
  // [ a_13 * (a_22 - r) - a_12 * a_23 ]
  x = coord3d(M(1, 2) * M(1, 2) - (M(1, 1) - lambda) * (M(2, 2) - lambda),
              M(0, 1) * (M(2, 2) - lambda) - M(0, 2) * M(1, 2),
              M(0, 2) * (M(1, 1) - lambda) - M(0, 1) * M(1, 2));
  if (x.norm() / (M(0, 0) + M(1, 1) + M(2, 2)) > 1.e-12) // not zero-ish
    return x / x.norm();

  cerr << "something is very wrong, possibly degenerate eigenvalues" << endl;
  return coord3d();
}


pair<coord3d, matrix3d> matrix3d::eigensystem() const {
  coord3d lambda(eigenvalues());

  // Sort eigenvalues by absolute value, smallest first
  if (fabs(lambda[0]) > fabs(lambda[1]))
    std::swap(lambda[0], lambda[1]);
  if (fabs(lambda[1]) > fabs(lambda[2]))
    std::swap(lambda[1], lambda[2]);
  if (fabs(lambda[0]) > fabs(lambda[1]))
    std::swap(lambda[0], lambda[1]);

  // Build eigenvector matrix
  matrix3d C;
  for (int i = 0; i < 3; i++) {
    coord3d c(eigenvector(lambda[i]));
    for (int j = 0; j < 3; j++)
      C(i, j) = c[j];
  }
  return make_pair(lambda, C);
}
