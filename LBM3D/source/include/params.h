#ifndef _PARAMS_H

#define _PARAMS_H

// struct Params {

//   int nx = 100;
//   int ny = 100;
//   int nz = 100;
//   int max_steps = 100;
//   int output_rate = 10;

//   double re = 100.;
//   double tau = 1.;
//   double lid_u = 0.01;
//   double lid_w = 0.0;
//   double lid_mag = 0.0;
//   double nu = 0.0;
//   double tol = 1e-2;
// };

struct Params {
  int nx;
  int ny;
  int nz;
  int max_steps;
  int output_rate;
  double re;
  double tau;
  double lid_u;
  double lid_w;
  double lid_mag;
  double nu;
  double tol;
};

#endif
