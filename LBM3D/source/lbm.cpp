#include "functors.hpp"
#include "params.h"
#include "output.h"

#include "athread.h"
#include <ostream>
#include <stdlib.h>
#include <memory>

#include "cstring"

#include <chrono>

typedef Kokkos::MDRangePolicy<Kokkos::Rank<3> > range_3d;
typedef Kokkos::MDRangePolicy<Kokkos::Rank<2> > range_2d;

extern void cpp_lbm_update (double* const fB, double* const fA, 
		double* const u, double* const v, double* const w, double* const rho,
				const Params &params, const int &step, int &converged, 
						const int &D1, const int &D2, const int &D3);

extern "C" void athread_lbm_update (double* fB, double* fA, 
		double* u, double* v, double* w, double* rho,
				struct Params *params, int step, int* converged, int D1, int D2, int D3);

void update(DistributionField fB, DistributionField fA, ScalarField u, ScalarField v, ScalarField w, ScalarField rho, Params &params, const int step,
    int &converged);

int main(int narg, char *arg[]) {

  Kokkos::initialize(narg, arg);

  struct Params params;

  params.nx = 100;
  params.ny = 100;
  params.nz = 100;
  params.max_steps = 100;
  params.output_rate = 10;
  params.re = 100.;
  params.tau = 1.;
  params.lid_u = 0.01;
  params.lid_w = 0.0;
  params.lid_mag = 0.0;
  params.nu = 0.0;
  params.tol = 1e-2;

  std::unique_ptr < Output > o(new Output());

  using elpased_time = std::chrono::duration<float, std::ratio<1, 1>>;
  auto time1 = std::chrono::high_resolution_clock::now();
  auto time2 = std::chrono::high_resolution_clock::now();
  double time_kok (0.0), time_cpp(0.0), time_ath (0.0);

  double* h_fA = nullptr;
  double* h_fB = nullptr;
  double* h_fC = nullptr;
  double* h_rho = nullptr;
  double* h_u = nullptr;
  double* h_v = nullptr;
  double* h_w = nullptr;

  double* d_fA = nullptr;
  double* d_fB = nullptr;
  double* d_fC = nullptr;
  double* d_rho = nullptr;
  double* d_u = nullptr;
  double* d_v = nullptr;
  double* d_w = nullptr;

  int converged = 0;
  int step = 0;
  // params
  params.nx = atoi(arg[1]);
  params.ny = atoi(arg[2]);
  params.nz = atoi(arg[3]);
  params.re = atof(arg[4]);
  params.tol = atof(arg[5]);
  params.max_steps = atoi(arg[6]);
  params.output_rate = atoi(arg[7]);

  params.lid_u = 0.05;
  params.lid_w = 0.0;
  params.lid_mag = sqrt(params.lid_u * params.lid_u + params.lid_w * params.lid_w);
  params.nu = params.lid_mag * double(params.ny - 2) / params.re;
  params.tau = 3. * params.nu + 0.5;

  const int nx = params.nx;
  const int ny = params.ny;
  const int nz = params.nz;
  int D1, D2, D3;

  {
    // allocate memory on device
    DistributionField v_fA("fA", ny, nx, nz);
    DistributionField v_fB("fB", ny, nx, nz);
    DistributionField v_fC = v_fA;

    D1 = v_fA.extent(0);
    D2 = v_fA.extent(1);
    D3 = v_fA.extent(2);

    ScalarField v_rho("rho", ny, nx, nz);
    ScalarField v_u("u", ny, nx, nz);
    ScalarField v_v("v", ny, nx, nz);
    ScalarField v_w("w", ny, nx, nz);

    // initialize values on device
    Kokkos::parallel_for("load_state", range_3d( { 0, 0, 0 }, { v_fA.extent(0), v_fA.extent(1), v_fA.extent(2) }), load_state(v_fA, v_fB, v_u, v_v, v_w, v_rho));

    DistributionField::HostMirror h_v_fA = Kokkos::create_mirror_view(v_fA);
    DistributionField::HostMirror h_v_fB = Kokkos::create_mirror_view(v_fB);
    DistributionField::HostMirror h_v_fC = Kokkos::create_mirror_view(v_fC);
    ScalarField::HostMirror h_v_u = Kokkos::create_mirror_view(v_u);
    ScalarField::HostMirror h_v_v = Kokkos::create_mirror_view(v_v);
    ScalarField::HostMirror h_v_w = Kokkos::create_mirror_view(v_w);
    ScalarField::HostMirror h_v_rho = Kokkos::create_mirror_view(v_rho);

    // write t = 0 to output file
    // o->write_view("output/u", h_v_u);
    // o->write_view("output/v", h_v_v);
    // o->write_view("output/w", h_v_w);
    // o->write_view("output/rho", h_v_rho);
    o->frame += 1;

    printf("Solving lid driven cavity (Re = %.2e, tau = %.2e, [%i x %i x %i] , RAM (MB) = %.1f...\n", params.re, params.tau, nx - 2, ny - 2, nz - 2,
        double(nx) * double(ny) * double(nz) * 8. * (2. * 19. + 4.) / (1024. * 1024.));

    {
      h_fA = (double*) malloc (sizeof(double) * ny * nx * nz * 19);
      h_fB = (double*) malloc (sizeof(double) * ny * nx * nz * 19);
      h_fC = (double*) malloc (sizeof(double) * ny * nx * nz * 19);
   
      h_u = (double*) malloc (sizeof(double) * ny * nx * nz);
      h_v = (double*) malloc (sizeof(double) * ny * nx * nz);
      h_w = (double*) malloc (sizeof(double) * ny * nx * nz);
      h_rho = (double*) malloc (sizeof(double) * ny * nx * nz);
   
      memcpy (h_fA, h_v_fA.data(), sizeof(double) * ny * nx * nz * 19);
      memcpy (h_fB, h_v_fB.data(), sizeof(double) * ny * nx * nz * 19);
      memcpy (h_fC, h_v_fC.data(), sizeof(double) * ny * nx * nz * 19);
      memcpy (h_u, h_v_u.data(), sizeof(double) * ny * nx * nz);
      memcpy (h_v, h_v_v.data(), sizeof(double) * ny * nx * nz);
      memcpy (h_w, h_v_w.data(), sizeof(double) * ny * nx * nz);
      memcpy (h_rho, h_v_rho.data(), sizeof(double) * ny * nx * nz);
    }
    {
      d_fA = (double*) malloc (sizeof(double) * ny * nx * nz * 19);
      d_fB = (double*) malloc (sizeof(double) * ny * nx * nz * 19);
      d_fC = (double*) malloc (sizeof(double) * ny * nx * nz * 19);
   
      d_u = (double*) malloc (sizeof(double) * ny * nx * nz);
      d_v = (double*) malloc (sizeof(double) * ny * nx * nz);
      d_w = (double*) malloc (sizeof(double) * ny * nx * nz);
      d_rho = (double*) malloc (sizeof(double) * ny * nx * nz);
   
      memcpy (d_fA, h_v_fA.data(), sizeof(double) * ny * nx * nz * 19);
      memcpy (d_fB, h_v_fB.data(), sizeof(double) * ny * nx * nz * 19);
      memcpy (d_fC, h_v_fC.data(), sizeof(double) * ny * nx * nz * 19);
      memcpy (d_u, h_v_u.data(), sizeof(double) * ny * nx * nz);
      memcpy (d_v, h_v_v.data(), sizeof(double) * ny * nx * nz);
      memcpy (d_w, h_v_w.data(), sizeof(double) * ny * nx * nz);
      memcpy (d_rho, h_v_rho.data(), sizeof(double) * ny * nx * nz);
    }

    // Kokkos::Timer timer;

    {
      time1 = std::chrono::high_resolution_clock::now();
      for (step = 0; step < params.max_steps; ++ step) {
        // collide stream
        cpp_lbm_update (h_fB, h_fA, h_u, h_v, h_w, h_rho,
				     params, step, converged, D1, D2, D3);
        h_fC = h_fA;
        h_fA = h_fB;
        h_fB = h_fC;
      }
      time2 = std::chrono::high_resolution_clock::now();
      time_cpp += elpased_time(time2 - time1).count();
    }
    {
      time1 = std::chrono::high_resolution_clock::now();
      for (step = 0; step < params.max_steps; ++ step) {
        // collide stream
        update(v_fB, v_fA, v_u, v_v, v_w, v_rho, params, step, converged);
        v_fC = v_fA;
        v_fA = v_fB;
        v_fB = v_fC;
      }
      time2 = std::chrono::high_resolution_clock::now();
      time_kok += elpased_time(time2 - time1).count();
    }
  }
  Kokkos::fence();
  Kokkos::finalize();
  {
    athread_init();
    time1 = std::chrono::high_resolution_clock::now();
    for (step = 0; step < params.max_steps; ++ step) {
      // collide stream
      athread_lbm_update (d_fB, d_fA, d_u, d_v, d_w, d_rho,
		    &params, step, &converged, D1, D2, D3);
      d_fC = d_fA;
      d_fA = d_fB;
      d_fB = d_fC;
    }
    time2 = std::chrono::high_resolution_clock::now();
    time_ath += elpased_time(time2 - time1).count();
    athread_halt();
  }

  // while (step < params.max_steps && !converged) {

  //   printf ("-------------------\n");

  //   time1 = std::chrono::high_resolution_clock::now();
  //   // collide stream
  //   // update(v_fB, v_fA, v_u, v_v, v_w, v_rho, params, step, converged);
  //   time2 = std::chrono::high_resolution_clock::now();
  //   time_kok += elpased_time(time2 - time1).count();

  //   printf ("Kokkos converged: %d;\t", converged);
  //   const int D1 = v_fA.extent(0);
  //   const int D2 = v_fA.extent(1);
  //   const int D3 = v_fA.extent(2);
  //   time1 = std::chrono::high_resolution_clock::now();
  //   // collide stream
  //   cpp_lbm_update (h_fB, h_fA, h_u, h_v, h_w, h_rho,
		 	    // params, step, converged, D1, D2, D3);
  //   time2 = std::chrono::high_resolution_clock::now();
  //   time_cpp += elpased_time(time2 - time1).count();
  //   printf ("C++ converged: %d;\t", converged);

  //   athread_init();
  //   time1 = std::chrono::high_resolution_clock::now();
  //   // collide stream
  //   athread_lbm_update (d_fB, d_fA, d_u, d_v, d_w, d_rho,
		 		  // &params, step, &converged, D1, D2, D3);
  //   time2 = std::chrono::high_resolution_clock::now();
  //   athread_halt();
  //   time_ath += elpased_time(time2 - time1).count();
  //   printf ("Athread converged: %d\n", converged);

  //   // swap pointers
  //   v_fC = v_fA;
  //   v_fA = v_fB;
  //   v_fB = v_fC;

  //   h_fC = h_fA;
  //   h_fA = h_fB;
  //   h_fB = h_fC;

  //   d_fC = d_fA;
  //   d_fA = d_fB;
  //   d_fB = d_fC;

  //   // write to file
  //   if ((step + 1) % params.output_rate == 0) {

  //     // printf("...output step = %i\n", step + 1);

  //     // deep copy from device to host
  //     Kokkos::deep_copy(h_v_u, v_u);
  //     Kokkos::deep_copy(h_v_v, v_v);
  //     Kokkos::deep_copy(h_v_w, v_w);
  //     Kokkos::deep_copy(h_v_rho, v_rho);

  //     // write to output file
  //     // o->write_view("output/u", h_v_u);
  //     // o->write_view("output/v", h_v_v);
  //     // o->write_view("output/w", h_v_w);
  //     // o->write_view("output/rho", h_v_rho);

  //     o->frame += 1;
  //   }
  //   step += 1;
  //   double abs_acc_err_u   = 0.0;
  //   double abs_acc_err_v   = 0.0;
  //   double abs_acc_err_w   = 0.0;
  //   double abs_acc_err_rho = 0.0;
  //   for (int i = 0; i < nx*ny*nz; ++i) {
  //     abs_acc_err_u   += std::abs (h_u[i] - h_v_u.data()[i]);
  //     abs_acc_err_v   += std::abs (h_v[i] - h_v_v.data()[i]);
  //     abs_acc_err_w   += std::abs (h_w[i] - h_v_w.data()[i]);
  //     abs_acc_err_rho += std::abs (h_rho[i] - h_v_rho.data()[i]);
  //   }
  //   printf ("Step: %d, absolute cumulative error per step, Kokkos code:\t%.3e (u), %.3e (v) %.3e (w), %.3e (rho)\n",
  //       step, abs_acc_err_u, abs_acc_err_v, abs_acc_err_w, abs_acc_err_rho);
  //   abs_acc_err_u   = 0.0;
  //   abs_acc_err_v   = 0.0;
  //   abs_acc_err_w   = 0.0;
  //   abs_acc_err_rho = 0.0;
  //   for (int i = 0; i < nx*ny*nz; ++i) {
  //     abs_acc_err_u   += std::abs (d_u[i] - h_u[i]);
  //     abs_acc_err_v   += std::abs (d_v[i] - h_v[i]);
  //     abs_acc_err_w   += std::abs (d_w[i] - h_w[i]);
  //     abs_acc_err_rho += std::abs (d_rho[i] - h_rho[i]);
  //   }
  //   printf ("Step: %d, absolute cumulative error per step, Athread code:\t%.3e (u), %.3e (v) %.3e (w), %.3e (rho)\n",
  //       step, abs_acc_err_u, abs_acc_err_v, abs_acc_err_w, abs_acc_err_rho);
  // }

  // double time = timer.seconds();
  double site_updates = double(nx - 2) * double(ny - 2) * double(nz - 2) * double(step) / (1000. * 1000.);
  // double msus = site_updates / time;
  // double bandwidth = msus * 1000. * 1000. * 2. * 19. * 8. / (1024. * 1024. * 1024.);
  // double cost = time / double(step);

  printf ("===================\n");
  if (converged) {
    printf("Solution converged to steady state tolerance of %.3e\n", params.tol);
  } else {
    printf("Solution did not converged within %i steps\n", params.max_steps);
  }

  printf ("Cost per step (s) without I/O = %.5f (C++), %.5f (Athread), %.5f (Kokkos)\n", 
      time_cpp / double (step), time_ath / double (step), time_kok / double(step));
  // printf("MLUPS = %.1f, GB/s = %.1f, cost per step (s) = %.5f\n", msus, bandwidth, cost);

  {
    free (h_fA);
    free (h_fB);
    free (h_fC);
    free (h_u);
    free (h_v);
    free (h_w);
    free (h_rho);
    free (d_fA);
    free (d_fB);
    free (d_fC);
    free (d_u);
    free (d_v);
    free (d_w);
    free (d_rho);
  }

  // Kokkos::finalize();

  exit(EXIT_SUCCESS);
}

void update(DistributionField fB, DistributionField fA, 
    ScalarField u, ScalarField v, ScalarField w, ScalarField rho, 
        Params &params, const int step, int &converged) {
  using koArr2D = Kokkos::Array<int64_t, 2>;
  using koArr3D = Kokkos::Array<int64_t, 3>;


  const int D1 = fA.extent(0);
  const int D2 = fA.extent(1);
  const int D3 = fA.extent(2);

  // const int t1 = (D1 + 64) / 65;
  // const int t2 = (D2 + 64) / 65;
  // const int t3 = (D3 + 64) / 65;

  // static auto tile2D = koArr2D{t1, t2};
  // static auto tile3D = koArr3D{t1, t2, t3};
  // static auto tile3D = koArr3D{1, t2, t3};
  // static auto tile3D = koArr3D{2, 2, 126};

  // static auto tile2D = koArr2D{2, 32};
  // static auto tile3D = koArr3D{2, 2, 32};

  static auto tile2D = koArr2D{2, 32};
  static auto tile3D = koArr3D{2, 2, 32};

  Kokkos::parallel_for("collide_stream", 
      range_3d ({1, 1, 1}, {D1 - 1, D2 - 1, D3 - 1}, tile3D), collide_stream(fB, fA, params));

  Kokkos::parallel_for("bb_left", 
      range_2d ({1, 1}, {D1 - 1, D3 - 1}, tile2D), bb_left(fB));

  Kokkos::parallel_for("bb_right", 
      range_2d ({1, 1}, {D1 - 1, D3 - 1}, tile2D), bb_right(fB));

  Kokkos::parallel_for("bb_front", 
      range_2d ({1, 1}, {D1 - 1, D2 - 1}, tile2D), bb_front(fB));

  Kokkos::parallel_for("bb_back", 
      range_2d ({1, 1}, {D1 - 1, D2 - 1}, tile2D), bb_back(fB));

  Kokkos::parallel_for ("bb_bottom", 
      range_2d({1, 1}, {D2 - 1, D3 - 1}, tile2D), bb_bottom(fB));

  Kokkos::parallel_for("bb_top", 
      range_2d({1, 1}, {D2 - 1, D3 - 1}, tile2D), bb_top(fB, params));

  if ((step + 1) % params.output_rate == 0) {

    Kokkos::parallel_for("compute_macroscopic", 
        range_3d ({1, 1, 1}, {D1 - 1, D2 - 1, D3 - 1}, tile3D), compute_macroscopic(fB, u, v, w, rho));

    converged = 1;
    Kokkos::parallel_reduce("check_if_steady_state", 
        range_3d({ 1, 1, 1}, {D1 - 1, D2 - 1, D3 - 1}, tile3D), is_steady_state(fA, u, v, w, rho, params.tol), converged);

    // Kokkos::parallel_reduce("check_if_steady_state", range_3d( { 1, 1, 1 }, { D1 - 1, D2 - 1, D3 - 1 }), is_steady_state(fA, u, v, w, rho, params.tol),
    //     Kokkos::BAnd<int>(converged));
  }
  return ;
}

  // printf ("start debug\n");
  // for (int i = 0; i < nx*ny*nz*19; ++i) {
  //   if (d_fB[i] != h_fB[i]) {
  //     printf ("i = %d, h_fB = %.32f, d_fb = %.32f, diff = %.32f\n", i, h_fB[i], d_fB[i], h_fB[i] - d_fB[i]);
  //   }
  // }
  // printf ("end debug\n");
  // return 0;