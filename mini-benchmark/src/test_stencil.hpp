#ifndef SRC_TEST_STENCIL_HPP_
#define SRC_TEST_STENCIL_HPP_

#include "kokkos_stencil.hpp"

#include "Kokkos_Core.hpp"
#include "athread.h"
// #include <simd.h>

#include <cstdio>

#include <math.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <type_traits>
#include <fstream>
#include <memory>


extern "C" void athread_5_point_stencil_float (int N, int I0, int I1, float* a, float* b);
extern "C" void athread_5_point_stencil_double (int N, int I0, int I1, double* a, double* b);

extern "C" void athread_5_point_stencil_double_1 (int N, int I0, int I1, double* a, double* b);
extern "C" void athread_5_point_stencil_double_2 (int N, int I0, int I1, double* a, double* b);
extern "C" void athread_5_point_stencil_double_3 (int N, int I0, int I1, double* a, double* b);
extern "C" void athread_5_point_stencil_double_4 (int N, int I0, int I1, double* a, double* b);
extern "C" void athread_5_point_stencil_double_5 (int N, int I0, int I1, double* a, double* b);
/*
template<typename T>
void test_mdrange_for_3d () {

	const int I0 = 32;
	const int I1 = 256;
	const int I2 = 10000;

	const int N = I0 * I1 * I2;
  size_t DATA_SIZE = N * sizeof(T);

	const bool PRINT_MATRIX = false;
	const bool DEBUG = true;

  // printf ("TIMES: %d, N: %d\n", TIMES, N);
  printf ("test 3d mdrange\n");
  std::cout<<"Data type: "<<typeid(T).name()<<", sizeof(T): "<<sizeof(T)<<std::endl;
  std::cout<<"Print matrix: "<<PRINT_MATRIX<<std::endl;

  double time_cpp;
  double time_ath;
  double time_kok;

  using elpased_time = std::chrono::duration<float, std::ratio<1, 1>>;
  auto time_cpp_1 = std::chrono::high_resolution_clock::now();
  auto time_cpp_2 = std::chrono::high_resolution_clock::now();
  auto time_kok_1 = std::chrono::high_resolution_clock::now();
  auto time_kok_2 = std::chrono::high_resolution_clock::now();
  auto time_ath_1 = std::chrono::high_resolution_clock::now();
  auto time_ath_2 = std::chrono::high_resolution_clock::now();

	T* h_a = nullptr;
	T* h_b = nullptr;
	T* h_c = nullptr;

	T* d_a = nullptr;
	T* d_b = nullptr;
	T* d_c = nullptr;

	T* t_c = nullptr;

  // 初始化主机端A和B数据
  {
    h_a = (T*) malloc (DATA_SIZE);
    h_b = (T*) malloc (DATA_SIZE);
    h_c = (T*) malloc (DATA_SIZE);
    d_a = (T*) malloc (DATA_SIZE);
    d_b = (T*) malloc (DATA_SIZE);
    d_c = (T*) malloc (DATA_SIZE);
    t_c = (T*) malloc (DATA_SIZE);

    for (int i = 0; i < N; ++i) {
      h_a[i] = rand() / static_cast<T>(RAND_MAX);
      h_b[i] = rand() / static_cast<T>(RAND_MAX);
      h_c[i] = rand() / static_cast<T>(RAND_MAX);
    }
  
    if (PRINT_MATRIX) {
      printf ("Matrix A: \n");
      for (int i = 0; i < N; ++i) {
        printf("%.2f\t", h_a[i]);
      }
      printf("\n");
      printf ("Matrix B: \n");
      for (int i = 0; i < N; ++i) {
        printf("%.2f\t", h_b[i]);
      }
      printf("\n");
    }
  }

  // Kokkos code
  auto kokkos_init_args = Kokkos::InitializationSettings()
    .set_print_configuration(true)
    .set_tune_internals(true);

  Kokkos::initialize(kokkos_init_args);

  using KokView = Kokkos::View<T ***>;
  KokView v_a("View Tensor A", I0, I1, I2);
  KokView v_b("View Tensor B", I0, I1, I2);
  KokView v_c("View Tensor C", I0, I1, I2);

  typename KokView::HostMirror h_v_a = Kokkos::create_mirror_view(v_a);
  typename KokView::HostMirror h_v_b = Kokkos::create_mirror_view(v_b);
  typename KokView::HostMirror h_v_c = Kokkos::create_mirror_view(v_c);

  // 初始化设备端A和B数据
  {
		for (int i0 = 0; i0 < I0; ++i0) {
		  for (int i1 = 0; i1 < I1; ++i1) {
		    for (int i2 = 0; i2 < I2; ++i2) {
      		h_v_a(i0, i1, i2) = h_a[i0*I1*I2 + i1*I2 + i2];
      		h_v_b(i0, i1, i2) = h_b[i0*I1*I2 + i1*I2 + i2];
      		d_a[i0*I1*I2 + i1*I2 + i2] = h_a[i0*I1*I2 + i1*I2 + i2];
      		d_b[i0*I1*I2 + i1*I2 + i2] = h_b[i0*I1*I2 + i1*I2 + i2];
        }
			}
		}
    Kokkos::deep_copy (v_a, h_v_a);
    Kokkos::deep_copy (v_b, h_v_b);
  }

	if (DEBUG) {
		// MPE, C++
		{
      time_cpp_1 = std::chrono::high_resolution_clock::now();
			for (int i0 = 0; i0 < I0; ++i0) {
			  for (int i1 = 0; i1 < I1; ++i1) {
			    for (int i2 = 0; i2 < I2; ++i2) {
						h_c[i0*I1*I2 + i1*I2 + i2] = h_a[i0*I1*I2 + i1*I2 + i2] * h_b[i0*I1*I2 + i1*I2 + i2];
          }
				}
			}
      time_cpp_2 = std::chrono::high_resolution_clock::now();
		}
		// MPE + 64 CPEs, Kokkos
		{
    	using Kokkos::parallel_for;
    	using Kokkos::MDRangePolicy;
    	using KoArr3D = Kokkos::Array<int64_t, 3>;
      int tile1 = (N + 64) / 65;
    	auto Tile3D   = KoArr3D{1, 1, tile1};
      time_kok_1 = std::chrono::high_resolution_clock::now();
    	parallel_for ("3d mult", MDRangePolicy<Kokkos::Rank<3>> (
        	KoArr3D{0, 0, 0}, KoArr3D{I0, I1, I2}, Tile3D), Functor3DMult<T>(v_a, v_b, v_c));
      time_kok_2 = std::chrono::high_resolution_clock::now();

    	Kokkos::deep_copy (h_v_c, v_c);

			for (int i0 = 0; i0 < I0; ++i0) {
			  for (int i1 = 0; i1 < I1; ++i1) {
			    for (int i2 = 0; i2 < I2; ++i2) {
					  t_c[i0*I1*I2 + i1*I2 + i2] = h_v_c(i0, i1, i2);
          }
				}
			}
      Kokkos::fence();
      Kokkos::finalize();
		}
		// 64 CPEs, Athread
		{
   		athread_init();
			if (std::is_same<typename std::decay<T>::type, float>::value) {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_3d_mult_float (I0, I1, I2, (float*)d_a, (float*)d_b, (float*)d_c);
        time_ath_2 = std::chrono::high_resolution_clock::now();
			} else if (std::is_same<typename std::decay<T>::type, double>::value) {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_3d_mult_double (I0, I1, I2, (double*)d_a, (double*)d_b, (double*)d_c);
        time_ath_2 = std::chrono::high_resolution_clock::now();
			}
	  	athread_halt();
		}
		double abs_errs = 0.0;
		for (int i = 0; i < N; ++i) {
			// printf ("i=%d, CPP: %.3f, Kokkos: %.3f, Athread: %.3f\t", 
			// 		i, h_c[i], t_c[i], d_c[i]);
			if (h_c[i] != t_c[i]) {
				abs_errs += std::abs (h_c[i] - t_c[i]);
				// abs_errs += 1;
				// printf ("err !");
			}
			// printf ("\n");
		}
		printf ("abs errors: %f\n", abs_errs);
    time_cpp = elpased_time(time_cpp_2 - time_cpp_1).count();
    time_ath = elpased_time(time_ath_2 - time_ath_1).count();
    time_kok = elpased_time(time_kok_2 - time_kok_1).count();
  	printf ("**********************************************\n");
  	printf ("Size, CPP, Athread, Kokkos, Speedup(Cpp/Kokkos), Efficiency(Athread/Kokkos)\n");
    printf ("%d, %.5f, %.5f, %.5f, %.3f, %.3f\n", N, time_cpp, time_ath, time_kok,
        time_cpp / time_kok, time_ath / time_kok);
	}
  // 释放资源
  {
    free (h_a);
    free (h_b);
    free (h_c);
    free (d_a);
    free (d_b);
    free (d_c);
    free (t_c);
#ifdef ENABLE_KOKKOS
    // Kokkos::fence();
    // Kokkos::finalize();
#endif // ENABLE_KOKKOS
  }
	return ;
}
*/

template<typename T>
void test_5_point_stencil () {

	// const int I0 = 267;
	// const int I1 = 420;

	// const int N = I0 * I1;
	const int TIMES = 200;
	const int START_NUM = 5e3;
  const int N = 5000;
  const int STRIDE = (N - START_NUM) / TIMES;

	const bool PRINT_MATRIX = false;
	const bool DEBUG = true;

  int I0 = N + 2;
  int I1 = N + 2;

  size_t DATA_SIZE = I0 * I1 * sizeof(T);

  // printf ("TIMES: %d, N: %d\n", TIMES, N);
  std::cout<<"Data type: "<<typeid(T).name()<<", sizeof(T): "<<sizeof(T)<<
      ", Matrix size (N*N): "<<N*N<<std::endl;
  std::cout<<"Print matrix: "<<PRINT_MATRIX<<std::endl;

  double time_cpp;
  double time_ath;
  double time_kok;

  using elpased_time = std::chrono::duration<float, std::ratio<1, 1>>;
  auto time_cpp_1 = std::chrono::high_resolution_clock::now();
  auto time_cpp_2 = std::chrono::high_resolution_clock::now();
  auto time_kok_1 = std::chrono::high_resolution_clock::now();
  auto time_kok_2 = std::chrono::high_resolution_clock::now();
  auto time_ath_1 = std::chrono::high_resolution_clock::now();
  auto time_ath_2 = std::chrono::high_resolution_clock::now();

	T* h_a = nullptr;
	T* h_b = nullptr;
	// T* h_c = nullptr;

	T* d_a = nullptr;
	T* d_b = nullptr;
	// T* d_c = nullptr;

  // 初始化主机端A和B数据
  {
    // h_a = (T*) malloc (DATA_SIZE);
    // h_b = (T*) malloc (DATA_SIZE);
    // // h_c = (T*) malloc (DATA_SIZE);
    // d_a = (T*) malloc (DATA_SIZE);
    // d_b = (T*) malloc (DATA_SIZE);
    // // d_c = (T*) malloc (DATA_SIZE);

    // h_a __attribute__ ((aligned(64))) = (T*) malloc (DATA_SIZE);
    // h_b __attribute__ ((aligned(64))) = (T*) malloc (DATA_SIZE);
    // d_a __attribute__ ((aligned(64))) = (T*) malloc (DATA_SIZE);
    // d_b __attribute__ ((aligned(64))) = (T*) malloc (DATA_SIZE);

    h_a = (T*) libc_aligned_malloc (DATA_SIZE);
    h_b = (T*) libc_aligned_malloc (DATA_SIZE);
    d_a = (T*) libc_aligned_malloc (DATA_SIZE);
    d_b = (T*) libc_aligned_malloc (DATA_SIZE);

    // libc_aligned_malloc
  
    for (int i = 0; i < I0 * I1; ++i) {
      h_a[i] = (rand()) / static_cast<T>(RAND_MAX);
      h_b[i] = (rand()) / static_cast<T>(RAND_MAX);
    }
    // memset (h_c, 0, DATA_SIZE);
    // memset (d_c, 0, DATA_SIZE);
    
    if (PRINT_MATRIX) {
      printf ("Matrix A: \n");
      for (int i = 0; i < I0 * I1; ++i) {
        printf("%.2f\t", h_a[i]);
      }
      printf("\n");
      printf ("Matrix B: \n");
      for (int i = 0; i < I0 * I1; ++i) {
        printf("%.2f\t", h_b[i]);
      }
      printf("\n");
    }
  }
  
  // Kokkos code
  auto kokkos_init_args = Kokkos::InitializationSettings()
    .set_print_configuration(true)
    .set_tune_internals(true);

  Kokkos::initialize(kokkos_init_args);
  using KokView = Kokkos::View<T **>;
  KokView v_a("View Matrix A", I0, I1);
  KokView v_b("View Matrix B", I0, I1);
  // KokView v_c("View Matrix C", I0, I1);
  
  typename KokView::HostMirror h_v_a = Kokkos::create_mirror_view(v_a);
  typename KokView::HostMirror h_v_b = Kokkos::create_mirror_view(v_b);
  // typename KokView::HostMirror h_v_c = Kokkos::create_mirror_view(v_c);
  
  // 初始化设备端A和B数据
  {
		for (int i0 = 0; i0 < I0; ++i0) {
		  for (int i1 = 0; i1 < I1; ++i1) {
    	  h_v_a(i0, i1) = h_a[i0*I1 + i1];
      	h_v_b(i0, i1) = h_b[i0*I1 + i1];
      	// h_v_c(i0, i1) = h_c[i0*I1 + i1];
      	d_a[i0*I1 + i1] = h_a[i0*I1 + i1];
      	d_b[i0*I1 + i1] = h_b[i0*I1 + i1];
			}
		}
    Kokkos::deep_copy (v_a, h_v_a);
    Kokkos::deep_copy (v_b, h_v_b);
    // Kokkos::deep_copy (v_c, h_v_c);
  }

	if (DEBUG) {
		// MPE, C++
		{
      time_cpp_1 = std::chrono::high_resolution_clock::now();
			for (int i0 = 1; i0 < N + 1; ++i0) {
			  for (int i1 = 1; i1 < N + 1; ++i1) {
					h_b[i0*I1 + i1] += h_a[(i0  )*I1 + (i1  )]
                           + h_a[(i0-1)*I1 + (i1  )]
                           + h_a[(i0+1)*I1 + (i1  )]
                           + h_a[(i0  )*I1 + (i1+1)]
                           + h_a[(i0  )*I1 + (i1-1)];
				}
			}
      time_cpp_2 = std::chrono::high_resolution_clock::now();
    }
		// MPE + 64 CPEs, Kokkos
    double t_1_32, t_1_64, t_1_128, t_1_256, t_1_512, t_1_1024;
    double t_2_512, t_4_512;
    double t_ath_1, t_ath_2, t_ath_3, t_ath_4, t_ath_5;
		{
    	using Kokkos::MDRangePolicy;
    	using KoArr2D = Kokkos::Array<int64_t, 2>;
      // auto Tile2D;
    	// Tile2D = KoArr2D{1, 32};
      time_kok_1 = std::chrono::high_resolution_clock::now();
    	Kokkos::parallel_for ("5-point-stencil", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{1, 1}, KoArr2D{N+1, N+1}, KoArr2D{1, 32}), Functor5PointStencil<T>(v_a, v_b));
      time_kok_2 = std::chrono::high_resolution_clock::now();
      t_1_32 = elpased_time(time_kok_2 - time_kok_1).count();

    	// Tile2D = KoArr2D{1, 64};
      time_kok_1 = std::chrono::high_resolution_clock::now();
    	Kokkos::parallel_for ("5-point-stencil", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{1, 1}, KoArr2D{N+1, N+1}, KoArr2D{1, 64}), Functor5PointStencil<T>(v_a, v_b));
      time_kok_2 = std::chrono::high_resolution_clock::now();
      t_1_64 = elpased_time(time_kok_2 - time_kok_1).count();

    	// Tile2D = KoArr2D{1, 128};
      time_kok_1 = std::chrono::high_resolution_clock::now();
    	Kokkos::parallel_for ("5-point-stencil", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{1, 1}, KoArr2D{N+1, N+1}, KoArr2D{1, 128}), Functor5PointStencil<T>(v_a, v_b));
      time_kok_2 = std::chrono::high_resolution_clock::now();
      t_1_128 = elpased_time(time_kok_2 - time_kok_1).count();

    	// Tile2D = KoArr2D{1, 256};
      time_kok_1 = std::chrono::high_resolution_clock::now();
    	Kokkos::parallel_for ("5-point-stencil", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{1, 1}, KoArr2D{N+1, N+1}, KoArr2D{1, 256}), Functor5PointStencil<T>(v_a, v_b));
      time_kok_2 = std::chrono::high_resolution_clock::now();
      t_1_256 = elpased_time(time_kok_2 - time_kok_1).count();

    	// Tile2D = KoArr2D{1, 512};
      time_kok_1 = std::chrono::high_resolution_clock::now();
    	Kokkos::parallel_for ("5-point-stencil", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{1, 1}, KoArr2D{N+1, N+1}, KoArr2D{1, 512}), Functor5PointStencil<T>(v_a, v_b));
      time_kok_2 = std::chrono::high_resolution_clock::now();
      t_1_512 = elpased_time(time_kok_2 - time_kok_1).count();

      time_kok_1 = std::chrono::high_resolution_clock::now();
    	Kokkos::parallel_for ("5-point-stencil", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{1, 1}, KoArr2D{N+1, N+1}, KoArr2D{1, 1024}), Functor5PointStencil<T>(v_a, v_b));
      time_kok_2 = std::chrono::high_resolution_clock::now();
      t_1_1024 = elpased_time(time_kok_2 - time_kok_1).count();

      time_kok_1 = std::chrono::high_resolution_clock::now();
    	Kokkos::parallel_for ("5-point-stencil", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{1, 1}, KoArr2D{N+1, N+1}, KoArr2D{2, 512}), Functor5PointStencil<T>(v_a, v_b));
      time_kok_2 = std::chrono::high_resolution_clock::now();
      t_2_512 = elpased_time(time_kok_2 - time_kok_1).count();

      time_kok_1 = std::chrono::high_resolution_clock::now();
    	Kokkos::parallel_for ("5-point-stencil", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{1, 1}, KoArr2D{N+1, N+1}, KoArr2D{4, 512}), Functor5PointStencil<T>(v_a, v_b));
      time_kok_2 = std::chrono::high_resolution_clock::now();
      t_4_512 = elpased_time(time_kok_2 - time_kok_1).count();

      Kokkos::fence();
      Kokkos::finalize();
		}
		// 64 CPEs, Athread
		{
   		athread_init();
			if (std::is_same<typename std::decay<T>::type, float>::value) {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_5_point_stencil_float (N, I0, I1, (float*)d_a, (float*)d_b);
        time_ath_2 = std::chrono::high_resolution_clock::now();
			} else if (std::is_same<typename std::decay<T>::type, double>::value) {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_5_point_stencil_double (N, I0, I1, (double*)d_a, (double*)d_b);
        time_ath_2 = std::chrono::high_resolution_clock::now();
			}
      {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_5_point_stencil_double_1 (N, I0, I1, (double*)d_a, (double*)d_b);
        time_ath_2 = std::chrono::high_resolution_clock::now();
        t_ath_1 = elpased_time(time_ath_2 - time_ath_1).count();

        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_5_point_stencil_double_2 (N, I0, I1, (double*)d_a, (double*)d_b);
        time_ath_2 = std::chrono::high_resolution_clock::now();
        t_ath_2 = elpased_time(time_ath_2 - time_ath_1).count();

        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_5_point_stencil_double_3 (N, I0, I1, (double*)d_a, (double*)d_b);
        time_ath_2 = std::chrono::high_resolution_clock::now();
        t_ath_3 = elpased_time(time_ath_2 - time_ath_1).count();

        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_5_point_stencil_double_4 (N, I0, I1, (double*)d_a, (double*)d_b);
        time_ath_2 = std::chrono::high_resolution_clock::now();
        t_ath_4 = elpased_time(time_ath_2 - time_ath_1).count();

        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_5_point_stencil_double_5 (N, I0, I1, (double*)d_a, (double*)d_b);
        time_ath_2 = std::chrono::high_resolution_clock::now();
        t_ath_5 = elpased_time(time_ath_2 - time_ath_1).count();

      }
	  	athread_halt();
		}
    time_cpp = elpased_time(time_cpp_2 - time_cpp_1).count();
    time_ath = elpased_time(time_ath_2 - time_ath_1).count();
    printf ("C++: %f, Athread: %f\n", time_cpp, time_ath);
    printf ("Kokkos:\n");
    printf ("Tile shape:\n");
    printf ("1x32: %f, 1x64: %f, 1x128: %f, 1x256: %f, 1x512: %f, 1x1024: %f\n",
        t_1_32, t_1_64, t_1_128, t_1_256, t_1_512, t_1_1024);
    printf ("1x32: %f, 1x64: %f, 1x128: %f, 2x512: %f, 4x512: %f\n",
        t_1_32, t_1_64, t_1_128, t_2_512, t_4_512);
    printf ("Athread, code1: %f, code2: %f, code3: %f, code4: %f, code5: %f\n",
        t_ath_1, t_ath_2, t_ath_3, t_ath_4, t_ath_5);
    // printf ("CPP: %.3f, Kokkos: %.3f, Raw kernel: %.3f\n", h_result, k_result, d_result);
    // time_cpp = elpased_time(time_cpp_2 - time_cpp_1).count();
    // time_ath = elpased_time(time_ath_2 - time_ath_1).count();
    // time_kok = elpased_time(time_kok_2 - time_kok_1).count();
  	// printf ("**********************************************\n");
  	// printf ("Size, CPP, Athread, Kokkos, Speedup(Cpp/Kokkos), Efficiency(Athread/Kokkos)\n");
    // printf ("%d, %.5f, %.5f, %.5f, %.3f, %.3f\n", N, time_cpp, time_ath, time_kok,
    //     time_cpp / time_kok, time_ath / time_kok);
	} else {
    /*
    std::unique_ptr<double []> time_cpp_arr = std::make_unique<double []>(TIMES + 1);
    std::unique_ptr<double []> time_ath_arr = std::make_unique<double []>(TIMES + 1);
    std::unique_ptr<double []> time_kok_arr = std::make_unique<double []>(TIMES + 1);
    std::ofstream file("dot_product_results.csv");
    file << std::fixed << std::setprecision(6);
    file << "Matrix size,CPP,Kokkos,Athread" << std::endl;
		// MPE, C++
		{
      int i = 0;
      for (int n = START_NUM; n <= N; n += STRIDE) {
        // h_result = 0.0;
        I0 = n;
        I1 = n;
        time_cpp_1 = std::chrono::high_resolution_clock::now();

    	// using Kokkos::MDRangePolicy;
    	// using KoArr2D = Kokkos::Array<int64_t, 2>;
      //   auto v_a_sub = Kokkos::subview(v_a, Kokkos::make_pair(0, I0), Kokkos::make_pair(0, I1));
      //   auto v_b_sub = Kokkos::subview(v_b, Kokkos::make_pair(0, I0), Kokkos::make_pair(0, I1));
      // 	// auto Tile2D = KoArr2D{2, 64};
      // 	auto Tile2D = KoArr2D{1, (I1+64)/65};
      	// Kokkos::parallel_reduce ("dot product 2d", MDRangePolicy<Kokkos::Rank<2>> (
        //   	KoArr2D{0, 0}, KoArr2D{I0, I1}, Tile2D), FunctorDotPro2D<T>(v_a_sub, v_b_sub), k_result);

				for (int i0 = 0; i0 < I0; ++i0) {
				  for (int i1 = 0; i1 < I1; ++i1) {
						h_result += h_a[i0*I1 + i1] * h_b[i0*I1 + i1];
					}
				}

        time_cpp_2 = std::chrono::high_resolution_clock::now();
        printf ("n = %d, h_result = %f, time: %f\n", n, 
            h_result, elpased_time(time_cpp_2 - time_cpp_1).count());
        time_cpp_arr[i++] = elpased_time(time_cpp_2 - time_cpp_1).count();
      }
      printf ("C++ ok\n");
		}
    // MPE + CPEs, Kokkos
		{
    	using Kokkos::MDRangePolicy;
    	using KoArr2D = Kokkos::Array<int64_t, 2>;
      int i = 0;
      for (int n = START_NUM; n <= N; n += STRIDE) {
        // k_result = 0.0;
        I0 = n;
        I1 = n;
        auto v_a_sub = Kokkos::subview(v_a, Kokkos::make_pair(0, I0), Kokkos::make_pair(0, I1));
        auto v_b_sub = Kokkos::subview(v_b, Kokkos::make_pair(0, I0), Kokkos::make_pair(0, I1));
      	auto Tile2D = KoArr2D{4, 1024};
        time_kok_1 = std::chrono::high_resolution_clock::now();
      	Kokkos::parallel_reduce ("dot product 2d", MDRangePolicy<Kokkos::Rank<2>> (
          	KoArr2D{0, 0}, KoArr2D{I0, I1}, Tile2D), FunctorDotPro2D<T>(v_a_sub, v_b_sub), k_result);
        time_kok_2 = std::chrono::high_resolution_clock::now();
        time_kok_arr[i++] = elpased_time(time_kok_2 - time_kok_1).count();
      }
      Kokkos::fence();
      Kokkos::finalize();
      printf ("Kokkos ok\n");
		}
		// CPEs, Athread
		{
   		athread_init();
      int i = 0;
      for (int n = START_NUM; n <= N; n += STRIDE) {
        // d_result = 0.0;
        I0 = n;
        I1 = n;
			  if (std::is_same<typename std::decay<T>::type, float>::value) {
          time_ath_1 = std::chrono::high_resolution_clock::now();
          athread_dot_product_2d_float (I0, I1, (float*)d_a, (float*)d_b, (float*)(&d_result));
          time_ath_2 = std::chrono::high_resolution_clock::now();
			  } else if (std::is_same<typename std::decay<T>::type, double>::value) {
          time_ath_1 = std::chrono::high_resolution_clock::now();
          athread_dot_product_2d_double (I0, I1, (double*)d_a, (double*)d_b, (double*)(&d_result));
          time_ath_2 = std::chrono::high_resolution_clock::now();
			  }
        time_ath_arr[i++] = elpased_time(time_ath_2 - time_ath_1).count();
      }
	  	athread_halt();
      printf ("Athread ok\n");
		}
    int i = 0;
    for (int n = START_NUM; n <= N; n += STRIDE) {
      file << n << "," << time_cpp_arr[i] << "," << time_kok_arr[i] << "," << time_ath_arr[i++] << std::endl;
    }

    file.close();
    */

  }
  // 释放内存
  {
    free (h_a);
    free (h_b);
    free (d_a);
    free (d_b);
  }
	return ;
}

#endif // SRC_TEST_STENCIL_HPP_