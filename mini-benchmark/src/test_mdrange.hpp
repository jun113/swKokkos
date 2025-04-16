#ifndef SRC_TEST_MDRANGE_HPP_
#define SRC_TEST_MDRANGE_HPP_

#include "kokkos_matrix_add.hpp"
#include "kokkos_3d_mult.hpp"
#include "kokkos_dot_product.hpp"

#include "Kokkos_Core.hpp"

#include <cstdio>

#include <math.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <type_traits>
#include <fstream>
#include <memory>

extern "C" void athread_matrix_add_double (int I0, int I1, double* a, double* b, double* c);
extern "C" void athread_matrix_add_float (int I0, int I1, float* a, float* b, float* c);

extern "C" void athread_3d_mult_double (int I0, int I1, int I2, double* a, double* b, double* c);
extern "C" void athread_3d_mult_float (int I0, int I1, int I2, float* a, float* b, float* c);
template<typename T>
void test_mdrange_for_2d () {

	// const int I0 = 267;
	// const int I1 = 420;

	const int I0 = 3;
	const int I1 = 153;

	const int N = I0 * I1;
  size_t DATA_SIZE = N * sizeof(T);

	const bool PRINT_MATRIX = false;
	const bool DEBUG = true;

  // printf ("TIMES: %d, N: %d\n", TIMES, N);
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

  using KokView = Kokkos::View<T **>;
  KokView v_a("View Matrix A", I0, I1);
  KokView v_b("View Matrix B", I0, I1);
  KokView v_c("View Matrix C", I0, I1);

  typename KokView::HostMirror h_v_a = Kokkos::create_mirror_view(v_a);
  typename KokView::HostMirror h_v_b = Kokkos::create_mirror_view(v_b);
  typename KokView::HostMirror h_v_c = Kokkos::create_mirror_view(v_c);

  // 初始化设备端A和B数据
  {
		for (int i0 = 0; i0 < I0; ++i0) {
		  for (int i1 = 0; i1 < I1; ++i1) {
      	h_v_a(i0, i1) = h_a[i0*I1 + i1];
      	h_v_b(i0, i1) = h_b[i0*I1 + i1];
      	d_a[i0*I1 + i1] = h_a[i0*I1 + i1];
      	d_b[i0*I1 + i1] = h_b[i0*I1 + i1];
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
					h_c[i0*I1 + i1] = h_a[i0*I1 + i1] + h_b[i0*I1 + i1];
				}
			}
      time_cpp_2 = std::chrono::high_resolution_clock::now();
		}
		// MPE + 64 CPEs, Kokkos
		{
    	using Kokkos::parallel_for;
    	using Kokkos::MDRangePolicy;
    	using KoArr2D = Kokkos::Array<int64_t, 2>;
    	auto Tile2D   = KoArr2D{3, 2};
      time_kok_1 = std::chrono::high_resolution_clock::now();
    	parallel_for ("Matrix addition", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{0, 0}, KoArr2D{I0, I1}, Tile2D), FunctorMatrixAdd<T>(v_a, v_b, v_c));
      time_kok_2 = std::chrono::high_resolution_clock::now();

    	Kokkos::deep_copy (h_v_c, v_c);

			for (int i0 = 0; i0 < I0; ++i0) {
			  for (int i1 = 0; i1 < I1; ++i1) {
					t_c[i0*I1 + i1] = h_v_c(i0, i1);
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
        athread_matrix_add_float (I0, I1, (float*)d_a, (float*)d_b, (float*)d_c);
        time_ath_2 = std::chrono::high_resolution_clock::now();
			} else if (std::is_same<typename std::decay<T>::type, double>::value) {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_matrix_add_double (I0, I1, (double*)d_a, (double*)d_b, (double*)d_c);
        time_ath_2 = std::chrono::high_resolution_clock::now();
			}
	  	athread_halt();
		}
		double abs_errs = 0.0;
		for (int i0 = 0; i0 < I0; ++i0) {
			for (int i1 = 0; i1 < I1; ++i1) {
        const int i = i0*I1 + i1;
				// printf ("i0=%d i1=%d, CPP: %.3f, Kokkos: %.3f, Athread: %.3f\t", 
				// 		i0, i1, h_c[i], t_c[i], d_c[i]);
				if (h_c[i] != t_c[i]) {
					abs_errs += std::abs (h_c[i] - t_c[i]);
					// printf ("err !");
				}
				// printf ("\n");
      }
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
  // 释放内存
  {
    free (h_a);
    free (h_b);
    free (h_c);
    free (d_a);
    free (d_b);
    free (d_c);
    free (t_c);
  }

	return ;
}

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

extern "C" void athread_dot_product_2d_double (int I0, int I1, double* a, double* b, double* result);
extern "C" void athread_dot_product_2d_float (int I0, int I1, float* a, float* b, float* result);

template<typename T>
void test_dot_product_2d () {

	// const int I0 = 267;
	// const int I1 = 420;


	// const int N = I0 * I1;
	const int TIMES = 200;
	const int START_NUM = 5e3;
  const int N = 1e4;
  const int STRIDE = (N - START_NUM) / TIMES;
  size_t DATA_SIZE = N * N * sizeof(T);

	const bool PRINT_MATRIX = false;
	const bool DEBUG = false;

  int I0 = N;
  int I1 = N;

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

	T* d_a = nullptr;
	T* d_b = nullptr;

  T h_result = 0.0;
  T k_result = 0.0;
  T d_result = 0.0;

  // Kokkos code
  auto kokkos_init_args = Kokkos::InitializationSettings()
    .set_print_configuration(true)
    .set_tune_internals(true);

  Kokkos::initialize(kokkos_init_args);
  using KokView = Kokkos::View<T **>;

  // 初始化主机端A和B数据
  {
    h_a = (T*) malloc (DATA_SIZE);
    h_b = (T*) malloc (DATA_SIZE);
    d_a = (T*) malloc (DATA_SIZE);
    d_b = (T*) malloc (DATA_SIZE);
  
    for (int i = 0; i < N*N; ++i) {
      h_a[i] = (rand()) / static_cast<T>(RAND_MAX);
      h_b[i] = (rand()) / static_cast<T>(RAND_MAX);
    }
    
    if (PRINT_MATRIX) {
      printf ("Matrix A: \n");
      for (int i = 0; i < N*N; ++i) {
        printf("%.2f\t", h_a[i]);
      }
      printf("\n");
      printf ("Matrix B: \n");
      for (int i = 0; i < N*N; ++i) {
        printf("%.2f\t", h_b[i]);
      }
      printf("\n");
    }
  }
  
  KokView v_a("View Matrix A", I0, I1);
  KokView v_b("View Matrix B", I0, I1);
  
  typename KokView::HostMirror h_v_a = Kokkos::create_mirror_view(v_a);
  typename KokView::HostMirror h_v_b = Kokkos::create_mirror_view(v_b);
  
  // 初始化设备端A和B数据
  {
		for (int i0 = 0; i0 < I0; ++i0) {
		  for (int i1 = 0; i1 < I1; ++i1) {
    	  h_v_a(i0, i1) = h_a[i0*I1 + i1];
      	h_v_b(i0, i1) = h_b[i0*I1 + i1];
      	d_a[i0*I1 + i1] = h_a[i0*I1 + i1];
      	d_b[i0*I1 + i1] = h_b[i0*I1 + i1];
			}
		}
    Kokkos::deep_copy (v_a, h_v_a);
    Kokkos::deep_copy (v_b, h_v_b);
  }

	if (DEBUG) {
		// MPE, C++
		{
      h_result = 0.0;
      time_cpp_1 = std::chrono::high_resolution_clock::now();
			for (int i0 = 0; i0 < I0; ++i0) {
			  for (int i1 = 0; i1 < I1; ++i1) {
					h_result += h_a[i0*I1 + i1] * h_b[i0*I1 + i1];
				}
			}
      time_cpp_2 = std::chrono::high_resolution_clock::now();
		}
		// MPE + 64 CPEs, Kokkos
		{
    	using Kokkos::MDRangePolicy;
    	using KoArr2D = Kokkos::Array<int64_t, 2>;
    	auto Tile2D   = KoArr2D{2, 64};
      time_kok_1 = std::chrono::high_resolution_clock::now();
    	Kokkos::parallel_reduce ("dot product 2d", MDRangePolicy<Kokkos::Rank<2>> (
        	KoArr2D{0, 0}, KoArr2D{I0, I1}, Tile2D), FunctorDotPro2D<T>(v_a, v_b), k_result);
      time_kok_2 = std::chrono::high_resolution_clock::now();

      Kokkos::fence();
      Kokkos::finalize();
		}
		// 64 CPEs, Athread
		{
   		athread_init();
      d_result = 0.0;
			if (std::is_same<typename std::decay<T>::type, float>::value) {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_dot_product_2d_float (I0, I1, (float*)d_a, (float*)d_b, (float*)(&d_result));
        time_ath_2 = std::chrono::high_resolution_clock::now();
			} else if (std::is_same<typename std::decay<T>::type, double>::value) {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_dot_product_2d_double (I0, I1, (double*)d_a, (double*)d_b, (double*)(&d_result));
        time_ath_2 = std::chrono::high_resolution_clock::now();
			}
	  	athread_halt();
		}
    printf ("CPP: %.3f, Kokkos: %.3f, Raw kernel: %.3f\n", h_result, k_result, d_result);
    time_cpp = elpased_time(time_cpp_2 - time_cpp_1).count();
    time_ath = elpased_time(time_ath_2 - time_ath_1).count();
    time_kok = elpased_time(time_kok_2 - time_kok_1).count();
  	printf ("**********************************************\n");
  	printf ("Size, CPP, Athread, Kokkos, Speedup(Cpp/Kokkos), Efficiency(Athread/Kokkos)\n");
    printf ("%d, %.5f, %.5f, %.5f, %.3f, %.3f\n", N, time_cpp, time_ath, time_kok,
        time_cpp / time_kok, time_ath / time_kok);
	} else {
    std::unique_ptr<double []> time_cpp_arr = std::make_unique<double []>(TIMES + 1);
    std::unique_ptr<double []> time_ath_arr = std::make_unique<double []>(TIMES + 1);
    std::unique_ptr<double []> time_kok_arr = std::make_unique<double []>(TIMES + 1);
    std::unique_ptr<double []> time_kok_sr_arr = std::make_unique<double []>(TIMES + 1);
    std::ofstream file("dot_product_results.csv");
    file << std::fixed << std::setprecision(6);
    file << "Matrix size,CPP,Kokkos,Athread,Kokkos-SR" << std::endl;
		// MPE, C++
		{
      int i = 0;
      for (int n = START_NUM; n <= N; n += STRIDE) {
        h_result = 0.0;
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
        k_result = 0.0;
        I0 = n;
        I1 = n;
        auto v_a_sub = Kokkos::subview(v_a, Kokkos::make_pair(0, I0), Kokkos::make_pair(0, I1));
        auto v_b_sub = Kokkos::subview(v_b, Kokkos::make_pair(0, I0), Kokkos::make_pair(0, I1));
      	// auto Tile2D = KoArr2D{2, 64};
        /*
        C++/Kokkos, max:  7.165649499457243 ; min:  5.835753026899688 ; mean:  6.500282354508347
        Athread/Kokkos, max:  0.8733265073508104 ; min:  0.7255465334166146 ; mean:  0.8095092475411221
        */
      	// auto Tile2D = KoArr2D{1, (I1+64)/65};
        /*
        C++/Kokkos, max:  7.240779755418331 ; min:  5.876027830487033 ; mean:  6.552180145548028
        Athread/Kokkos, max:  0.8818196029388776 ; min:  0.7277237048665619 ; mean:  0.8152050536854023
        */
      	// auto Tile2D = KoArr2D{2, (I1+64)/65};
        /*
        C++/Kokkos, max:  7.236814064997337 ; min:  5.940151977369793 ; mean:  6.583136882044338
        Athread/Kokkos, max:  0.8852742077896595 ; min:  0.7317057784506686 ; mean:  0.8176878067376913
        */
      	// auto Tile2D = KoArr2D{4, (I1+64)/65};
        /*
        C++/Kokkos, max:  7.298149766899766 ; min:  5.899489965628119 ; mean:  6.5775620030032735
        Athread/Kokkos, max:  0.888510051241624 ; min:  0.733341754563254 ; mean:  0.8184620431756622
        */
      	// auto Tile2D = KoArr2D{8, (I1+64)/65};
        /*
        C++/Kokkos, max:  7.822340530415064 ; min:  7.253524025105464 ; mean:  7.508374378758875
        Athread/Kokkos, max:  0.9458006085049177 ; min:  0.9139314744315258 ; mean:  0.9331485848433432
        */
      	// auto Tile2D = KoArr2D{8, 256};
        /*
        C++/Kokkos, max:  7.299687010954616 ; min:  6.846000740603592 ; mean:  7.0070103587737
        Athread/Kokkos, max:  0.8796187582015332 ; min:  0.854358527845685 ; mean:  0.8708860080864806
        */
      	// auto Tile2D = KoArr2D{8, 128};
        /*
        C++/Kokkos, max:  8.088552054248105 ; min:  7.39142091152815 ; mean:  7.744563217913966
        Athread/Kokkos, max:  0.9742387649760856 ; min:  0.9257756655461534 ; mean:  0.9622044601499248
        */
      	// auto Tile2D = KoArr2D{8, 512};
        /*
        C++/Kokkos, max:  8.19208633093525 ; min:  7.295465094689936 ; mean:  7.837051997888109
        Athread/Kokkos, max:  0.9886748992257982 ; min:  0.9162424270346226 ; mean:  0.9744892403637877
        */
      	// auto Tile2D = KoArr2D{8, 1024};
        /*
        C++/Kokkos, max:  8.206660006660007 ; min:  6.76500447894894 ; mean:  7.7003066692911
        Athread/Kokkos, max:  0.9925810230378758 ; min:  0.8505819182861348 ; mean:  0.9586783588278902
        */
      	// auto Tile2D = KoArr2D{8, 2048};
        /*
        C++/Kokkos, max:  7.299045521292217 ; min:  5.909333333333334 ; mean:  6.608965433100767
        Athread/Kokkos, max:  0.8903135586543764 ; min:  0.7343355229914098 ; mean:  0.8211541807064793
        */
      	// auto Tile2D = KoArr2D{(I1+64)/65, (I1+64)/65};
        /*
        C++/Kokkos, max:  8.193599477097894 ; min:  7.298170901982686 ; mean:  7.839369087941249
        Athread/Kokkos, max:  0.9885363768917492 ; min:  0.9163919688342822 ; mean:  0.9759811369437982 
        */
      	auto Tile2D = KoArr2D{1, N};
        time_kok_1 = std::chrono::high_resolution_clock::now();
      	Kokkos::parallel_reduce ("dot product 2d", MDRangePolicy<Kokkos::Rank<2>> (
          	KoArr2D{0, 0}, KoArr2D{I0, I1}, Tile2D), FunctorDotPro2D<T>(v_a_sub, v_b_sub), k_result);
        time_kok_2 = std::chrono::high_resolution_clock::now();
        time_kok_arr[i++] = elpased_time(time_kok_2 - time_kok_1).count();
      }
      i = 0;
      for (int n = START_NUM; n <= N; n += STRIDE) {
        k_result = 0.0;
        I0 = n;
        I1 = n;
        auto v_a_sub = Kokkos::subview(v_a, Kokkos::make_pair(0, I0), Kokkos::make_pair(0, I1));
        auto v_b_sub = Kokkos::subview(v_b, Kokkos::make_pair(0, I0), Kokkos::make_pair(0, I1));
        time_kok_1 = std::chrono::high_resolution_clock::now();
      	Kokkos::parallel_reduce ("dot product 1d", I0*I1, FunctorDotPro1D<T>(v_a_sub, v_b_sub), k_result);
        time_kok_2 = std::chrono::high_resolution_clock::now();
        time_kok_sr_arr[i++] = elpased_time(time_kok_2 - time_kok_1).count();
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
        d_result = 0.0;
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
      // file << n << "," << time_cpp_arr[i] << "," << time_kok_arr[i] << "," << time_ath_arr[i++] << std::endl;
      file << n << "," << time_cpp_arr[i] << "," << time_kok_arr[i] << "," << time_ath_arr[i] << "," << time_kok_sr_arr[i++] << std::endl;
    }

    file.close();

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

#endif // SRC_TEST_MDRANGE_HPP_