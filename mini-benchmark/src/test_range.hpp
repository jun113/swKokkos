#ifndef SRC_TEST_RANGE_HPP_
#define SRC_TEST_RANGE_HPP_

#include "kokkos_axpy.hpp"

#include "Kokkos_Core.hpp"

#include <cstdio>

#include <math.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <type_traits>
#include <fstream>
#include <memory>

extern "C" void athread_axpy_float (int, float, float*, float*);
extern "C" void athread_axpy_double (int, double, double*, double*);

template<typename T>
void test_range_axpy () {

	const int TIMES = 200;
	const int START_NUM = 1e7;
  const int N = 1e8;
  const int STRIDE = (N - START_NUM) / TIMES;
	const bool PRINT_VECTOR = false;
	const bool DEBUG = false;

  size_t DATA_SIZE = N * sizeof(T);

  printf ("TIMES: %d, STRIDE = %d, N: %d\n", TIMES, STRIDE, N);
  std::cout<<"Data type: "<<typeid(T).name()<<", sizeof(T): "<<sizeof(T)<<std::endl;
  std::cout<<"Print vector: "<<PRINT_VECTOR<<std::endl;

  // ==========================
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

  T alpha = 0.3;
	// x and y in host
  T* h_x = nullptr;
  T* h_y = nullptr;

  // x and y in device
  T* d_x = nullptr;
  T* d_y = nullptr;

	// tmp y
	T* t_y = nullptr;


  // 初始化主机端A和B数据
  {
    h_x = (T*) malloc (DATA_SIZE);
    h_y = (T*) malloc (DATA_SIZE);
    d_x = (T*) malloc (DATA_SIZE);
    d_y = (T*) malloc (DATA_SIZE);
    t_y = (T*) malloc (DATA_SIZE);

    for (int i = 0; i < N; ++i) {
      h_x[i] = rand() / static_cast<T>(RAND_MAX);
      h_y[i] = rand() / static_cast<T>(RAND_MAX);
    }
  
    if (PRINT_VECTOR) {
      printf ("Vector A: \n");
      for (int i = 0; i < N; ++i) {
        printf("%.2f\t", h_x[i]);
      }
      printf("\n");
      printf ("Vector B: \n");
      for (int i = 0; i < N; ++i) {
        printf("%.2f\t", h_y[i]);
      }
      printf("\n");
    }
  }

  printf("=================================\n");

  // Kokkos code
  auto kokkos_init_args = Kokkos::InitializationSettings()
    .set_print_configuration(true)
    .set_tune_internals(true);

  Kokkos::initialize(kokkos_init_args);

  using KokView1D = Kokkos::View<T *>;
  KokView1D v_x("View Vector X", N);
  KokView1D v_y("View Vector Y", N);

  typename KokView1D::HostMirror h_v_x = Kokkos::create_mirror_view(v_x);
  typename KokView1D::HostMirror h_v_y = Kokkos::create_mirror_view(v_y);

  // 初始化设备端A和B数据
  {
    for (int i = 0; i < N; ++i) {
      h_v_x(i) = h_x[i];
      h_v_y(i) = h_y[i];
      d_x[i] = h_x[i];
      d_y[i] = h_y[i];
    }
    Kokkos::deep_copy (v_x, h_v_x);
    Kokkos::deep_copy (v_y, h_v_y);
  }
	if (DEBUG) {
		// MPE, C++
		{
      time_cpp_1 = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < N; ++i) {
        h_y[i] = alpha * h_x[i] + h_y[i];
      }
      time_cpp_2 = std::chrono::high_resolution_clock::now();
		}
		// MPE + 64 CPEs, Kokkos
		{
      time_kok_1 = std::chrono::high_resolution_clock::now();
      Kokkos::parallel_for ("AXPY", N, FunctorAXPY<T>(alpha, v_x, v_y));
      time_kok_2 = std::chrono::high_resolution_clock::now();

    	Kokkos::deep_copy (h_v_y, v_y);

			for (int i = 0; i < N; ++i) {
				t_y[i] = h_v_y(i);
			}

      Kokkos::fence();
      Kokkos::finalize();
		}
		// 64 CPEs, Athread
		{
   		athread_init();
			if (std::is_same<typename std::decay<T>::type, float>::value) {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_axpy_float (N, (float)alpha, (float*)d_x, (float*)d_y);
        time_ath_2 = std::chrono::high_resolution_clock::now();
			} else if (std::is_same<typename std::decay<T>::type, double>::value) {
        time_ath_1 = std::chrono::high_resolution_clock::now();
        athread_axpy_double (N, (double)(alpha), (double*)(d_x), (double*)(d_y));
        time_ath_2 = std::chrono::high_resolution_clock::now();
			}
	  	athread_halt();
		}
		double abs_errs = 0.0;
		for (int i = 0; i < N; ++i) {
			printf ("i=%d, CPP: %.3f, Kokkos: %.3f, Athread: %.3f\t", 
					i, h_y[i], t_y[i], d_y[i]);
			if (h_y[i] != t_y[i]) {
				abs_errs += std::abs (h_y[i] - t_y[i]);
				printf ("err !");
			}
			printf ("\n");
		}
		printf ("abs errors: %f\n", abs_errs);
    time_cpp = elpased_time(time_cpp_2 - time_cpp_1).count();
    time_ath = elpased_time(time_ath_2 - time_ath_1).count();
    time_kok = elpased_time(time_kok_2 - time_kok_1).count();
  	printf ("**********************************************\n");
  	printf ("Size, CPP, Athread, Kokkos, Speedup(Cpp/Kokkos), Efficiency(Athread/Kokkos)\n");
    printf ("%d, %.5f, %.5f, %.5f & %.3f & %.3f\n", N, time_cpp, time_ath, time_kok,
        time_cpp / time_kok, time_ath / time_kok);
	} else {
    // std::unique_ptr<double []> time_cpp_arr = std::make_unique<double []>(TIMES);
    // std::unique_ptr<double []> time_ath_arr = std::make_unique<double []>(TIMES);
    // std::unique_ptr<double []> time_kok_arr = std::make_unique<double []>(TIMES);
    std::unique_ptr<double []> time_cpp_arr = std::make_unique<double []>(TIMES + 1);
    std::unique_ptr<double []> time_ath_arr = std::make_unique<double []>(TIMES + 1);
    std::unique_ptr<double []> time_kok_arr = std::make_unique<double []>(TIMES + 1);
    std::ofstream file("axpy_results.csv");
    file << std::fixed << std::setprecision(6);
    file << "Vector size,CPP,Kokkos,Athread" << std::endl;
		// MPE, C++
		{
      // for (int i = 1; i <= TIMES; ++i) {
      //   const int n = i * START_NUM;
      int i = 0;
      for (int n = START_NUM; n <= N; n += STRIDE) {
        time_cpp_1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < n; ++i) {
          h_y[i] = alpha * h_x[i] + h_y[i];
        }
        time_cpp_2 = std::chrono::high_resolution_clock::now();
        // time_cpp_arr[i-1] = elpased_time(time_cpp_2 - time_cpp_1).count();
        time_cpp_arr[i++] = elpased_time(time_cpp_2 - time_cpp_1).count();
      }
		}
    // MPE + CPEs, Kokkos
    {
      // for (int i = 1; i <= TIMES; ++i) {
      //   const int n = i * START_NUM;
      int i = 0;
      for (int n = START_NUM; n <= N; n += STRIDE) {
        time_kok_1 = std::chrono::high_resolution_clock::now();
        Kokkos::parallel_for ("AXPY", n, FunctorAXPY<T>(alpha, v_x, v_y));
        time_kok_2 = std::chrono::high_resolution_clock::now();
        // time_kok_arr[i-1] = elpased_time(time_kok_2 - time_kok_1).count();
        time_kok_arr[i++] = elpased_time(time_kok_2 - time_kok_1).count();
      }
      Kokkos::fence();
    }
    Kokkos::finalize();
		// CPEs, Athread
		{
   		athread_init();
      // for (int i = 1; i <= TIMES; ++i) {
      //   const int n = i * START_NUM;
      int i = 0;
      for (int n = START_NUM; n <= N; n += STRIDE) {
			  if (std::is_same<typename std::decay<T>::type, float>::value) {
          time_ath_1 = std::chrono::high_resolution_clock::now();
          athread_axpy_float (n, (float)alpha, (float*)d_x, (float*)d_y);
          time_ath_2 = std::chrono::high_resolution_clock::now();
			  } else if (std::is_same<typename std::decay<T>::type, double>::value) {
          time_ath_1 = std::chrono::high_resolution_clock::now();
          athread_axpy_double (n, (double)(alpha), (double*)(d_x), (double*)(d_y));
          time_ath_2 = std::chrono::high_resolution_clock::now();
			  }
        // time_ath_arr[i-1] = elpased_time(time_ath_2 - time_ath_1).count();
        time_ath_arr[i++] = elpased_time(time_ath_2 - time_ath_1).count();
      }
	  	athread_halt();
		}
    int i = 0;
    for (int n = START_NUM; n <= N; n += STRIDE) {
    // for (int i = 1; i <= TIMES; ++i) {
    //   const int n = i * START_NUM;
      // file << n << "," << time_cpp_arr[i-1] << "," << time_kok_arr[i-1] << "," << time_ath_arr[i-1] << std::endl;
      file << n << "," << time_cpp_arr[i] << "," << time_kok_arr[i] << "," << time_ath_arr[i++] << std::endl;
    }

    file.close();
    // if (PRINTMATRIX) {
    //   printf ("Result:\n");
    //   for (int i = 0; i < N; ++i) {
    //     printf("CPP (MPE) = %.2f, Kokkos (CPEs) = %.2f\n", 
    //         h_y[i], v_y(i));
    //   }
    // }
	}

  // 释放资源
  {
    free (h_x);
    free (h_y);
    free (d_x);
    free (d_y);
    free (t_y);
#ifdef ENABLE_KOKKOS
    // Kokkos::fence();
    // Kokkos::finalize();
#endif // ENABLE_KOKKOS
  }

	return ;
}



#endif // SRC_TEST_RANGE_HPP_
