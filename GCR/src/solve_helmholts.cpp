//#define CUDA
#include <iostream>
#ifdef CUDA
#include <hip/hip_runtime.h>
#include <hipblas.h>
#endif

#include "athread.h"
#include "Kokkos_Core.hpp"
#include <mpi.h>
//#include "gptl.h"
#include<sys/time.h>
#include <chrono>
//#include "inc/matrixpro.hpp"

#define TIMERSTART(tag)  auto tag##_start = std::chrono::steady_clock::now()
#define TIMEREND(tag)  auto tag##_end =  std::chrono::steady_clock::now()
#define DURATION_s(tag) printf("%s costs %d s\n",#tag,std::chrono::duration_cast<std::chrono::seconds>(tag##_end - tag##_start).count())
#define DURATION_ms(tag) \
   int c1,c2,numprocs; \
c1=std::chrono::duration_cast<std::chrono::milliseconds>(tag##_end - tag##_start).count();\
        MPI_Allreduce(&c1,&c2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);\
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);\
if(mytid[0]==0)printf("%s costs %d ms\n",#tag,c2/numprocs);

extern "C" void __module_halo_MOD_glob_updatehalo_real_3d(float* h_p, int *,int *); 
extern "C" void glob_updatehalo_p_(float *h_p, int *m, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte){
     __module_halo_MOD_glob_updatehalo_real_3d(h_p,m,iter_max);
        }

extern "C" void cu_psolve_main(float ep, float* a0, float* f0, float* cm_helm,int iter_max, double *x0, int* ijk_index, int *mytid);
extern "C" void psolve_main(float ep, float* a0, float* f0, float* cm_helm,int iter_max, double *x0, int* ijk_index, int *mytid);
extern "C" void kokkos_psolve_main(float ep, float* a0, float* f0, float* cm_helm,int iter_max, double *x0, int* ijk_index, int *mytid);
extern "C" void athread_psolve_main(float ep, float* a0, float* f0, float* cm_helm,int iter_max, double *x0, int* ijk_index, int *mytid);
extern "C" void cu_solve_helmholts_(float* a_helm, float* b_helm, float* cm_helm, float *threshold, double *pi, int *ijk_index,  int *mytid, size_t size_pi){

        int iter_max,max_iteration;
        iter_max= 5;//restart number is 5
        //max_iteration=20;
        float ep = 1e-11;
        //static int  init_flag=0;
        #ifdef CUDA
        //if ( init_flag==0) {
            int deviceCount;
            hipSetDevice(mytid[0]%4);
        hipblasHandle_t h;
        hipblasCreate(&h);
        hipblasSetPointerMode(h, HIPBLAS_POINTER_MODE_DEVICE);
    	hipDeviceSynchronize(); //test gptl time
        //}
        #endif	
 
  using elpased_time = std::chrono::duration<float, std::ratio<1, 1>>;
  auto time1 = std::chrono::high_resolution_clock::now();
  auto time2 = std::chrono::high_resolution_clock::now();
  double time_kok (0.0), time_cpp(0.0), time_ath (0.0);

        //GPTLsetoption (//GPTLoverhead, 0); // Don't print overhead estimate
        //GPTLinitialize ();               // Initialize //GPTL
    TIMERSTART(gcr_main);
        //GPTLstart ("gcr_main");  
if(mytid[0]==0) printf("begin psolve_main\n");
        // #ifdef CUDA
        // cu_psolve_main(ep, a_helm, b_helm, cm_helm, iter_max, pi, ijk_index, mytid);
        // #else
      time1 = std::chrono::high_resolution_clock::now();
      psolve_main(ep, a_helm, b_helm, cm_helm, iter_max, pi, ijk_index, mytid);
      time2 = std::chrono::high_resolution_clock::now();
      time_cpp = elpased_time(time2 - time1).count();

      athread_init();
      for (int i = 0; i < size_pi; ++i) {
        pi[i] = 1.0;
      }
      time1 = std::chrono::high_resolution_clock::now();
      athread_psolve_main(ep, a_helm, b_helm, cm_helm, iter_max, pi, ijk_index, mytid);
      time2 = std::chrono::high_resolution_clock::now();
      time_ath = elpased_time(time2 - time1).count();
      athread_halt();

      Kokkos::initialize();
      for (int i = 0; i < size_pi; ++i) {
        pi[i] = 1.0;
      }
      time1 = std::chrono::high_resolution_clock::now();
      kokkos_psolve_main(ep, a_helm, b_helm, cm_helm, iter_max, pi, ijk_index, mytid);
      time2 = std::chrono::high_resolution_clock::now();
      time_kok = elpased_time(time2 - time1).count();

        // #endif
        // iter_max = 5;
        // athread_psolve_main(ep, a_helm, b_helm, cm_helm, iter_max, pi, ijk_index, mytid);
if(mytid[0]==0) printf("end psolve_main\n");
      if (mytid[0] == 0) {
        printf ("C++: %.3f, Athread: %.3f, Kokkos: %.3f\n",
            time_cpp, time_ath, time_kok);
        printf ("C++/Athread: %.2f, C++/Kokkos: %.2f, Athread/Kokkos: %.2f\n",
            time_cpp/time_ath, time_cpp/time_kok, time_ath/time_kok);
      }
        //GPTLstop ("gcr_main");               // Stop the manual timer
    TIMEREND(gcr_main);   //GPTLfinalize();
    // DURATION_ms(gcr_main);
  Kokkos::finalize();
  return ;
}

