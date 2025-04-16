#include "Kokkos_Athread_ParamWrap.h"
#include "Kokkos_Athread_FastThreadSpawn.h"
#include "simd.h"
#include "athread.h"

extern void SLAVE_FUN(register_kernel)();
void launch_register_kernel () {
  // athread_spawn (register_kernel, NULL);
  // athread_join ();
  spawn_proxy_run(register_kernel, NULL);
  spawn_proxy_join();
  return ;
}

// parallel_for
extern void SLAVE_FUN (kokkos_athread_parallel_for) (struct AthreadParamWrap*);
void kokkos_athread_parallel_for_launch (struct AthreadParamWrap* param) {
  // athread_spawn (kokkos_athread_parallel_for, param);
  spawn_proxy_run (kokkos_athread_parallel_for, param);
  return ;
}

// parallel_reduce
extern void SLAVE_FUN (kokkos_athread_parallel_reduce) (struct AthreadParamWrap*);
void kokkos_athread_parallel_reduce_launch (struct AthreadParamWrap* param) {
  // athread_spawn (kokkos_athread_parallel_reduce, param);
  spawn_proxy_run (kokkos_athread_parallel_reduce, param);
  return ;
}

// extern void SLAVE_FUN(parallel_for_1D)(struct AthreadParamWrap*);
// void athread_parallel_for_launch_1D (struct AthreadParamWrap* param) {
//   spawn_proxy_run(parallel_for_1D, param);
//   return ;
// }

// extern void SLAVE_FUN(parallel_for_2D)(struct AthreadParamWrap*);
// void athread_parallel_for_launch_2D (struct AthreadParamWrap* param) {
//   spawn_proxy_run(parallel_for_2D,param);
//   spawn_proxy_join();
//   return ;
// }

// extern void SLAVE_FUN(parallel_for_3D)(struct AthreadParamWrap*);
// void athread_parallel_for_launch_3D (struct AthreadParamWrap* param) {

// //#if 1 //Kokkos_ATHREAD_FAST
//     spawn_proxy_run(parallel_for_3D,param);
//     spawn_proxy_join();
// //#else
// //	athread_spawn(parallel_for_3D, param);
// //	athread_join();
// //#endif    
// 	return ;
// }
// extern void SLAVE_FUN(parallel_for_4D)(struct AthreadParamWrap*);
// void athread_parallel_for_launch_4D (struct AthreadParamWrap* param) {

// //#if 1 //Kokkos_ATHREAD_FAST
//     spawn_proxy_run(parallel_for_4D,param);
//     spawn_proxy_join();
// //#else
// //	athread_spawn(parallel_for_4D, param);
// //	athread_join();
// //#endif
// 	return ;
// }

// // parallel_reduce
// extern void SLAVE_FUN(parallel_reduce_1D)(struct AthreadParamWrap*);
// void athread_parallel_reduce_launch_1D (struct AthreadParamWrap* param) {

// //#if 1 //Kokkos_ATHREAD_FAST
//     spawn_proxy_run(parallel_reduce_1D,param);
//     spawn_proxy_join();
// //#else
// //	athread_spawn(parallel_reduce_1D, param);
// //	athread_join();
// //#endif
// 	return ;
// }
// extern void SLAVE_FUN(parallel_reduce_2D)(struct AthreadParamWrap*);
// void athread_parallel_reduce_launch_2D (struct AthreadParamWrap* param) {

// //#if 1 //Kokkos_ATHREAD_FAST
//     spawn_proxy_run(parallel_reduce_2D,param);
//     spawn_proxy_join();
// //#else
// //	athread_spawn(parallel_reduce_2D, param);
// //	athread_join();
// //#endif
// 	return ;
// }
