// __INCLUDE_FUNCTOR_HPP_START__
#include "/home/export/base/shisuan/swiap/online/wjl/swKokkos-paper/dev/PCG/SWCG-WJL/src/inc/functor.hpp"
// __INCLUDE_FUNCTOR_HPP_END__

#include "Athread/Kokkos_Athread_Utils_Slave.h"
#include "Athread/Kokkos_Athread_HashTable.hpp"

#include "Athread/Kokkos_Athread_RegisterFunction.hpp"
#include "Athread/Kokkos_Athread_ParamWrap.h"
#include "Athread/Kokkos_Athread_FastThreadSpawn.h"

#include "simd.h"
#include "slave.h"

#include <string>

#include <time.h>

//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT  = 0;
__thread_local crts_rply_t  dma_rply = 0;
__thread_local crts_rply_t  l_rply   = 0;
__thread_local crts_rply_t  r_rply   = 0;

__thread_local_fix AthreadParamWrap kokkos_athread_local_param __attribute__ ((aligned(64)));

// Buffer for the operator of a reduction between slave cores
__thread_local double buf_reduce __attribute__ ((aligned(64)));
__thread_local int buf_reduce_int __attribute__ ((aligned(64)));


extern "C" void register_kernel() {

  kokkos_athread_hash_table_init ();
  // __REGISTER_START__
  AthreadHashNode* func0 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func0->key = arr_char_to_arr_int (
      "FunctorSpmv1");
  func0->fp = spmv_1D;
  kokkos_athread_hash_table_insert (func0);
  AthreadHashNode* func1 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func1->key = arr_char_to_arr_int (
      "FunctorAXPY<float>1");
  func1->fp = axpy_f_1D;
  kokkos_athread_hash_table_insert (func1);
  AthreadHashNode* func2 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func2->key = arr_char_to_arr_int (
      "FunctorAXPY<double>1");
  func2->fp = axpy_d_1D;
  kokkos_athread_hash_table_insert (func2);
  AthreadHashNode* func3 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func3->key = arr_char_to_arr_int (
      "MyReduce1");
  func3->fp = my_reduce_1D;
  kokkos_athread_hash_table_insert (func3);
  // __REGISTER_END__
  
  // if (athread_tid == 0) {
  //   kokkos_athread_hash_table_print (); 
  // }
  return ;
}

inline void check_function_pointer (const void* fp, const char* func_name) {
  if (fp == nullptr) {
    if (athread_tid == 0) {
      // printf ("Error in \"%s\" -> \"", __PRETTY_FUNCTION__);
      printf ("Error in \"%s\" -> \"", func_name);
      for (int i = 0; i < KOKKOS_ATHREAD_KEY_LEN; ++i) {
        printf ("%c", g_athread_functor_key[i]);
      }
      printf ("\". Hash value: %d\n", kokkos_athread_local_param.hash_value);
    }
    athread_ssync_array(); 
    exit (EXIT_FAILURE);
  }
  return ;
}

extern "C" void kokkos_athread_parallel_for (AthreadParamWrap* host_param) {
  athread_dma_iget(&kokkos_athread_local_param, host_param, 
      sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);
  const FunctionPointer fp __attribute__ ((aligned(64))) = 
      kokkos_athread_hash_table_lookup(&kokkos_athread_local_param)->fp;

  check_function_pointer ((void*)fp, __PRETTY_FUNCTION__);

  fp(&kokkos_athread_local_param, athread_tid);
  // if (athread_tid == 0) {
  //   printf("launch fp: %p\n", fp);
  // }
  return ;
}

extern "C" void kokkos_athread_parallel_reduce (AthreadParamWrap* host_param) {
  athread_dma_iget(&kokkos_athread_local_param, host_param, 
      sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);
  const FunctionPointer fp __attribute__ ((aligned(64))) = 
      kokkos_athread_hash_table_lookup(&kokkos_athread_local_param)->fp;
  check_function_pointer ((void*)fp, __PRETTY_FUNCTION__);
  fp(&kokkos_athread_local_param, athread_tid);
  // reduce from all slave cores
  slave_cores_reduce(&(kokkos_athread_local_param.reduce_double),
      &(kokkos_athread_local_param.reduce_double), 1, athread_double, OP_add, &buf_reduce, 1);
  // slave_cores_reduce(&(kokkos_athread_local_param.reduce_int),
  //     &(kokkos_athread_local_param.reduce_int), 1, athread_int, OP_and, &buf_reduce_int, 1);
  if (athread_tid == 0) {
    host_param->reduce_int = kokkos_athread_local_param.reduce_int;
    host_param->reduce_double = kokkos_athread_local_param.reduce_double;
  }
  return ;
}

// extern "C" void parallel_reduce_1D(AthreadParamWrap* host_param) {

//   athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
//   D_COUNT++;
//   athread_dma_wait_value(&dma_rply, D_COUNT);

//   const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

//   const int ATHREAD_SLAVE_CORES = 64;

//   const int start = kokkos_athread_local_param.range[0][0];
//   const int end   = kokkos_athread_local_param.range[0][1];
//   const int len   = end - start;
//   const int times = (len + ATHREAD_SLAVE_CORES - 1) / ATHREAD_SLAVE_CORES;

//   for (int i = 0; i < times; ++i) {
//     int index = start + athread_tid * times + i;
//     if (index < end) {
//       kokkos_athread_local_param.index[0] = index;
//       fp(&kokkos_athread_local_param);
//     }
//   }
//   // reduce from all slave cores
//   slave_cores_reduce(&(kokkos_athread_local_param.reduce_result),
//       &(kokkos_athread_local_param.reduce_result), 1, athread_double, OP_add, &buf_reduce, 1);
//   if (athread_tid == 0) {
//     host_param->reduce_result = kokkos_athread_local_param.reduce_result;
//   }
//   return ;
// }

// extern "C" void parallel_reduce_2D(AthreadParamWrap* host_param) {

//   athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
//   D_COUNT++;
//   athread_dma_wait_value(&dma_rply, D_COUNT);

//   const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

//   const int ATHREAD_SLAVE_CORES = 64;

//   const int start0 = kokkos_athread_local_param.range[0][0];
//   const int end0   = kokkos_athread_local_param.range[0][1];
//   const int start1 = kokkos_athread_local_param.range[1][0];
//   const int end1   = kokkos_athread_local_param.range[1][1];
//   const int tile0  = kokkos_athread_local_param.tile[0];
//   const int tile1  = kokkos_athread_local_param.tile[1];

//   const int num_tiles0 = (end0 - start0 + tile0 - 1) / tile0;
//   const int num_tiles1 = (end1 - start1 + tile1 - 1) / tile1;

//   const int num_tiles = num_tiles0 * num_tiles1;
  
//   for (int index_tile = athread_tid; index_tile < num_tiles; 
//       index_tile += ATHREAD_SLAVE_CORES) {
//     const int index_tile0 = index_tile / num_tiles1;
//     const int index_tile1 = index_tile % num_tiles1;

//     for (int index00 = 0; index00 < tile0; ++ index00) {
//       const int index0 = start0 + index_tile0 * tile0 + index00;
//       if (index0 >= end0) { break; }
//       for (int index_11 = 0; index_11 < tile1; ++ index_11) {
//         const int index1 = start1 + index_tile1 * tile1 + index_11;
//         if (index1 >= end1) { break; }
//         kokkos_athread_local_param.index[0] = index0;
//         kokkos_athread_local_param.index[1] = index1;
//         fp(&kokkos_athread_local_param);
//       }
//     }
//   }
//   // reduce from all slave cores
//   slave_cores_reduce(&(kokkos_athread_local_param.reduce_result),
//       &(kokkos_athread_local_param.reduce_result), 1, athread_double, OP_add, &buf_reduce, 1);
//   if (athread_tid == 0) {
//     host_param->reduce_result = kokkos_athread_local_param.reduce_result;
//   }
//   return ;
// }

// extern "C" void parallel_reduce_3D(AthreadParamWrap* host_param) {

//   athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
//   D_COUNT++;
//   athread_dma_wait_value(&dma_rply, D_COUNT);

//   const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

//   const int ATHREAD_SLAVE_CORES = 64;

//   const int start[3] = {kokkos_athread_local_param.range[0][0], kokkos_athread_local_param.range[1][0], kokkos_athread_local_param.range[2][0]};
//   const int end[3]   = {kokkos_athread_local_param.range[0][1], kokkos_athread_local_param.range[1][1], kokkos_athread_local_param.range[2][1]};
//   const int tile[3]  = {kokkos_athread_local_param.tile[0],     kokkos_athread_local_param.tile[1],     kokkos_athread_local_param.tile[2]};

//   const int num_tiles0 = (end[0] - start[0] + tile[0] - 1) / tile[0];
//   const int num_tiles1 = (end[1] - start[1] + tile[1] - 1) / tile[1];
//   const int num_tiles2 = (end[2] - start[2] + tile[2] - 1) / tile[2];

//   const int tmp0 = num_tiles1 * num_tiles2;
//   const int num_tiles = num_tiles0 * tmp0;
  
//   for (int index_tile = athread_tid; index_tile < num_tiles; index_tile += ATHREAD_SLAVE_CORES) {

//     const int index_tile0 = index_tile / tmp0;
//     const int tmp1        = index_tile % tmp0;
//     const int index_tile1 = tmp1      / num_tiles2;
//     const int index_tile2 = tmp1      % num_tiles2;

//     for (int index_00 = 0; index_00 < tile[0]; ++ index_00) {
//       const int index0 = start[0] + index_tile0 * tile[0] + index_00;
//       if (index0 >= end[0]) { break; }
//       for (int index_11 = 0; index_11 < tile[1]; ++ index_11) {
//         const int index1 = start[1] + index_tile1 * tile[1] + index_11;
//         if (index1 >= end[1]) { break; }
//         for (int index_22 = 0; index_22 < tile[2]; ++ index_22) {
//           const int index2 = start[2] + index_tile2 * tile[2] + index_22;
//           if (index2 >= end[2]) { break; }

//           kokkos_athread_local_param.index[0] = index0;
//           kokkos_athread_local_param.index[1] = index1;
//           kokkos_athread_local_param.index[2] = index2;
//           fp(&kokkos_athread_local_param);
//         }
//       }
//     }
//   }
//   // reduce from all slave cores
//   slave_cores_reduce(&(kokkos_athread_local_param.reduce_result),
//       &(kokkos_athread_local_param.reduce_result), 1, athread_double, OP_add, &buf_reduce, 1);
//   if (athread_tid == 0) {
//     host_param->reduce_result = kokkos_athread_local_param.reduce_result;
//   }
//   return ;
// }
