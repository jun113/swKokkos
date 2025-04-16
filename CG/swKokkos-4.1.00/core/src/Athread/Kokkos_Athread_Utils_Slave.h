#ifndef KOKKOS_ATHREAD_UTILS_SLAVE_H_
#define KOKKOS_ATHREAD_UTILS_SLAVE_H_

#include "Athread/Kokkos_Athread_RegisterFunction.hpp"

#include "simd.h"

extern "C" __thread_local_fix intv16 result_simd;
extern "C" __thread_local_fix AthreadParamWrap 
    kokkos_athread_local_param __attribute__ ((aligned(64)));

extern int* g_athread_functor_key __attribute__ ((aligned(64)));

typedef void (*FunctionPointer)(AthreadParamWrap *, const int &);

extern int* arr_char_to_arr_int(const char *str);

extern bool key_cmp (int* const key1, int* const key2);
extern bool key_cmp (int* const key1, int* const key2, const int &num_intv16);

extern FunctionPointer lookup_fp (AthreadParamWrap*);

extern "C" void slave_cores_reduce(void* src_addr, void* dest_addr, 
		int units, int dtype, int optype, void *buf, int buf_item);

#endif  // KOKKOS_ATHREAD_UTILS_SLAVE_H_
