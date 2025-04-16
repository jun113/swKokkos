#include "Athread/Kokkos_Athread_Utils_Slave.h"
#include "Athread/Kokkos_Athread_ParamWrap.h"

#include "simd.h"
#include "slave.h"

int* arr_char_to_arr_int(const char *str) {
  // This function used SIMD for str cmp
  int key_len = 0;
  while (str[key_len] != '\0') {
    key_len += 1;
  }
  if (key_len > KOKKOS_ATHREAD_KEY_LEN) {
    printf("len=%d, %s\n", key_len, str);
    printf("Current version the length of Functor's name must be less than 32.\n");
    exit(EXIT_FAILURE);
  }

  int *arr_int __attribute__ ((aligned(64))) = nullptr;

  arr_int = (int*)libc_aligned_malloc(sizeof(int) * KOKKOS_ATHREAD_KEY_LEN);

  for (int i = 0; i < key_len; ++i) {
    arr_int[i] = str[i];
  }

  for (int i = key_len; i < KOKKOS_ATHREAD_KEY_LEN; ++i) {
    arr_int[i] = '0';
  }
  return arr_int;
}

int* arr_char_to_arr_int(const char *str, int &num_intv16) {
  // This function used SIMD for str cmp
  int key_len = 0;
  while (str[key_len] != '\0') {
    key_len += 1;
  }
  if (key_len > 32) {
    printf("len=%d, %s\n", key_len, str);
    printf("Current version the length of Functor's name must be less than 32.\n");
    exit(EXIT_FAILURE);
  }

  int *arr_int __attribute__ ((aligned(64))) = nullptr;

  if (key_len > 16) {
    arr_int = (int*)libc_aligned_malloc(sizeof(int) * 32);
  } else {
    arr_int = (int*)libc_aligned_malloc(sizeof(int) * 16);
  }

  for (int i = 0; i < key_len; ++i) {
    arr_int[i] = str[i];
  }

  if (key_len > 16) {
    for (int i = key_len; i < 32; ++i) {
      arr_int[i] = '0';
    }
    num_intv16 = 2;
  } else {
    for (int i = key_len; i < 16; ++i) {
      arr_int[i] = '0';
    }
    num_intv16 = 1;
  }
  return arr_int;
}

bool str_cmp(char* const str1, char* const str2) {
  int i1 = 0;
  int i2 = 0;
  while (str1[i1] != '\0' && str2[i2] != '\0') {
    if (str1[i1] != str2[i2]) {
      return true;
    }
    ++i1;
    ++i2;
  }
  return false;
}

__thread_local_fix intv16 result_simd;

inline bool str_cmp_intv16 (const intv16 &str1_simd, const intv16 &str2_simd) {
  // about 1 cycle
  result_simd = simd_vcmpeqw(str1_simd, str2_simd);
  // about 29 cycles
  if (simd_reduc_plusw(result_simd) == 16) {
    return true;
  } else {
    return false;
  }
}

__thread_local_fix intv16 key1_simd;
__thread_local_fix intv16 key2_simd;

bool key_cmp (int* const key1, int* const key2) {
  simd_load(key1_simd, key1);
  simd_load(key2_simd, key2);
  if (!str_cmp_intv16(key1_simd, key2_simd)) {
    // compare first intv16
    return false;
  }
  simd_load(key1_simd, &(key1[16]));
  simd_load(key2_simd, &(key2[16]));
  if (str_cmp_intv16(key1_simd, key2_simd)) {
    // compare second intv16
    return true;
  } else {
    return false;
  }
}

bool key_cmp (int* const key1, int* const key2, const int &num_intv16) {
  if (num_intv16 == 1) {
    simd_load(key1_simd, key1);
    simd_load(key2_simd, key2);
    if (str_cmp_intv16(key1_simd, key2_simd)) {
      return true;
    } else {
      return false;
    }
  } else if (num_intv16 == 2) {
    simd_load(key1_simd, key1);
    simd_load(key2_simd, key2);
    if (!str_cmp_intv16(key1_simd, key2_simd)) {
      // compare first intv16
      return false;
    }
    simd_load(key1_simd, &(key1[16]));
    simd_load(key2_simd, &(key2[16]));
    if (str_cmp_intv16(key1_simd, key2_simd)) {
      // compare second intv16
      return true;
    } else {
      return false;
    }
  }
}

FunctionPointer lookup_fp (AthreadParamWrap* param) {
  return kokkos_athread_hash_table_lookup (param)->fp;
}