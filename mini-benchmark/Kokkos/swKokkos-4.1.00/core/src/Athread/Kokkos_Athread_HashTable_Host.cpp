#include "Athread/Kokkos_Athread_ParamWrap.h"
#include "Athread/Kokkos_Athread_HashTable.hpp"
#include "athread.h"
#include "simd.h"

unsigned int kokkos_athread_get_hash_host (int* const key) {
  unsigned int hash_value = 0;

  // printf ("key char: ");
  // for (int i = 0; i < 32; ++i) {
  //   printf ("%c ", key[i]);
  // }
  // printf ("\n");
  // printf ("key int: ");
  // for (int i = 0; i < 32; ++i) {
  //   printf ("%d ", key[i]);
  // }
  // printf ("\n");
  // printf ("================================\n");

#if defined (_SIMD_HOST)
  const int len_simd = 8;
  intv8 key_simd;
	for (int i = 0; i < KOKKOS_ATHREAD_KEY_LEN; i += len_simd) {
    simd_load (key_simd, &(key[i]));
    hash_value += simd_reduc_plusw (key_simd);
    // for (int j = 0; j < 8; ++j) {
    //   printf ("%d ", key[i + j]);
    // }
    // printf ("\n");
    // simd_print_intv8 (key_simd);
    // printf ("i = %d, hash_value = %d\n", i, hash_value);
    // printf ("------------------------\n");
	}
#else
	for (int i = 0; i < KOKKOS_ATHREAD_KEY_LEN; ++i) {
    hash_value += key[i];
	}
#endif

  // hash_value = 0;
	// for (int i = 0; i < KOKKOS_ATHREAD_KEY_LEN; ++i) {
  //   hash_value += key[i];
  //   // hash_value  = (hash_value * key[i]) % KOKKOS_ATHREAD_TABLE_SIZE;
	// }

  hash_value %= KOKKOS_ATHREAD_TABLE_SIZE;
  return hash_value;
}