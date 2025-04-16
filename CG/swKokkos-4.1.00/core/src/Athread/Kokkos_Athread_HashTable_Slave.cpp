#include "Athread/Kokkos_Athread_HashTable.hpp"
#include "Athread/Kokkos_Athread_Utils_Slave.h"

#include "simd.h"
#include "slave.h"

#include <stdio.h>
#include <string.h>
#include <stdbool.h>

__thread_local_fix AthreadHashNode* g_kokkos_athread_hash_table[KOKKOS_ATHREAD_TABLE_SIZE] __attribute__((aligned(64)));

void kokkos_athread_hash_table_init () {
  for (int i = 0; i < KOKKOS_ATHREAD_TABLE_SIZE; ++i) {
    g_kokkos_athread_hash_table[i] = nullptr;
  }
  return ;
}

inline unsigned int kokkos_athread_get_hash_slave (int* const key) {
  unsigned int hash_value = 0;
#if defined (_SIMD_SLAVE)
  const int len_simd = 16;
  intv16 key_simd;
	for (int i = 0; i < KOKKOS_ATHREAD_KEY_LEN; i += len_simd) {
    simd_load (key_simd, &(key[i]));
    hash_value += simd_reduc_plusw (key_simd);
	}
#else
	for (int i = 0; i < KOKKOS_ATHREAD_KEY_LEN; ++i) {
    hash_value += key[i];
    // hash_value  = (hash_value * key[i]) % KOKKOS_ATHREAD_TABLE_SIZE;
	}
#endif
  hash_value %= KOKKOS_ATHREAD_TABLE_SIZE;
  return hash_value;
}

bool kokkos_athread_hash_table_insert (AthreadHashNode* const hash_node) {
  if (hash_node == nullptr) return false;
  // int index = kokkos_athread_hash (hash_node->key);
  int index = kokkos_athread_get_hash_slave (hash_node->key);
	hash_node->next = g_kokkos_athread_hash_table[index];
	g_kokkos_athread_hash_table[index] = hash_node;
  return true;
}

// AthreadHashNode *hash_table_lookup (int* key) {
//   unsigned int index = kokkos_athread_hash (key);
// 	AthreadHashNode *tmp = g_kokkos_athread_hash_table[index];
// 	// while (tmp != NULL && strncmp (tmp->key, key, MAX_NAME) != 0) {
// 	// 	tmp = tmp->next;
// 	// }
// 	return tmp;
// }
AthreadHashNode *kokkos_athread_hash_table_lookup (AthreadParamWrap* const param) {
	AthreadHashNode *cur = g_kokkos_athread_hash_table[param->hash_value];
  // if (athread_tid == 0) {
  //   printf ("cur: %p, %p\n", cur, cur->fp);
  // }
	while (cur != nullptr && (! key_cmp (cur->key, g_athread_functor_key))) {
		cur = cur->next;
	}
	return cur;
}

void kokkos_athread_hash_table_print () {
  printf ("Start\n");
  for (int i = 0; i < KOKKOS_ATHREAD_TABLE_SIZE; ++i) {
    if (g_kokkos_athread_hash_table[i] == nullptr) {
      printf ("\t%i\t----\n", i);
    } else {
      printf ("\t%i\t", i);
			AthreadHashNode* tmp = g_kokkos_athread_hash_table[i];
			while (tmp != nullptr) {
        printf ("-> ");
        for (int i = 0; i < KOKKOS_ATHREAD_KEY_LEN; ++i) {
          printf ("%c", tmp->key[i]);
          if (tmp->key[i] == '>') {
            printf ("%c, ", tmp->key[i+1]);
            break;
          }
        }
				printf ("\"%p\" ", tmp->fp);
				tmp = tmp->next;
			}
			printf ("\n");
    }
  }
  printf ("End\n");
  return ;
}
