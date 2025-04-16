#ifndef KOKKOS_ATHREAD_HASHTABLE_HPP_
#define KOKKOS_ATHREAD_HASHTABLE_HPP_

#include "Athread/Kokkos_Athread_ParamWrap.h"

#if defined (__sw_host__)
#include "athread.h"
#elif defined (__sw_slave__)
#include "slave.h"
#endif

#define KOKKOS_ATHREAD_TABLE_SIZE 47

typedef struct AthreadHashNode {
  int* key;
  void (*fp)(AthreadParamWrap *, const int&);
	struct AthreadHashNode* next;
} AthreadHashNode;

extern void kokkos_athread_hash_table_init ();

extern unsigned int kokkos_athread_get_hash_host (int* const key);

extern bool kokkos_athread_hash_table_insert (AthreadHashNode* const hash_node);

extern AthreadHashNode *kokkos_athread_hash_table_lookup (AthreadParamWrap* const param);

extern void kokkos_athread_hash_table_print ();

inline unsigned int kokkos_athread_hash (int* const key) {
  unsigned int hash_value = 0;

	for (int i = 0; i < KOKKOS_ATHREAD_KEY_LEN; ++i) {
    hash_value += key[i];
    // hash_value  = (hash_value * key[i]) % KOKKOS_ATHREAD_TABLE_SIZE;
	}

  hash_value %= KOKKOS_ATHREAD_TABLE_SIZE;
  return hash_value;
}



// AthreadHashNode *hash_table_delete (int* key) {
//   int index = kokkos_athread_hash (key);
// 	AthreadHashNode *tmp = g_kokkos_athread_hash_table[index];
// 	AthreadHashNode *prev = nullptr;
// 	// while (tmp != nullptr && strncmp (tmp->key, key, MAX_NAME) != 0) {
// 	// 	prev = tmp;
// 	// 	tmp = tmp->next;
// 	// }
// 	// if (tmp == nullptr) {
// 	// 	return nullptr;
// 	// }
// 	// if (prev == nullptr) {
// 	// 	hash_table[index] = tmp->next;
// 	// } else {
// 	// 	prev->next = tmp->next;
// 	// }
// 	return tmp;
// }

#endif // KOKKOS_ATHREAD_HASHTABLE_HPP_