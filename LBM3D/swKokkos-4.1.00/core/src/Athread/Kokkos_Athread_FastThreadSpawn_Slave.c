#include "Athread/Kokkos_Athread_FastThreadSpawn.h"

#include "slave.h"

extern "C" void spawn_proxy_stub (void *_) {
  while (1) {
    struct spawn_args l_spawn_arg;
    l_spawn_arg.flag = 0;
    if (CRTS_tid == 0) {
      while (l_spawn_arg.flag == 0) {
        l_spawn_arg.flag = g_spawn_args[0].flag;
      }
      asm volatile ("memb");
      l_spawn_arg.fnptr = g_spawn_args[0].fnptr;
      l_spawn_arg.arg = g_spawn_args[0].arg;
    }
    CRTS_rma_bcast_coll(&l_spawn_arg, &l_spawn_arg, sizeof l_spawn_arg, 0);
    if (l_spawn_arg.flag == -1) {
      break;
    }
    ((void (*)(void *))l_spawn_arg.fnptr)(l_spawn_arg.arg);
    flush_slave_cache();
    CRTS_ssync_array();
    asm volatile ("memb");
    if (CRTS_tid == 0) {
      g_spawn_args[0].flag = 0;

    }
  }
  return ;
}