#ifdef __sw_host__

#include "../../kokkos/swKokkos-4.1.00-opt/core/src/Athread/Kokkos_Athread_FastThreadSpawn.h"
#include "athread_param.h"

#include "athread.h"

extern void SLAVE_FUN(athread_get_halo_transpose_double_slave)(struct HaloTransposeDouble3D*);
extern void SLAVE_FUN(athread_put_halo_transpose_double_slave)(struct HaloTransposeDouble3D*);

extern void SLAVE_FUN(athread_get_halo_transpose_double2float_slave)(struct HaloTransposeDouble2Float3D*);
extern void SLAVE_FUN(athread_put_halo_transpose_float2double_slave)(struct HaloTransposeFloat2Double3D*);

void athread_get_halo_transpose_double_host (double* arrSrc, double* arrObj,
    int startLayer, int lenLayer, int lenA, int lenB, int lenC) {
  // int startB, endB, startC, endC;
  // for (int a = 0; a < lenA; ++a) {

    // startB = startLayer;
    // endB   = lenB - startLayer;
    // startC = startLayer;
    // endC   = startLayer + lenLayer;
    // for (int b = startB; b < endB; ++b) {
    //   for (int c = startC; c < endC; ++c) {
    //     arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
    //   }
    // }

  //   startB = startLayer;
  //   endB   = startLayer + lenLayer;
  //   startC = startLayer + lenLayer;
  //   endC   = lenC - startLayer - lenLayer;
  //   for (int b = startB; b < endB; ++b) {
  //     for (int c = startC; c < endC; ++c) {
  //       arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
  //     }
  //   }

  //   startB = lenB - startLayer - lenLayer;
  //   endB   = lenB - startLayer;
  //   startC = startLayer + lenLayer;
  //   endC   = lenC - startLayer - lenLayer;
  //   for (int b = startB; b < endB; ++b) {
  //     for (int c = startC; c < endC; ++c) {
  //       arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
  //     }
  //   }

  //   startB = startLayer;
  //   endB   = lenB - startLayer;
  //   startC = lenC - startLayer - lenLayer;
  //   endC   = lenC - startLayer;
  //   for (int b = startB; b < endB; ++b) {
  //     for (int c = startC; c < endC; ++c) {
  //       arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
  //     }
  //   }
  // }

  // for (int a = 0; a < lenA; ++a) {

  //   startB = startLayer;
  //   endB   = startLayer + lenLayer;
  //   startC = startLayer;
  //   endC   = lenC - startLayer;
  //   for (int b = startB; b < endB; ++b) {
  //     for (int c = startC; c < endC; ++c) {
  //       arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
  //     }
  //   }

  //   startB = startLayer + lenLayer;
  //   endB   = lenB - startLayer - lenLayer;
  //   startC = startLayer;
  //   endC   = startLayer + lenLayer;
  //   for (int b = startB; b < endB; ++b) {
  //     for (int c = startC; c < endC; ++c) {
  //       arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
  //     }
  //   }

  //   startB = startLayer + lenLayer;
  //   endB   = lenB - startLayer - lenLayer;
  //   startC = lenC - startLayer - lenLayer;
  //   endC   = lenC - startLayer;
  //   for (int b = startB; b < endB; ++b) {
  //     for (int c = startC; c < endC; ++c) {
  //       arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
  //     }
  //   }

  //   startB = lenB - startLayer - lenLayer;
  //   endB   = lenB - startLayer;
  //   startC = startLayer;
  //   endC   = lenC - startLayer;
  //   for (int b = startB; b < endB; ++b) {
  //     for (int c = startC; c < endC; ++c) {
  //       arrObj[b*lenC*lenA + c*lenA + a] = arrSrc[a*lenB*lenC + b*lenC + c];
  //     }
  //   }
  // }

  struct HaloTransposeDouble3D param;
  // param.startB[0]  = startLayer;
  // param.startB[1]  = startLayer;
  // param.startB[2]  = lenB - startLayer - lenLayer;
  // param.startB[3]  = startLayer;
  // param.endB[0]    = lenB - startLayer;
  // param.endB[1]    = startLayer + lenLayer;
  // param.endB[2]    = lenB - startLayer;
  // param.endB[3]    = lenB - startLayer;
  // param.startC[0]  = startLayer;
  // param.startC[1]  = startLayer + lenLayer;
  // param.startC[2]  = startLayer + lenLayer;
  // param.startC[3]  = lenC - startLayer - lenLayer;
  // param.endC[0]    = startLayer + lenLayer;
  // param.endC[1]    = lenC - startLayer - lenLayer;
  // param.endC[2]    = lenC - startLayer - lenLayer;
  // param.endC[3]    = lenC - startLayer;
  // Priority C
  param.startB[0] = startLayer;
  param.startB[1] = startLayer + lenLayer;
  param.startB[2] = startLayer + lenLayer;
  param.startB[3] = lenB - startLayer - lenLayer;
  param.endB[0]   = startLayer + lenLayer;
  param.endB[1]   = lenB - startLayer - lenLayer;
  param.endB[2]   = lenB - startLayer - lenLayer;
  param.endB[3]   = lenB - startLayer;
  param.startC[0] = startLayer;
  param.startC[1] = startLayer;
  param.startC[2] = lenC - startLayer - lenLayer;
  param.startC[3] = startLayer;
  param.endC[0]   = lenC - startLayer;
  param.endC[1]   = startLayer + lenLayer;
  param.endC[2]   = lenC - startLayer;
  param.endC[3]   = lenC - startLayer;

  param.arrSrc = arrSrc;
  param.arrObj = arrObj;
  param.lenA   = lenA;
  param.lenB   = lenB;
  param.lenC   = lenC;

	// athread_spawn(athread_get_halo_transpose_double_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_get_halo_transpose_double_slave, &param);
	spawn_proxy_join();

  return ;
}

void athread_put_halo_transpose_double_host (double* arrSrc, double* arrObj,
    int startLayer, int lenLayer, int lenA, int lenB, int lenC) {

  // Priority C
  struct HaloTransposeDouble3D param;
  param.startB[0] = startLayer;
  param.startB[1] = startLayer + lenLayer;
  param.startB[2] = startLayer + lenLayer;
  param.startB[3] = lenB - startLayer - lenLayer;
  param.endB[0]   = startLayer + lenLayer;
  param.endB[1]   = lenB - startLayer - lenLayer;
  param.endB[2]   = lenB - startLayer - lenLayer;
  param.endB[3]   = lenB - startLayer;
  param.startC[0] = startLayer;
  param.startC[1] = startLayer;
  param.startC[2] = lenC - startLayer - lenLayer;
  param.startC[3] = startLayer;
  param.endC[0]   = lenC - startLayer;
  param.endC[1]   = startLayer + lenLayer;
  param.endC[2]   = lenC - startLayer;
  param.endC[3]   = lenC - startLayer;

  param.arrSrc = arrSrc;
  param.arrObj = arrObj;
  param.lenA   = lenA;
  param.lenB   = lenB;
  param.lenC   = lenC;

	// athread_spawn(athread_put_halo_transpose_double_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_put_halo_transpose_double_slave, &param);
	spawn_proxy_join();

  return ;
}
void athread_get_halo_transpose_double2float_host (double* arrSrc, float* arrObj,
    int startLayer, int lenLayer, int lenA, int lenB, int lenC) {

  struct HaloTransposeDouble2Float3D param;

  param.startB[0] = startLayer;
  param.startB[1] = startLayer + lenLayer;
  param.startB[2] = startLayer + lenLayer;
  param.startB[3] = lenB - startLayer - lenLayer;
  param.endB[0]   = startLayer + lenLayer;
  param.endB[1]   = lenB - startLayer - lenLayer;
  param.endB[2]   = lenB - startLayer - lenLayer;
  param.endB[3]   = lenB - startLayer;
  param.startC[0] = startLayer;
  param.startC[1] = startLayer;
  param.startC[2] = lenC - startLayer - lenLayer;
  param.startC[3] = startLayer;
  param.endC[0]   = lenC - startLayer;
  param.endC[1]   = startLayer + lenLayer;
  param.endC[2]   = lenC - startLayer;
  param.endC[3]   = lenC - startLayer;

  param.arrSrc = arrSrc;
  param.arrObj = arrObj;
  param.lenA   = lenA;
  param.lenB   = lenB;
  param.lenC   = lenC;

	// athread_spawn(athread_get_halo_transpose_double2float_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_get_halo_transpose_double2float_slave, &param);
	spawn_proxy_join();

  return ;
}
void athread_put_halo_transpose_float2double_host (float* arrSrc, double* arrObj,
    int startLayer, int lenLayer, int lenA, int lenB, int lenC) {

  // Priority C
  struct HaloTransposeFloat2Double3D param;
  param.startB[0] = startLayer;
  param.startB[1] = startLayer + lenLayer;
  param.startB[2] = startLayer + lenLayer;
  param.startB[3] = lenB - startLayer - lenLayer;
  param.endB[0]   = startLayer + lenLayer;
  param.endB[1]   = lenB - startLayer - lenLayer;
  param.endB[2]   = lenB - startLayer - lenLayer;
  param.endB[3]   = lenB - startLayer;
  param.startC[0] = startLayer;
  param.startC[1] = startLayer;
  param.startC[2] = lenC - startLayer - lenLayer;
  param.startC[3] = startLayer;
  param.endC[0]   = lenC - startLayer;
  param.endC[1]   = startLayer + lenLayer;
  param.endC[2]   = lenC - startLayer;
  param.endC[3]   = lenC - startLayer;

  param.arrSrc = arrSrc;
  param.arrObj = arrObj;
  param.lenA   = lenA;
  param.lenB   = lenB;
  param.lenC   = lenC;

	// athread_spawn(athread_put_halo_transpose_float2double_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_put_halo_transpose_float2double_slave, &param);
	spawn_proxy_join();

  return ;
}


#endif // __sw_host__
