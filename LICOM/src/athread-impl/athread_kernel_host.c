#ifdef __sw_host__

#include "../../kokkos/swKokkos-4.1.00-opt/core/src/Athread/Kokkos_Athread_FastThreadSpawn.h"
#include "athread_param.h"

#include "athread.h"

extern void SLAVE_FUN(athread_readyt_4_1_slave)(struct paramReadyt4*);
void athread_readyt_4_1 (int KM, int JMT, int IMT, 
    double *dzp, double *h0, double *h0l, double *h0f, double *psa, 
        double *gg, double *pp, double *ppa, double *vit) {
  
  // for (int j = 0; j < JMT; ++j) {
  //   for (int i = 0; i < IMT; ++i) {
  //     h0l[j*IMT + i] = h0f[j*IMT + i];
  //     h0f[j*IMT + i] =  h0[j*IMT + i];
  //     pp [j*IMT + i] =  gg[j*IMT + i] * vit[j*IMT + i] * 0.5 * dzp[0];
  //     ppa[j*IMT + i] = psa[j*IMT + i] * vit[j*IMT + i];
  //   }
  // }
  struct paramReadyt4 param;
  param.IMT = IMT;
  param.JMT = JMT;
  param.KM  = KM;
  param.dzp = dzp;
  param.h0  = h0;
  param.h0l = h0l;
  param.h0f = h0f;
  param.psa = psa;
  param.gg  = gg;
  param.pp  = pp;
  param.ppa = ppa;
  param.vit = vit;
	// athread_spawn(athread_readyt_4_1_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_readyt_4_1_slave, &param);
	spawn_proxy_join();
  return ;
}
extern void SLAVE_FUN(athread_readyt_4_2_slave)(struct paramReadyt4*);
void athread_readyt_4_2 (int KM, int JMT, int IMT, 
    double *dzp, double *h0, double *h0l, double *h0f, double *psa, 
        double *gg, double *pp, double *ppa, double *vit) {
  
  // TODO SIMD
  // int strideK = JMT * IMT;
  // for (int k = 1; k < KM; ++k) {
  //   for (int j = 0; j < JMT; ++j) {
  //     for (int i = 0; i < IMT; ++i) {
  //       pp[k*strideK + j*IMT + i] = vit[k*strideK + j*IMT + i] *
  //         (pp[(k-1)*strideK + j*IMT + i] + 0.5 * 
  //         (gg[(k)  *strideK + j*IMT + i] * dzp[k] +
  //          gg[(k-1)*strideK + j*IMT + i] * dzp[k-1]));

  //       ppa[k*strideK + j*IMT + i] = vit[k*strideK + j*IMT + i] *
  //         (ppa[(k-1)*strideK + j*IMT + i] + gg[(k-1)*strideK + j*IMT + i] * dzp[k-1]);
  //     }
  //   }
  // }
  struct paramReadyt4 param;
  param.IMT = IMT;
  param.JMT = JMT;
  param.KM  = KM;
  param.dzp = dzp;
  param.h0  = h0;
  param.h0l = h0l;
  param.h0f = h0f;
  param.psa = psa;
  param.gg  = gg;
  param.pp  = pp;
  param.ppa = ppa;
  param.vit = vit;
	// athread_spawn(athread_readyt_4_2_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_readyt_4_2_slave, &param);
  spawn_proxy_join ();
  return ;
}

extern void SLAVE_FUN(athread_upwell_slave)(struct paramUpwell*);
void athread_upwell(int KM, int JMT, int IMT,
    double* uwk, double* vwk, double* h0wk, double* ws, double* ohbt,
    double* vit, double* dzp, double* work, double* wka) {

  struct paramUpwell param;
  param.KM = KM;
  param.JMT = JMT;
  param.IMT = IMT;
  param.uwk = uwk;
  param.vwk = vwk;
  param.h0wk = h0wk;
  param.ws   = ws;
  param.ohbt = ohbt;
  param.vit  = vit;
  param.dzp  = dzp;
  param.work = work;
  param.wka  = wka;
	// athread_spawn(athread_upwell_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_upwell_slave, &param);
	spawn_proxy_join();

  return ;
}

extern void SLAVE_FUN(athread_bclinc_4_slave)(struct paramBclinc4*);
void athread_bclinc_4_host (int KM, int JMT, int IMT, double G,
    double P5, double onbb, double od0, double aa,
    double* h0bf, double* h0bl, double* work, double* psa,
    double* vit,  double* gg,   double* dzp,  double* wka) {

  struct paramBclinc4 param;
  param.KM = KM;
  param.JMT = JMT;
  param.IMT = IMT;
  param.G = G;
  param.P5 = P5;
  param.onbb = onbb;
  param.od0   = od0;
  param.aa = aa;
  param.h0bf  = h0bf;
  param.h0bl  = h0bl;
  param.work = work;
  param.psa  = psa;
  param.vit  = vit;
  param.gg  = gg;
  param.dzp  = dzp;
  param.wka  = wka;
	// athread_spawn(athread_bclinc_4_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_bclinc_4_slave, &param);
	spawn_proxy_join();

  return;
}
extern void SLAVE_FUN(athread_bclinc_14_slave)(struct paramBclinc14*);
void athread_bclinc_14_host (int KM, int JMT, int IMT,
    double dtc2, double aidif, int* kmu,
    double* odzt, double* odzp, double* sbcy, double* bbcy,
    double* vp, double* dlv, double* wka, double* viv, double* akmu) {
  struct paramBclinc14 param;
  param.KM = KM;
  param.JMT = JMT;
  param.IMT = IMT;
  param.dtc2 = dtc2;
  param.aidif = aidif;
  param.kmu = kmu;
  param.odzt = odzt;
  param.odzp = odzp;
  param.sbcy = sbcy;
  param.bbcy = bbcy;
  param.vp = vp;
  param.dlv = dlv;
  param.wka = wka;
  param.viv = viv;
  param.akmu = akmu;
	// athread_spawn(athread_bclinc_14_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_bclinc_14_slave, &param);
	spawn_proxy_join();
  return ;
}

extern void SLAVE_FUN(athread_invtriu_slave)(struct paramInvtriu*);
void athread_invtriu_host (int KM, int JMT, int IMT,
    double c2dtc, double aidif, int* kmu,
    double* wk, double* topbc, double* bombc, double* dcb,
    double* odzt, double* odzp, double* viv) {
  struct paramInvtriu param;
  param.KM    = KM;
  param.JMT   = JMT;
  param.IMT   = IMT;
  param.c2dtc = c2dtc;
  param.aidif = aidif;
  param.kmu   = kmu;
  param.wk    = wk;
  param.topbc = topbc;
  param.bombc = bombc;
  param.dcb   = dcb;
  param.odzt  = odzt;
  param.odzp  = odzp;
  param.viv   = viv;
	// athread_spawn (athread_invtriu_slave, &param);
	// athread_join ();
	spawn_proxy_run (athread_invtriu_slave, &param);
	spawn_proxy_join ();
  return;
}

extern void SLAVE_FUN(athread_vinteg_slave)(struct paramVinteg*);
void athread_vinteg_host (int KM, int JMT, int IMT,
    double* wk3, double* wk2, double* dzp, double* viv, double* ohbu) {
  struct paramVinteg param;
  param.KM   = KM;
  param.JMT  = JMT;
  param.IMT  = IMT;
  param.wk3  = wk3;
  param.wk2  = wk2;
  param.dzp  = dzp;
  param.viv  = viv;
  param.ohbu = ohbu;
	// athread_spawn (athread_vinteg_slave, &param);
	// athread_join ();
	spawn_proxy_run (athread_vinteg_slave, &param);
	spawn_proxy_join ();

  return ;
}
extern void SLAVE_FUN(athread_tracer_15_slave)(struct paramTracer15*);
void athread_tracer_15_host (int KM, int JMT, int IMT,
    double* tf, double* adv_tt, double* vit) {

  struct paramTracer15 param;
  param.KM = KM;
  param.JMT = JMT;
  param.IMT = IMT;
  param.adv_tt = adv_tt;
  param.tf = tf;
  param.vit = vit;
	// athread_spawn(athread_tracer_15_slave, &param);
	// athread_join();
	spawn_proxy_run(athread_tracer_15_slave, &param);
	spawn_proxy_join();

  return ;
}

#endif // __sw_host__
