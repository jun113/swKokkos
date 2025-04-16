#ifdef __sw_host__

#include "../../kokkos/swKokkos-4.1.00-opt/core/src/Athread/Kokkos_Athread_FastThreadSpawn.h"
#include "athread_param.h"

#include "athread.h"

extern void SLAVE_FUN(advection_tracer_1_0)(advection_Para*);
extern void SLAVE_FUN(advection_tracer_2)(advection2_Para*);
extern void SLAVE_FUN(advection_tracer_3)(advection2_Para*);


void athread_advection_1_0_host (int imt,int jmt,int km,double *wkb,double *wkd,
	double *u_wface,double *v_sface,double *hts,double *htw,double *dxu,double *dyu,int tid) {

    advection_Para param;

    param.cuBlock[0] = imt;
    param.cuBlock[1] = 1;
    param.cuBlock[2] = 1;
    param.cuGrid[0]  = (imt + param.cuBlock[0] - 1) / param.cuBlock[0];
    param.cuGrid[1]  = (jmt + param.cuBlock[1] - 1) / param.cuBlock[1];
    param.cuGrid[2]  = (km + param.cuBlock[2] - 1) / param.cuBlock[2];

    param.imt        = imt;
    param.jmt        = jmt;
    param.km         = km;
    param.rank       = tid;
    param.d_wkb      = wkb;
    param.d_wkd      = wkd;
    param.d_u_wface  = u_wface;
    param.d_v_sface  = v_sface;
    param.d_hts      = hts;
    param.d_htw      = htw;
    param.d_dxu      = dxu;
    param.d_dyu      = dyu;

	// athread_spawn(advection_tracer_1_0, &param);
	// athread_join();
	spawn_proxy_run(advection_tracer_1_0, &param);
	spawn_proxy_join();

  return ;
}

void athread_advection_2_3_host (int imt,int jmt,int km,double *wkb,
	double *wkd,double *u_wface,double *v_sface,double *hts,double *htw,double *dxu,
	double *dyu,double *ws,double *at,double *adv_tt,double *tarea_r,double *odzp,
	double *hun,double *hue,double *at00,double *atmax,double *atmin,double *odzt,double *vit,double dts,int tid) {

    advection2_Para param;

    param.cuBlock[0] = imt;
    param.cuBlock[1] = 1;
    param.cuBlock[2] = 1;
    param.cuGrid[0]  = (imt + param.cuBlock[0] - 1) / param.cuBlock[0];
    param.cuGrid[1]  = (jmt + param.cuBlock[1] - 1) / param.cuBlock[1];
    param.cuGrid[2]  = (km + param.cuBlock[2] - 1) / param.cuBlock[2];
    param.imt        = imt;
    param.jmt        = jmt;
    param.km         = km;
    param.rank       = tid;
    param.d_wkb      = wkb;
    param.d_wkd      = wkd;
    param.d_u_wface  = u_wface;
    param.d_v_sface  = v_sface;
    param.d_hts      = hts;
    param.d_htw      = htw;
    param.d_dxu      = dxu;
    param.d_dyu      = dyu;
    param.d_ws       = ws;
    param.d_at       = at;
    param.d_adv_tt   = adv_tt;
    param.d_tarea_r  = tarea_r;
    param.d_odzp     = odzp;
    param.d_hun      = hun;
    param.d_hue      = hue;
    param.d_at00     = at00;
    param.d_atmax    = atmax;
    param.d_atmin    = atmin;
    param.d_odzt     = odzt;
    param.d_vit      = vit;
    param.d_dts      = dts;


	// athread_spawn(advection_tracer_2, &param);
	// athread_join();
	spawn_proxy_run(advection_tracer_2, &param);
	spawn_proxy_join();

  // athread_spawn(advection_tracer_3, &param);
	// athread_join();
  spawn_proxy_run(advection_tracer_3, &param);
	spawn_proxy_join();

  return ;
}


#endif // __sw_host__
