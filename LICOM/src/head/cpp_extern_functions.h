#ifndef LICOM3_KOKKOS_SRC_HEAD_CPP_EXTERN_FUNCTION_H_
#define LICOM3_KOKKOS_SRC_HEAD_CPP_EXTERN_FUNCTION_H_

#include "def-undef.h"

#include "cpp_param_mod.h"
#include "cpp_pop_halo_mod.hpp"
#include "cpp_isopyc_mod.h"

#include "licom_test_time.hpp"

#if (defined LICOM_ENABLE_KOKKOS)
#include "kokkos_config.hpp"
#endif

using CppParamMod::MAX_BLOCKS_CLINIC;
using CppParamMod::KM;
using CppParamMod::KMP1;
using CppParamMod::JMT;
using CppParamMod::IMT;
using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;

#ifdef ISO
using CppIsopycMod::NRPL;
#endif // ISO

#ifdef LICOM_ENABLE_TEST_TIME
namespace MyTest {
extern TestTime::MyTime my_time;
} // namespace MyTest
#endif // LICOM_ENABLE_TEST_TIME

extern void cpp_readyt();
extern void cpp_readyc();
extern void cpp_barotr();
extern void cpp_bclinc();
extern void cpp_tracer();
extern void cpp_icesnow();
extern void cpp_convadj();

extern void kokkos_init_view();

extern void stepon(const int &);

extern void kokkos_readyt();
extern void kokkos_readyc();
extern void kokkos_barotr();
extern void kokkos_bclinc();
extern void kokkos_tracer();
extern void kokkos_icesnow();
extern void kokkos_convadj();

extern void read_jra (const int &start_day, const int &num_day, 
		const char* fname, float* const val);
extern void read_jra (const char* fname, double* const buffer, 
    double* const s_lon, double* const s_lat);
extern void read_jra (const int &start_day, const int &num_day, const char* fname,
		double* const buffer, double* const s_lon, double* const s_lat, float* const val);

extern void kokkos_init_jra_daily_low();
extern void kokkos_init_jra_daily_high();
extern void kokkos_jra_daily_low  (const int &iday);
extern void kokkos_jra_daily_high (const int &iday);

extern void daily_update_h2d();
extern void daily_update_d2h();

extern void kokkos_nextstep();

extern void pop_haloupdate_barotr_1(const int &, const int &);
extern void pop_haloupdate_barotr_2(const int &, const int &);

extern void pop_haloupdate_bclinc_2(const int &, const int &);
extern void pop_haloupdate_bclinc_3(const int &, const int &);
extern void pop_haloupdate_bclinc_33(const int &, const int &);

extern void pop_haloupdate_tracer_1(const int &, const int &);
extern void pop_haloupdate_tracer_2(const int &);

#if (defined LICOM_ENABLE_KOKKOS)

extern void gpu_get_halo_transpose_double (const ViewDouble4D &viewSrc, double* const arrObj,
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC);
extern void gpu_put_halo_transpose_double (double* const arrSrc, const ViewDouble4D &viewDst, 
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC);
extern void gpu_get_halo_transpose_bclinc (const ViewDouble4D &viewSrc, double* const arrObj,
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC);
extern void gpu_put_halo_transpose_bclinc (double* const arrSrc, const ViewDouble4D &viewDst, 
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC);
extern void gpu_get_halo_transpose_tracer (const ViewDouble4D &viewSrc, double* const arrObj,
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC);
extern void gpu_put_halo_transpose_tracer (double* const arrSrc, const ViewDouble4D &viewDst, 
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC);
#endif // (defined LICOM_ENABLE_KOKKOS)

#ifdef LOWRES 
extern void cpp_init_jra_daily_low();
extern void cpp_jra_daily_low(const int &iday);
#endif // LOWRES
extern void cpp_init_jra_daily_high();
extern void cpp_jra_daily_high(const int &iday);

extern void upwell (
    const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
    const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
    const double (&)[MAX_BLOCKS_CLINIC][JMT][IMT]);

extern double dens (const double &, const double &, const int &);

extern void vinteg(const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
    double (&)[MAX_BLOCKS_CLINIC][JMT][IMT]);

extern void div(const int &, 
    double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK]);

extern void grad(const int &,
          double (&)[NY_BLOCK][NX_BLOCK], 
          double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK]);

extern void zcurl(const int &, 
    double (&)[NY_BLOCK][NX_BLOCK], 
    const double (&)[NY_BLOCK][NX_BLOCK], 
    const double (&)[NY_BLOCK][NX_BLOCK]);

void invtrit(double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
    const double (&)[MAX_BLOCKS_CLINIC][JMT][IMT], 
				const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT], 
        		const double &, const double &);

extern void invtriu (double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
    const double (&)[MAX_BLOCKS_CLINIC][JMT][IMT],
    const double (&)[MAX_BLOCKS_CLINIC][JMT][IMT],
    const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
    const double &, const double &);

extern void tgrid_to_ugrid (double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK], const int &);

extern void ugrid_to_tgrid(double (&)[NY_BLOCK][NX_BLOCK],
    const double (&)[NY_BLOCK][NX_BLOCK],
        const int &, const int &);

extern void density ();

extern void thermal (
    const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
    const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
	const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
		  double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
		  double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT],
	const double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]);


extern void turb_2 (
    const double (&)[KM],   // wp8
    const double (&)[KM],   // wp1
    const double (&)[KM],   // wp2
    const double (&)[KM],   // wp3
          double (&)[KM],   // wp4
    const double (&)[KM],   // wp5
          double (&)[KM],   // wp6
    const double  &,        // DFRICMX * 1.0e+4
    const double  &,        // DWNDMIX * 1.0e+4
          double (&)[KM-1], // akm_back
          double (&)[KM-1], // akt_back
          double (&)[KM-1], // aks_back
    const double (&)[KM],   // wp7
    const double  &,        // wp9
    const double  &,        // wp10
    const double  &,        // wp11
    const double  &,        // fcort[iblock][j][i]
          double  &,        // mldtmp
          double (&)[KM-1], // wk1
          double (&)[KM-1], // wk2
          double (&)[KM-1], // wk3
    const int     &,        // iwk
    const int     &,        // iwk1
    const int     &,        // km - 1
    const int     &,        // 1
    const int     &,        // 0
    const int     &,        // 0
    const int     &,        // i
    const int     &         /* j */);

#ifdef ISO

void isoflux(const int &, 
    const double (&)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
    const double (&)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
    const double (&)[MAX_BLOCKS_CLINIC][3][JMT][KM+1][IMT],
    const double (&)[MAX_BLOCKS_CLINIC][JMT][KM][IMT],
        double (&)[MAX_BLOCKS_CLINIC][KM][JMT][IMT]);

extern void isopyc(double (&)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
  					double (&)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
  					double (&)[MAX_BLOCKS_CLINIC][3][JMT][KM+1][IMT],
						double (&)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]);

extern void isoadv (const double (&)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
		const double (&)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT],
  			double (&)[MAX_BLOCKS_CLINIC][JMT][KM][IMT]);

extern void k2_3 (const double (&)[MAX_BLOCKS_CLINIC][NRPL][JMT][KM+1][IMT],
		double (&)[MAX_BLOCKS_CLINIC][3][JMT][KMP1][IMT],
				double (&)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT]);
extern void k1_3 (const double (&)[MAX_BLOCKS_CLINIC][NRPL][JMT][KM+1][IMT],
		double (&)[MAX_BLOCKS_CLINIC][3][JMT][KMP1][IMT],
				double (&)[MAX_BLOCKS_CLINIC][1][JMT][KM+1][IMT]);
extern void k3_123 (const double (&)[MAX_BLOCKS_CLINIC][NRPL][JMT][KM+1][IMT],
		double (&)[MAX_BLOCKS_CLINIC][3][JMT][KMP1][IMT],
				double (&)[MAX_BLOCKS_CLINIC][3][JMT][KM+1][IMT]);

#endif // ISO
extern int str_trim_cmp(char const *cmp_str, char const *tar_str);

#ifdef __sw_host__
extern "C" void athread_get_halo_transpose_double_host (double* arrSrc, double* arrObj,
    int startLayer, int lenLayer, int lenA, int lenB, int lenC);
extern "C" void athread_put_halo_transpose_double_host (double* arrSrc, double* arrObj,
    int startLayer, int lenLayer, int lenA, int lenB, int lenC);

extern "C" void athread_readyt_4_1 (int KM, int JMT, int IMT, 
    double *dzp, double *h0, double *h0l, double *h0f, double *psa, 
        double *gg, double *pp, double *ppa, double *vit);
extern "C" void athread_readyt_4_2 (int KM, int JMT, int IMT, 
    double *dzp, double *h0, double *h0l, double *h0f, double *psa, 
        double *gg, double *pp, double *ppa, double *vit);
extern "C" void athread_upwell(int KM, int JMT, int IMT,
    double* uwk, double* vwk, double* h0wk, double* ws, double* ohbt,
    double* vit, double* dzp, double* work, double* wka);

extern "C" void athread_advection_1_0_host (int imt,int jmt,int km,double *wkb,double *wkd,
	double *u_wface,double *v_sface,double *hts,double *htw,double *dxu,double *dyu,int tid);

extern "C" void athread_advection_2_3_host (int imt,int jmt,int km,double *wkb,
	double *wkd,double *u_wface,double *v_sface,double *hts,double *htw,double *dxu,
	double *dyu,double *ws,double *at,double *adv_tt,double *tarea_r,double *odzp,
	double *hun,double *hue,double *at00,double *atmax,double *atmin,double *odzt,double *vit,double dts,int tid);

extern "C" void athread_bclinc_4_host (int KM, int JMT, int IMT, double G,
    double P5, double onbb, double od0, double aa,
    double* h0bf, double* h0bl, double* work, double* psa,
    double* vit,  double* gg,   double* dzp,  double* wka);

extern "C" void athread_bclinc_14_host (int KM, int JMT, int IMT,
    double dtc2, double aidif, int* kmu,
    double* odzt, double* odzp, double* sbcy, double* bbcy,
    double* vp, double* dlv, double* wka, double* viv, double* akmu);

extern "C" void athread_invtriu_host (int KM, int JMT, int IMT,
    double c2dtc, double aidif, int* kmu,
    double* wk, double* topbc, double* bombc, double* dcb,
    double* odzt, double* odzp, double* viv);
extern "C" void athread_vinteg_host (int KM, int JMT, int IMT,
    double* wk3, double* wk2, double* dzp, double* viv, double* ohbu);
extern "C" void athread_tracer_15_host (int KM, int JMT, int IMT,
    double* tf, double* adv_tt, double* vit);
#endif // __sw_host__


#endif // LICOM3_KOKKOS_SRC_HEAD_CPP_EXTERN_FUNCTION_H_