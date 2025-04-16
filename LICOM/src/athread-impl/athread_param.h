#ifndef LICOM3_KOKKOS_SRC_SW_IMPL_SW_PARAM_H_
#define LICOM3_KOKKOS_SRC_SW_IMPL_SW_PARAM_H_
#if (defined __sw_host__) || (defined __sw_slave__)

struct HaloTransposeDouble3D {
  int     startB[4];
  int     endB[4];
  int     startC[4];
  int     endC[4];
  double* arrSrc;
  double* arrObj;
  int     lenA;
  int     lenB;
  int     lenC;
};

struct HaloTransposeDouble2Float3D {
  int     startB[4];
  int     endB[4];
  int     startC[4];
  int     endC[4];
  double* arrSrc;
  float*  arrObj;
  int     lenA;
  int     lenB;
  int     lenC;
};
struct HaloTransposeFloat2Double3D {
  int     startB[4];
  int     endB[4];
  int     startC[4];
  int     endC[4];
  float* arrSrc;
  double* arrObj;
  int     lenA;
  int     lenB;
  int     lenC;
};

struct paramReadyt4 {
  int IMT;
  int JMT;
  int KM;
  double* dzp;
  double* h0;
  double* h0l;
  double* h0f;
  double* psa;
  double* gg;
  double* pp;
  double* ppa;
  double* vit;
};

struct paramUpwell {
  int KM;
  int JMT;
  int IMT;
  double* uwk;
  double* vwk;
  double* h0wk;
  double* ws;
  double* ohbt;
  double* vit;
  double* dzp;
  double* work;
  double* wka;
};

struct paramBclinc4 {
  int KM; 
  int JMT; 
  int IMT; 
  double G;
  double P5; 
  double onbb; 
  double od0;
  double aa;
  double* h0bf; 
  double* h0bl; 
  double* work; 
  double* psa;
  double* vit; 
  double* gg;   
  double* dzp;  
  double* wka;
};

typedef struct kernelPara{
    int cuGrid[3];
	int cuBlock[3];
    int imt;
    int jmt;
    int km;
	int rank;
	double *d_wkb;
	double *d_wkd;
	double *d_u_wface;
	double *d_v_sface;
	double *d_hts;
	double *d_htw;
	double *d_dxu;
	double *d_dyu;

}advection_Para;


typedef struct advection2Para{
    int cuGrid[3];
	int cuBlock[3];
    int imt;
    int jmt;
    int km;
	int rank;
	double *d_wkb;
	double *d_wkd;
	double *d_u_wface;
	double *d_v_sface;
	double *d_hts;
	double *d_htw;
	double *d_dxu;
	double *d_dyu;
	double *d_ws;
	double *d_at;
	double *d_adv_tt;
	double *d_tarea_r;
	double *d_odzp;
	double *d_hun;
	double *d_hue;
	double *d_at00;
	double *d_atmax;
	double *d_atmin;
	double *d_odzt;
	double *d_vit;
	double d_dts;

}advection2_Para;

struct paramBclinc14 {
  int KM; 
  int JMT; 
  int IMT;
  double dtc2; 
  double aidif; 
  int* kmu;
  double* odzt; 
  double* odzp; 
  double* sbcy; 
  double* bbcy;
  double* vp; 
  double* dlv; 
  double* wka; 
  double* viv; 
  double* akmu;
};

struct paramInvtriu {
  int KM; 
  int JMT; 
  int IMT;
  double c2dtc; 
  double aidif; 
  int* kmu;
  double* wk; 
  double* topbc; 
  double* bombc; 
  double* dcb;
  double* odzt; 
  double* odzp; 
  double* viv; 
};
struct paramVinteg {
  int KM;
  int JMT; 
  int IMT;
  double* wk3;
  double* wk2; 
  double* dzp; 
  double* viv; 
  double* ohbu;
};

struct paramReadyt10 {
  int KM;
  int JMT;
  int IMT;
  double  G;
  double* dlu;
  double* dlv;
  double* dzp;
  double* zkt;
  double* ohbt;
  double* gg;
};

struct paramTracer15 {
  int KM; 
  int JMT; 
  int IMT;
  double* tf; 
  double* vit; 
  double* adv_tt; 
};

#endif // (defined __sw_host__) || (defined __sw_slave__)
#endif // LICOM3_KOKKOS_SRC_SW_IMPL_SW_PARAM_H_