#include "functor.hpp"
#include "Kokkos_Core.hpp"
#include "norm.hpp"
#include "gcr.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <mpi.h>
//#include <mkl.h>
//#include <mkl_cblas.h>
//#include "gptl.h"
//#include "pop_halo_c.hpp"
//#include "matrixpro.hpp"

extern "C" void matrixpro_c(float *a,float *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend);
extern "C" void svrasr( float*  yl,int m, float*  b, int ni, int nk, int nj);
extern "C" void glob_updatehalo_p_(float *h_p, int *m, int *iter_max,int *its, int *ite, int *kts, int *kte, int *jts, int *jte);
extern "C" void matrixpro_x(float *a,double *b,float *c,int its, int ite, int jts, int jte, int kts, int kte, int jend, int ims, int ime,  int jms, int jme);


inline double norm2(const float *v, const  int n)
{
    double tmp_sum = 0;
    #pragma omp parallel for  reduction(+:tmp_sum) schedule(runtime)
    for(int i=0; i<n; i++) 
        tmp_sum += v[i] * v[i];
    return tmp_sum;//sqrt(tmp);
}

extern "C" void  kokkos_psolve_main (float ep, float *a0, float *f0, float *cm_helm,int iter_max, double *x0, int* ijk_index, int *mytid){

  int idep, jdep, ids,ide,jds,jde,kds,kde,ims,ime,
      jms,jme,kms,kme, its,ite,jts,jte,kts,kte;
  assign_indices(ijk_index, idep, jdep, ids, ide, jds, jde, kds, kde,
      ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte);

  int i,j,k,m,l;
  double d;
  // float d_r4;
  float c11_tmp;
  float c12_tmp;
 
  int max_iteration = 50;
  int ibegin = its;
  int iend   = ite;
  int jbegin = jts;
  int jend   = std::min (jde - 1, jte);
 
  int NX = ite - its + 3;
  int NY = kte - kts + 3;
  int NZ = jte - jts + 3;
  int NG1 = NX * NY * NZ;
  int NG = (ite - its + 1 ) * NY * (jte - jts + 1);
  int NG2 = (ime - ims + 1 ) * NY * (jme - jms + 1);

  float *h_p  = new float[NG1*iter_max];
  float *h_ap = new float[NG*(iter_max+1)];
  float *h_ar = new float[NG];
  float *h_r  = new float[NG];

  double aps[iter_max+1], b[iter_max+1];
  double c1[iter_max+2], c2[iter_max+2];

  ViewDouble1D v_aps (aps, iter_max + 1);
  ViewDouble1D v_b   (b,   iter_max + 1);
  ViewDouble1D v_c1  (c1,  iter_max + 2);
  ViewDouble1D v_c2  (c2,  iter_max + 2);

  memset (h_p,  0, sizeof(float) * NG1 * iter_max);
  memset (h_ap, 0, sizeof(float) * NG * (iter_max+1));
  memset (h_ar, 0, sizeof(float) * NG);
  memset (h_r,  0, sizeof(float) * NG);

  ViewFloat1D v_p  (h_p,  NG1 * (iter_max+1));
  ViewFloat1D v_ap (h_ap, NG * (iter_max+1));
  ViewFloat1D v_ar (h_ar, NG);
  ViewFloat1D v_r  (h_r,  NG);
  size_t size_pi = (ime-ims+1) * (kte-kts+3) * (jme-jms+1);
  ViewDouble1D v_x0 (x0, size_pi);

  size_t size_a_helm = 19 * (ite-its+1) * (kte-kts+3) * (jte-jts+1);
  size_t size_b_helm = (ite-its+1) * (kte-kts+3) * (jte-jts+1);
  size_t size_cm_helm = 7 * (ite-its+3) * (kte-kts+3) * (jte-jts+3);
  ViewFloat1D v_a0 (a0, size_a_helm);
  ViewFloat1D v_f0 (f0, size_b_helm);
  ViewFloat1D v_cm_helm (cm_helm, size_cm_helm);

  using KoArr2D = Kokkos::Array<int64_t, 2>;
  using KoArr3D = Kokkos::Array<int64_t, 3>;

  // auto tile2D = KoArr2D{1, (ite+2 - (its-1)) >> 1};
  // auto tile3D = KoArr3D{2, (kte+2 - (kts-1)) >> 1, ite+2 - (its-1)};
  auto tile2D = KoArr2D{1, (ite+2 - (its-1))};
  auto tile3D = KoArr3D{1, (kte+2 - (kts-1)), ite+2 - (its-1)};

  // para_jte1(h_p[index4b(i,k,j,0)]=x0[index_x(i,k,j)];)
  Kokkos::parallel_for ("h_p=x0", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
      KoArr3D{jts-1, kts-1, its-1}, KoArr3D{jte+2, kte+2, ite+2}, tile3D), 
          Functor1 (its, kts, jts, kte, ims, ime, jms, NX, NY, NZ, v_p, v_x0));

  // matrixpro_c(a0, h_p,h_r,its, ite, jts, jte, kts, kte, jend);
  Kokkos::parallel_for ("matrixpro,kts-1", Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
      KoArr2D{jts, its}, KoArr2D{jte+1, ite+1}, tile2D), 
          MatrixProKts (its, kts, jts, ite, NX, NY, 0, 0, 0, v_a0, v_p, v_r));
  Kokkos::parallel_for ("matrixpro,[kts,kte]", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
      KoArr3D{jts, kts, its}, KoArr3D{jte+1, kte+1, ite+1}, tile3D), 
          MatrixProKtsKte (its, kts, jts, ite, NX, NY, 0, 0, 0, v_a0, v_p, v_r));

  Kokkos::parallel_for ("matrixpro,kte+1", Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
      KoArr2D{jts, its}, KoArr2D{jte+1, ite+1}, tile2D), 
          MatrixProKte (its, kts, jts, ite, kte, NX, NY, 0, 0, 0, v_a0, v_p, v_r));
  d = 0.0;
  // gcr_init
  m = 0;
  // para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]-h_r[index3(i,k,j)];)
  Kokkos::parallel_for ("r=f0-r", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
      KoArr3D{jts, kts-1, its}, KoArr3D{jte+1, kte+2, ite+1}, tile3D), 
          Functor2 (its, kts, jts, ite, NX, NY, v_r, v_f0));
  printr();
  // para_jte(h_p[index4b(i,k,j,m)] = h_r[index3(i,k,j)];)
  Kokkos::parallel_for ("p=r", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
      KoArr3D{jts, kts-1, its}, KoArr3D{jte+1, kte+2, ite+1}, tile3D), 
          Functor3 (m, its, kts, jts, ite, NX, NY, NZ, v_p, v_r));

  glob_updatehalo_p_(h_p, &m, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    
  // svrasr(h_p,m, cm_helm, NX, NY, NZ);
  Kokkos::parallel_for ("svrasr", Kokkos::RangePolicy(1, NZ+1), 
      Svrasr (m, NX, NY, NZ, v_p, v_cm_helm));
  glob_updatehalo_p_(h_p, &m, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
  d = 0.0;
  Kokkos::parallel_reduce ("d+=h_r^2", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
      KoArr3D{jbegin, kts-1, ibegin}, KoArr3D{jend+1, kte+2, iend+1}, tile3D), 
          FunctorReduceD(its, kts, jts, ite, NY, v_r), d);
  //#pragma omp parallel for private(j,k,i) reduction(+:d) schedule(runtime) 
  // for(j=jbegin;j<=jend;j++){
  //   for(k=kts-1;k<=kte+1;k++){
  //     for(i=ibegin;i<=iend;i++){
  //       d = d + h_r[index3(i,k,j)]*h_r[index3(i,k,j)];
  //     }
  //   }
  // }
  c1[0] = d;
  MPI_Allreduce (c1, c2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  d = c2[0];
  int s;
  for (s = 0; s < max_iteration; s++) {
    //GPTLstart ("gcr_iteration");
    m = s % (iter_max - 1);
    if ((m == 0) && (s != 0)) {
      // para_jte1(h_p[index4b(i,k,j,0)] = h_p[index4b(i,k,j,iter_max-1)];)
      Kokkos::parallel_for ("h_p(0)=h_p(iter_max-1)", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
          KoArr3D{jts-1, kts-1, its-1}, KoArr3D{jte+2, kte+2, ite+2}, tile3D), 
              Functor4 (its, kts, jts, NX, NY, NZ, iter_max,v_p));
    }
    //GPTLstart ("maxtrixpro");  
    // matrixpro_c(a0, h_p+m*NG1,h_ap+m*NG, its, ite, jts, jte, kts, kte, jend);
    Kokkos::parallel_for ("matrixpro,kts-1", Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
        KoArr2D{jts, its}, KoArr2D{jte+1, ite+1}, tile2D), 
            MatrixProKts (its, kts, jts, ite, NX, NY, 0, m*NG1, m*NG, v_a0, v_p, v_ap));
    Kokkos::parallel_for ("matrixpro,[kts,kte]", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        KoArr3D{jts, kts, its}, KoArr3D{jte+1, kte+1, ite+1}, tile3D), 
            MatrixProKtsKte (its, kts, jts, ite, NX, NY, 0, m*NG1, m*NG, v_a0, v_p, v_ap));
    Kokkos::parallel_for ("matrixpro,kte+1", Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
        KoArr2D{jts, its}, KoArr2D{jte+1, ite+1}, tile2D), 
            MatrixProKte (its, kts, jts, ite, kte, NX, NY, 0, m*NG1, m*NG, v_a0, v_p, v_ap));
    //GPTLstop ("maxtrixpro");  
    double c11 (0.0), c12 (0.0);
        
    //#pragma omp parallel for private(j,k,i) reduction(+:c11,c12) schedule(runtime) 
    // para_jend(
    //   c11 = c11 +  h_r[index3(i,k,j)]* h_ap[index4(i,k,j,m)];
    //   c12 = c12 + h_ap[index4(i,k,j,m)]* h_ap[index4(i,k,j,m)];
    // )
    Kokkos::parallel_reduce ("c11+=r*ap", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        KoArr3D{jts, kts-1, its}, KoArr3D{jend+1, kte+2, ite+1}, tile3D), 
            FunctorReduceC11 (m, its, kts, jts, ite, jte, NY, v_r, v_ap), c11);
    Kokkos::parallel_reduce ("c12+=ap^2", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        KoArr3D{jts, kts-1, its}, KoArr3D{jend+1, kte+2, ite+1}, tile3D), 
            FunctorReduceC12 (m, its, kts, jts, ite, jte, NY, v_ap), c12);
    //GPTLstop ("reduction");  
            
    c1[0] = c11;
    c1[1] = c12;
    MPI_Allreduce (c1, c2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double ac = c2[0] / c2[1];
    aps[m] = c2[1];
  
        //if(mytid[0]==0) printf("ac=%.3e\n",ac);
        //#pragma omp parallel for private(j,k,i) schedule(runtime) 
    // para_jte1(
    //     x0[index_x(i,k,j)]= x0[index_x(i,k,j)]+ac * h_p[index4b(i,k,j,m)];
    // )
    Kokkos::parallel_for ("h_p(0)=h_p(iter_max-1)", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        KoArr3D{jts-1, kts-1, its-1}, KoArr3D{jte+2, kte+2, ite+2}, tile3D), 
            Functor5 (its, kts, jts, kte, ims, ime, jms, NX, NY, NZ, m, ac, v_p, v_x0));
    // para_jend(
    //     h_r[index3(i,k,j)] = h_r[index3(i,k,j)]-ac * h_ap[index4(i,k,j,m)];
    // )
    Kokkos::parallel_for ("r=r-ac*ap", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        KoArr3D{jts, kts-1, its}, KoArr3D{jend+1, kte+2, ite+1}, tile3D), 
            Functor6 (its, kts, jts, ite, jte, NY, m, ac, v_r, v_ap));
    d = 0.0;
    //#pragma omp parallel for private(j,k,i) reduction(+:d) schedule(runtime) 
    // para_jend(
    //     d = d  + h_r[index3(i,k,j)]*h_r[index3(i,k,j)];
    // )
    Kokkos::parallel_reduce ("d+=r^2", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        KoArr3D{jts, kts-1, its}, KoArr3D{jend+1, kte+2, ite+1}, tile3D), 
            FunctorReduceD(its, kts, jts, ite, NY, v_r), d);

    //#pragma omp parallel for private(j,k,i) schedule(runtime) 
    // para_jte(h_p[index4b(i,k,j,m+1)] = h_r[index3(i,k,j)];)
    Kokkos::parallel_for ("p=r", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        KoArr3D{jts, kts-1, its}, KoArr3D{jte+1, kte+2, ite+1}, tile3D), 
            Functor3 (m+1, its, kts, jts, ite, NX, NY, NZ, v_p, v_r));

    int m_tmp = m + 1;
    //GPTLstart ("glob_updatehalo");  
    glob_updatehalo_p_(h_p, &m_tmp, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    //GPTLstop ("glob_updatehalo");
        
    //GPTLstart ("precondition");  
    // svrasr(h_p,m_tmp, cm_helm, NX, NY, NZ);
    Kokkos::parallel_for ("svrasr", Kokkos::RangePolicy(1, NZ+1), 
        Svrasr (m_tmp, NX, NY, NZ, v_p, v_cm_helm));
    //GPTLstop ("precondition");  
    //GPTLstart ("glob_updatehalo");  
    glob_updatehalo_p_ (h_p, &m_tmp, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    //GPTLstop ("glob_updatehalo");  
    //GPTLstart ("maxtrixpro");  
    // matrixpro_c(a0, h_p+m_tmp*NG1,h_ar, its, ite, jts, jte, kts, kte, jend);
    Kokkos::parallel_for ("matrixpro,kts-1", Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
        KoArr2D{jts, its}, KoArr2D{jte+1, ite+1}, tile2D), 
            MatrixProKts (its, kts, jts, ite, NX, NY, 0, m_tmp*NG1, 0, v_a0, v_p, v_ar));
    Kokkos::parallel_for ("matrixpro,[kts,kte]", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
        KoArr3D{jts, kts, its}, KoArr3D{jte+1, kte+1, ite+1}, tile3D), 
            MatrixProKtsKte (its, kts, jts, ite, NX, NY, 0, m_tmp*NG1, 0, v_a0, v_p, v_ar));
    Kokkos::parallel_for ("matrixpro,kte+1", Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
        KoArr2D{jts, its}, KoArr2D{jte+1, ite+1}, tile2D), 
            MatrixProKte (its, kts, jts, ite, kte, NX, NY, 0, m_tmp*NG1, 0, v_a0, v_p, v_ar));
    //GPTLstop ("maxtrixpro");  
    //GPTLstart ("reduction_cl");  
    for (l = 0; l <= m; l++) {
      double cl = 0.0;
      //#pragma omp parallel for private(j,k,i) reduction(+:cl) schedule(runtime) 
      // para_jend(cl = cl + h_ar[index3(i,k,j)]*h_ap[index4(i,k,j,l)];)

      Kokkos::parallel_reduce ("cl+=ar*ap", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
          KoArr3D{jts, kts-1, its}, KoArr3D{jend+1, kte+2, ite+1}, tile3D), 
              FunctorReduceCl (l, its, kts, jts, ite, jte, NY, v_ar, v_ap), cl);
      c1[l] = cl;
    }
    //GPTLstop ("reduction_cl");  

    c1[m+1] = d;
    MPI_Allreduce (c1, c2, m + 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    d = c2[m+1];

    // for (l = 0; l <= m; l++){
    //   b[l] = - c2[l] / aps[l];
    // }
    Kokkos::parallel_for ("b=-c2/aps", m+1, Functor7 (v_b, v_c2, v_aps));

    //#pragma omp parallel for private(j,k,i,l) schedule(runtime)   
    // for(j=jts-1;j<=jte+1;j++){
    //   for(l=0;l<=m;l++){
    //     for(k=kts-1;k<=kte+1;k++){
    //       for(i=its-1;i<=ite+1;i++){
    //         h_p[index4b(i,k,j,m+1)]+=b[l]*h_p[index4b(i,k,j,l)];//)
    //       }
    //     }
    //   }
    // } 
    for(l = 0; l <= m; ++l) {
      Kokkos::parallel_for ("p+=b*p", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
          KoArr3D{jts-1, kts-1, its-1}, KoArr3D{jte+2, kte+2, ite+2}, tile3D), 
              Functor8 (m, l, its, kts, jts, NX, NY, NZ, v_p, v_b));
    }

    //GPTLstop ("reductionp");  
    //GPTLstop ("gcr_iteration");  
    if (d <= ep || s == max_iteration-1 ) { 
      if(mytid[0]==0)printf("RES of gcr is %e in %d iterations\n", sqrt(d) ,s);
      break;
    }
    if (mytid[0]==0) printf("RES of gcr is %e in %d iterations\n", sqrt(d) ,s);
    if (s == 10) return;
  }
  // End iter
  // matrixpro_x(a0, x0,h_ap,its, ite, jts, jte, kts, kte, jend, ims, ime, jms, jme);
  Kokkos::parallel_for ("matrixpro,kts-1", Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
      KoArr2D{jts, its}, KoArr2D{jte+1, ite+1}, tile2D), 
          MatrixProKtsD (its, kts, jts, ite, NX, NY, 0, 0, 0, v_a0, v_x0, v_ap));
  Kokkos::parallel_for ("matrixpro,[kts,kte]", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
      KoArr3D{jts, kts, its}, KoArr3D{jte+1, kte+1, ite+1}, tile3D), 
          MatrixProKtsKteD (its, kts, jts, ite, NX, NY, 0, 0, 0, v_a0, v_x0, v_ap));
  Kokkos::parallel_for ("matrixpro,kte+1", Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
      KoArr2D{jts, its}, KoArr2D{jte+1, ite+1}, tile2D), 
          MatrixProKteD (its, kts, jts, ite, kte, NX, NY, 0, 0, 0, v_a0, v_x0, v_ap));
  //initialize(h_r,NG);
  memset (v_r.data(), 0, sizeof(float) * NG);
  // para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]- h_ap[index3(i,k,j)];)
  Kokkos::parallel_for ("r=f0-ap", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
      KoArr3D{jts, kts-1, its}, KoArr3D{jend+1, kte+2, ite+1}, tile3D), 
          Functor9 (its, kts, jts, ite, NX, NY, v_r, v_f0, v_ap));

  // c1[0] = norm2(h_r, NG);
  Kokkos::parallel_reduce ("norm2", NG, FunctorNorm2(v_r), c1[0]);

  MPI_Allreduce(c1, c2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  c2[0] = sqrt (c2[0]);

  if(mytid[0] == 5 || mytid[0] == 12) {
    printf("RES of gmres is %e \n",c2[0]);
  }

  delete [] h_p;
  delete [] h_ap; 
  delete [] h_ar; 
  delete [] h_r;
  return ;
}
