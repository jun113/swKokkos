extern "C" void glob_updatehalo_x0_(double *x0, int *ims, int *ime, int *jms, int *jme, int *kts, int *kte);
extern "C" void glob_updatehalo_p_(float* p, int* m, int* iter_max, int *its, int *ite, int *kts, int *kte, int *jts, int *jte);


void  psolve_main(float ep, float *a0, float *f0, float *cm_helm,int iter_max, double *x0, 
    int dep, int  jdep, int  ids, int  ide, int  jds, int  jde, int  kds, int  kde, 
    int ims, int  ime, int  jms, int  jme, int  kms, int  kme, 
    int its, int  ite, int  jts, int  jte, int  kts, int  kte, int *mytid){

    #include "inc/defvar.hpp"

    //glob_updatehalo_x0_(x0, &ims, &ime, &jms, &jme, &kts, &kte);

    initialize(h_p,NG1*(iter_max+1));
    initialize(h_ap,NG*(iter_max+1));

    para_jte1(h_p[index4b(i,k,j,0)]=x0[index_x(i,k,j)];)

    matrixpro_c(a0, h_p,h_r,its, ite, jts, jte, kts, kte, jde);
    m = 0;
    d = 0.0;
    // gcr_init
    para_jend(h_r[index3(i,k,j)]  = f0[index3(i,k,j)]-h_r[index3(i,k,j)];)
    para_jte(h_p[index4b(i,k,j,m)] = h_r[index3(i,k,j)];)

    glob_updatehalo_p_(h_p, &m, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);

    svrasr(h_p,m, cm_helm, NX, NY, NZ);
    glob_updatehalo_p_(h_p, &m, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
    d = 0.0;
    for(j=jbegin;j<=jend;j++){
        for(k=kts-1;k<=kte+1;k++){
            for(i=ibegin;i<=iend;i++){
                d = d + h_r[index3(i,k,j)]*h_r[index3(i,k,j)];
            }
        }
    }
    c1[0] = d;
    MPI_Allreduce(c1,c2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    d = c2[0];
    int s;
    for(s=0;s<max_iteration;s++){
        GPTLstart ("iteration");  
        m = s%(iter_max-1);
        if ((m == 0) &&( s!=0) ) {
            para_jte1(h_p[index4b(i,k,j,0)] = h_p[index4b(i,k,j,iter_max-1)];)
        }
        matrixpro_c(a0, h_p+m*NG1,h_ap+m*NG, its, ite, jts, jte, kts, kte, jde);
        double c11=0.0,c12=0.0;
        para_jend(
            c11 = c11 +  h_r[index3(i,k,j)]* h_ap[index4(i,k,j,m)];
            c12 = c12 + h_ap[index4(i,k,j,m)]* h_ap[index4(i,k,j,m)];
        )
            
        c1[0]=c11;
        c1[1]=c12;
        MPI_Allreduce(c1,c2,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        double ac=c2[0]/c2[1];
        aps[m]=c2[1];
  
        d = 0.0;
        para_jend(
            x0[index_x(i,k,j)]= x0[index_x(i,k,j)]+ac * (double)h_p[index4b(i,k,j,m)];
            h_r[index3(i,k,j)] = h_r[index3(i,k,j)]-ac * h_ap[index4(i,k,j,m)];
            d = d  + h_r[index3(i,k,j)]*h_r[index3(i,k,j)];
        )

        para_jte(h_p[index4b(i,k,j,m+1)] = h_r[index3(i,k,j)];)

        int m_tmp=m+1;
        glob_updatehalo_p_(h_p, &m_tmp, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);
        svrasr(h_p,m_tmp, cm_helm, NX, NY, NZ);
        glob_updatehalo_p_(h_p, &m_tmp, &iter_max, &its, &ite, &kts, &kte, &jts, &jte);

        matrixpro_c(a0, h_p+m_tmp*NG1,h_ar, its, ite, jts, jte, kts, kte, jde);

        for(l=0;l<=m;l++){
            double cl=0.0;
            para_jend(cl = cl + h_ar[index3(i,k,j)]*h_ap[index4(i,k,j,l)];)
            c1[l]=cl;
        }  

        c1[m+1] = d;
        MPI_Allreduce(c1,c2,m+2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        d = c2[m+1];
        if (d <= ep || s == max_iteration-1 ) { 
            if(mytid[0]==0)printf("RES of gcr is %e in %d iterations\n", d ,s);
            break;}
        for(l=0;l<=m;l++){
            double bl=-c2[l]/aps[l];
            para_jte1( h_p[index4b(i,k,j,m+1)]=h_p[index4b(i,k,j,m+1)]+bl*h_p[index4b(i,k,j,l)];)
        } 
        if(mytid[0]==0)printf("RES of gcr is %e in %d iterations\n", d ,s);
        GPTLstop ("iteration");  
    }

	delete[] h_p;
    delete[] h_ap; 
    delete[] h_ar; delete[] h_r;
}

