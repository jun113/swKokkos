#ifndef NORM_HPP
#define NORM_HPP
#include <stdio.h>
static float cblas_sdot(int N, const float *X, int incX, const float *Y, int incY) {
    float result = 0.0f;
    int ix = 0;
    int iy = 0;
    for (int i = 0; i < N; i++) {
        result += X[ix] * Y[iy];
        ix += incX;
        iy += incY;
    }
    return result;
}
static void initialize_ijk_index(int* mytid, int* ijk_index) {
     int imt=1000;
     int jmt=501;
     int kmt=88;
    // 固定值
    ijk_index[0] = 2;
    ijk_index[1] = 2;
    ijk_index[2] = 1;
    ijk_index[3] = imt;
    ijk_index[4] = 1;
    ijk_index[5] = jmt;
    ijk_index[6] = 1;
    ijk_index[7] = kmt;
    ijk_index[12] = 0;
    ijk_index[13] = kmt;
    ijk_index[18] = 1;
    ijk_index[19] = kmt-1;
    int rank=mytid[0];
    int imt_sub=int(imt/mytid[8]);
    int jmt_sub=int(jmt/mytid[7]);
    // mytid[5]=int(rank / mytid[8])
    // mytid[6] = (rank % mytid[8]) 
    //printf("mytid[5]=%d,%d\n",mytid[5],int(rank / mytid[8]));
    // 根据rank计算的值
    int base = mytid[6] * imt_sub;
    int offset = mytid[5] * jmt_sub;
    if(mytid[5]==mytid[7]-1) jmt_sub++; 
    //printf("imt_sub=%d,%d\n",imt_sub,jmt_sub); 
    //printf("tid=%d,base=%d,%d\n",mytid[0],base,offset); 
    ijk_index[8] = base - 1;
    ijk_index[9] = base + imt_sub+2;
    ijk_index[10] = offset - 1;
    ijk_index[11] = offset + jmt_sub+2;
    ijk_index[14] = base + 1;
    ijk_index[15] = base + imt_sub;
    ijk_index[16] = offset + 1;
    ijk_index[17] = offset + jmt_sub;
    //ijk_index[17] = 1;
    //if(mytid[0]==0)printf("ijk_index[13]=%d\n",ijk_index[13]);
}

static void assign_indices(int* ind, int &idep, int &jdep, int &ids, int &ide, int &jds, int &jde, int &kds, int &kde,
                    int &ims, int &ime, int &jms, int &jme, int &kms, int &kme,
                    int &its, int &ite, int &jts, int &jte, int &kts, int &kte) {
    idep = ind[0];jdep = ind[1];ids = ind[2];ide = ind[3];
    jds = ind[4];jde = ind[5];kds = ind[6];kde = ind[7];
    ims = ind[8];ime = ind[9];jms = ind[10];jme = ind[11];
    kms = ind[12];kme = ind[13];its = ind[14];ite = ind[15];
    jts = ind[16];jte = ind[17];kts = ind[18];kte = ind[19];
    
//if(mytid[0]==0){
//    printf("kms: %d, kme: %d, its: %d, ite: %d\n", kms, kme, its, ite);
//}
}

#ifdef CUDA

#define cuprintr(s)\
		hipblasSdot(h, NG, d_r, 1, d_r, 1, dd);\
		CHECK(hipMemcpy(&c12_tmp, dd, sizeof(float), hipMemcpyDeviceToHost));\
		if(mytid[0]==12)printf("iter=%d,r_norm=%E\n",s,c12_tmp);

#define cuprintar()\
		hipblasSdot(h, NG, d_ar, 1, d_ar, 1, dd);\
		CHECK(hipMemcpy(&c12_tmp, dd, sizeof(float), hipMemcpyDeviceToHost));\
		if(mytid[0]==12)printf("line=%d,ar_norm=%E\n",__LINE__,c12_tmp);

#define cuprintp(m)\
		hipblasSdot(h, NG1, d_p+(m)*NG1, 1, d_p+(m)*NG1, 1, dd);\
		CHECK(hipMemcpy(&c12_tmp, dd, sizeof(float), hipMemcpyDeviceToHost));\
		if(mytid[0]==12)printf("mytid=%d,m=%d,p_norm=%E\n",mytid[0],m,c12_tmp);
/*
#define cuprintp1(m)\
		hipblasSdot(h, NG1, d_p1+(m)*NG1, 1, d_p1+(m)*NG1, 1, dd);\
		CHECK(hipMemcpy(&c12_tmp, dd, sizeof(float), hipMemcpyDeviceToHost));\
		if(mytid[0]==12||mytid[0]==1)printf("mytid=%d,p1_norm=%E\n",mytid[0],c12_tmp);
#define cuprintap(m)\
		hipblasSdot(h, NG, d_p+(m)*NG, 1, d_p+(m)*NG, 1, dd);\
		CHECK(hipMemcpy(&c12_tmp, dd, sizeof(float), hipMemcpyDeviceToHost));\
		if(mytid[0]==12)printf("line=%d,ap_norm=%E\n",__LINE__,c12_tmp);
*/
#endif
#define print_norm(a,b,n) \
                if(mytid[0]==0){\
                    float a_norm=cblas_sdot(n,a,1,a,1);\
                    printf("line=%d,norm of %s =%E\n",__LINE__,b,a_norm);\
                }


#define printres() \
                if(mytid[0]==0){\
                    para_jte1(h_p1[index_x(i,k,j)]= x0[index_x(i,k,j)]-x1[index_x(i,k,j)];)\
                    c12_tmp=cblas_sdot(NG1,h_p1,1,h_p1,1);\
                    printf("line=%d,res=%E\n",__LINE__,c12_tmp);\
                }

#define printr()\
                c11_tmp=cblas_sdot(NG,h_r,1,h_r,1);\
       MPI_Allreduce(&c11_tmp,&c12_tmp,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD);\
                if(mytid[0]==12) printf("mytid=%d,r_norm=%E,tol=%E\n",mytid[0],sqrt(c12_tmp),sqrt(c12_tmp));
#define printar()\
                if(mytid[0]==12) c12_tmp=cblas_sdot(NG,h_ar,1,h_ar,1);\
                if(mytid[0]==12) printf("line=%d,ar_norm=%E\n",__LINE__,c12_tmp);
#define printx()\
                 if(mytid[0]==12) c12_tmp=(float)cblas_ddot(NG2,x0,1,x0,1);\
                 if(mytid[0]==12) printf("m=%d,x_norm=%E\n",m,c12_tmp);

#define printy()\
                if(mytid[0]==12) c12_tmp=cblas_sdot(NG,y,1,y,1);\
                if(mytid[0]==12) printf("line=%d,m=%d,y_norm=%E\n",__LINE__,m,c12_tmp);
#define printp(m)\
                if(mytid[0]==12) c12_tmp=cblas_sdot(NG1,h_p+(m)*NG1,1,h_p+(m)*NG1,1);\
                if(mytid[0]==12) printf("line=%d,m=%d,p_norm=%E\n",__LINE__,m,c12_tmp);
#define printap(m)\
                if(mytid[0]==12) c12_tmp=cblas_sdot(NG,h_ap+(m)*NG,1,h_ap+(m)*NG,1);\
                if(mytid[0]==12) printf("line=%d,m=%d,ap_norm=%E\n",__LINE__,m,c12_tmp);




#ifdef CUDA

#define cupara_jend(status)\
        if(j>=jts && j<=jend){\
            if(k>=kts-1 && k<=kte+1){\
                if(i>=its && i<=ite){\
                    status\
                }\
            }\
        }

#define cupara_jte(status)\
        if(j>=jts && j<=jte){\
            if(k>=kts-1 && k<=kte+1){\
                if(i>=its && i<=ite){\
                    status\
                }\
            }\
        }
#define cupara_jte1(status)\
        if(j>=jts-1 && j<=jte+1){\
            if(k>=kts-1 && k<=kte+1){\
                if(i>=its-1 && i<=ite+1){\
                    status\
                }\
            }\
        }
#endif
#define para_jend(status)\
        for(j=jts;j<=jend;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its;i<=ite;i++){\
                    status\
                }\
            }\
        }

#define para_jte(status)\
        for(j=jts;j<=jte;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its;i<=ite;i++){\
                    status\
                }\
            }\
        }
#define para_jte1(status)\
        for(j=jts-1;j<=jte+1;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its-1;i<=ite+1;i++){\
                    status\
                }\
            }\
        }

#define para_step(status) for(int l=0;l<2*step+1;l++){ status }
#define para_m(status) for(int l=0;l<iter_max;l++){ status }
#define print_line()\
    if(tid==0)printf("line=%d\n",__LINE__)
#endif
#define boundary(status)\
        for(j=jts-2;j<=jts-1;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its-2;i<=ite+2;i++){\
                    status\
        }}}\
        for(j=jte+1;j<=jte+2;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its-2;i<=ite+2;i++){\
                    status\
        }}}\
        for(i=its-1;i<=its;i++){\
            for(j=jts-2;j<=jte+2;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                    status\
                }\
            }}\
        for(i=ite+1;i<=ite+2;i++){\
            for(j=jts-2;j<=jte+2;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                    status\
                }\
        }}

#define boundary1(status)\
        for(j=jts-2;j<=jts;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its-2;i<=ite+2;i++){\
                    status\
        }}}\
        for(j=jte;j<=jte+2;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                for(i=its-2;i<=ite+2;i++){\
                    status\
        }}}\
        for(i=its-2;i<=its;i++){\
            for(j=jts-2;j<=jte+2;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                    status\
                }\
            }}\
        for(i=ite;i<=ite+2;i++){\
            for(j=jts-2;j<=jte+2;j++){\
            for(k=kts-1;k<=kte+1;k++){\
                    status\
                }\
        }}


#define printdiff()\
         {int cnt=0;\
            fprintf(fp,"diff=\n");\
                para_jte1(if(j==jts&& fabs(h_p[index4b(i,k,j,0)]-h_p[index4b(i,k,j,2)])>0.1*(fabs(h_p[index4b(i,k,j,0)])+fabs(h_p[index4b(i,k,j,2)]))){\
                        cnt++;fprintf(fp,"  %3d,%3d,%2d  ",i-(its-1),k-(kts-1),mytid[0]);\
                        if(cnt %11==0) fprintf(fp,"\n");\
                })}
#define printdiff_r()\
         {int cnt=0;\
            fprintf(fp,"diff_r=\n");\
                para_jte(if(j==jts&& fabs(h_r[index3(i,k,j)]-h_r1[index3(i,k,j)])>0.1*(fabs(h_r1[index3(i,k,j)])+fabs(h_r[index3(i,k,j)]))){\
                        cnt++;fprintf(fp,"  %3d,%3d,%3d  ",i-(its),k-(kts),j-jts);\
                        if(cnt %11==0) fprintf(fp,"\n");\
                })}

