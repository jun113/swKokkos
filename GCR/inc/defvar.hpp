//#define share_mem
#ifdef share_mem
        #define cu_matrixpro matrixpro_sme
#endif
#define halfwrp
#ifdef halfwrp
        #define cu_svrasr pacu_svrasr
#endif


        int i,j,k,m,l,ierr;
        double d;
        float d_r4;
        float c11_tmp;
        float c12_tmp;
 
        //iter_max=10;
        int max_iteration = 50;
        int ibegin = its;
        int iend   = ite;
        int jbegin = jts;
        int jend   = std::min (jde-1,jte);
 
 
        int NX = ite - its + 3;
        int NY = kte - kts + 3;
        int NZ = jte - jts + 3;
        int NG1 = NX * NY * NZ;
        int NG = (ite - its + 1 ) * NY * (jte - jts + 1);
        int NG2 = (ime - ims + 1 ) * NY * (jme - jms + 1);

        float* d_p=NULL;
        float* d_p1=NULL;
        float* d_zp=NULL;
	double *d_x0=NULL;
	float* d_a0=NULL;
	float* d_a_helm=NULL;
	float* d_f0=NULL;
	float* d_c=NULL;
	float* d_r=NULL;
	float* d_ar=NULL;
	float* d_ap=NULL;
	
	double *d_aps=NULL; 
	float *result=NULL;
	float *dd=NULL;
	double *d_c1=NULL;
	float *d_xc=NULL; 
	float *d_Gm=NULL; 
	float *d_buff=NULL; 
	float *d_cm_helm=NULL; 
	float *d_ac=NULL; 

        static int init_flag=0; 
        #ifdef cagcr
        float *h_p=new float[NG1*(iter_max+2)];
        float *h_ap=new float[NG*(iter_max+1)];
        float *h_r1=new float[NG1];
        #else
        float *h_p=new float[NG1*iter_max];
        float *h_ap=new float[NG*(iter_max+1)];
        #endif
        #ifdef CUDA
        float *buff=new float[2*(NX*NY+NY*NZ)];
        float *h_r=new float[NG];
        #else
        float *h_ar=new float[NG];
        float *h_r=new float[NG];
        #endif

        #ifdef CUDA
        #ifdef share_mem
        dim3 threadsPerBlock_mat(16, 4, 6);
        dim3 numBlocks_mat(ceil(((float)NX)/(threadsPerBlock_mat.x-2)),ceil(((float) NY-1)/(threadsPerBlock_mat.y)),ceil(((float)NZ)/(threadsPerBlock_mat.z-2)));
        #else
        dim3 threadsPerBlock_mat(16, 4, 4);
        dim3 numBlocks_mat(ceil(((float)NX)/16),ceil(((float) NY)/4),ceil(((float)NZ)/4));
        #endif
        hipblasHandle_t h;
        hipblasCreate(&h);
        hipblasSetPointerMode(h, HIPBLAS_POINTER_MODE_DEVICE);

        //GPTLstart ("hipmalloc");
        if ( init_flag==0) {
	    hipMalloc( (void **)(&dd), sizeof(float) );
            #ifdef CA
            CHECK(hipMalloc( (void **)(&d_xc), (2*iter_max+1)*sizeof(float) )); 
            CHECK(hipMalloc( (void **)(&result), (iter_max+1)*(iter_max+1)*sizeof(float) )); 
            CHECK(hipMalloc( (void **)(&d_Gm), (iter_max+1)*sizeof(float) )); 
            CHECK(hipMalloc((void **) &d_p, NG1*(iter_max+2)*sizeof(float)));
            CHECK(hipMalloc((void **) &d_p1, NG1*sizeof(float)));
	    CHECK(hipMemset(result,0,(iter_max+1)*(iter_max+1)*sizeof(float)));
            #else
            CHECK(hipMalloc( (void **)(&result), iter_max*sizeof(float) )); 
            CHECK(hipMalloc((void **) &d_p, NG1*(iter_max+1)*sizeof(float)));
            CHECK(hipMalloc((void **) &d_zp, NG1*sizeof(float)));
            CHECK(hipMalloc((void **) &d_ar, NG*sizeof(float))); //ar(its:ite,kts-1:kte+1,jts:jte)
            CHECK(hipMalloc((void **) &d_aps, (iter_max+1)*sizeof(double)));
	    CHECK(hipMalloc( (void **)(&d_c1), (iter_max+2)*sizeof(double)));
	    CHECK(hipMemset(d_ar,0,NG*sizeof(float)));
            #endif

            CHECK(hipMalloc((void **) &d_ap, NG*(iter_max+1)*sizeof(float)));
            CHECK(hipMemset(d_ap,0,NG*(iter_max+1)*sizeof(float)));
	    CHECK(hipMalloc( (void **)(&d_buff), 2*(NY*NX+NY*NZ)*sizeof(float) ));
            CHECK(hipMalloc((void **) &d_x0, NG2*sizeof(double)));
            CHECK(hipMalloc((void **) &d_a0, NG*19*sizeof(float)));
            CHECK(hipMalloc((void **) &d_f0, NG*sizeof(float)));
            CHECK(hipMalloc((void **) &d_cm_helm, NG1*7*sizeof(float)));
            CHECK(hipMalloc((void **) &d_r, NG*sizeof(float)));
	    CHECK(hipMemset(d_r,0,NG*sizeof(float)));
            
#ifdef trans
            CHECK(hipMalloc((void **) &d_a_helm, NG*19*sizeof(float)));
            CHECK(hipMemcpy(d_a_helm,a0,19*NG*sizeof(float),hipMemcpyHostToDevice));
#else
        CHECK(hipMemcpy(d_a0,a0,19*NG*sizeof(float),hipMemcpyHostToDevice));
#endif
            }
        //GPTLstop ("hipmalloc");

        //GPTLstart ("hipmemcpy");
        //CHECK(hipMemcpy(d_a_helm,a0,19*NG*sizeof(float),hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_f0,f0,NG*sizeof(float),hipMemcpyHostToDevice));
	CHECK(hipMemcpy(d_cm_helm,cm_helm,7*NG1*sizeof(float),hipMemcpyHostToDevice));
        //CHECK(hipMemcpy(d_x0,x0,NG2*sizeof(double),hipMemcpyHostToDevice));
        //GPTLstop ("hipmemcpy");

        #endif
        #ifndef CA
        double aps[iter_max+1],b[iter_max+1];
        double c1[iter_max+2],c2[iter_max+2];
        float c1_r4[iter_max+2];
        #endif
        init_flag=1;
        // float ac[iter_max+1],aps[iter_max+1],b[iter_max+1];
        // float c1[iter_max+4],c2[iter_max+4];
