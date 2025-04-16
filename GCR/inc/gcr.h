#include<string.h> 
#define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
#define index3b(i,k,j) (i-its+1 + (k-kts+1)*NX + (j-jts+1)*NX*NY)
#define index4(i,k,j,m) (i-its + ((k)+(1-kts))*(ite-its+1) + (j-jts)*NY*(ite-its+1) + (m) * NY*(ite-its+1)*(jte-jts+1))
#define index4b(i,k,j,m) (i+(1-its) + ((k)+(1-kts))*NX + ((j)+(1-jts))*NX*NY + (m) * NY*NX*NZ)
#ifdef CUDA
#define trans
#endif
#ifdef trans
        #define index_a(m, i,k,j) ((i-its)  + (k-kts+1) * (ite-its+1) + (j-jts) * (ite-its+1) * (kte-kts+3) +(m-1)* (jte-jts+1) * (ite-its+1) * (kte-kts+3))
#else
        #define index_a(m,i,k,j) (m-1+(i-its)*19+ (k-kts+1)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)*NY)
#endif
        #define index_a0(m,i,k,j) (m-1+(i-its)*19+ (k-kts+1)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)*NY)
//trans

#define index_x(i,k,j) (i-ims + (k-kts+1)*(ime-ims+1) + (j-jms)*(ime-ims+1)*(kte-kts+3))
#define index_x0(i,k,j) (i+1 + (k)*(ime-ims+1) + (j+1)*(ime-ims+1)*NY)
//#define smeIndex(i,k,j) (i + (k)*(blockDim.x)+ (j)*((blockDim.x)*(blockDim.y)))
	#define cuinitialize(i,j)\
		CHECK(hipMemset(i,0,(j)*sizeof(float)));
	#define initialize(i,j) for(int ii=0;ii<j;ii++) i[ii]=0.0;
#define CHECK(call)                                                            \
{                                                                              \
        const hipError_t error = call;                                            \
        if (error != hipSuccess )                       \
        {                                                                          \
                    fprintf(stderr,"Error: %s:%d, ", __FILE__, __LINE__);                 \
                    fprintf(stderr,"code: %d, reason: %s\n", error,                       \
                                            hipGetErrorString(error));                                    \
                    exit(1);                                                               \
                }                                                                          \
}

#ifdef __HIP_PLATFORM_HCC__
const int warplen = 64;
#else
const int warplen = 32;
#endif


