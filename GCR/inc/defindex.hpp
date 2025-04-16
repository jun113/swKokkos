	int NX = ite - its + 3;
	int NY = kte - kts + 3;
        int NZ = jte - jts + 3;
	int ims=its-2;
	int ime=ite+2;
	int jms=jts-2;
	int jme=jte+2;
        int ibegin = its;
        int iend   = ite;
        int jbegin = jts;
        int jend   = min(jde-1,jte);
        int i= blockDim.x * blockIdx.x + threadIdx.x +its-1;
        int k= blockDim.y * blockIdx.y + threadIdx.y +kts-1;
        int j= blockDim.z * blockIdx.z + threadIdx.z +jts-1;

