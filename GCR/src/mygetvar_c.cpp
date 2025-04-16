#include <iostream>
#ifdef netcdf
#include <netcdf.h>
#include "getVar.hpp"
#endif
#include <cstring>
#include <math.h>
#include "gcr.h"
#include "norm.hpp"
//#include <mkl.h>
//#include <mkl_cblas.h>
//#define ERR {if(status!=NC_NOERR)printf("Error at line=%d: %s %s\n", __LINE__, varname, nc_strerror(status));}
//#define ERR1 {if(status!=NC_NOERR)printf("Error at line=%d:  %s\n", __LINE__, nc_strerror(status));}
//#define ISNETCDF4ATT "_IsNetcdf4"
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std ;
//#include <filesystem>

#ifdef netcdf
void check(int line, int status){
      if (status != NC_NOERR) {
          std::cout << "line "<<line<<", error code " <<status<<": "<<nc_strerror(status) << std::endl;
      }
}
void mygetVar_int(int ncid,const char* varname, int * data)
{
    int varid, status;
    check(23, nc_inq_varid(ncid, varname, &varid));
    check(38,nc_get_var_int(ncid,varid,data));
}
#endif

// 定义读取文件中数字的函数
std::vector<int> readNumbersFromFile(const std::string& filePath) {
    std::vector<int> numbers;
    std::ifstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return numbers; // 返回空的向量
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        int number;
        char comma;

        while (ss >> number) {
            numbers.push_back(number);
            ss >> comma; // 跳过逗号
        }
    }

    file.close();
    return numbers;
}


// Function to initialize a float array with random values
void mygetVar_rand(int ncid, const char* varName, float* a, size_t size, unsigned int seed, int its, int ite, int jts, int jte, int kts, int kte) {
    // Create a random number generator with a fixed seed
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dis(0.8f, 1.2f); // Random numbers between 0 and 1
    //std::uniform_real_distribution<float> dis2(0.9f, 1.1f); // Random numbers between 0 and 1
    //std::normal_distribution<float> dis(1.0f, 0.5f);
    //std::normal_distribution<float> dis2(1.0f, 0.1f);

    int NX = ite - its + 3;
    int NY = kte - kts + 3;
    int NZ = jte - jts + 3;
    // Fill the array with random numbers
    //for (size_t i = 0; i < size; ++i) {
    //    array[i] = dis(gen);
    //}
    int i,j,k;
    for(j=jts;j<=jte;j++){
        for(k=kts-1;k<=kte+1;k++){
            for(i=its;i<=ite;i++){
                 a[index_a(1,i,k,j)]= 1;
                 a[index_a(2,i,k,j)]=-0.17 * dis(gen);
                 a[index_a(3,i,k,j)]= a[index_a(2,i,k,j)]; // * dis2(gen);
                 a[index_a(4,i,k,j)]=-0.04 * dis(gen);
                 a[index_a(5,i,k,j)]= a[index_a(4,i,k,j)];// * dis2(gen);
                 a[index_a(6,i,k,j)]=3.5e-6 * dis(gen);
                 a[index_a(7,i,k,j)]=-3.5e-6 * dis(gen);
                 a[index_a(8,i,k,j)]= a[index_a(6,i,k,j)];// * dis2(gen);
                 a[index_a(9,i,k,j)]= a[index_a(7,i,k,j)];// * dis2(gen);
                 a[index_a(10,i,k,j)]=0.45 * dis(gen);
                 a[index_a(15,i,k,j)]= a[index_a(10,i,k,j)];// * dis2(gen);
                 a[index_a(11,i,k,j)]= 2.5e-3 * dis(gen);
                 a[index_a(12,i,k,j)]= a[index_a(11,i,k,j)];// * dis2(gen);
                 a[index_a(13,i,k,j)]= -3e-4 * dis(gen);
                 a[index_a(14,i,k,j)]= a[index_a(13,i,k,j)];// * dis2(gen);
                 a[index_a(16,i,k,j)]=-1.4e-3 * dis(gen);
                 a[index_a(17,i,k,j)]= a[index_a(16,i,k,j)];// * dis2(gen);
                 a[index_a(18,i,k,j)]=2.8e-4 * dis(gen);
                 a[index_a(19,i,k,j)]= a[index_a(18,i,k,j)];// * dis2(gen);
            }
        }
    }
 //   std::cout << "Initialized " << varName << " with random float values.\n";
}
/*
// Function to initialize a double array with random values
void mygetVar_rand_d(int ncid, const char* varName, double* array, size_t size, unsigned int seed) {
    // Create a random number generator with a fixed seed
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0); // Random numbers between 0 and 1

    // Fill the array with random numbers
    for (size_t i = 0; i < size; ++i) {
        array[i] = dis(gen);
    }

 //   std::cout << "Initialized " << varName << " with random double values.\n";
}
*/
//extern "C" void module_gcr_mp_solve_helmholts_(float *a_helm,float *b_helm,float *cm_helm,float *threshold, double *pi,int* ijk_index, int *myprcid);
extern "C" void cu_solve_helmholts_(float *a_helm,float *b_helm,float *cm_helm,float *threshold, double *pi, int *ijk_index,  int *myprcid, size_t size_pi);
extern "C" void assign_indices_(int * ijk_index);

//extern "C" void module_halo_mp_getijk_(int *);
extern "C" void ILU_5(float *a_helm, float *cm_helm,
	int jend, int its, int ite, int jts, int jte, int kts, int kte);



void compareAndPrintDifferences(int rank, const std::vector<int>& vec1, const std::vector<int>& vec2) {
    if (vec1.size() != vec2.size()) {
        std::cerr << "Vectors are of different sizes!" << std::endl;
        return;
    }

    for (size_t i = 0; i < vec1.size(); ++i) {
        if (vec1[i] != vec2[i]) {
            std::cout <<"Rank:"<< rank<< ", Difference at index " << i << ": " << vec1[i] << " vs " << vec2[i] << std::endl;
        }
    }
}
extern "C" void mygetvar_c_(int *mytid){
    int ncid,ncid1, dim_len;
    int rank=mytid[0];
    int  idep, jdep, ids,ide,jds,jde,kds,kde,ims,ime,
         jms,jme,kms,kme,its,ite,jts,jte,kts,kte;
    std::string dirname = std::string("../data/postvar") + (rank < 10 ? "0" : "") + std::to_string(rank);
#ifdef netcdf
    check(75, nc_open(dirname.c_str(), NC_NOWRITE, &ncid));
#endif

    std::vector<int> ijk_index(20);// = readNumbersFromFile(filename);
    initialize_ijk_index(mytid, ijk_index.data());
    
   // int  idep, jdep, ids,ide,jds,jde,kds,kde,ims,ime,
    //     jms,jme,kms,kme,its,ite,jts,jte,kts,kte;
    assign_indices(ijk_index.data(), idep, jdep, ids, ide, jds, jde, 
            kds, kde, ims, ime, jms, jme, kms, kme,
                   its, ite, jts, jte, kts, kte);
    assign_indices_(ijk_index.data());

    size_t size_a_helm = 19 * (ite-its+1) * (kte-kts+3) * (jte-jts+1);
    size_t size_b_helm = (ite-its+1) * (kte-kts+3) * (jte-jts+1);
    size_t size_cm_helm = 7 * (ite-its+3) * (kte-kts+3) * (jte-jts+3);
    size_t size_pi = (ime-ims+1) * (kte-kts+3) * (jme-jms+1);

    float* a_helm = new float[size_a_helm];
    float* b_helm = new float[size_b_helm];
    float* cm_helm = new float[size_cm_helm];
    float* cm_helm1 = new float[size_cm_helm];
    float* pi = new float[size_pi];
    double* dbl_pi = new double[size_pi];

    memset (a_helm,   0, sizeof(float) * size_a_helm);
    memset (b_helm,   0, sizeof(float) *size_b_helm);
    memset (cm_helm,  0, sizeof(float) *size_cm_helm);
    memset (cm_helm1, 0, sizeof(float) *size_cm_helm);
    memset (pi,       0, sizeof(float) *size_pi);
    memset (dbl_pi,   0, sizeof(double) *size_pi);
#ifdef netcdf
    mygetVar(ncid, "a_helm",a_helm);
    mygetVar(ncid, "b_helm",b_helm);
    mygetVar(ncid, "cm_helm",cm_helm);
    mygetVar(ncid, "pi",pi);
#endif
#define rand
#ifdef rand
    unsigned int seed = 42; // Fixed seed for reproducibility
    mygetVar_rand(ncid, "a_helm", a_helm, size_a_helm, seed, its, ite, jts, jte, kts, kte);
    //mygetVar_rand(ncid, "b_helm", b_helm, size_b_helm, seed);
    for(int i=0;i<size_b_helm;i++){
        b_helm[i]=1.0;
    }
#endif
    print_norm(a_helm, "a_helm", size_a_helm);
    print_norm(cm_helm, "cm_helm", size_cm_helm);
    float threshold=1e-11;
    int i,j,k;
    int NX = ite - its + 3;
    int NY = kte - kts + 3;
    int NZ = jte - jts + 3;
    int NG1 = NX * NY * NZ;
    int NG2= (ime - ims + 1 ) * NY * (jme - jms + 1);
    int jend   = min(jde-1,jte);
    //jend=jte;
    ILU_5(a_helm,cm_helm1, jend,its, ite, jts, jte, kts, kte);
    print_norm(cm_helm1, "cm_helm", size_cm_helm);

    for(i=0;i<size_pi;i++){
        dbl_pi[i]=1.0;  //pi[i];
    }
    //print_norm(b_helm, "b_helm", size_b_helm);
    cu_solve_helmholts_(a_helm,b_helm,cm_helm1,&threshold, dbl_pi, 
        ijk_index.data(), mytid, size_pi);

    delete [] a_helm;
    delete [] b_helm;
    delete [] cm_helm;
    delete [] cm_helm1;
    delete [] pi;
    delete [] dbl_pi;
return;
}
