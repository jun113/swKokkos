// Automatically translated using m2cpp 2.0 on 2022-06-20 00:24:05

#ifndef MYGETVAR_M_HPP
#define MYGETVAR_M_HPP
//#include <cstdio>

//#include <armadillo>
#include <netcdf.h>
//using namespace arma ;
#define ERR {if(status!=NC_NOERR)printf("Error at line=%d: %s %s\n", __LINE__, varname, nc_strerror(status));}

int mygetVar(int ncid,const char* varname, float * data)
{
    int varid, status;
    int  rh_id;
    size_t dim_len, max_len=0;
    nc_type rh_type;
    int rh_ndims;
    int rh_dimids[NC_MAX_VAR_DIMS];
    int rh_natts;
    status = nc_inq_varid(ncid, varname, &varid);ERR

    status = nc_inq_var(ncid, varid, 0, &rh_type, &rh_ndims, rh_dimids, &rh_natts);ERR
    for (int i = 0; i < rh_ndims; i++) {
        nc_inq_dimlen(ncid, rh_dimids[i], &dim_len);
	if (dim_len > max_len) {
    	    max_len = dim_len;
	}
    }
    //if(rh_ndims>1) printf("rh_ndims=%d, len=%ld\n",rh_ndims, max_len);
    status = nc_get_var_float(ncid,varid,data);ERR
    return max_len;
}

int mygetVar(int ncid,const char* varname, double * data)
{
    int varid, status;
    int  rh_id;
    size_t dim_len, max_len=0;
    nc_type rh_type;
    int rh_ndims;
    int rh_dimids[NC_MAX_VAR_DIMS];
    int rh_natts;
    status = nc_inq_varid(ncid, varname, &varid);ERR

    status = nc_inq_var(ncid, varid, 0, &rh_type, &rh_ndims, rh_dimids, &rh_natts);ERR
    for (int i = 0; i < rh_ndims; i++) {
        nc_inq_dimlen(ncid, rh_dimids[i], &dim_len);
	if (dim_len > max_len) {
    	    max_len = dim_len;
	}
    }
    //if(rh_ndims>1) printf("rh_ndims=%d, len=%ld\n",rh_ndims, max_len);
    status = nc_get_var_double(ncid,varid,data);ERR
    return max_len;
}
#endif
