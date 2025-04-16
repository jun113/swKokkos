#ifdef __sw_slave__

#include "athread_param.h"
#include <simd.h>

#include "slave.h"
//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT  = 0;
__thread_local crts_rply_t  dma_rply = 0;

void athread_readyt_4_1_slave (struct paramReadyt4* arg_param) {
	int total_len = arg_param->IMT * arg_param->JMT;
	double* h0l = arg_param->h0l;
	double* h0f = arg_param->h0f;
	double* h0  = arg_param->h0;
	double* pp  = arg_param->pp;
	double* gg  = arg_param->gg;
	double* vit = arg_param->vit;
	double* dzp = arg_param->dzp;
	double* ppa = arg_param->ppa;
	double* psa = arg_param->psa;

	int NUM_CORES = 64;
	int len = (total_len + NUM_CORES - 1) / NUM_CORES;

  double dzp_zero_half = dzp[0] * 0.5;

	for (int i = 0; i < len; ++i) {
		int idx = athread_tid*len + i;
		if (idx < total_len) {
      h0l[idx] = h0f[idx];
      h0f[idx] =  h0[idx];
      pp [idx] =  gg[idx] * vit[idx] * dzp_zero_half;
      ppa[idx] = psa[idx] * vit[idx];
		}
	}
	return ;
}

void athread_readyt_4_2_slave (struct paramReadyt4* arg_param) {
	int IMT = arg_param-> IMT;
	int JMT = arg_param-> JMT;
	int KM  = arg_param-> KM;
	double* h0l = arg_param->h0l;
	double* h0f = arg_param->h0f;
	double* h0  = arg_param->h0;
	double* pp  = arg_param->pp;
	double* gg  = arg_param->gg;
	double* vit = arg_param->vit;
	double* dzp = arg_param->dzp;
	double* ppa = arg_param->ppa;
	double* psa = arg_param->psa;

	// 4: gg pp ppa vit
	int allocatable_num_blk = get_allocatable_size() / (4 * KM * sizeof(double));
	if (allocatable_num_blk == 0) {
		printf ("Err in athread readyt 4 2 slave\n");
		exit(0);
	}

	int NUM_CORES = 64;
	int total_len = JMT * IMT;
	int blk_per_cpe = (total_len + NUM_CORES - 1) / NUM_CORES;

	int num_blk = (blk_per_cpe + allocatable_num_blk - 1) / allocatable_num_blk;

	int local_blk = blk_per_cpe;
	if (blk_per_cpe > allocatable_num_blk) {
		local_blk = allocatable_num_blk;
	}
	int remainder_blk = blk_per_cpe % local_blk;

  size_t size_tile = local_blk * sizeof(double);
	size_t size_blk = size_tile * KM;
	double* ldm_gg  = (double*) ldm_malloc (size_blk);
	double* ldm_pp  = (double*) ldm_malloc (size_blk);
	double* ldm_ppa = (double*) ldm_malloc (size_blk);
	double* ldm_vit = (double*) ldm_malloc (size_blk);

	int strideK = IMT * JMT;

	for (int idx_blk = 0; idx_blk < num_blk; ++idx_blk) {
		if ((idx_blk != num_blk - 1) || (idx_blk == 0) 
        || (idx_blk == num_blk - 1 && remainder_blk == 0)) {
			// dst, src, len, bsize, stride
      // Blocked
			// athread_dma_get_stride (ldm_gg, 
			// 		&gg[athread_tid*blk_per_cpe + idx_blk*local_blk], 
			// 		size_blk, local_blk * sizeof(double),
			// 		(strideK - local_blk) * sizeof(double));
			// athread_dma_get_stride (ldm_pp, 
			// 		&pp[athread_tid*blk_per_cpe + idx_blk*local_blk], 
			// 		size_blk, local_blk * sizeof(double),
			// 		(strideK - local_blk)  * sizeof(double));
			// athread_dma_get_stride (ldm_ppa, 
			// 		&ppa[athread_tid*blk_per_cpe + idx_blk*local_blk], 
			// 		size_blk, local_blk * sizeof(double),
			// 		(strideK - local_blk) * sizeof(double));
			// athread_dma_get_stride (ldm_vit, 
			// 		&vit[athread_tid*blk_per_cpe + idx_blk*local_blk], 
			// 		size_blk, local_blk * sizeof(double),
			// 		(strideK - local_blk)  * sizeof(double));

      // No blocked
			athread_dma_iget_stride (ldm_gg, 
					&gg[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					size_blk, size_tile,
					(strideK - local_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_pp, 
					&pp[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					size_blk, size_tile,
					(strideK - local_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_ppa, 
					&ppa[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					size_blk, size_tile,
					(strideK - local_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_vit, 
					&vit[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					size_blk, size_tile,
					(strideK - local_blk) * sizeof(double), &dma_rply);
      D_COUNT++;

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

			for (int idx = 0; idx < local_blk; ++idx) {
				int ii = idx_blk * local_blk + idx;
				if (ii >= blk_per_cpe) {break;}

        double pre_pp  = ldm_pp[idx];
        double pre_ppa = ldm_ppa[idx];
        double pre_gg  = ldm_gg[idx];
        double pre_dzp = dzp[0];
 
				for (int k = 1; k < KM; ++k) {
     	    ldm_pp[k*local_blk + idx] = ldm_vit[k*local_blk + idx] *
      		    (pre_pp + 0.5 * (ldm_gg[k*local_blk + idx] * dzp[k] 
                  + pre_gg * pre_dzp));
		      ldm_ppa[k*local_blk + idx] = ldm_vit[k*local_blk + idx] *
  		        (pre_ppa + pre_gg * pre_dzp);

          pre_gg  = ldm_gg [k*local_blk + idx];
          pre_pp  = ldm_pp [k*local_blk + idx];
          pre_ppa = ldm_ppa[k*local_blk + idx];
          pre_dzp = dzp[k];
 
				}
			}
			athread_dma_iput_stride (&pp[athread_tid*blk_per_cpe + idx_blk*local_blk],
          ldm_pp,  size_blk, size_tile, 
          (strideK - local_blk)  * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iput_stride (&ppa[athread_tid*blk_per_cpe + idx_blk*local_blk],
          ldm_ppa, size_blk, size_tile, 
          (strideK - local_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
		} else {
      // No blocked
			athread_dma_iget_stride (ldm_gg, 
					&gg[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					remainder_blk * KM * sizeof(double), remainder_blk * sizeof(double),
					(strideK - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_pp, 
					&pp[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					remainder_blk * KM * sizeof(double), remainder_blk * sizeof(double),
					(strideK - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_ppa, 
					&ppa[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					remainder_blk * KM * sizeof(double), remainder_blk * sizeof(double),
					(strideK - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_vit, 
					&vit[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					remainder_blk * KM * sizeof(double), remainder_blk * sizeof(double),
					(strideK - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

			for (int idx = 0; idx < remainder_blk; ++idx) {
				int ii = idx_blk * local_blk + idx;
				if (ii >= blk_per_cpe) {break;}

        double pre_pp  = ldm_pp[idx];
        double pre_ppa = ldm_ppa[idx];
        double pre_gg  = ldm_gg[idx];
        double pre_dzp = dzp[0];
 
				for (int k = 1; k < KM; ++k) {
     	    ldm_pp[k*remainder_blk + idx] = ldm_vit[k*remainder_blk + idx] *
      		    (pre_pp + 0.5 * (ldm_gg[k*remainder_blk + idx] * dzp[k] 
                  + pre_gg * pre_dzp));
		      ldm_ppa[k*remainder_blk + idx] = ldm_vit[k*remainder_blk + idx] *
  		        (pre_ppa + pre_gg * pre_dzp);

          pre_gg  = ldm_gg [k*remainder_blk + idx];
          pre_pp  = ldm_pp [k*remainder_blk + idx];
          pre_ppa = ldm_ppa[k*remainder_blk + idx];
          pre_dzp = dzp[k];
 
				}
			}
			athread_dma_iput_stride (
          &pp[athread_tid*blk_per_cpe + idx_blk*local_blk], ldm_pp, 
					remainder_blk * KM * sizeof(double), remainder_blk * sizeof(double),
					(strideK - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iput_stride (
          &ppa[athread_tid*blk_per_cpe + idx_blk*local_blk], ldm_ppa, 
					remainder_blk * KM * sizeof(double), remainder_blk * sizeof(double),
					(strideK - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
		}
	}

	ldm_free (ldm_gg,  size_blk);
	ldm_free (ldm_pp,  size_blk);
	ldm_free (ldm_ppa, size_blk);
	ldm_free (ldm_vit, size_blk);

	return ;
}
void athread_upwell_slave (struct paramUpwell *arg_param) {
  int KM  = arg_param->KM;
  int JMT = arg_param->JMT;
  int IMT = arg_param->IMT;

  double* work = arg_param->work;
  double* ws   = arg_param->ws;
  double* wka  = arg_param->wka;
  double* vit  = arg_param->vit;
  double* dzp  = arg_param->dzp;
  double* ohbt  = arg_param->ohbt;
  double* h0wk  = arg_param->h0wk;
  int total_len = JMT * IMT;
  int blk_per_cpe = (total_len + 63) / 64;

	// 3: ws, vit, wka
  // get_allocatable_size() - blk_per_cpe), Considering work
	int allocatable_num_blk = (get_allocatable_size() - blk_per_cpe * sizeof(double)) 
      / (3 * KM * sizeof(double));
	if (allocatable_num_blk == 0) {
		printf ("Err in athread upwell slave\n");
		exit(0);
	}
	int num_blk = (blk_per_cpe + allocatable_num_blk - 1) 
      / allocatable_num_blk;

	int local_blk = blk_per_cpe;
	if (blk_per_cpe > allocatable_num_blk) {
		local_blk = allocatable_num_blk;
	}
	int remainder_blk = blk_per_cpe % local_blk;

  size_t size_tile = local_blk * sizeof(double);
	size_t size_blk = size_tile * KM;

  double* ldm_work = (double*) ldm_malloc (size_tile);
	double* ldm_ws   = (double*) ldm_malloc (size_blk);
	double* ldm_wka  = (double*) ldm_malloc (size_blk);
	double* ldm_vit  = (double*) ldm_malloc (size_blk);

  int start_cpe = athread_tid * blk_per_cpe;
  // printf("localblk %d, blk_per_cpe %d, size tile %d, blk %d, remainder %d\n", 
  // local_blk, blk_per_cpe,
  // size_tile, size_blk, remainder_blk);


	for (int idx_blk = 0; idx_blk < num_blk; ++idx_blk) {
    int stride_blk = idx_blk * local_blk;
		if ((idx_blk != num_blk - 1) || (idx_blk == 0) 
        || (idx_blk == num_blk - 1 && remainder_blk == 0)) {
      // No blocked
			athread_dma_iget_stride (ldm_ws, 
					&ws[start_cpe + stride_blk], 
					size_blk, size_tile,
					(total_len - local_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_wka, 
					&wka[start_cpe + stride_blk], 
					size_blk, size_tile,
					(total_len - local_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_vit, 
					&vit[start_cpe + stride_blk], 
					size_blk, size_tile,
					(total_len - local_blk) * sizeof(double), &dma_rply);
      D_COUNT++;

			for (int idx = 0; idx < local_blk; ++idx) {
				int ii = stride_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        ldm_work[idx] = 0.0;
      }

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = stride_blk + idx;
			  if (ii >= blk_per_cpe) {break;}

        int index = start_cpe + stride_blk + idx;
        if (index > total_len) {break;}
        int j = index / IMT;
        // if (j < 1 || j >= JMT-1) {break;}
        int i = index % IMT;
        // if (i < 1 || i >= IMT-1) {break;}
        if (i >= 1 && i < IMT-1 && j >= 1 && j < JMT-1) {
          for (int k = 0; k < KM; ++k) {
            ldm_work[idx] = ldm_work[idx] - dzp[k]
                * ldm_wka[k*local_blk + idx] * ldm_vit[k*local_blk + idx];
          }
        }

      }
      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = stride_blk + idx;
			  if (ii >= blk_per_cpe) {break;}

        int index = start_cpe + stride_blk + idx;
        if (index >= total_len) {break;}
        int j = index / IMT;
        int i = index % IMT;
        if (i >= 1 && i < IMT-1 && j >= 1 && j < JMT-1) {

          double pre_ws  = ldm_ws[idx];

          for (int k = 1; k < KM; ++k) {
            ldm_ws[k*local_blk + idx] = ldm_vit[k*local_blk + idx]
                * (pre_ws + dzp[k-1] * (ldm_work[idx] 
                * ohbt[j*IMT + i] + ldm_wka[(k-1)*local_blk + idx]));
            pre_ws = ldm_ws[k*local_blk + idx];
          }
        }
      }

			for (int idx = 0; idx < local_blk; ++idx) {
        int ii = stride_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + stride_blk + idx;
        if (index >= total_len) {break;}
        int j = index / IMT;
        // if (j < 1 || j >= JMT-1) {break;}
        int i = index % IMT;
        // if (i < 1 || i >= IMT-1) {break;}
        if (i >= 1 && i < IMT-1 && j >= 1 && j < JMT-1) {
          ldm_work[idx] = 1.0 / (1.0 + h0wk[j*IMT + i] * ohbt[j*IMT + i]);
        }
      }
			for (int idx = 0; idx < local_blk; ++idx) {
        int ii = stride_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + stride_blk + idx;
        if (index >= total_len) {break;}
        int j = index / IMT;
        // if (j < 1 || j >= JMT-1) {break;}
        int i = index % IMT;
        // if (i < 1 || i >= IMT-1) {break;}
        if (i >= 1 && i < IMT-1 && j >= 1 && j < JMT-1) {
          for (int k = 1; k < KM; ++k) {
            ldm_ws[k*local_blk + idx] *= ldm_work[idx];
          }
        }
      }

			athread_dma_iput_stride (&ws[athread_tid*blk_per_cpe + idx_blk*local_blk],
          ldm_ws, size_blk, size_tile, 
          (total_len - local_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);

    } else {
      // No blocked
			athread_dma_iget_stride (ldm_ws, 
					&ws[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_wka, 
					&wka[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_vit, 
					&vit[athread_tid*blk_per_cpe + idx_blk*local_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;

			for (int idx = 0; idx < remainder_blk; ++idx) {
				int ii = stride_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        ldm_work[idx] = 0.0;
      }

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      for (int idx = 0; idx < remainder_blk; ++idx) {
        int ii = stride_blk + idx;
			  if (ii >= blk_per_cpe) {break;}

        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int j = index / IMT;
        int i = index % IMT;
        if (i >= 1 && i < IMT-1 && j >= 1 && j < JMT-1) {

          for (int k = 0; k < KM; ++k) {
            ldm_work[idx] = ldm_work[idx] - dzp[k]
                * ldm_wka[k*remainder_blk + idx] 
                * ldm_vit[k*remainder_blk + idx];
          }
        }
      }
      for (int idx = 0; idx < remainder_blk; ++idx) {
        int ii = idx_blk * local_blk + idx;
			  if (ii >= blk_per_cpe) {break;}

        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int j = index / IMT;
        int i = index % IMT;
        if (i >= 1 && i < IMT-1 && j >= 1 && j < JMT-1) {

          double pre_ws  = ldm_ws[idx];

          for (int k = 1; k < KM; ++k) {
            ldm_ws[k*remainder_blk + idx] = ldm_vit[k*remainder_blk + idx]
                * (pre_ws + dzp[k-1] * (ldm_work[idx] 
                * ohbt[j*IMT + i] + ldm_wka[(k-1)*remainder_blk + idx]));
            pre_ws = ldm_ws[k*remainder_blk + idx];
          }
        }
      }

			for (int idx = 0; idx < remainder_blk; ++idx) {
				int ii = stride_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int j = index / IMT;
        int i = index % IMT;
        if (i >= 1 && i < IMT-1 && j >= 1 && j < JMT-1) {
          ldm_work[idx] = 1.0 / (1.0 + h0wk[j*IMT + i] * ohbt[j*IMT + i]);
        }
      }
			for (int idx = 0; idx < remainder_blk; ++idx) {
				int ii = stride_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int j = index / IMT;
        int i = index % IMT;
        if (i >= 1 && i < IMT-1 && j >= 1 && j < JMT-1) {
          for (int k = 1; k < KM; ++k) {
            ldm_ws[k*remainder_blk + idx] *= ldm_work[idx];
          }
        }
      }

			athread_dma_iput_stride (&ws[start_cpe + stride_blk], ldm_ws,
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);

    }
  }

  ldm_free (ldm_work, size_tile);
  ldm_free (ldm_ws,   size_blk);
  ldm_free (ldm_wka,  size_blk);
  ldm_free (ldm_vit,  size_blk);

  return ;
}

void athread_bclinc_4_slave (struct paramBclinc4 *arg_param) {
  int IMT = arg_param->IMT;
  int JMT = arg_param->JMT;
  int KM  = arg_param->KM;
  double G  = arg_param->G;
  double P5 = arg_param->P5;
  double onbb = arg_param->onbb;
  double od0 = arg_param->od0;
  double aa = arg_param->aa;
  double* h0bf = arg_param->h0bf;
  double* h0bl = arg_param->h0bl;
  double* work = arg_param->work;
  double* psa = arg_param->psa;
  double* vit = arg_param->vit;
  double* gg = arg_param->gg;
  double* dzp = arg_param->dzp;
  double* wka = arg_param->wka;

  int total_len = JMT * IMT;
  int blk_per_cpe = (total_len + 63) / 64;
	// 3: gg, vit, wka
	int allocatable_num_blk = get_allocatable_size() 
      / (3 * KM * sizeof(double));
	if (allocatable_num_blk == 0) {
		printf ("Err in athread upwell slave\n");
		exit(0);
	}

	int num_blk = (blk_per_cpe + allocatable_num_blk - 1) 
      / allocatable_num_blk;

	int local_blk = blk_per_cpe;
	if (blk_per_cpe > allocatable_num_blk) {
		local_blk = allocatable_num_blk;
	}
	int remainder_blk = blk_per_cpe % local_blk;

  size_t size_tile = local_blk * sizeof(double);
	size_t size_blk = size_tile * KM;
  size_t size_stride_layer = (total_len - local_blk) * sizeof(double);

	double* ldm_gg   = (double*) ldm_malloc (size_blk);
	double* ldm_wka  = (double*) ldm_malloc (size_blk);
	double* ldm_vit  = (double*) ldm_malloc (size_blk);

  int start_cpe = athread_tid * blk_per_cpe;

	for (int idx_blk = 0; idx_blk < num_blk; ++idx_blk) {
    int start_blk = idx_blk * local_blk;
		if ((idx_blk != num_blk - 1) || (idx_blk == 0) 
        || (idx_blk == num_blk - 1 && remainder_blk == 0)) {
      // No blocked
			athread_dma_iget_stride (ldm_gg, &gg[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_vit, &vit[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;

			for (int idx = 0; idx < local_blk; ++idx) {
				int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        h0bf[index] = h0bf[index] * onbb;
        work[index] = aa * h0bf[index] + (1.0 - aa) * h0bl[index];
      }

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        double pre_wkk;
        double cur_wkk;
        pre_wkk = (psa[index] * od0 + work[index] * G) * ldm_vit[idx];
        for (int k = 0; k < KM; ++k) {
          cur_wkk = pre_wkk - ldm_gg[k*local_blk + idx] 
              * dzp[k] * ldm_vit[k*local_blk + idx];
          ldm_wka[k*local_blk + idx] = P5 * (pre_wkk + cur_wkk);
          pre_wkk = cur_wkk;
        }
      }
			athread_dma_iput_stride (&wka[start_cpe + start_blk],
          ldm_wka, size_blk, size_tile, 
          (total_len - local_blk)  * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
    } else {
      // No blocked
			athread_dma_iget_stride (ldm_gg, 
					&gg[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_vit, 
					&vit[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;

			for (int idx = 0; idx < remainder_blk; ++idx) {
				int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        h0bf[index] = h0bf[index] * onbb;
        work[index] = aa * h0bf[index] + (1.0 - aa) * h0bl[index];
      }

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      for (int idx = 0; idx < remainder_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        double pre_wkk;
        double cur_wkk;
        pre_wkk = (psa[index] * od0 + work[index] * G) * ldm_vit[idx];
        for (int k = 0; k < KM; ++k) {
          cur_wkk = pre_wkk - ldm_gg[k*remainder_blk + idx] 
              * dzp[k] * ldm_vit[k*remainder_blk + idx];
          ldm_wka[k*remainder_blk + idx] = P5 * (pre_wkk + cur_wkk);
          pre_wkk = cur_wkk;
        }
      }
			athread_dma_iput_stride (&wka[start_cpe + start_blk], ldm_wka, 
          remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
    }
  }

  ldm_free (ldm_gg,  size_blk);
  ldm_free (ldm_wka, size_blk);
  ldm_free (ldm_vit, size_blk);
  return ;
}

void athread_bclinc_14_slave (
    struct paramBclinc14 *arg_param) {
  int KM = arg_param->KM;
  int JMT = arg_param->JMT;
  int IMT = arg_param->IMT;
  double dtc2 = arg_param->dtc2;
  double aidif = arg_param->aidif;
  int* kmu = arg_param->kmu;
  double* odzt = arg_param->odzt;
  double* odzp = arg_param->odzp;
  double* sbcy = arg_param->sbcy;
  double* bbcy = arg_param->bbcy;
  double* vp = arg_param->vp;
  double* dlv = arg_param->dlv;
  double* wka = arg_param->wka;
  double* viv = arg_param->viv;
  double* akmu = arg_param->akmu;

  int total_len = JMT * IMT;
  double a8[KM], b8[KM], c8[KM], d8[KM];
  double e8[KM+1], f8[KM+1];
  int blk_per_cpe = (total_len + 63) / 64;
	// 5: wka, vp, dlv, akmu, viv
	int allocatable_num_blk = (get_allocatable_size() - 7 * KM * sizeof(double)) 
      / (5 * KM * sizeof(double));
	if (allocatable_num_blk == 0) {
		printf ("Err in athread upwell slave\n");
		exit(0);
	}

	int num_blk = (blk_per_cpe + allocatable_num_blk - 1) 
      / allocatable_num_blk;

	int local_blk = blk_per_cpe;
	if (blk_per_cpe > allocatable_num_blk) {
		local_blk = allocatable_num_blk;
	}
	int remainder_blk = blk_per_cpe % local_blk;

  size_t size_tile = local_blk * sizeof(double);
	size_t size_blk = size_tile * KM;
  size_t size_stride_layer = (total_len - local_blk) * sizeof(double);

	double* ldm_wka  = (double*) ldm_malloc (size_blk);
	double* ldm_vp   = (double*) ldm_malloc (size_blk);
	double* ldm_dlv  = (double*) ldm_malloc (size_blk);
	double* ldm_akmu = (double*) ldm_malloc (size_blk);
	double* ldm_viv  = (double*) ldm_malloc (size_blk);

  int start_cpe = athread_tid * blk_per_cpe;

	for (int idx_blk = 0; idx_blk < num_blk; ++idx_blk) {
    int start_blk = idx_blk * local_blk;
		if ((idx_blk != num_blk - 1) || (idx_blk == 0) 
        || (idx_blk == num_blk - 1 && remainder_blk == 0)) {
      // No blocked
			athread_dma_iget_stride (ldm_vp, &vp[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_dlv, &dlv[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_akmu, &akmu[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_viv, &viv[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;

      // vp, dlv, wka
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int i = index % IMT;
        int j = index / IMT;
        if (i >= 2 && i < IMT-2 && j >= 2 && j < JMT-2) {
          for (int k = 0; k < KM; ++k) {
            ldm_wka[k*local_blk + idx] = ldm_vp[k*local_blk + idx]
                + ldm_dlv[k*local_blk + idx] * dtc2;
          }
        }
      }
      // akmu, vlv
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      // invtriu
      double* wk    = ldm_wka;
      double* topbc = sbcy;
      double* bombc = bbcy;
      double* dcb   = ldm_akmu;
      double  c2dtc = dtc2;
      int k;
      double c2dtc_times_aidif = c2dtc * aidif;

      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int i = index % IMT;
        int j = index / IMT;
        if (i >= 2 && i < IMT-2 && j >= 2 && j < JMT-2) {
          int kz = kmu[j*IMT + i];
          if (kz > 0 && kz <= KM) {
            for (k = 1; k < kz; ++k) {
              a8[k] = dcb[(k-1)*local_blk + idx] * odzt[k] 
                  * odzp[k] * c2dtc_times_aidif;
              d8[k] = wk[k*local_blk + idx];
            }
            for (k = 1; k < kz - 1 ; ++k) { 
              c8[k] = dcb[k*local_blk + idx] * odzt[k + 1] 
                  * odzp[k] * c2dtc_times_aidif;
              b8[k] = 1.0 + a8[k] + c8[k];
              e8[k] = 0.0;
              f8[k] = 0.0;
            }
            //B.C. AT TOP
            k = 0;
            a8[k] = odzp[k] * c2dtc_times_aidif;
            c8[k] = dcb[k*local_blk + idx] * odzt[k+1] * odzp[k] 
                * c2dtc_times_aidif;
            b8[k] = 1.0 + c8[k];
            d8[k] = wk[k*local_blk + idx];
            e8[k] = 0.0;
            f8[k] = 0.0;
            //B.C. AT BOTTOM
            b8[kz-1] = 1.0 + a8[kz-1];
            c8[kz-1] = odzp[kz-1] * c2dtc_times_aidif;
            e8[kz] = 0.0;
            f8[kz] = 0.0;
            d8[kz-1] = wk[(kz-1)*local_blk + idx]
                - bombc[j*IMT + i] * odzp[kz-1] * c2dtc_times_aidif;
            //NOW INVERT
            for (k = kz - 1; k >= 0; --k) {
              double g0 = 1.0 / (b8[k] - c8[k] * e8[k+1]);
              e8[k] = a8[k] * g0;
              f8[k] = (d8[k] + c8[k] * f8[k+1]) * g0;
            }
            //B.C. AT SURFACE 
            double pre_wk = (e8[0] * topbc[j*IMT + i] + f8[0]) * ldm_viv[idx];
            wk[idx] = pre_wk;
            for (k = 1; k < kz; ++k) {
              pre_wk = (e8[k] * pre_wk + f8[k]) * ldm_viv[k*local_blk + idx];
              wk[k*local_blk + idx] = pre_wk;
            } 
          }
        }
      }

			athread_dma_iput_stride (&wka[start_cpe + start_blk],
          ldm_wka, size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
			// athread_dma_put_stride (&wka[start_cpe + start_blk],
      //     ldm_wka, size_blk, size_tile, size_stride_layer);
    } else {
      // No blocked
			athread_dma_iget_stride (ldm_vp, &vp[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_dlv, &dlv[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;

			athread_dma_iget_stride (ldm_akmu, &akmu[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_viv, &viv[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;

      // vp, dlv
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

			// athread_dma_get_stride (ldm_vp, &vp[start_cpe + start_blk], 
			// 		remainder_blk * sizeof(double) * KM, 
			// 		remainder_blk * sizeof(double), 
			// 		(total_len - remainder_blk) * sizeof(double));
			// athread_dma_get_stride (ldm_dlv, &dlv[start_cpe + start_blk], 
			// 		remainder_blk * sizeof(double) * KM, 
			// 		remainder_blk * sizeof(double), 
			// 		(total_len - remainder_blk) * sizeof(double));
			// athread_dma_get_stride (ldm_akmu, &akmu[start_cpe + start_blk], 
			// 		remainder_blk * sizeof(double) * KM, 
			// 		remainder_blk * sizeof(double), 
			// 		(total_len - remainder_blk) * sizeof(double));
			// athread_dma_get_stride (ldm_viv, &viv[start_cpe + start_blk], 
			// 		remainder_blk * sizeof(double) * KM, 
			// 		remainder_blk * sizeof(double), 
			// 		(total_len - remainder_blk) * sizeof(double));

      for (int idx = 0; idx < remainder_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int i = index % IMT;
        int j = index / IMT;
        if (i >= 2 && i < IMT-2 && j >= 2 && j < JMT-2) {
          for (int k = 0; k < KM; ++k) {
            ldm_wka[k*remainder_blk + idx] = ldm_vp[k*remainder_blk + idx]
                + ldm_dlv[k*remainder_blk + idx] * dtc2;
          }
        }
      }
      // akmu, vlv
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      // invtriu
      double* wk    = ldm_wka;
      double* topbc = sbcy;
      double* bombc = bbcy;
      double* dcb   = ldm_akmu;
      double  c2dtc = dtc2;
      int k;
      double c2dtc_times_aidif = c2dtc * aidif;

      for (int idx = 0; idx < remainder_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int i = index % IMT;
        int j = index / IMT;
        if (i >= 2 && i < IMT-2 && j >= 2 && j < JMT-2) {
          int kz = kmu[j*IMT + i];
          if (kz > 0 && kz <= KM) {
            for (k = 1; k < kz; ++k) {
              a8[k] = dcb[(k-1)*remainder_blk + idx] * odzt[k] 
                  * odzp[k] * c2dtc_times_aidif;
              d8[k] = wk[k*remainder_blk + idx];
            }
            for (k = 1; k < kz - 1 ; ++k) { 
              c8[k] = dcb[k*remainder_blk + idx] * odzt[k + 1] 
                  * odzp[k] * c2dtc_times_aidif;
              b8[k] = 1.0 + a8[k] + c8[k];
              e8[k] = 0.0;
              f8[k] = 0.0;
            }
            //B.C. AT TOP
            k = 0;
            a8[k] = odzp[k] * c2dtc_times_aidif;
            c8[k] = dcb[k*remainder_blk + idx] * odzt[k+1] * odzp[k] 
                * c2dtc_times_aidif;
            b8[k] = 1.0 + c8[k];
            d8[k] = wk[k*remainder_blk + idx];
            e8[k] = 0.0;
            f8[k] = 0.0;
            //B.C. AT BOTTOM
            b8[kz-1] = 1.0 + a8[kz-1];
            c8[kz-1] = odzp[kz-1] * c2dtc_times_aidif;
            e8[kz] = 0.0;
            f8[kz] = 0.0;
            d8[kz-1] = wk[(kz-1)*remainder_blk + idx]
                - bombc[j*IMT + i] * odzp[kz-1] * c2dtc_times_aidif;
            //NOW INVERT
            for (k = kz - 1; k >= 0; --k) {
              double g0 = 1.0 / (b8[k] - c8[k] * e8[k+1]);
              e8[k] = a8[k] * g0;
              f8[k] = (d8[k] + c8[k] * f8[k+1]) * g0;
            }
            //B.C. AT SURFACE 
            double pre_wk = (e8[0] * topbc[j*IMT + i] + f8[0]) * ldm_viv[idx]; 
            wk[idx] = pre_wk;
            for (k = 1; k < kz; ++k) {
              pre_wk = (e8[k] * pre_wk + f8[k]) * ldm_viv[k*remainder_blk + idx];
              wk[k*remainder_blk + idx] = pre_wk;
            } 
          }
        }
      }

			athread_dma_iput_stride (&wka[start_cpe + start_blk], ldm_wka, 
          remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
			// athread_dma_put_stride (&wka[start_cpe + start_blk], ldm_wka, 
      //     remainder_blk * sizeof(double) * KM, 
			// 		remainder_blk * sizeof(double), 
			// 		(total_len - remainder_blk) * sizeof(double));
    }
  }
  ldm_free (ldm_dlv,  size_blk);
  ldm_free (ldm_vp,   size_blk);
  ldm_free (ldm_wka,  size_blk);
  ldm_free (ldm_viv,  size_blk);
  ldm_free (ldm_akmu, size_blk);

  return ;
}

void athread_invtriu_slave (struct paramInvtriu* arg_param) {
  int KM        = arg_param->KM; 
  int JMT       = arg_param->JMT; 
  int IMT       = arg_param->IMT;
  double c2dtc  = arg_param->c2dtc; 
  double aidif  = arg_param->aidif; 
  int* kmu      = arg_param->kmu;
  double* wk    = arg_param->wk; 
  double* topbc = arg_param->topbc; 
  double* bombc = arg_param->bombc; 
  double* dcb   = arg_param->dcb;
  double* odzt  = arg_param->odzt; 
  double* odzp  = arg_param->odzp; 
  double* viv   = arg_param->viv; 

  int total_len = JMT * IMT;
  double a8[KM], b8[KM], c8[KM], d8[KM];
  double e8[KM+1], f8[KM+1];
  int blk_per_cpe = (total_len + 63) / 64;
	// 3: wk, dcb, viv
	int allocatable_num_blk = (get_allocatable_size() - (4 * KM + 2*(KM+1)) * sizeof(double)) 
      / (3 * KM * sizeof(double));
	if (allocatable_num_blk == 0) {
		printf ("Err in athread upwell slave\n");
		exit(0);
	}

	int num_blk = (blk_per_cpe + allocatable_num_blk - 1) 
      / allocatable_num_blk;

	int local_blk = blk_per_cpe;
	if (blk_per_cpe > allocatable_num_blk) {
		local_blk = allocatable_num_blk;
	}
	int remainder_blk = blk_per_cpe % local_blk;


  size_t size_tile = local_blk * sizeof(double);
	size_t size_blk = size_tile * KM;
  size_t size_stride_layer = (total_len - local_blk) * sizeof(double);

	double* ldm_wk  = (double*) ldm_malloc (size_blk);
	double* ldm_dcb = (double*) ldm_malloc (size_blk);
	double* ldm_viv = (double*) ldm_malloc (size_blk);

  // printf("localblk %d, blk_per_cpe %d, size tile %d, blk %d, remainder %d, numblk %d, alloc %d\n", 
  // local_blk, 
  // blk_per_cpe,
  // size_tile, 
  // size_blk, 
  // remainder_blk, 
  // num_blk, 
  // allocatable_num_blk);

  int start_cpe = athread_tid * blk_per_cpe;
  int k;
  double c2dtc_times_aidif = c2dtc * aidif;

	for (int idx_blk = 0; idx_blk < num_blk; ++idx_blk) {
    int start_blk = idx_blk * local_blk;
		if ((idx_blk != num_blk - 1) || (idx_blk == 0) 
        || (idx_blk == num_blk - 1 && remainder_blk == 0)) {
  // printf("1 idx_blk %d, localblk %d, blk_per_cpe %d, size tile %d, blk %d, remainder %d, numblk %d, alloc %d\n", 
  // idx_blk,
  // local_blk, 
  // blk_per_cpe,
  // size_tile, 
  // size_blk, 
  // remainder_blk, 
  // num_blk, 
  // allocatable_num_blk);
      // No blocked
			athread_dma_iget_stride (ldm_wk, &wk[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_dcb, &dcb[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_viv, &viv[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      // invtriu
      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int i = index % IMT;
        int j = index / IMT;
        if (i >= 2 && i < IMT-2 && j >= 2 && j < JMT-2) {
          int kz = kmu[j*IMT + i];
          if (kz > 0 && kz <= KM) {
            for (k = 1; k < kz; ++k) {
              a8[k] = ldm_dcb[(k-1)*local_blk + idx] * odzt[k] 
                  * odzp[k] * c2dtc_times_aidif;
              d8[k] = ldm_wk[k*local_blk + idx];
            }
            for (k = 1; k < kz - 1 ; ++k) { 
              c8[k] = ldm_dcb[k*local_blk + idx] * odzt[k + 1] 
                  * odzp[k] * c2dtc_times_aidif;
              b8[k] = 1.0 + a8[k] + c8[k];
              e8[k] = 0.0;
              f8[k] = 0.0;
            }
            //B.C. AT TOP
            k = 0;
            a8[k] = odzp[k] * c2dtc_times_aidif;
            c8[k] = ldm_dcb[k*local_blk + idx] * odzt[k+1] * odzp[k] 
                * c2dtc_times_aidif;
            b8[k] = 1.0 + c8[k];
            d8[k] = ldm_wk[k*local_blk + idx];
            e8[k] = 0.0;
            f8[k] = 0.0;
            //B.C. AT BOTTOM
            b8[kz-1] = 1.0 + a8[kz-1];
            c8[kz-1] = odzp[kz-1] * c2dtc_times_aidif;
            e8[kz] = 0.0;
            f8[kz] = 0.0;
            d8[kz-1] = ldm_wk[(kz-1)*local_blk + idx]
                - bombc[j*IMT + i] * odzp[kz-1] * c2dtc_times_aidif;
            //NOW INVERT
            for (k = kz - 1; k >= 0; --k) {
              double g0 = 1.0 / (b8[k] - c8[k] * e8[k+1]);
              e8[k] = a8[k] * g0;
              f8[k] = (d8[k] + c8[k] * f8[k+1]) * g0;
            }
            //B.C. AT SURFACE 
            double pre_wk = (e8[0] * topbc[j*IMT + i] + f8[0]) * ldm_viv[idx];
            ldm_wk[idx] = pre_wk; 
            for (k = 1; k < kz; ++k) {
              pre_wk = (e8[k] * pre_wk + f8[k]) * ldm_viv[k*local_blk + idx];
              ldm_wk[k*local_blk + idx] = pre_wk;
            } 
          }
        }
      }

			athread_dma_iput_stride (&wk[start_cpe + start_blk],
          ldm_wk, size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
    } else {
  // printf("2 idx_blk %d, localblk %d, blk_per_cpe %d, size tile %d, blk %d, remainder %d, numblk %d, alloc %d\n", 
  // idx_blk,
  // local_blk, 
  // blk_per_cpe,
  // size_tile, 
  // size_blk, 
  // remainder_blk, 
  // num_blk, 
  // allocatable_num_blk);
      // No blocked
			athread_dma_iget_stride (ldm_wk, &wk[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_dcb, &dcb[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_viv, &viv[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

			// athread_dma_get_stride (ldm_wk, &wk[start_cpe + start_blk], 
			// 		remainder_blk * sizeof(double) * KM, 
			// 		remainder_blk * sizeof(double), 
			// 		(total_len - remainder_blk) * sizeof(double));
			// athread_dma_get_stride (ldm_dcb, &dcb[start_cpe + start_blk], 
			// 		remainder_blk * sizeof(double) * KM, 
			// 		remainder_blk * sizeof(double), 
			// 		(total_len - remainder_blk) * sizeof(double));
			// athread_dma_get_stride (ldm_viv, &viv[start_cpe + start_blk], 
			// 		remainder_blk * sizeof(double) * KM, 
			// 		remainder_blk * sizeof(double), 
			// 		(total_len - remainder_blk) * sizeof(double));

      // invtriu

      for (int idx = 0; idx < remainder_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        int i = index % IMT;
        int j = index / IMT;
        if (i >= 2 && i < IMT-2 && j >= 2 && j < JMT-2) {
          int kz = kmu[j*IMT + i];
          if (kz > 0 && kz <= KM) {
            for (k = 1; k < kz; ++k) {
              a8[k] = ldm_dcb[(k-1)*remainder_blk + idx] * odzt[k] 
                  * odzp[k] * c2dtc_times_aidif;
              d8[k] = ldm_wk[k*remainder_blk + idx];
            }
            for (k = 1; k < kz - 1 ; ++k) { 
              c8[k] = ldm_dcb[k*remainder_blk + idx] * odzt[k + 1] 
                  * odzp[k] * c2dtc_times_aidif;
              b8[k] = 1.0 + a8[k] + c8[k];
              e8[k] = 0.0;
              f8[k] = 0.0;
            }
            //B.C. AT TOP
            k = 0;
            a8[k] = odzp[k] * c2dtc_times_aidif;
            c8[k] = ldm_dcb[k*remainder_blk + idx] * odzt[k+1] * odzp[k] 
                * c2dtc_times_aidif;
            b8[k] = 1.0 + c8[k];
            d8[k] = ldm_wk[k*remainder_blk + idx];
            e8[k] = 0.0;
            f8[k] = 0.0;
            //B.C. AT BOTTOM
            b8[kz-1] = 1.0 + a8[kz-1];
            c8[kz-1] = odzp[kz-1] * c2dtc_times_aidif;
            e8[kz] = 0.0;
            f8[kz] = 0.0;
            d8[kz-1] = ldm_wk[(kz-1)*remainder_blk + idx]
                - bombc[j*IMT + i] * odzp[kz-1] * c2dtc_times_aidif;
            //NOW INVERT
            for (k = kz - 1; k >= 0; --k) {
              double g0 = 1.0 / (b8[k] - c8[k] * e8[k+1]);
              e8[k] = a8[k] * g0;
              f8[k] = (d8[k] + c8[k] * f8[k+1]) * g0;
            }
            //B.C. AT SURFACE 
            double pre_wk = (e8[0] * topbc[j*IMT + i] + f8[0]) * ldm_viv[idx];
            ldm_wk[idx] = pre_wk; 
            for (k = 1; k < kz; ++k) {
              pre_wk = (e8[k] * pre_wk + f8[k]) * ldm_viv[k*remainder_blk + idx];
              ldm_wk[k*remainder_blk + idx] = pre_wk;
            } 
          }
        }
      }

			athread_dma_iput_stride (&wk[start_cpe + start_blk], ldm_wk, 
          remainder_blk * sizeof(double) * KM, 
          remainder_blk * sizeof(double), 
          (total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
			// athread_dma_put_stride (&wk[start_cpe + start_blk], ldm_wk, 
      //     remainder_blk * sizeof(double) * KM, 
      //     remainder_blk * sizeof(double), 
      //     (total_len - remainder_blk) * sizeof(double));
    }
  }
  ldm_free (ldm_wk,  size_blk);
  ldm_free (ldm_dcb, size_blk);
  ldm_free (ldm_viv, size_blk);

  return ;
}

void athread_vinteg_slave (struct paramVinteg* arg_param) {
  int KM = arg_param->KM;
  int JMT = arg_param->JMT;
  int IMT = arg_param->IMT;
  double* wk3 = arg_param->wk3;
  double* wk2 = arg_param->wk2;
  double* dzp = arg_param->dzp;
  double* viv = arg_param->viv;
  double* ohbu = arg_param->ohbu;

  int total_len = JMT * IMT;
  int blk_per_cpe = (total_len + 63) / 64;
  // 2: viv, wk3
	int allocatable_num_blk = get_allocatable_size() / (2 * KM * sizeof(double));
	if (allocatable_num_blk == 0) {
		printf ("Err in athread upwell slave\n");
		exit(0);
	}
	int num_blk = (blk_per_cpe + allocatable_num_blk - 1) 
      / allocatable_num_blk;

	int local_blk = blk_per_cpe;
	if (blk_per_cpe > allocatable_num_blk) {
		local_blk = allocatable_num_blk;
	}
	int remainder_blk = blk_per_cpe % local_blk;

  size_t size_tile = local_blk * sizeof(double);
	size_t size_blk = size_tile * KM;
  size_t size_stride_layer = (total_len - local_blk) * sizeof(double);

	double* ldm_wk3 = (double*) ldm_malloc (size_blk);
	double* ldm_viv = (double*) ldm_malloc (size_blk);

  int start_cpe = athread_tid * blk_per_cpe;

	for (int idx_blk = 0; idx_blk < num_blk; ++idx_blk) {
    int start_blk = idx_blk * local_blk;
		if ((idx_blk != num_blk - 1) || (idx_blk == 0) 
        || (idx_blk == num_blk - 1 && remainder_blk == 0)) {
      // No blocked
			athread_dma_iget_stride (ldm_wk3, &wk3[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_viv, &viv[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;

      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        wk2[index] = 0.0;
      }

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        for (int k = 0; k < KM; ++k) {
          wk2[index] += dzp[k] * ohbu[index] *
              ldm_wk3[k*local_blk + idx] * ldm_viv[k*local_blk + idx];
        }
      }
    } else {
      // No blocked
			athread_dma_iget_stride (ldm_wk3, &wk3[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
          remainder_blk * sizeof(double), 
          (total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iget_stride (ldm_viv, &viv[start_cpe + start_blk], 
					remainder_blk * sizeof(double) * KM, 
          remainder_blk * sizeof(double), 
          (total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;

      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        wk2[index] = 0.0;
      }

      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        for (int k = 0; k < KM; ++k) {
          wk2[index] += dzp[k] * ohbu[index] *
              ldm_wk3[k*remainder_blk + idx] * ldm_viv[k*remainder_blk + idx];
        }
      }
    }
  }
  ldm_free (ldm_viv, size_blk);
  ldm_free (ldm_wk3, size_blk);

  return ;
}

void athread_tracer_15_slave (struct paramTracer15* arg_param) {
  int KM = arg_param->KM;
  int JMT = arg_param->JMT;
  int IMT = arg_param->IMT;

  double* adv_tt = arg_param->adv_tt;
  double* tf = arg_param->tf;
  double* vit = arg_param->vit;

  int total_len = KM * JMT * IMT;
  int blk_per_cpe = (total_len + 63) / 64;
	// 3: tf, adv_tt, vit
	int allocatable_num_blk = get_allocatable_size() / (3 * sizeof(double));
	allocatable_num_blk = (allocatable_num_blk / 8) * 8;
	if (allocatable_num_blk == 0) {
		printf ("Err in athread tracer 15 slave\n");
		exit(0);
	}

	int num_blk = (blk_per_cpe + allocatable_num_blk - 1)  / allocatable_num_blk;
 	int local_blk = blk_per_cpe;
	if (blk_per_cpe > allocatable_num_blk) {
        	local_blk = allocatable_num_blk;
	}
	int remainder_blk = blk_per_cpe % local_blk;

        size_t size_tile = local_blk * sizeof(double);
	size_t size_blk = size_tile;
        size_t size_stride_layer = (total_len - local_blk) * sizeof(double);

	double* ldm_tf  = (double*) ldm_malloc (size_blk);
	double* ldm_adv_tt   = (double*) ldm_malloc (size_blk);
	double* ldm_vit  = (double*) ldm_malloc (size_blk);

        int start_cpe = athread_tid * blk_per_cpe;

	for (int idx_blk = 0; idx_blk < num_blk; ++idx_blk) {
          int start_blk = idx_blk * local_blk;
		if ((idx_blk != num_blk - 1) || (idx_blk == 0)) {
      // No blocked
		athread_dma_iget_stride (ldm_tf, &tf[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
		athread_dma_iget_stride (ldm_adv_tt, &adv_tt[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
		athread_dma_iget_stride (ldm_vit, &vit[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

      doublev8 tf_v,adv_tt_v,vit_v;
      for (int idx = 0; idx < local_blk; idx+=8) {
        int ii = start_blk + idx + 7;
	  if (ii >= blk_per_cpe) {
             for (int j = ii - 7; j < blk_per_cpe; j++) {
               ldm_tf[j] = ldm_adv_tt[j] * ldm_vit[j];
	     }    
		  break;
	  }
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
//        int i = index % IMT;
//        int j = (index % (JMT * IMT)) / IMT;
//        if (i >= 2 && i < IMT-2 && j >= 2 && j < JMT-2) {
//            ldm_tf[idx] = ldm_adv_tt[idx] * ldm_vit[idx];
            simd_load(adv_tt_v,&(ldm_adv_tt[idx]));
            simd_load(vit_v,&(ldm_vit[idx]));
            tf_v = adv_tt_v * vit_v;
            simd_store(tf_v,&(ldm_tf[idx]));
//        }
      }
     
		athread_dma_iput_stride (&tf[start_cpe + start_blk],
          ldm_tf, size_blk, size_tile, 
          (total_len - local_blk)  * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
    } else{
		athread_dma_iget_stride (ldm_tf, &tf[start_cpe + start_blk], 
					remainder_blk * sizeof(double), 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
		athread_dma_iget_stride (ldm_adv_tt, &tf[start_cpe + start_blk], 
					remainder_blk * sizeof(double), 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
		athread_dma_iget_stride (ldm_vit, &tf[start_cpe + start_blk], 
					remainder_blk * sizeof(double), 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      doublev8 tf_v,adv_tt_v,vit_v;
      for (int idx = 0; idx < local_blk; idx+=8) {
        int ii = start_blk + idx + 7;
	  if (ii >= blk_per_cpe) {
             for (int j = ii - 7; j < blk_per_cpe; j++) {
               ldm_tf[j] = ldm_adv_tt[j] * ldm_vit[j];
	     }    
		  break;
	  }
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
//        int i = index % IMT;
//        int j = (index % (JMT * IMT)) / IMT;
//        if (i >= 2 && i < IMT-2 && j >= 2 && j < JMT-2) {
//            ldm_tf[idx] = ldm_adv_tt[idx] * ldm_vit[idx];
            simd_load(adv_tt_v,&(ldm_adv_tt[idx]));
            simd_load(vit_v,&(ldm_vit[idx]));
            tf_v = adv_tt_v * vit_v;
            simd_store(tf_v,&(ldm_tf[idx]));
//        }
      }
			athread_dma_iput_stride (&tf[start_cpe + start_blk], ldm_tf, 
                                        remainder_blk * sizeof(double) , 
					remainder_blk * sizeof(double), 
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
    }
  }
  ldm_free (ldm_tf,  size_blk);
  ldm_free (ldm_adv_tt,   size_blk);
  ldm_free (ldm_vit,  size_blk);
  
}

void athread_readyt_10_slave (struct paramReadyt10* arg_param) {

  int KM       = arg_param->KM;
  int JMT      = arg_param->JMT;
  int IMT      = arg_param->JMT;
  double G     = arg_param->G;
  double* dlu  = arg_param->dlu;
  double* dlv  = arg_param->dlv;
  double* dzp  = arg_param->dzp;
  double* zkt  = arg_param->zkt;
  double* ohbt = arg_param->ohbt;
  double* gg   = arg_param->gg;

  int total_len = JMT * IMT;
  int blk_per_cpe = (total_len + 63) / 64;
  // 1: gg
	int allocatable_num_blk = (get_allocatable_size() 
      - 4 * blk_per_cpe * sizeof(double)) 
          / (1 * KM * sizeof(double));
	if (allocatable_num_blk == 0) {
		printf ("Err in athread upwell slave\n");
		exit(0);
	}

	int num_blk = (blk_per_cpe + allocatable_num_blk - 1) 
      / allocatable_num_blk;

	int local_blk = blk_per_cpe;
	if (blk_per_cpe > allocatable_num_blk) {
		local_blk = allocatable_num_blk;
	}
	int remainder_blk = blk_per_cpe % local_blk;

  size_t size_tile = local_blk * sizeof(double);
	size_t size_blk = size_tile * KM;
  size_t size_stride_layer = (total_len - local_blk) * sizeof(double);

	double* ldm_gg   = (double*) ldm_malloc (size_blk);
	double* ldm_dlu0 = (double*) ldm_malloc (size_tile);
	double* ldm_dlu1 = (double*) ldm_malloc (size_tile);
	double* ldm_dlv0 = (double*) ldm_malloc (size_tile);
	double* ldm_dlv1 = (double*) ldm_malloc (size_tile);
  
  double abcd;
  int start_cpe = athread_tid * blk_per_cpe;
	for (int idx_blk = 0; idx_blk < num_blk; ++idx_blk) {
    int start_blk = idx_blk * local_blk;
		if ((idx_blk != num_blk - 1) || (idx_blk == 0) 
        || (idx_blk == num_blk - 1 && remainder_blk == 0)) {
      // No blocked
			athread_dma_iget_stride (ldm_gg, &gg[start_cpe + start_blk], 
					size_blk, size_tile, size_stride_layer, &dma_rply);
      D_COUNT++;
      for (int idx = 0; idx < local_blk; ++idx) {
        ldm_dlu0[idx] = 0.0;
        ldm_dlu1[idx] = 0.0;
      }
      // gg
      athread_dma_wait_value(&dma_rply, D_COUNT);
      for (int idx = 0; idx < local_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        for (int k = 0; k < KM; ++k) {
          abcd = ldm_gg[k*local_blk + idx] * ohbt[index] * dzp[k];
          ldm_dlu0[idx] += abcd;
          ldm_dlu1[idx] += abcd * zkt[k];
        }
        // for (int k = 0; k < KM; ++k) {
          // ldm_dlv0[idx] = (ldm_dlu0[idx] + ldm_dlu1[idx]
          //     * ohbt[index]) / G;
          // ldm_dlv1[idx] = ldm_dlu1[idx] * ohbt[index] * ohbt[index];
        // }
      }
			athread_dma_iput (&dlu[start_cpe + start_blk], 
          ldm_dlu0, size_tile, &dma_rply);
      D_COUNT++;
			athread_dma_iput (&dlu[total_len + start_cpe + start_blk], 
          ldm_dlu1, size_tile, &dma_rply);
      D_COUNT++;
			// athread_dma_iput (&dlv[start_cpe + start_blk], 
      //     ldm_dlv0, size_tile, &dma_rply);
      // D_COUNT++;
			// athread_dma_iput (&dlv[total_len + start_cpe + start_blk], 
      //     ldm_dlv1, size_tile, &dma_rply);
      // D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      // athread_dma_wait_value(&dma_rply, D_COUNT);
      // athread_dma_wait_value(&dma_rply, D_COUNT);
    } else {
      // No blocked
			athread_dma_iget_stride (ldm_gg, &gg[start_cpe + start_blk], 
					remainder_blk * KM * sizeof(double), remainder_blk * sizeof(double),
					(total_len - remainder_blk) * sizeof(double), &dma_rply);
      D_COUNT++;
      for (int idx = 0; idx < remainder_blk; ++idx) {
        ldm_dlu0[idx] = 0.0;
        ldm_dlu1[idx] = 0.0;
      }
      // gg
      athread_dma_wait_value(&dma_rply, D_COUNT);
      for (int idx = 0; idx < remainder_blk; ++idx) {
        int ii = start_blk + idx;
				if (ii >= blk_per_cpe) {break;}
        int index = start_cpe + ii;
        if (index >= total_len) {break;}
        for (int k = 0; k < KM; ++k) {
          abcd = ldm_gg[k*remainder_blk + idx] * ohbt[index] * dzp[k];
          ldm_dlu0[idx] += abcd;
          ldm_dlu1[idx] += abcd * zkt[k];
        }
        // for (int k = 0; k < KM; ++k) {
          ldm_dlv0[idx] = (ldm_dlu0[idx] + ldm_dlu1[idx]
              * ohbt[index]) / G;
          ldm_dlv1[idx] = ldm_dlu1[idx] * ohbt[index] * ohbt[index];
        // }
      }
			athread_dma_iput (&dlu[start_cpe + start_blk], 
          ldm_dlu0, remainder_blk * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iput (&dlu[total_len + start_cpe + start_blk], 
          ldm_dlu1, remainder_blk * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iput (&dlv[start_cpe + start_blk], 
          ldm_dlv0, remainder_blk * sizeof(double), &dma_rply);
      D_COUNT++;
			athread_dma_iput (&dlv[total_len + start_cpe + start_blk], 
          ldm_dlv1, remainder_blk * sizeof(double), &dma_rply);
      D_COUNT++;
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);
      athread_dma_wait_value(&dma_rply, D_COUNT);

    }
  }

  ldm_free (ldm_gg,   size_blk);
  ldm_free (ldm_dlu0, size_tile);
  ldm_free (ldm_dlu1, size_tile);
  ldm_free (ldm_dlv0, size_tile);
  ldm_free (ldm_dlv1, size_tile);
  return ;
}

#endif // __sw_slave__