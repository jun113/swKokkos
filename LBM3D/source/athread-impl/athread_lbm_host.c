#include "../include/params.h"
#include "../include/athread_lbm_params.h"

#include "athread.h"

#include <math.h>

extern void SLAVE_FUN (collide_stream_slave) (struct AthreadParams*);
extern void SLAVE_FUN (bb_left_slave) (struct AthreadParams*);
extern void SLAVE_FUN (bb_right_slave) (struct AthreadParams*);
extern void SLAVE_FUN (bb_front_slave) (struct AthreadParams*);
extern void SLAVE_FUN (bb_back_slave) (struct AthreadParams*);
extern void SLAVE_FUN (bb_bottom_slave) (struct AthreadParams*);
extern void SLAVE_FUN (bb_top_slave) (struct AthreadParams*);
extern void SLAVE_FUN (compute_macroscopic_slave) (struct AthreadParams*);
extern void SLAVE_FUN (check_if_steady_state_slave) (struct AthreadParams*);

void athread_lbm_update (double* fB, double* fA, 
		double* u, double* v, double* w, double* rho,
				struct Params *params, int step, int* converged, int D1, int D2, int D3) {
	struct AthreadParams athread_params;

	athread_params.fB = fB;
	athread_params.fA = fA;
	athread_params.u = u;
	athread_params.v = v;
	athread_params.w = w;
	athread_params.rho = rho;
	athread_params.params = params;
	athread_params.D1 = D1;
	athread_params.D2 = D2;
	athread_params.D3 = D3;
	athread_params.converged = converged;
	// collide_stream
	athread_spawn (collide_stream_slave, &athread_params);
	athread_join ();
	// bb_left
	athread_spawn (bb_left_slave, &athread_params);
	athread_join ();
	// bb_right
	athread_spawn (bb_right_slave, &athread_params);
	athread_join ();
	// bb_front
	athread_spawn (bb_front_slave, &athread_params);
	athread_join ();
	// bb_back
	athread_spawn (bb_back_slave, &athread_params);
	athread_join ();
	// bb_bottom
	athread_spawn (bb_bottom_slave, &athread_params);
	athread_join ();
	// bb_top
	athread_spawn (bb_top_slave, &athread_params);
	athread_join ();
  if ((step + 1) % params->output_rate == 0) {
		// compute_macroscopic
		athread_spawn (compute_macroscopic_slave, &athread_params);
		athread_join ();
		// check_if_steady_state
		athread_spawn (check_if_steady_state_slave, &athread_params);
		athread_join ();
	}
	return ;
}