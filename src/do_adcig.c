/* Extracting angle domain common image gather (ADCIG) from RTM
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *----------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
#include "acq.h"
#include "fwi.h"
#include <mpi.h>

void read_data(sim_t *sim, acq_t *acq);
void write_data(sim_t *sim, acq_t *acq);
void setup_data_weight(acq_t *acq, sim_t *sim);

void check_cfl(sim_t *sim);

void fdtd_init(sim_t *sim, int flag);
void fdtd_null(sim_t *sim, int flag);
void fdtd_free(sim_t *sim, int flag);
void fdtd_update_v(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy);
void fdtd_update_p(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy);

void extend_model_init(sim_t *sim);
void extend_model(sim_t *sim, float ***vp, float ***rho, float ***kappa, float ***buz, float ***bux, float ***buy);
void extend_model_free(sim_t *sim);

void computing_box_init(acq_t *acq, sim_t *sim, int adj);
void computing_box_free(sim_t *sim, int adj);

void cpml_init(sim_t *sim);
void cpml_free(sim_t *sim);

void decimate_interp_init(sim_t *sim, int flag);
void decimate_interp_free(sim_t *sim, int flag);
void decimate_interp_bndr(sim_t *sim, int flag, int it, int interp, float **face1, float **face2, float **face3);

void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it);
void inject_adjoint_source(sim_t *sim, acq_t *acq, float ***rp, float **dres, int it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);

void do_adcig(sim_t *sim, acq_t *acq)
{

}
