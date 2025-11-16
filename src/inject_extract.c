/* inject sources, extract seismic data from receivers
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------*/
#include "sim.h"
#include "acq.h"

//inject point source based on delta distribution
void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it)
{
  int isrc, i1, i2, i3, i1_, i2_, i3_;
  float w1, w2, w3, s;

  float dt = sim->sign_dt*sim->dt;
  s = stf_it*dt/sim->volume;
  for(isrc=0; isrc<acq->nsrc; isrc++){
    for(i3=-sim->ri; i3<=sim->ri; i3++){
      if(sim->n3>1) {
	w3 = acq->src_w3[isrc][i3+sim->ri];
	i3_ = acq->src_i3[isrc] + i3;
      }else{
	w3 = 1.;
	i3_ = 0;
      }
      for(i2=-sim->ri; i2<=sim->ri; i2++){
	w2 = acq->src_w2[isrc][i2+sim->ri];
	i2_ = acq->src_i2[isrc] + i2;
	for(i1=-sim->ri; i1<=sim->ri; i1++){
	  w1 = acq->src_w1[isrc][i1+sim->ri];
	  i1_ = acq->src_i1[isrc] + i1;
	  sp[i3_][i2_][i1_] += s*w1*w2*w3;

          //for points above free surface, subtract response at mirror location
	  if(sim->freesurf && acq->src_nm[isrc]>0){
	    w1 = acq->src_w1m[isrc][i1+sim->ri];
	    i1_ = acq->src_i1m[isrc] + i1;
	    sp[i3_][i2_][i1_] -= s*w1*w2*w3;
	  }
	}
      }
    }
    
  }//end for isrc
}

void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it)
/*< extract data from wavefield using Kaiser windowed sinc interpolation >*/
{   
  int irec, i1, i2, i3, i1_, i2_, i3_;
  float w1, w2, w3, s;

  for(irec=0; irec<acq->nrec; irec++) {
    s = 0;
    for(i3=-sim->ri; i3<=sim->ri; i3++){
      if(sim->n3>1) {
	w3 = acq->rec_w3[irec][i3+sim->ri];
	i3_ = acq->rec_i3[irec] + i3;
      }else{
	w3 = 1.;
	i3_ = 0;
      }
      for(i2=-sim->ri; i2<=sim->ri; i2++){
	w2 = acq->rec_w2[irec][i2+sim->ri];
	i2_ = acq->rec_i2[irec] + i2;
	for(i1=-sim->ri; i1<=sim->ri; i1++){
	  w1 = acq->rec_w1[irec][i1+sim->ri];
	  i1_ = acq->rec_i1[irec] + i1;
	  s += sp[i3_][i2_][i1_]*w1*w2*w3;

	  if(sim->freesurf && acq->rec_nm[irec]>0){//subtract response at mirror location
	    w1 = acq->rec_w1m[irec][i1+sim->ri];
	    i1_ = acq->rec_i1m[irec] + i1;
	    s -= sp[i3_][i2_][i1_]*w1*w2*w3;
	  }
	}
      }
    }
    dat[irec][it] = s;
  }//end for irec
}


void inject_adjoint_source(sim_t *sim, acq_t *acq, float ***rp, float **dres, int it)
/*< inject adjoint source using Kaiser windowed sinc interpolation >*/
{   
  int irec, i1, i2, i3, i1_, i2_, i3_;
  float w1, w2, w3, s;

  float dt = -sim->sign_dt*sim->dt;
  for(irec=0; irec<acq->nrec; irec++){
    for(i3=-sim->ri; i3<=sim->ri; i3++){
      if(sim->n3>1) {
	w3 = acq->rec_w3[irec][i3+sim->ri];
	i3_ = acq->rec_i3[irec] + i3;
      }else{
	w3 = 1.;
	i3_ = 0;
      }
      for(i2=-sim->ri; i2<=sim->ri; i2++){
	w2 = acq->rec_w2[irec][i2+sim->ri];
	i2_ = acq->rec_i2[irec] + i2;
	for(i1=-sim->ri; i1<=sim->ri; i1++){
	  w1 = acq->rec_w1[irec][i1+sim->ri];
	  i1_ = acq->rec_i1[irec] + i1;
	  s = sim->kappa[i3_][i2_][i1_]*dres[irec][it]*dt/sim->volume;
	  rp[i3_][i2_][i1_] += s*w1*w2*w3;

	  if(sim->freesurf && acq->rec_nm[irec]>0){//subtract response at mirror location
	    w1 = acq->rec_w1m[irec][i1+sim->ri];
	    i1_ = acq->rec_i1m[irec] + i1;
	    s = sim->kappa[i3_][i2_][i1_]*dres[irec][it]*dt/sim->volume;
	    rp[i3_][i2_][i1_] -= s*w1*w2*w3;
	  }
	}
      }
    }
  }//end for irec

}

