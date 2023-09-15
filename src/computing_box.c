/* 2D/3D seismic modelling, RTM and FWI code
 *-----------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
#include "acq.h"

void computing_box_init(acq_t *acq, sim_t *sim, int adj)
/*< specify computing box to do simulation efficiently >*/
{
  int irec,isrc,it,nxe,nze,nye,shift;
  float t;

  shift = sim->ri + 2;//1 due to NINT + 2 due to 4th order FD
  if(adj){
    sim->i1min_adj=alloc1int(sim->nt);
    sim->i1max_adj=alloc1int(sim->nt);
    sim->i2min_adj=alloc1int(sim->nt);
    sim->i2max_adj=alloc1int(sim->nt);
    sim->i3min_adj=alloc1int(sim->nt);
    sim->i3max_adj=alloc1int(sim->nt);
    sim->i1min_adj[sim->nt-1]=sim->i1max_adj[sim->nt-1]=acq->rec_i1[0];
    sim->i2min_adj[sim->nt-1]=sim->i2max_adj[sim->nt-1]=acq->rec_i2[0];
    if(sim->n3>1) 
      sim->i3min_adj[sim->nt-1]=sim->i3max_adj[sim->nt-1]=acq->rec_i3[0];
    else
      sim->i3min_adj[sim->nt-1]=sim->i3max_adj[sim->nt-1]=0;
    for(irec=0; irec<acq->nrec; irec++){
      sim->i1min_adj[sim->nt-1]=MIN(sim->i1min_adj[sim->nt-1],acq->rec_i1[irec]);
      sim->i1max_adj[sim->nt-1]=MAX(sim->i1max_adj[sim->nt-1],acq->rec_i1[irec]);
      sim->i2min_adj[sim->nt-1]=MIN(sim->i2min_adj[sim->nt-1],acq->rec_i2[irec]);
      sim->i2max_adj[sim->nt-1]=MAX(sim->i2max_adj[sim->nt-1],acq->rec_i2[irec]);
      if(sim->n3>1){
	sim->i3min_adj[sim->nt-1]=MIN(sim->i3min_adj[sim->nt-1],acq->rec_i3[irec]);
	sim->i3max_adj[sim->nt-1]=MAX(sim->i3max_adj[sim->nt-1],acq->rec_i3[irec]);
      }
    }
    sim->i1min_adj[sim->nt-1] -= shift;
    sim->i1max_adj[sim->nt-1] += shift;
    sim->i2min_adj[sim->nt-1] -= shift;
    sim->i2max_adj[sim->nt-1] += shift;
    if(sim->n3>1){
      sim->i3min_adj[sim->nt-1] -= shift;
      sim->i3max_adj[sim->nt-1] += shift;
    }
	
    for(it=sim->nt-2; it>=0; it--){
      t=(sim->nt-1-it)*sim->dt;
      nze=NINT(sim->vmax*t/sim->d1)+shift;
      nxe=NINT(sim->vmax*t/sim->d2)+shift;
      sim->i1min_adj[it]=MAX(sim->i1min_adj[sim->nt-1]-nze,0);
      sim->i1max_adj[it]=MIN(sim->i1max_adj[sim->nt-1]+nze,sim->n1pad-1);
      sim->i2min_adj[it]=MAX(sim->i2min_adj[sim->nt-1]-nxe,0);
      sim->i2max_adj[it]=MIN(sim->i2max_adj[sim->nt-1]+nxe,sim->n2pad-1);
      if(sim->n3>1){
	nye=NINT(sim->vmax*t/sim->d3)+shift;
	sim->i3min_adj[it]=MAX(sim->i3min_adj[sim->nt-1]-nye,0);
	sim->i3max_adj[it]=MIN(sim->i3max_adj[sim->nt-1]+nye,sim->n3pad-1);
      }else
	sim->i3min_adj[it]=sim->i3max_adj[it]=0;
    }
  }else{
    sim->i1min_fwd=alloc1int(sim->nt);
    sim->i1max_fwd=alloc1int(sim->nt);
    sim->i2min_fwd=alloc1int(sim->nt);
    sim->i2max_fwd=alloc1int(sim->nt);
    sim->i3min_fwd=alloc1int(sim->nt);
    sim->i3max_fwd=alloc1int(sim->nt);
    sim->i1min_fwd[0]=acq->src_i1[0];
    sim->i1max_fwd[0]=acq->src_i1[0];
    sim->i2min_fwd[0]=acq->src_i2[0];
    sim->i2max_fwd[0]=acq->src_i2[0];
    if(sim->n3>1) {
      sim->i3min_fwd[0]=acq->src_i3[0];
      sim->i3max_fwd[0]=acq->src_i3[0];
    }else{
      sim->i3min_fwd[0]=0;
      sim->i3max_fwd[0]=0;	
    }
    for(isrc=0; isrc<acq->nsrc; isrc++){
      sim->i1min_fwd[0]=MIN(sim->i1min_fwd[0],acq->src_i1[isrc]);
      sim->i1max_fwd[0]=MAX(sim->i1max_fwd[0],acq->src_i1[isrc]);
      sim->i2min_fwd[0]=MIN(sim->i2min_fwd[0],acq->src_i2[isrc]);
      sim->i2max_fwd[0]=MAX(sim->i2max_fwd[0],acq->src_i2[isrc]);
      if(sim->n3>1){
	sim->i3min_fwd[0]=MIN(sim->i3min_fwd[0],acq->src_i3[isrc]);
	sim->i3max_fwd[0]=MAX(sim->i3max_fwd[0],acq->src_i3[isrc]);
      }
    }
    sim->i1min_fwd[sim->nt-1] -= shift;
    sim->i1max_fwd[sim->nt-1] += shift;
    sim->i2min_fwd[sim->nt-1] -= shift;
    sim->i2max_fwd[sim->nt-1] += shift;
    if(sim->n3>1){
      sim->i3min_fwd[sim->nt-1] -= shift;
      sim->i3max_fwd[sim->nt-1] += shift;
    }

    for(it=1; it<sim->nt; it++){
      t=it*sim->dt;
      nze=NINT(sim->vmax*t/sim->d1)+shift;
      nxe=NINT(sim->vmax*t/sim->d2)+shift;
      sim->i1min_fwd[it]=MAX(sim->i1min_fwd[0]-nze,0);
      sim->i1max_fwd[it]=MIN(sim->i1max_fwd[0]+nze,sim->n1pad-1);
      sim->i2min_fwd[it]=MAX(sim->i2min_fwd[0]-nxe,0);
      sim->i2max_fwd[it]=MIN(sim->i2max_fwd[0]+nxe,sim->n2pad-1);
      
      if(sim->n3>1){
	nye=NINT(sim->vmax*t/sim->d3)+shift;
	sim->i3min_fwd[it]=MAX(sim->i3min_fwd[0]-nye,0);
	sim->i3max_fwd[it]=MIN(sim->i3max_fwd[0]+nye,sim->n3pad-1);
      }else
	sim->i3min_fwd[it]=sim->i3max_fwd[it]=0;
    }
  }
}

void computing_box_close(sim_t *sim, int adj)
{
  if(adj){
    free1int(sim->i1min_adj);
    free1int(sim->i1max_adj);
    free1int(sim->i2min_adj);
    free1int(sim->i2max_adj);
    free1int(sim->i3min_adj);
    free1int(sim->i3max_adj);
  }else{
    free1int(sim->i1min_fwd);
    free1int(sim->i1max_fwd);
    free1int(sim->i2min_fwd);
    free1int(sim->i2max_fwd);
    free1int(sim->i3min_fwd);
    free1int(sim->i3max_fwd);
  }
}
