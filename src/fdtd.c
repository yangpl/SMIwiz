/* 2D/3D seismic modelling, RTM and FWI code
 *-----------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------
 */
#include "cstd.h"
#include "sim.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void fdtd_init(sim_t *sim, int flag)
{ 
  if(flag==0){//scattering field after Born approximation
    sim->p0  = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->vz0 = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->vx0 = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->memD1p0 = alloc3float(2*sim->nb, sim->n2pad, sim->n3pad);
    sim->memD2p0 = alloc3float(sim->n1pad, 2*sim->nb, sim->n3pad);
    sim->memD1vz0 = alloc3float(2*sim->nb, sim->n2pad, sim->n3pad);
    sim->memD2vx0 = alloc3float(sim->n1pad, 2*sim->nb, sim->n3pad);
    if(sim->n3>1){
      sim->vy0 = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
      sim->memD3p0 = alloc3float(sim->n1pad, sim->n2pad, 2*sim->nb);
      sim->memD3vy0 = alloc3float(sim->n1pad, sim->n2pad, 2*sim->nb);
    }

  }else if(flag==1){//forward field
    sim->p1  = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->vz1 = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->vx1 = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->memD1p1 = alloc3float(2*sim->nb, sim->n2pad, sim->n3pad);
    sim->memD2p1 = alloc3float(sim->n1pad, 2*sim->nb, sim->n3pad);
    sim->memD1vz1 = alloc3float(2*sim->nb, sim->n2pad, sim->n3pad);
    sim->memD2vx1 = alloc3float(sim->n1pad, 2*sim->nb, sim->n3pad);
    if(sim->n3>1){
      sim->vy1 = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
      sim->memD3p1 = alloc3float(sim->n1pad, sim->n2pad, 2*sim->nb);
      sim->memD3vy1 = alloc3float(sim->n1pad, sim->n2pad, 2*sim->nb);
    }
    //backup divergence and time derivative of particle velocity
    sim->divv = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->dvxdt = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->dvydt = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->dvzdt = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    
  }else if(flag==2){//adjoint field
    sim->p2  = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->vz2 = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->vx2 = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
    sim->memD1p2 = alloc3float(2*sim->nb, sim->n2pad, sim->n3pad);
    sim->memD2p2 = alloc3float(sim->n1pad, 2*sim->nb, sim->n3pad);
    sim->memD1vz2 = alloc3float(2*sim->nb, sim->n2pad, sim->n3pad);
    sim->memD2vx2 = alloc3float(sim->n1pad, 2*sim->nb, sim->n3pad);
    if(sim->n3>1){
      sim->vy2 = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
      sim->memD3p2 = alloc3float(sim->n1pad, sim->n2pad, 2*sim->nb);
      sim->memD3vy2 = alloc3float(sim->n1pad, sim->n2pad, 2*sim->nb);
    }
  }

}

void fdtd_null(sim_t *sim, int flag)
{
  if(flag==0){//scattering field after Born approximation
    memset(sim->p0 [0][0], 0, sim->n123pad*sizeof(float));
    memset(sim->vz0[0][0], 0, sim->n123pad*sizeof(float));
    memset(sim->vx0[0][0], 0, sim->n123pad*sizeof(float));
    memset(sim->memD1p0[0][0], 0, 2*sim->nb*sim->n2pad*sim->n3pad*sizeof(float));
    memset(sim->memD2p0[0][0], 0, sim->n1pad*2*sim->nb*sim->n3pad*sizeof(float));
    memset(sim->memD1vz0[0][0], 0, 2*sim->nb*sim->n2pad*sim->n3pad*sizeof(float));
    memset(sim->memD2vx0[0][0], 0, sim->n1pad*2*sim->nb*sim->n3pad*sizeof(float));
    if(sim->n3>1){
      memset(sim->vy0[0][0], 0, sim->n123pad*sizeof(float));
      memset(sim->memD3p0[0][0], 0, sim->n1pad*sim->n2pad*2*sim->nb*sizeof(float));
      memset(sim->memD3vy0[0][0], 0, sim->n1pad*sim->n2pad*2*sim->nb*sizeof(float));
    }

  }else if(flag==1){//forward field
    memset(sim->p1 [0][0], 0, sim->n123pad*sizeof(float));
    memset(sim->vz1[0][0], 0, sim->n123pad*sizeof(float));
    memset(sim->vx1[0][0], 0, sim->n123pad*sizeof(float));
    memset(sim->memD1p1[0][0], 0, 2*sim->nb*sim->n2pad*sim->n3pad*sizeof(float));
    memset(sim->memD2p1[0][0], 0, sim->n1pad*2*sim->nb*sim->n3pad*sizeof(float));
    memset(sim->memD1vz1[0][0], 0, 2*sim->nb*sim->n2pad*sim->n3pad*sizeof(float));
    memset(sim->memD2vx1[0][0], 0, sim->n1pad*2*sim->nb*sim->n3pad*sizeof(float));
    if(sim->n3>1){
      memset(sim->vy1[0][0], 0, sim->n123pad*sizeof(float));
      memset(sim->memD3p1[0][0], 0, sim->n1pad*sim->n2pad*2*sim->nb*sizeof(float));
      memset(sim->memD3vy1[0][0], 0, sim->n1pad*sim->n2pad*2*sim->nb*sizeof(float));
    }

    memset(&sim->divv[0][0][0], 0, sim->n123pad*sizeof(float));
    memset(&sim->dvzdt[0][0][0], 0, sim->n123pad*sizeof(float));
    memset(&sim->dvxdt[0][0][0], 0, sim->n123pad*sizeof(float));
    memset(&sim->dvydt[0][0][0], 0, sim->n123pad*sizeof(float));

  }else if(flag==2){//adjoint field
    memset(sim->p2 [0][0], 0, sim->n123pad*sizeof(float));
    memset(sim->vz2[0][0], 0, sim->n123pad*sizeof(float));
    memset(sim->vx2[0][0], 0, sim->n123pad*sizeof(float));
    memset(sim->memD1p2[0][0], 0, 2*sim->nb*sim->n2pad*sim->n3pad*sizeof(float));
    memset(sim->memD2p2[0][0], 0, sim->n1pad*2*sim->nb*sim->n3pad*sizeof(float));
    memset(sim->memD1vz2[0][0], 0, 2*sim->nb*sim->n2pad*sim->n3pad*sizeof(float));
    memset(sim->memD2vx2[0][0], 0, sim->n1pad*2*sim->nb*sim->n3pad*sizeof(float));
    if(sim->n3>1){
      memset(sim->vy2[0][0], 0, sim->n123pad*sizeof(float));
      memset(sim->memD3p2[0][0], 0, sim->n1pad*sim->n2pad*2*sim->nb*sizeof(float));
      memset(sim->memD3vy2[0][0], 0, sim->n1pad*sim->n2pad*2*sim->nb*sizeof(float));
    }
  }

}

void fdtd_close(sim_t *sim, int flag)
{  
  if(flag==0){//scattering field for Born modelling
    free3float(sim->p0);
    free3float(sim->vx0);
    free3float(sim->vz0);
    free3float(sim->memD1p0);
    free3float(sim->memD2p0);
    free3float(sim->memD1vz0);
    free3float(sim->memD2vx0);
    if(sim->n3>1){
      free3float(sim->vy0);
      free3float(sim->memD3p0);
      free3float(sim->memD3vy0);
    }

  }else if(flag==1){//forward field
    free3float(sim->p1);
    free3float(sim->vx1);
    free3float(sim->vz1);
    free3float(sim->memD1p1);
    free3float(sim->memD2p1);
    free3float(sim->memD1vz1);
    free3float(sim->memD2vx1);
    if(sim->n3>1){
      free3float(sim->vy1);
      free3float(sim->memD3p1);
      free3float(sim->memD3vy1);
    }
    free3float(sim->divv);
    free3float(sim->dvxdt);
    free3float(sim->dvydt);
    free3float(sim->dvzdt);

  }else if(flag==2){//adjoint field
    free3float(sim->p2);
    free3float(sim->vx2);
    free3float(sim->vz2);
    free3float(sim->memD1p2);
    free3float(sim->memD2p2);
    free3float(sim->memD1vz2);
    free3float(sim->memD2vx2);
    if(sim->n3>1){
      free3float(sim->vy2);
      free3float(sim->memD3p2);
      free3float(sim->memD3vy2);
    }

  }
    
}


void fdtd_update_v(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy)
{
  int i1, i2, i3, j1, j2, j3, k1, k2, k3;
  float D1p, D2p, D3p;
  int i1min, i1max, i2min, i2max, i3min, i3max;
  float ***p, ***vz, ***vx, ***vy, ***memD1p, ***memD2p, ***memD3p;
  float dt = sim->sign_dt*sim->dt;
  float _d1 = 1./sim->d1;
  float _d2 = 1./sim->d2;
  float _d3 = 1./sim->d3;

  if(sim->order==4){
    i1min=1;
    i1max=sim->n1pad-3;
    i2min=1;
    i2max=sim->n2pad-3;
    i3min=(sim->n3>1)?1:0;
    i3max=(sim->n3>1)?(sim->n3pad-3):0;
  }else if(sim->order==8){
    i1min=3;
    i1max=sim->n1pad-5;
    i2min=3;
    i2max=sim->n2pad-5;
    i3min=(sim->n3>1)?3:0;
    i3max=(sim->n3>1)?(sim->n3pad-5):0;
  }

  if(sim->ibox){
    if(adj){
      i1min=MAX(sim->i1min_adj[it], i1min);
      i1max=MIN(sim->i1max_adj[it], i1max);
      i2min=MAX(sim->i2min_adj[it], i2min);
      i2max=MIN(sim->i2max_adj[it], i2max);
      i3min=(sim->n3>1)?MAX(sim->i3min_adj[it], i3min):0;
      i3max=(sim->n3>1)?MIN(sim->i3max_adj[it], i3max):0;
    }else{
      i1min=MAX(sim->i1min_fwd[it], i1min);
      i1max=MIN(sim->i1max_fwd[it], i1max);
      i2min=MAX(sim->i2min_fwd[it], i2min);
      i2max=MIN(sim->i2max_fwd[it], i2max);
      i3min=(sim->n3>1)?MAX(sim->i3min_fwd[it], i3min):0;
      i3max=(sim->n3>1)?MIN(sim->i3max_fwd[it], i3max):0;
    }
    if(sim->sign_dt<0 && flag==1){
      i1min = sim->nb;
      i1max = sim->n1+sim->nb-1;
      i2min = sim->nb;
      i2max = sim->n2+sim->nb-1;
      if(sim->n3>1){
	i3min = sim->nb;
	i3max = sim->n3pad-1-sim->nb;
      }
    }
    if(sim->freesurf) i1min = sim->nb;
  }


  if(flag==0){
    p = sim->p0;
    vz = sim->vz0;
    vx = sim->vx0;
    vy = sim->vy0;
    memD1p = sim->memD1p0;
    memD2p = sim->memD2p0;
    memD3p = sim->memD3p0;
  }else if(flag==1){
    p = sim->p1;
    vz = sim->vz1;
    vx = sim->vx1;
    vy = sim->vy1;
    memD1p = sim->memD1p1;
    memD2p = sim->memD2p1;
    memD3p = sim->memD3p1;
  }else if(flag==2){
    p = sim->p2;
    vz = sim->vz2;
    vx = sim->vx2;
    vy = sim->vy2;
    memD1p = sim->memD1p2;
    memD2p = sim->memD2p2;
    memD3p = sim->memD3p2;
  }
  
#ifdef _OPENMP
#pragma omp parallel for default(none)					\
  schedule(static)							\
  private(i1, i2, i3, j1, j2, j3, k1, k2, k3, D1p, D2p, D3p)		\
  shared(i1min, i1max, i2min, i2max, i3min, i3max, _d1, _d2, _d3, dt,	\
	 flag, p, vz, vx, vy, memD1p, memD2p, memD3p, sim,		\
	 kappa, buz, bux, buy)
#endif  
  for(i3=i3min; i3<=i3max; i3++){
    for(i2=i2min; i2<=i2max; i2++){
      for(i1=i1min; i1<=i1max; i1++){
	if(sim->order==4){
	  D1p = 1.125*(p[i3][i2][i1+1]-p[i3][i2][i1])
	    -0.041666666666666664*(p[i3][i2][i1+2]-p[i3][i2][i1-1]);
	  D2p = 1.125*(p[i3][i2+1][i1]-p[i3][i2][i1])
	    -0.041666666666666664*(p[i3][i2+2][i1]-p[i3][i2-1][i1]);
	}else if(sim->order==8){
	  D1p = 1.196289062500000*(p[i3][i2][i1+1]-p[i3][i2][i1])
	    -0.079752604166667*(p[i3][i2][i1+2]-p[i3][i2][i1-1])
	    +0.009570312500000*(p[i3][i2][i1+3]-p[i3][i2][i1-2])
	    -0.000697544642857*(p[i3][i2][i1+4]-p[i3][i2][i1-3]);
	  D2p = 1.196289062500000*(p[i3][i2+1][i1]-p[i3][i2][i1])
	    -0.079752604166667*(p[i3][i2+2][i1]-p[i3][i2-1][i1])
	    +0.009570312500000*(p[i3][i2+3][i1]-p[i3][i2-2][i1])
	    -0.000697544642857*(p[i3][i2+4][i1]-p[i3][i2-3][i1]);
	}
	D1p *= _d1;
	D2p *= _d2;
      
	if(i1<sim->nb) {
	  memD1p[i3][i2][i1] = sim->pmlb[i1]*memD1p[i3][i2][i1] + sim->pmla[i1]*D1p;
	  D1p += memD1p[i3][i2][i1];
	}else if(i1>=sim->n1pad-sim->nb){
	  j1 = sim->n1pad -1 - i1;
	  k1 = j1 + sim->nb;
	  memD1p[i3][i2][k1] = sim->pmlb[j1]*memD1p[i3][i2][k1] + sim->pmla[j1]*D1p;
	  D1p += memD1p[i3][i2][k1];
	}
	if(i2<sim->nb) {
	  memD2p[i3][i2][i1] = sim->pmlb[i2]*memD2p[i3][i2][i1] + sim->pmla[i2]*D2p;
	  D2p += memD2p[i3][i2][i1];
	}else if(i2>=sim->n2pad-sim->nb){
	  j2 = sim->n2pad-1-i2;
	  k2 = j2 + sim->nb;
	  memD2p[i3][k2][i1] = sim->pmlb[j2]*memD2p[i3][k2][i1] + sim->pmla[j2]*D2p;
	  D2p += memD2p[i3][k2][i1];
	}
      
	vz[i3][i2][i1] -= dt*buz[i3][i2][i1]*D1p;
	vx[i3][i2][i1] -= dt*bux[i3][i2][i1]*D2p;
	if(flag==1){
	  sim->dvzdt[i3][i2][i1] = -buz[i3][i2][i1]*D1p;
	  sim->dvxdt[i3][i2][i1] = -bux[i3][i2][i1]*D2p;
	}
	
	if(sim->n3>1){
	  if(sim->order==4){
	    D3p = 1.125*(p[i3+1][i2][i1]-p[i3][i2][i1])
	      -0.041666666666666664*(p[i3+2][i2][i1]-p[i3-1][i2][i1]);
	  }else if(sim->order==8){
	    D3p = 1.196289062500000*(p[i3+1][i2][i1]-p[i3][i2][i1])
	      -0.079752604166667*(p[i3+2][i2][i1]-p[i3-1][i2][i1])
	      +0.009570312500000*(p[i3+3][i2][i1]-p[i3-2][i2][i1])
	      -0.000697544642857*(p[i3+4][i2][i1]-p[i3-3][i2][i1]);
	  }
	  D3p *= _d3;
	  if(i3<sim->nb){
	    memD3p[i3][i2][i1] = sim->pmlb[i3]*memD3p[i3][i2][i1] + sim->pmla[i3]*D3p;
	    D3p += memD3p[i3][i2][i1];
	  }else if(i3>=sim->n3pad-sim->nb){
	    j3 = sim->n3pad-1-i3;
	    k3 = j3 + sim->nb;
	    memD3p[k3][i2][i1] = sim->pmlb[j3]*memD3p[k3][i2][i1] + sim->pmla[j3]*D3p;
	    D3p += memD3p[k3][i2][i1];
	  }
	  
	  vy[i3][i2][i1] -= dt*buy[i3][i2][i1]*D3p;
	  if(flag==1) sim->dvydt[i3][i2][i1] = -buy[i3][i2][i1]*D3p;
	}
	
      }//end for i3
    }//end for i2
  }//end for i1
  
  if(sim->freesurf){
    if(sim->order==4){
#ifdef _OPENMP
#pragma omp parallel for default(none)		\
  schedule(static)				\
  private(i2, i3)				\
  shared(i2min, i2max, i3min, i3max, vz, sim)
#endif  
      for(i3=i3min; i3<=i3max; i3++){
	for(i2=i2min; i2<=i2max; i2++){
	  vz[i3][i2][sim->nb-1] = vz[i3][i2][sim->nb];
	  vz[i3][i2][sim->nb-2] = vz[i3][i2][sim->nb+1];
	}
      }

    }else if(sim->order==8){
#ifdef _OPENMP
#pragma omp parallel for default(none)		\
  schedule(static)				\
  private(i2, i3)				\
  shared(i2min, i2max, i3min, i3max, vz, sim)
#endif  
      for(i3=i3min; i3<=i3max; i3++){
	for(i2=i2min; i2<=i2max; i2++){
	  vz[i3][i2][sim->nb-1] = vz[i3][i2][sim->nb];
	  vz[i3][i2][sim->nb-2] = vz[i3][i2][sim->nb+1];
	  vz[i3][i2][sim->nb-3] = vz[i3][i2][sim->nb+2];
	  vz[i3][i2][sim->nb-4] = vz[i3][i2][sim->nb+3];
	}
      }

    }//end if order
  }//end if freesurf
  
}


void fdtd_update_p(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy)
{
  int i1, i2, i3, j1, j2, j3, k1, k2, k3;
  float D1vz, D2vx, D3vy, divv;
  int i1min, i1max, i2min, i2max, i3min, i3max;
  float ***p, ***vz, ***vx, ***vy, ***memD1vz, ***memD2vx, ***memD3vy;
  float dt = sim->sign_dt*sim->dt;
  float _d1 = 1./sim->d1;
  float _d2 = 1./sim->d2;
  float _d3 = 1./sim->d3;

  if(sim->order==4){
    i1min=2;
    i1max=sim->n1pad-2;
    i2min=2;
    i2max=sim->n2pad-2;
    i3min=(sim->n3>1)?2:0;
    i3max=(sim->n3>1)?(sim->n3pad-2):0;
  }else if(sim->order==8){
    i1min=4;
    i1max=sim->n1pad-4;
    i2min=4;
    i2max=sim->n2pad-4;
    i3min=(sim->n3>1)?4:0;
    i3max=(sim->n3>1)?(sim->n3pad-4):0;
  }
  
  if(sim->ibox){
    if(adj){
      i1min=MAX(sim->i1min_adj[it], i1min);
      i1max=MIN(sim->i1max_adj[it], i1max);
      i2min=MAX(sim->i2min_adj[it], i2min);
      i2max=MIN(sim->i2max_adj[it], i2max);
      i3min=(sim->n3>1)?MAX(sim->i3min_adj[it], i3min):0;
      i3max=(sim->n3>1)?MIN(sim->i3max_adj[it], i3max):0;
    }else{
      i1min=MAX(sim->i1min_fwd[it], i1min);
      i1max=MIN(sim->i1max_fwd[it], i1max);
      i2min=MAX(sim->i2min_fwd[it], i2min);
      i2max=MIN(sim->i2max_fwd[it], i2max);
      i3min=(sim->n3>1)?MAX(sim->i3min_fwd[it], i3min):0;
      i3max=(sim->n3>1)?MIN(sim->i3max_fwd[it], i3max):0;
    }
    if(sim->sign_dt<0 && flag==1){
      i1min = sim->nb;
      i1max = sim->n1pad-1-sim->nb;
      i2min = sim->nb;
      i2max = sim->n2pad-1-sim->nb;
      if(sim->n3>1){
	i3min = sim->nb;
	i3max = sim->n3pad-1-sim->nb;
      }
    }
    if(sim->freesurf) i1min = sim->nb;
  }
  
  if(flag==0){
    p = sim->p0;
    vz = sim->vz0;
    vx = sim->vx0;
    vy = sim->vy0;
    memD1vz = sim->memD1vz0;
    memD2vx = sim->memD2vx0;
    memD3vy = sim->memD3vy0;
  }else if(flag==1){
    p = sim->p1;
    vz = sim->vz1;
    vx = sim->vx1;
    vy = sim->vy1;
    memD1vz = sim->memD1vz1;
    memD2vx = sim->memD2vx1;
    memD3vy = sim->memD3vy1;
  }else if(flag==2){
    p = sim->p2;
    vz = sim->vz2;
    vx = sim->vx2;
    vy = sim->vy2;
    memD1vz = sim->memD1vz2;
    memD2vx = sim->memD2vx2;
    memD3vy = sim->memD3vy2;
  }

#ifdef _OPENMP
#pragma omp parallel for default(none)					\
  schedule(static)							\
  private(i1, i2, i3, j1, j2, j3, k1, k2, k3, D1vz, D2vx, D3vy, divv)	\
  shared(i1min, i1max, i2min, i2max, i3min, i3max, _d1, _d2, _d3, dt,	\
	 flag, p, vz, vx, vy, memD1vz, memD2vx, memD3vy, sim,		\
	 kappa, buz, bux, buy)
#endif  
  for(i3=i3min; i3<=i3max; i3++) {
    for(i2=i2min; i2<=i2max; i2++) {
      for(i1=i1min; i1<=i1max; i1++) {
	if(sim->order==4){
	  D1vz = 1.125*(vz[i3][i2][i1]-vz[i3][i2][i1-1])
	    -0.041666666666666664*(vz[i3][i2][i1+1]-vz[i3][i2][i1-2]);
	  D2vx = 1.125*(vx[i3][i2][i1]-vx[i3][i2-1][i1])
	    -0.041666666666666664*(vx[i3][i2+1][i1]-vx[i3][i2-2][i1]);
	}else if(sim->order==8){
	  D1vz = 1.196289062500000*(vz[i3][i2][i1]-vz[i3][i2][i1-1])
	    -0.079752604166667*(vz[i3][i2][i1+1]-vz[i3][i2][i1-2])
	    +0.009570312500000*(vz[i3][i2][i1+2]-vz[i3][i2][i1-3])
	    -0.000697544642857*(vz[i3][i2][i1+3]-vz[i3][i2][i1-4]);
	  D2vx = 1.196289062500000*(vx[i3][i2][i1]-vx[i3][i2-1][i1])
	    -0.079752604166667*(vx[i3][i2+1][i1]-vx[i3][i2-2][i1])
	    +0.009570312500000*(vx[i3][i2+2][i1]-vx[i3][i2-3][i1])
	    -0.000697544642857*(vx[i3][i2+3][i1]-vx[i3][i2-4][i1]);
	}
	D1vz *= _d1;
	D2vx *= _d2;
      
	if(i1<sim->nb) {
	  memD1vz[i3][i2][i1] = sim->pmlb[i1]*memD1vz[i3][i2][i1] + sim->pmla[i1]*D1vz;
	  D1vz += memD1vz[i3][i2][i1];
	}else if(i1>=sim->n1pad-sim->nb){
	  j1 = sim->n1pad -1 - i1;
	  k1 = j1 + sim->nb;
	  memD1vz[i3][i2][k1] = sim->pmlb[j1]*memD1vz[i3][i2][k1] + sim->pmla[j1]*D1vz;
	  D1vz += memD1vz[i3][i2][k1];
	}
	if(i2<sim->nb) {
	  memD2vx[i3][i2][i1] = sim->pmlb[i2]*memD2vx[i3][i2][i1] + sim->pmla[i2]*D2vx;
	  D2vx += memD2vx[i3][i2][i1];
	}else if(i2>=sim->n2pad-sim->nb){
	  j2 = sim->n2pad-1-i2;
	  k2 = j2 + sim->nb;
	  memD2vx[i3][k2][i1] = sim->pmlb[j2]*memD2vx[i3][k2][i1] + sim->pmla[j2]*D2vx;
	  D2vx += memD2vx[i3][k2][i1];
	}

	if(sim->n3>1){
	  if(sim->order==4){
	    D3vy = 1.125*(vy[i3][i2][i1]-vy[i3-1][i2][i1])
	      -0.041666666666666664*(vy[i3+1][i2][i1]-vy[i3-2][i2][i1]);	  
	  }else if(sim->order==8){
	    D3vy = 1.196289062500000*(vy[i3][i2][i1]-vy[i3-1][i2][i1])
	      -0.079752604166667*(vy[i3+1][i2][i1]-vy[i3-2][i2][i1])
	      +0.009570312500000*(vy[i3+2][i2][i1]-vy[i3-3][i2][i1])
	      -0.000697544642857*(vy[i3+3][i2][i1]-vy[i3-4][i2][i1]);
	  }
	  D3vy *= _d3;
      
	  if(i3<sim->nb) {
	    memD3vy[i3][i2][i1] = sim->pmlb[i3]*memD3vy[i3][i2][i1] + sim->pmla[i3]*D3vy;
	    D3vy += memD3vy[i3][i2][i1];
	  }else if(i3>=sim->n3pad-sim->nb){
	    j3 = sim->n3pad-1-i3;
	    k3 = j3 + sim->nb;
	    memD3vy[k3][i2][i1] = sim->pmlb[j3]*memD3vy[k3][i2][i1] + sim->pmla[j3]*D3vy;
	    D3vy += memD3vy[k3][i2][i1];
	  }
	}else
	  D3vy = 0.;

	divv = (D1vz + D2vx + D3vy);
	p[i3][i2][i1] -= dt*kappa[i3][i2][i1]*divv;
	if(flag==1) sim->divv[i3][i2][i1] = divv;
      }
    }
  }

  if(sim->freesurf){
    if(sim->order==4){
#ifdef _OPENMP
#pragma omp parallel for default(none)		\
  schedule(static)				\
  private(i2, i3)				\
  shared(i2min, i2max, i3min, i3max, p, sim)
#endif  
      for(i3=i3min; i3<=i3max; i3++){
	for(i2=i2min; i2<=i2max; i2++){
	  p[i3][i2][sim->nb] = 0;
	  p[i3][i2][sim->nb-1] = -p[i3][i2][sim->nb+1];
	}
      }
    }else if(sim->order==8){
#ifdef _OPENMP
#pragma omp parallel for default(none)		\
  schedule(static)				\
  private(i2, i3)				\
  shared(i2min, i2max, i3min, i3max, p, sim)
#endif  
      for(i3=i3min; i3<=i3max; i3++){
	for(i2=i2min; i2<=i2max; i2++){
	  p[i3][i2][sim->nb] = 0;
	  p[i3][i2][sim->nb-1] = -p[i3][i2][sim->nb+1];
	  p[i3][i2][sim->nb-2] = -p[i3][i2][sim->nb+2];
	  p[i3][i2][sim->nb-3] = -p[i3][i2][sim->nb+3];
	}
      }

    }//end if order
  }//end if freesurf

}
