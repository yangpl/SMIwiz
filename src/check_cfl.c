/* 2D/3D seismic modelling, RTM and FWI code
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
 

void check_cfl(sim_t *sim)
{
  float tmp;// freqmax, lambda_min, ppw1, ppw2, ppw3;
  int i1, i2, i3;

  sim->vmin = sim->vp[0][0][0];
  sim->vmax = sim->vp[0][0][0];
  sim->rhomin = sim->rho[0][0][0];
  sim->rhomax = sim->rho[0][0][0];
  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      for(i1=0; i1<sim->n1; i1++){
	sim->vmin = MIN(sim->vmin, sim->vp[i3][i2][i1]);
	sim->vmax = MAX(sim->vmax, sim->vp[i3][i2][i1]);
	sim->rhomin = MIN(sim->rhomin, sim->rho[i3][i2][i1]);
	sim->rhomax = MAX(sim->rhomax, sim->rho[i3][i2][i1]);
      }
    }
  }
  if(iproc==0) {
    printf("--------- check CFL ----------------\n");
    printf("[rhomin, rhomax]=[%g, %g]\n", sim->rhomin, sim->rhomax);
    printf("[vmin, vmax]=[%g, %g]\n", sim->vmin, sim->vmax);
  }

  /* CFL = dt*vmax* \sum_i |c_i|sqrt(1/dx^2 + 1/dy^2 + 1/dz^2), where c_i are finite difference coefficients */
  if(sim->n3>1)
    tmp = 1./(sim->d1*sim->d1) + 1./(sim->d2*sim->d2) + 1./(sim->d3*sim->d3);
  else
    tmp = 1./(sim->d1*sim->d1) + 1./(sim->d2*sim->d2);
  sim->cfl = sim->dt*sim->vmax*sqrt(tmp);
  if(sim->order==4) tmp = 1.125 + 0.041666666666666664;
  else tmp = (1.196289062500000 + 0.079752604166667 + 0.009570312500000 + 0.000697544642857);
  sim->cfl *= tmp;
  if(iproc==0) printf("cfl=%g\n", sim->cfl);
  if(sim->cfl>=1.0) err("CFL stability condition unsatisifed!");

  //printf("[ppw1, ppw2, ppw3]=[%g, %g, %g] (points per wavelength)\n", ppw1, ppw2, ppw3);

  
  /* freqmax = sim->fm;//estimate maximum frequency=2*fpeak frequency */
  /* lambda_min = sim->vmin/freqmax;//minimum wavelength */
  /* ppw1 = lambda_min/sim->d1; */
  /* ppw2 = lambda_min/sim->d2; */
  /* ppw3 = (sim->n3>1)?lambda_min/sim->d3:ppw2; */
    
  /* if(sim->order==4){ */
  /*   if(ppw1<5.) err("ppw1<5"); */
  /*   if(ppw2<5.) err("ppw2<5"); */
  /*   if(ppw3<5.) err("ppw3<5"); */
  /* }else{//if(sim->order==8){ */
  /*   if(ppw1<4.) err("ppw1<4"); */
  /*   if(ppw2<4.) err("ppw2<4"); */
  /*   if(ppw3<4.) err("ppw3<4"); */
  /* } */

}
