/* 2D/3D seismic modelling, RTM and FWI code
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
 
void cpml_init(sim_t *sim)
/*< initialize PML abosorbing coefficients >*/
{
  int ib;
  float x, lx, damp, damp0;
  float Rc = 1e-5;
  float alpha = 3.1415926*sim->freq;

  sim->pmla = alloc1float(sim->nb);
  sim->pmlb = alloc1float(sim->nb);
  sim->pmla_mh = alloc1float(sim->nb);
  sim->pmlb_mh = alloc1float(sim->nb);
  sim->pmla_ph = alloc1float(sim->nb);
  sim->pmlb_ph = alloc1float(sim->nb);

  lx = sim->nb*sim->d1;
  damp0 = -3.*sim->vmax*logf(Rc)/(2.*lx);
  for(ib=0; ib<sim->nb; ib++)   {
    x = (sim->nb-ib)*sim->d1;//should not allow x=0
    x /= lx;
    damp = damp0*x*x;    
    sim->pmlb[ib] = exp(-(damp+alpha)*sim->dt);
    sim->pmla[ib] = damp*(sim->pmlb[ib] - 1.0)/(damp + alpha);    
    
    x = (sim->nb-ib-0.5)*sim->d1;//should not allow x=0, half grid shifted
    x /= lx;
    damp = damp0*x*x;    
    sim->pmlb_ph[ib] = exp(-(damp+alpha)*sim->dt);
    sim->pmla_ph[ib] = damp*(sim->pmlb_ph[ib] - 1.0)/(damp + alpha);

    x = (sim->nb-ib+0.5)*sim->d1;//should not allow x=0, half grid shifted
    x /= lx;
    damp = damp0*x*x;    
    sim->pmlb_mh[ib] = exp(-(damp+alpha)*sim->dt);
    sim->pmla_mh[ib] = damp*(sim->pmlb_mh[ib] - 1.0)/(damp + alpha);
  }

}

void cpml_free(sim_t *sim)
{
  free1float(sim->pmla);
  free1float(sim->pmlb);
  free1float(sim->pmla_ph);
  free1float(sim->pmlb_ph);
  free1float(sim->pmla_mh);
  free1float(sim->pmlb_mh);
}
