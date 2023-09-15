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
#include "mpi_info.h"

void cpml_init(sim_t *sim)
/*< initialize PML abosorbing coefficients >*/
{
  int ib;
  float x, lx, damp, damp0;

  float freq = 10.;
  float Rc = 1e-4;
  float alpha = 3.1415926*freq;

  sim->pmla = alloc1float(sim->nb);
  sim->pmlb = alloc1float(sim->nb);

  lx = sim->nb*sim->d1;
  damp0 = -3.*sim->vmax*logf(Rc)/(2.*lx);
  for(ib=0; ib<sim->nb; ib++)   {
    x = (sim->nb-ib)*sim->d1;//should not allow x=0
    x /= lx;
    damp = damp0*x*x;
    
    sim->pmlb[ib] = exp(-(damp+alpha)*sim->dt);
    sim->pmla[ib] = damp*(sim->pmlb[ib] - 1.0)/(damp + alpha);
  }

}

void cpml_close(sim_t *sim)
{
  free1float(sim->pmla);
  free1float(sim->pmlb);
}
