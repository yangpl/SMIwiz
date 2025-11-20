/* extend model domain with absorbing boundary layers
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *----------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"

void extend_model_init(sim_t *sim)
{
  sim->buz = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  sim->bux = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  sim->buy = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  sim->kappa = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  sim->coef1p2epsil = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  sim->sqrt1p2delta = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);

  int i1,i2,i3;
  int i1_,i2_,i3_;
  float ***vpmod, ***rhomod;
  float ***epsilmod, ***deltamod;
  
  vpmod = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  rhomod = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  epsilmod = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  deltamod = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  
  for(i3=0; i3<sim->n3; i3++){
    i3_ = (sim->n3>1)?i3+sim->nb:0;
    for(i2=0; i2<sim->n2; i2++){
      i2_ = i2+sim->nb;
      for(i1=0; i1<sim->n1; i1++){
	i1_ = i1+sim->nb;
	vpmod[i3_][i2_][i1_] = sim->vp[i3][i2][i1];
	rhomod[i3_][i2_][i1_] = sim->rho[i3][i2][i1];
	epsilmod[i3_][i2_][i1_] = sim->epsil[i3][i2][i1];
	deltamod[i3_][i2_][i1_] = sim->delta[i3][i2][i1];
      }
    }
  }
  for(i3=0; i3<sim->n3pad; i3++){
    for(i2=0; i2<sim->n2pad; i2++){
      for(i1=0; i1<sim->nb; i1++){
	vpmod[i3][i2][i1] = vpmod[i3][i2][sim->nb];
	rhomod[i3][i2][i1] = rhomod[i3][i2][sim->nb];
	epsilmod[i3][i2][i1] = epsilmod[i3][i2][sim->nb];
	deltamod[i3][i2][i1] = deltamod[i3][i2][sim->nb];

	i1_ = sim->n1pad-1-i1;
	vpmod[i3][i2][i1_] = vpmod[i3][i2][sim->n1pad-1-sim->nb];
	rhomod[i3][i2][i1_] = rhomod[i3][i2][sim->n1pad-1-sim->nb];
	epsilmod[i3][i2][i1_] = epsilmod[i3][i2][sim->n1pad-1-sim->nb];
	deltamod[i3][i2][i1_] = deltamod[i3][i2][sim->n1pad-1-sim->nb];
      }
    }
  }
  for(i3=0; i3<sim->n3pad; i3++){
    for(i2=0; i2<sim->nb; i2++){
      for(i1=0; i1<sim->n1pad; i1++){
	vpmod[i3][i2][i1] = vpmod[i3][sim->nb][i1];
	rhomod[i3][i2][i1] = rhomod[i3][sim->nb][i1];
	epsilmod[i3][i2][i1] = epsilmod[i3][sim->nb][i1];
	deltamod[i3][i2][i1] = deltamod[i3][sim->nb][i1];

	i2_ = sim->n2pad-1-i2;
	vpmod[i3][i2_][i1] = vpmod[i3][sim->n2pad-1-sim->nb][i1];
	rhomod[i3][i2_][i1] = rhomod[i3][sim->n2pad-1-sim->nb][i1];
	epsilmod[i3][i2_][i1] = epsilmod[i3][sim->n2pad-1-sim->nb][i1];
	deltamod[i3][i2_][i1] = deltamod[i3][sim->n2pad-1-sim->nb][i1];
      }
    }
  }
  if(sim->n3>1){
    for(i3=0; i3<sim->nb; i3++){
      for(i2=0; i2<sim->n2pad; i2++){
	for(i1=0; i1<sim->n1pad; i1++){
	  vpmod[i3][i2][i1] = vpmod[sim->nb][i2][i1];
	  rhomod[i3][i2][i1] = rhomod[sim->nb][i2][i1];
	  epsilmod[i3][i2][i1] = epsilmod[sim->nb][i2][i1];
	  deltamod[i3][i2][i1] = deltamod[sim->nb][i2][i1];

	  i3_ = sim->n3pad-1-i3;
	  vpmod[i3_][i2][i1] = vpmod[sim->n3pad-1-sim->nb][i2][i1];
	  rhomod[i3_][i2][i1] = rhomod[sim->n3pad-1-sim->nb][i2][i1];
	  epsilmod[i3_][i2][i1] = epsilmod[sim->n3pad-1-sim->nb][i2][i1];
	  deltamod[i3_][i2][i1] = deltamod[sim->n3pad-1-sim->nb][i2][i1];
	}
      }
    }
  }

  /* build kappa and buoyancy from rho and vp */
  for(i3=0; i3<sim->n3pad; i3++){
    for(i2=0; i2<sim->n2pad; i2++){
      for(i1=0; i1<sim->n1pad; i1++){
	sim->kappa[i3][i2][i1] = rhomod[i3][i2][i1]*vpmod[i3][i2][i1]*vpmod[i3][i2][i1];
	sim->coef1p2epsil[i3][i2][i1] = 1. + 2.*epsilmod[i3][i2][i1];
	sim->sqrt1p2delta[i3][i2][i1] = sqrt(1. + 2.*deltamod[i3][i2][i1]);
      }
    }
  }
  for(i3=0; i3<sim->n3pad; i3++){
    for(i2=0; i2<sim->n2pad; i2++){
      sim->buz[i3][i2][0] = 1./rhomod[i3][i2][0];
      for(i1=1; i1<sim->n1pad; i1++){
	sim->buz[i3][i2][i1] = 0.5*(1./rhomod[i3][i2][i1] +1./rhomod[i3][i2][i1-1]);
      }
    }
  }
  for(i3=0; i3<sim->n3pad; i3++){
    for(i1=0; i1<sim->n1pad; i1++){
      sim->bux[i3][0][i1] = 1./rhomod[i3][0][i1];
      for(i2=1; i2<sim->n2pad; i2++){
	sim->bux[i3][i2][i1] = 0.5*(1./rhomod[i3][i2][i1] +1./rhomod[i3][i2-1][i1]);
      }
    }
  }

  if(sim->n3>1){
    for(i2=0; i2<sim->n2pad; i2++){
      for(i1=0; i1<sim->n1pad; i1++){
	sim->buy[0][i2][i1] = 1./rhomod[0][i2][i1];
	for(i3=1; i3<sim->n3pad; i3++){
	  sim->buy[i3][i2][i1] = 0.5*(1./rhomod[i3][i2][i1] +1./rhomod[i3-1][i2][i1]);
	}
      }
    }
  }  
  
  free3float(vpmod);
  free3float(rhomod);
  free3float(epsilmod);
  free3float(deltamod);
}

void extend_model_free(sim_t *sim)
{
  free3float(sim->kappa);
  free3float(sim->buz);
  free3float(sim->bux);
  free3float(sim->buy);
  free3float(sim->coef1p2epsil);
  free3float(sim->sqrt1p2delta);
}


void rwi_extend_model_init(sim_t *sim, float ***vp, float ***rho, float ***kappa, float ***buz, float ***bux, float ***buy)
{
  int i1,i2,i3;
  int i1_,i2_,i3_;
  float ***vpmod, ***rhomod;
  
  vpmod = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
  rhomod = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);
   
  for(i3=0; i3<sim->n3; i3++){
    i3_ = (sim->n3>1)?i3+sim->nb:0;
    for(i2=0; i2<sim->n2; i2++){
      i2_ = i2+sim->nb;
      for(i1=0; i1<sim->n1; i1++){
	i1_ = i1+sim->nb;
	vpmod[i3_][i2_][i1_] = vp[i3][i2][i1];
	rhomod[i3_][i2_][i1_] = rho[i3][i2][i1];
      }
    }
  }
  for(i3=0; i3<sim->n3pad; i3++){
    for(i2=0; i2<sim->n2pad; i2++){
      for(i1=0; i1<sim->nb; i1++){
	vpmod[i3][i2][i1] = vpmod[i3][i2][sim->nb];
	rhomod[i3][i2][i1] = rhomod[i3][i2][sim->nb];

	i1_ = sim->n1pad-1-i1;
	vpmod[i3][i2][i1_] = vpmod[i3][i2][sim->n1pad-1-sim->nb];
	rhomod[i3][i2][i1_] = rhomod[i3][i2][sim->n1pad-1-sim->nb];
      }
    }
  }
  for(i3=0; i3<sim->n3pad; i3++){
    for(i2=0; i2<sim->nb; i2++){
      for(i1=0; i1<sim->n1pad; i1++){
	vpmod[i3][i2][i1] = vpmod[i3][sim->nb][i1];
	rhomod[i3][i2][i1] = rhomod[i3][sim->nb][i1];

	i2_ = sim->n2pad-1-i2;
	vpmod[i3][i2_][i1] = vpmod[i3][sim->n2pad-1-sim->nb][i1];
	rhomod[i3][i2_][i1] = rhomod[i3][sim->n2pad-1-sim->nb][i1];
      }
    }
  }
  if(sim->n3>1){
    for(i3=0; i3<sim->nb; i3++){
      for(i2=0; i2<sim->n2pad; i2++){
	for(i1=0; i1<sim->n1pad; i1++){
	  vpmod[i3][i2][i1] = vpmod[sim->nb][i2][i1];
	  rhomod[i3][i2][i1] = rhomod[sim->nb][i2][i1];

	  i3_ = sim->n3pad-1-i3;
	  vpmod[i3_][i2][i1] = vpmod[sim->n3pad-1-sim->nb][i2][i1];
	  rhomod[i3_][i2][i1] = rhomod[sim->n3pad-1-sim->nb][i2][i1];
	}
      }
    }
  }

  /* build kappa and buoyancy from rho and vp */
  for(i3=0; i3<sim->n3pad; i3++){
    for(i2=0; i2<sim->n2pad; i2++){
      for(i1=0; i1<sim->n1pad; i1++){
	kappa[i3][i2][i1] = rhomod[i3][i2][i1]*vpmod[i3][i2][i1]*vpmod[i3][i2][i1];
      }
    }
  }
  for(i3=0; i3<sim->n3pad; i3++){
    for(i2=0; i2<sim->n2pad; i2++){
      buz[i3][i2][0] = 1./rhomod[i3][i2][0];
      for(i1=1; i1<sim->n1pad; i1++){
	buz[i3][i2][i1] = 0.5*(1./rhomod[i3][i2][i1] +1./rhomod[i3][i2][i1-1]);
      }
    }
  }
  for(i3=0; i3<sim->n3pad; i3++){
    for(i1=0; i1<sim->n1pad; i1++){
      bux[i3][0][i1] = 1./rhomod[i3][0][i1];
      for(i2=1; i2<sim->n2pad; i2++){
	bux[i3][i2][i1] = 0.5*(1./rhomod[i3][i2][i1] +1./rhomod[i3][i2-1][i1]);
      }
    }
  }

  if(sim->n3>1){
    for(i2=0; i2<sim->n2pad; i2++){
      for(i1=0; i1<sim->n1pad; i1++){
	buy[0][i2][i1] = 1./rhomod[0][i2][i1];
	for(i3=1; i3<sim->n3pad; i3++){
	  buy[i3][i2][i1] = 0.5*(1./rhomod[i3][i2][i1] +1./rhomod[i3-1][i2][i1]);
	}
      }
    }
  }  
  
  free3float(vpmod);
  free3float(rhomod);
}
