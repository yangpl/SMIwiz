#include "cstd.h"
#include "sim.h"
#include "fwi.h"

void precondition(sim_t *sim, fwi_t *fwi, float *x)
/*< depth precondition to compensate geometrical spreading >*/
{
  int i1, i2, i3, ipar, k;
  float z, s1, s2;

  if(fwi->preco==2){//pseudo-Hessian precondition
    s1 = 0.;
    s2 = 0.;
    for(ipar=0; ipar<fwi->npar; ipar++){
      z = 0;
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    k = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	    z = MAX(z, fwi->hess[k]);//find the maximum value of pseudo-Hessian
	  }
	}
      }
      
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    k = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));

	    s1 += x[k]*x[k];
	    x[k] /= (fwi->hess[k] + 1e-5*z);//regularize to avoid division by 0
	    s2 += x[k]*x[k];
	  }
	}
      }
    }//end for ipar

    z = sqrt(s1/s2);
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    k = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	    x[k] *= z;
	  }
	}
      }
    }//end for ipar

  }else{//depth preconditioning
  
    s1 = 0.;
    s2 = 0.;
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    z = i1*sim->d1;
	    if(sim->n3>1) z = z*z;

	    k = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	    s1 += x[k]*x[k];
	    x[k] *= z;
	    s2 += x[k]*x[k];
	  }
	}
      }
    }//end for ipar

    z = sqrt(s1/s2);
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    k = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	    x[k] *= z;
	  }
	}
      }
    }//end for ipar

  }//end if
  
}

