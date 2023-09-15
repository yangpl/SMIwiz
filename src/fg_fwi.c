/*----------------------------------------------------------------------
 Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology
 Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes
 Homepage: https://yangpl.wordpress.com
 E-mail: ypl.2100@gmail.com
 -------------------------------------------------------------------*/
#include <mpi.h>
#include "cstd.h"
#include "sim.h"
#include "acq.h"
#include "fwi.h"
#include "mpi_info.h"

sim_t *sim;
acq_t *acq;
fwi_t *fwi;


void check_cfl(sim_t *sim);

void fdtd_init(sim_t *sim, int flag);
void fdtd_null(sim_t *sim, int flag);
void fdtd_close(sim_t *sim, int flag);
void fdtd_update_v(sim_t *sim, int flag, int it, int adj);
void fdtd_update_p(sim_t *sim, int flag, int it, int adj);

void decimate_interp_init(sim_t *sim);
void decimate_interp_close(sim_t *sim);
void decimate_interp_bndr(sim_t *sim, int interp, int it);

void check_cfl(sim_t *sim);

void extend_model_init(sim_t *sim);
void extend_model(sim_t *sim);
void extend_model_close(sim_t *sim);

void computing_box_init(acq_t *acq, sim_t *sim, int adj);
void computing_box_close(sim_t *sim, int adj);

void cpml_init(sim_t *sim);
void cpml_close(sim_t *sim);

void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);
void inject_adjoint_source(sim_t *sim, acq_t *acq, float ***rp, float **dres, int it);

void read_data(sim_t *sim, acq_t *acq);
void write_data(sim_t *sim, acq_t *acq);
void setup_data_weight(acq_t *acq, sim_t *sim);

float regularization_tikhonov(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3);
float regularization_tv(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3);
void triangle_smoothing(float ***mod, int n1, int n2, int n3, int r1, int r2, int r3, int repeat);



float ***h1, ***h2;

/*--------------------------------------------------------------*/
void fg_fwi_init(sim_t *sim_, acq_t *acq_, fwi_t *fwi_)
{
  sim = sim_;
  acq = acq_;
  fwi = fwi_;
  
  sim->dobs = alloc2float(sim->nt,acq->nrec);
  sim->dcal = alloc2float(sim->nt,acq->nrec);
  sim->dres = alloc2float(sim->nt,acq->nrec);
  memset(sim->dobs[0], 0, sim->nt*acq->nrec*sizeof(float));
  memset(sim->dcal[0], 0, sim->nt*acq->nrec*sizeof(float));
  memset(sim->dres[0], 0, sim->nt*acq->nrec*sizeof(float));

  read_data(sim, acq);//read observed data
  setup_data_weight(acq, sim);

  if(!getparint("check", &sim->check)) sim->check = 0;
  if(!getparint("itcheck", &sim->itcheck)) sim->itcheck = sim->nt/2;
  
  fwi->iter = 0;
  fwi->alpha = 1;
  fwi->firstgrad = 1;
  if(fwi->preco==2){
    if(iproc==0) printf("pseudo-Hessian activated\n");
    fwi->hess = alloc1float(fwi->n);//pseudo-Hessian preconditioning
    h1 = alloc3float(sim->n1, sim->n2, sim->n3);
    h2 = alloc3float(sim->n1, sim->n2, sim->n3);
  }
  
  if(!getparint("objopt", &fwi->objopt)) fwi->objopt = 0;//0=L2; 1=AWI;
  if(!getparfloat("gamma1", &fwi->gamma1)) fwi->gamma1 = 0;//Tikhonov reg
  if(!getparfloat("gamma2", &fwi->gamma2)) fwi->gamma2 = 0;//TV reg
  
  if(!getparint("r1", &fwi->r1)) fwi->r1 = 1;
  if(!getparint("r2", &fwi->r2)) fwi->r2 = 1;
  if(!getparint("r3", &fwi->r3)) fwi->r3 = 1;
  if(!getparint("repeat", &fwi->repeat)) fwi->repeat = 3;

  if(iproc==0) printf("-----------------fwi init done----------------------\n"); 
}


/*--------------------------------------------------------------*/
void fg_fwi_close(sim_t *sim)
{
  free2float(sim->dobs);
  free2float(sim->dcal);
  free2float(sim->dres);
  if(fwi->preco==2) free1float(fwi->hess);
  if(fwi->preco==2){
    free3float(h1);
    free3float(h2);
  }
}

/*--------------------------------------------------------------*/
void fg_mod_reg(sim_t *sim, fwi_t *fwi, float *x, float *g)
/*< model regularization (assume fwi gradient has been computed and stored in grad) >*/
{
  int i1,i2,i3,i, j, k, ipar;
  
  float tmp = pow(0.85,fwi->iter);
  float gamma1 =  fwi->gamma1*tmp;
  float gamma2 =  fwi->gamma2*tmp;
  float fcost_mod = 0;
  float ***grad;

  grad = alloc3float(sim->n1, sim->n2, sim->n3);
  for(ipar=0; ipar<fwi->npar; ipar++){
    if(fwi->gamma1>0){//Tikhonov regularizaiton
      fcost_mod = regularization_tikhonov(&x[ipar*sim->n123], &grad[0][0][0], sim->n1, sim->n2, sim->n3, sim->d1, sim->d2, sim->d3);
      fwi->fcost += fcost_mod*gamma1;
      for(i3 = 0; i3<sim->n3; i3++){
	for(i2 = 0; i2<sim->n2; i2++){
	  for(i1 = 0; i1<sim->n1; i1++){
	    i = i1 + sim->n1*(i2 + sim->n2*i3);
	    g[i + sim->n123*ipar] += gamma1*grad[i3][i2][i1];
	  }
	}
      }      
    }//end if
    
    if(fwi->gamma2>0){//TV regularization
      fcost_mod = regularization_tv(&x[ipar*sim->n123], &grad[0][0][0], sim->n1, sim->n2, sim->n3, sim->d1, sim->d2, sim->d3);
      fwi->fcost += fcost_mod*gamma2; 
      for(i3 = 0; i3<sim->n3; i3++){
	for(i2 = 0; i2<sim->n2; i2++){
	  for(i1 = 0; i1<sim->n1; i1++){
	    i = i1 + sim->n1*(i2 + sim->n2*i3);
	    g[i + sim->n123*ipar] += gamma2*grad[i3][i2][i1];
	  }
	}
      }      
    }//end if

    /* mirror around bathymetry before smoothing, keep the average around bathymetry */
    for(i3 = 0; i3<sim->n3; i3++){
      for(i2 = 0; i2<sim->n2; i2++){
	for(j = 0,i1 = fwi->ibathy[i3][i2]+fwi->itransition-1; i1>= 0; i1--,j++){
	  i = i1 + sim->n1*(i2 + sim->n2*i3);
	  k = fwi->ibathy[i3][i2]+fwi->itransition+j +sim->n1*(i2 + sim->n2*i3);
	  g[i+sim->n123*ipar] = g[k+sim->n123*ipar];
	}
      }
    }
    memcpy(&grad[0][0][0], &g[ipar*sim->n123], sim->n123*sizeof(float));
    triangle_smoothing(grad, sim->n1, sim->n2, sim->n3, fwi->r1, fwi->r2, fwi->r3, fwi->repeat);
    for(i3 = 0; i3<sim->n3; i3++){
      for(i2 = 0; i2<sim->n2; i2++){
	for(i1 = 0; i1<sim->n1; i1++){
	  i = i1 + sim->n1*(i2 + sim->n2*i3);
	  g[i + sim->n123*ipar] = (i1<= fwi->ibathy[i3][i2])?0:grad[i3][i2][i1];
	}
      }
    }
  }//end for ipar
  
  free3float(grad);
}

/*--------------------------------------------------------------*/
float fg_fwi(float *x, float *g)
/*< misfit function and gradient evaluation of FWI >*/
{
  int j;
  int it, irec, i1, i2, i3, i1_, i2_, i3_, ipar;
  float s1, s2, tmp;
  float ***g1, ***g2;
  char fname[sizeof("dsyn_0000")];
  FILE *fp;

  g1 = alloc3float(sim->n1, sim->n2, sim->n3);
  g2 = alloc3float(sim->n1, sim->n2, sim->n3);
  memset(&g1[0][0][0], 0, sim->n123*sizeof(float));
  memset(&g2[0][0][0], 0, sim->n123*sizeof(float));
  if(fwi->preco==2){//allocate memory for pseudo-Hessian
    memset(&h1[0][0][0], 0, sim->n123*sizeof(float));
    memset(&h2[0][0][0], 0, sim->n123*sizeof(float));
  }
  
  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(i1=0; i1<sim->n1; i1++){
	  j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));

	  tmp = exp(x[j]);
	  if(fwi->family==1){//vp-rho
	    if(fwi->idxpar[ipar]==1) sim->vp[i3][i2][i1] = tmp;
	    if(fwi->idxpar[ipar]==2) sim->rho[i3][i2][i1] = tmp;
	  }
	  if(fwi->family==2){//vp-ip
	    if(fwi->idxpar[ipar]==1) sim->vp[i3][i2][i1] = tmp;
	    if(fwi->idxpar[ipar]==2) sim->rho[i3][i2][i1] = tmp/sim->vp[i3][i2][i1];
	  }
	}
      }
    }
  }

  check_cfl(sim);
  cpml_init(sim);  
  extend_model_init(sim);
  fdtd_init(sim, 1);//flag=1, incident field
  fdtd_init(sim, 2);//flag=2, adjoint field
  fdtd_null(sim, 1);//flag=1, incident field
  fdtd_null(sim, 2);//flag=2, adjoint field
  decimate_interp_init(sim);
  extend_model(sim);
  computing_box_init(acq, sim, 0);
  computing_box_init(acq, sim, 1);
  
  /*--------------------------------------------------------------*/
  if(iproc==0) printf("----stage 1: forward modelling!--------\n");
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    decimate_interp_bndr(sim, 0, it);/* interp=0 */
    fdtd_update_v(sim, 1, it, 0);
    fdtd_update_p(sim, 1, it, 0);
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);
    
    if(iproc==0 && sim->check==1 && it==sim->itcheck){
      fp = fopen("wave1.bin", "wb");
      for(i3=0; i3<sim->n3; i3++){
	i3_ = (sim->n3>1)?i3 + sim->nb:0;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2 + sim->nb;
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1 + sim->nb;
	    
	    fwrite(&sim->p1[i3_][i2_][i1_], sizeof(float), 1, fp);
	  }
	}
      }
      fclose(fp);
    }
  }
  
  /*--------------------------------------------------------------*/
  fwi->fcost = 0;
  for(irec=0; irec<acq->nrec; irec++){
    for(it=0; it<sim->nt; it++){
      sim->dres[irec][it] = (sim->dobs[irec][it]-sim->dcal[irec][it])*acq->wdat[irec][it];
      fwi->fcost += sim->dres[irec][it]*sim->dres[irec][it];
      sim->dres[irec][it] *= acq->wdat[irec][it];      
    }
  }
  MPI_Allreduce(&fwi->fcost, &tmp, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  fwi->fcost = 0.5*tmp*sim->dt;
  if(iproc==0) printf("fcost=%g\n", fwi->fcost);
  
  sprintf(fname, "dsyn_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"error opening file\n"); exit(1);}
  fwrite(&sim->dcal[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);
  sprintf(fname, "dres_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"error opening file\n"); exit(1);}
  fwrite(&sim->dres[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);
  
  /*--------------------------------------------------------------*/
  if(iproc==0) printf("----stage 2: adjoint modelling!-----------\n");
  sim->sign_dt = -1;
  for(it=sim->nt-1; it>=0; it--){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);
  
    inject_adjoint_source(sim, acq, sim->p2, sim->dres, it);
    fdtd_update_v(sim, 2, it, 1);
    fdtd_update_p(sim, 2, it, 1);

    inject_source(sim, acq, sim->p1, sim->stf[it]);
    fdtd_update_p(sim, 1, it, 0);
    fdtd_update_v(sim, 1, it, 0);
    decimate_interp_bndr(sim, 1, it);/* interp=1 */

    if(iproc==0 && sim->check==1 && it==sim->itcheck){
      fp = fopen("wave2.bin", "wb");
      for(i3=0; i3<sim->n3; i3++){
	i3_ = (sim->n3>1)?i3 + sim->nb:0;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2 + sim->nb;
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1 + sim->nb;
	    
	    fwrite(&sim->p1[i3_][i2_][i1_], sizeof(float), 1, fp);
	  }
	}
      }
      fclose(fp);
    }
    
    for(i3=0; i3<sim->n3; i3++){
      i3_ = (sim->n3>1)?i3+sim->nb:0;
      for(i2=0; i2<sim->n2; i2++){
	i2_ = i2+sim->nb;
	for(i1=0; i1<sim->n1; i1++){
	  i1_ = i1+sim->nb;
	  g1[i3][i2][i1] += sim->p2[i3_][i2_][i1_]*sim->divv[i3_][i2_][i1_];
	  g2[i3][i2][i1] += sim->vz2[i3_][i2_][i1_]*sim->dvzdt[i3_][i2_][i1_] + sim->vz2[i3_][i2_][i1_-1]*sim->dvzdt[i3_][i2_][i1_-1];
	  g2[i3][i2][i1] += sim->vx2[i3_][i2_][i1_]*sim->dvxdt[i3_][i2_][i1_] + sim->vx2[i3_][i2_-1][i1_]*sim->dvxdt[i3_][i2_-1][i1_];
	  if(sim->n3>1) g2[i3][i2][i1] += sim->vz2[i3_][i2_][i1_]*sim->dvzdt[i3_][i2_][i1_] + sim->vz2[i3_-1][i2_][i1_]*sim->dvzdt[i3_-1][i2_][i1_];
	}
      }
    }
    if(fwi->preco==2){
      for(i3=0; i3<sim->n3; i3++){
	i3_ = (sim->n3>1)?i3+sim->nb:0;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2+sim->nb;
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1+sim->nb;
	    h1[i3][i2][i1] += sim->divv[i3_][i2_][i1_]*sim->divv[i3_][i2_][i1_];
	    h2[i3][i2][i1] += sim->dvzdt[i3_][i2_][i1_]*sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1]*sim->dvzdt[i3_][i2_][i1_-1];
	    h2[i3][i2][i1] += sim->dvxdt[i3_][i2_][i1_]*sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_]*sim->dvxdt[i3_][i2_-1][i1_];
	    if(sim->n3>1) h2[i3][i2][i1] += sim->dvzdt[i3_][i2_][i1_]*sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_-1][i2_][i1_]*sim->dvzdt[i3_-1][i2_][i1_];
	  }
	}
      }
    }//end if    

  }  
  cpml_close(sim);
  extend_model_close(sim);
  fdtd_close(sim, 1);
  fdtd_close(sim, 2);
  decimate_interp_close(sim);
  computing_box_close(sim, 0);
  computing_box_close(sim, 1);

  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      for(i1=0; i1<sim->n1; i1++){
	tmp = sim->volume*sim->dt;
	g1[i3][i2][i1] *= tmp;//dJ/dln(kappa)
	g2[i3][i2][i1] *= sim->rho[i3][i2][i1]*0.5*tmp;//dJ/dln(rho)
	if(i1<= fwi->ibathy[i3][i2]){//reset again to avoid leakage
	  g1[i3][i2][i1] = 0.;
	  g2[i3][i2][i1] = 0.;
	}
      }
    }
  }
  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(i1=0; i1<sim->n1; i1++){
	  j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	  if(fwi->family==1){//vp-rho
	    if(fwi->idxpar[ipar]==1) g[j] = 2.0*g1[i3][i2][i1];//dJ/dln(vp)
	    if(fwi->idxpar[ipar]==2) g[j] = g1[i3][i2][i1] + g2[i3][i2][i1];//dJ/dln(rho)
	  }
	  if(fwi->family==2){//vp-ip
	    if(fwi->idxpar[ipar]==1) g[j] = g1[i3][i2][i1] - g2[i3][i2][i1];//dJ/dln(vp)
	    if(fwi->idxpar[ipar]==2) g[j] = g1[i3][i2][i1] + g2[i3][i2][i1];//dJ/dln(ip)
	  }
	}
      }
    }
  }
  if(fwi->preco==2){
    MPI_Allreduce(&h1[0][0][0], fwi->hess, sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&h1[0][0][0], fwi->hess, sim->n123*sizeof(float));
    MPI_Allreduce(&h2[0][0][0], fwi->hess, sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&h2[0][0][0], fwi->hess, sim->n123*sizeof(float));
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(i1=0; i1<sim->n1; i1++){
	  tmp = sim->volume*sim->dt;
	  s1 = sim->rho[i3][i2][i1]*sim->vp[i3][i2][i1]*sim->vp[i3][i2][i1];
	  s2 = sim->rho[i3][i2][i1];
	  h1[i3][i2][i1] *= s1*s1*tmp;//H_{ln(kappa),ln(kappa)}
	  h2[i3][i2][i1] *= s2*s2*0.5*tmp;//H_{ln(rho),ln(rho)}
	}
      }
    }
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	    if(fwi->family==1){//vp-rho
	      if(fwi->idxpar[ipar]==1) fwi->hess[j] = 4.0*h1[i3][i2][i1];//H_{ln(vp),ln(vp)}
	      if(fwi->idxpar[ipar]==2) fwi->hess[j] = h1[i3][i2][i1] + h2[i3][i2][i1];//H_{ln(rho),ln(rho)}
	    }
	    if(fwi->family==2){//vp-ip
	      if(fwi->idxpar[ipar]==1) fwi->hess[j] = h1[i3][i2][i1] + h2[i3][i2][i1];//H_{ln(vp),ln(vp)}
	      if(fwi->idxpar[ipar]==2) fwi->hess[j] = h1[i3][i2][i1] + h2[i3][i2][i1];//H_{ln(ip),ln(ip)}
	    }
	  }
	}
      }
    }

  }//end if
  

  for(ipar=0; ipar<fwi->npar; ipar++){
    MPI_Allreduce(&g[ipar*sim->n123], &g1[0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&g[ipar*sim->n123], &g1[0][0][0], sim->n123*sizeof(float));
  }
  free3float(g1);
  free3float(g2);
  fg_mod_reg(sim, fwi, x, g);
  if(iproc==0 && fwi->firstgrad){
    fp=fopen("gradient_fwi","wb");
    fwrite(g, fwi->n*sizeof(float), 1, fp);
    fclose(fp);
  }
  
  if(sim->mode==1 && fwi->firstgrad){
    s1 = fabs(x[0]);
    s2 = fabs(g[0]);
    for(j=0; j<fwi->n; j++){
      s1 = MAX(s1, fabs(x[j]));
      s2 = MAX(s2, fabs(g[j]));
    }
    fwi->alpha = 0.02*s1/s2;
    /*
    s1 = 0;
    s2 = 0;
    for(j=0; j<fwi->n; j++){
      s1 += fabs(x[j]);
      s2 += fabs(g[j]);
    }
    fwi->alpha = 5e-4*s1/s2;
    */
    fwi->firstgrad = 0;
    if(iproc==0) printf("scaling=%e\n", fwi->alpha);
  }
  if(!fwi->firstgrad){
    for(j=0; j<fwi->n; j++) g[j] *= fwi->alpha;
  }
  fwi->fcost *= fwi->alpha;
  if(iproc==0) printf("scaled fcost=%g\n", fwi->fcost);
  
  return fwi->fcost;

}

