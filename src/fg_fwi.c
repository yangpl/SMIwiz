/* function and gradient evaluation for full waveform inversion (FWI) 
 * and reflection waveform inversion (RWI)
 *--------------------------------------------------------------------
 * Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology
 * Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes
 * Homepage: https://yangpl.wordpress.com
 * E-mail: ypl.2100@gmail.com
 *-------------------------------------------------------------------*/
#include <mpi.h>
#include "cstd.h"
#include "sim.h"
#include "acq.h"
#include "fwi.h"
 

sim_t *sim;
acq_t *acq;
fwi_t *fwi;
float ***h1, ***h2;

void check_cfl(sim_t *sim);

void fdtd_init(sim_t *sim, int flag);
void fdtd_null(sim_t *sim, int flag);
void fdtd_free(sim_t *sim, int flag);
void fdtd_update_v(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy);
void fdtd_update_p(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy);

void decimate_interp_init(sim_t *sim, int flag);
void decimate_interp_free(sim_t *sim, int flag);
void decimate_interp_bndr(sim_t *sim, int flag, int it, int interp, float **face1, float **face2, float **face3);

void check_cfl(sim_t *sim);

void extend_model_init(sim_t *sim);
void extend_model(sim_t *sim, float ***vp, float ***rho, float ***kappa, float ***buz, float ***bux, float ***buy);
void extend_model_free(sim_t *sim);

void computing_box_init(acq_t *acq, sim_t *sim, int adj);
void computing_box_free(sim_t *sim, int adj);

void cpml_init(sim_t *sim);
void cpml_free(sim_t *sim);

void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);
void inject_adjoint_source(sim_t *sim, acq_t *acq, float ***rp, float **dres, int it);

float awi_adjoint_source(acq_t *acqui, sim_t *sim, fwi_t *fwi);

void read_data(sim_t *sim, acq_t *acq);
void write_data(sim_t *sim, acq_t *acq);
void setup_data_weight(acq_t *acq, sim_t *sim);

float regularization_tikhonov(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3);
float regularization_tv(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3);
void triangle_smoothing(float ***mod, int n1, int n2, int n3, int r1, int r2, int r3, int repeat);

/*--------------------------------------------------------------*/
void fg_fwi_init(sim_t *sim_, acq_t *acq_, fwi_t *fwi_)
{
  sim = sim_;
  acq = acq_;
  fwi = fwi_;
  
  read_data(sim, acq);//read observed data
  setup_data_weight(acq, sim);

  if(!getparint("itcheck", &sim->itcheck)) sim->itcheck = sim->nt/2;
  
  fwi->iter = 0;
  fwi->alpha = 1;
  fwi->firstgrad = 1;
  if(fwi->rwi){
    sim->ip = alloc3float(sim->n1, sim->n2, sim->n3);
    sim->dm = alloc3float(sim->n1, sim->n2, sim->n3);

    char *dmfile;
    if(!getparstring("dmfile", &dmfile)) err("must give dmfile= ");

    for(int i3=0; i3<sim->n3; i3++){
      for(int i2=0; i2<sim->n2; i2++){
	for(int i1=0; i1<sim->n1; i1++){
	  sim->ip[i3][i2][i1] = sim->vp[i3][i2][i1]*sim->rho[i3][i2][i1];
	}
      }
    }
    
    FILE *fp = fopen(dmfile, "rb");
    fread(&sim->dm[0][0][0], sim->n123*sizeof(float), 1, fp);
    fclose(fp);
  }
  if(fwi->preco==2){
    if(iproc==0) printf("pseudo-Hessian activated\n");
    fwi->hess = alloc1float(fwi->n);//pseudo-Hessian preconditioning
    h1 = alloc3float(sim->n1, sim->n2, sim->n3);
    h2 = alloc3float(sim->n1, sim->n2, sim->n3);
  }

  if(!getparint("objopt", &fwi->objopt)) fwi->objopt = 0;//0=L2-FWI; 1=AWI
  if(!getparfloat("gamma1", &fwi->gamma1)) fwi->gamma1 = 0;//Tikhonov reg
  if(!getparfloat("gamma2", &fwi->gamma2)) fwi->gamma2 = 0;//TV reg  
  if(!getparint("r1", &fwi->r1)) fwi->r1 = 1;
  if(!getparint("r2", &fwi->r2)) fwi->r2 = 1;
  if(!getparint("r3", &fwi->r3)) fwi->r3 = 1;
  if(!getparint("repeat", &fwi->repeat)) fwi->repeat = 3;

  if(iproc==0) printf("-----------------fwi init done----------------------\n"); 
}


/*--------------------------------------------------------------*/
void fg_fwi_free(sim_t *sim)
{

  if(fwi->rwi){
    free3float(sim->ip);
    free3float(sim->dm);
  }
  if(fwi->preco==2){
    free1float(fwi->hess);
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

float fg_rwi(float *x, float *g);

/*--------------------------------------------------------------*/
float fg_fwi(float *x, float *g)
/*< misfit function and gradient evaluation of FWI >*/
{
  int j;
  int it, irec, i1, i2, i3, i1_, i2_, i3_, ipar;
  //int in1, in2, in3;
  float s1, s2, tmp;
  float ***g1, ***g2;
  char fname[sizeof("dsyn_0000")];
  FILE *fp;

  if(fwi->rwi) fwi->fcost = fg_rwi(x, g);
  else{
    g1 = alloc3float(sim->n1, sim->n2, sim->n3);
    g2 = alloc3float(sim->n1, sim->n2, sim->n3);
    memset(&g1[0][0][0], 0, sim->n123*sizeof(float));
    memset(&g2[0][0][0], 0, sim->n123*sizeof(float));
    if(fwi->preco==2){//allocate memory for pseudo-Hessian
      memset(&h1[0][0][0], 0, sim->n123*sizeof(float));
      memset(&h2[0][0][0], 0, sim->n123*sizeof(float));
    }
  
    for(ipar=0; ipar<fwi->npar; ipar++){
      if(fwi->idxpar[ipar]==1){//1st parameter - vp
	for(i3=0; i3<sim->n3; i3++){
	  for(i2=0; i2<sim->n2; i2++){
	    for(i1=0; i1<sim->n1; i1++){
	      j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	      tmp = exp(x[j]);
	      if(fwi->family==1) sim->vp[i3][i2][i1] = tmp;//vp-rho
	      if(fwi->family==2) sim->vp[i3][i2][i1] = tmp;//vp-ip
	    }
	  }
	}
      }//end if
    }//end for ipar
    for(ipar=0; ipar<fwi->npar; ipar++){
      if(fwi->idxpar[ipar]==2){//2nd parameter - rho/ip
	for(i3=0; i3<sim->n3; i3++){
	  for(i2=0; i2<sim->n2; i2++){
	    for(i1=0; i1<sim->n1; i1++){
	      j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	      tmp = exp(x[j]);
	      if(fwi->family==1) sim->rho[i3][i2][i1] = tmp;//vp-rho
	      if(fwi->family==2) sim->rho[i3][i2][i1] = tmp/sim->vp[i3][i2][i1];//vp-ip
	    }
	  }
	}
      }//end if
    }//end for ipar
    
    check_cfl(sim);
    cpml_init(sim);  
    extend_model_init(sim);
    fdtd_init(sim, 1);//flag=1, incident field
    fdtd_init(sim, 2);//flag=2, adjoint field
    fdtd_null(sim, 1);//flag=1, incident field
    fdtd_null(sim, 2);//flag=2, adjoint field
    decimate_interp_init(sim, 1);
    extend_model(sim, sim->vp, sim->rho, sim->kappa, sim->buz, sim->bux, sim->buy);
    computing_box_init(acq, sim, 0);
    computing_box_init(acq, sim, 1);
  
    /*--------------------------------------------------------------*/
    if(iproc==0) printf("----stage 1: forward modelling!--------\n");
    sim->sign_dt = 1;
    for(it=0; it<sim->nt; it++){
      if(iproc==0 && it%100==0) printf("it-----%d\n", it);

      decimate_interp_bndr(sim, 1, it, 0, sim->face1, sim->face2, sim->face3);/* interp=0 */
      fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
      fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
      inject_source(sim, acq, sim->p1, sim->stf[it]);
      extract_wavefield(sim, acq, sim->p1, sim->dcal, it);
    
      if(iproc==0 && it==sim->itcheck){
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
    if(fwi->objopt == 1)//AWI by Michael Warner
      fwi->fcost = awi_adjoint_source(acq, sim, fwi);
    else{//classic FWI using l2 norm
      for(irec=0; irec<acq->nrec; irec++){
	for(it=0; it<sim->nt; it++){
	  sim->dres[irec][it] = (sim->dobs[irec][it]-sim->dcal[irec][it])*acq->wdat[irec][it];
	  fwi->fcost += sim->dres[irec][it]*sim->dres[irec][it];
	  sim->dres[irec][it] *= acq->wdat[irec][it];      
	}
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
      fdtd_update_v(sim, 2, it, 1, sim->kappa, sim->buz, sim->bux, sim->buy);
      fdtd_update_p(sim, 2, it, 1, sim->kappa, sim->buz, sim->bux, sim->buy);

      inject_source(sim, acq, sim->p1, sim->stf[it]);
      fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
      fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
      decimate_interp_bndr(sim, 1, it, 1, sim->face1, sim->face2, sim->face3);/* interp=1 */

      if(iproc==0 && it==sim->itcheck){
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
	//in3 = (sim->n3>1)?(i3>=sim->order/2 && i3<sim->n3-sim->order):1;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2+sim->nb;
	  //in2 =  (i2>=sim->order/2 && i2<sim->n2-sim->order);
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1+sim->nb;
	    //in1 = (i1>=sim->order/2 && i1<sim->n1-sim->order);
	    //in1 = in1 && (i1>fwi->ibathy[i3][i2]);

	    if(i1>fwi->ibathy[i3][i2]){//reset again to avoid leakage
	      g1[i3][i2][i1] += sim->p2[i3_][i2_][i1_]*sim->divv[i3_][i2_][i1_];
	      g2[i3][i2][i1] += (sim->vz2[i3_][i2_][i1_] + sim->vz2[i3_][i2_][i1_-1])*(sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1]);
	      g2[i3][i2][i1] += (sim->vx2[i3_][i2_][i1_] + sim->vx2[i3_][i2_-1][i1_])*(sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_]);
	      if(sim->n3>1) g2[i3][i2][i1] += (sim->vz2[i3_][i2_][i1_] + sim->vz2[i3_-1][i2_][i1_])*(sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_-1][i2_][i1_]);
	    }//end if
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
	      h2[i3][i2][i1] += (sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1])*(sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1]);
	      h2[i3][i2][i1] += (sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_])*(sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_]);
	      if(sim->n3>1) h2[i3][i2][i1] += (sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_-1][i2_][i1_])*(sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_-1][i2_][i1_]);
	    }
	  }
	}
      }//end if    

    }  
    cpml_free(sim);
    extend_model_free(sim);
    fdtd_free(sim, 1);
    fdtd_free(sim, 2);
    decimate_interp_free(sim, 1);
    computing_box_free(sim, 0);
    computing_box_free(sim, 1);

    tmp = sim->volume*sim->dt;
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(i1=0; i1<sim->n1; i1++){
	  g1[i3][i2][i1] *= tmp;//dJ/dln(kappa)
	  g2[i3][i2][i1] *= sim->rho[i3][i2][i1]*0.25*tmp;//dJ/dln(rho)
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
	    h2[i3][i2][i1] *= s2*s2*0.25*tmp;//H_{ln(rho),ln(rho)}
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
      fwi->alpha = 0.01*s1/s2;
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
  }

  return fwi->fcost;
}

//===================================================================
float fg_rwi(float *x, float *g)
/*< misfit function and gradient evaluation of RWI >*/
{
  int j;
  int it, irec, i1, i2, i3, i1_, i2_, i3_, ipar;
  int in1, in2, in3;
  float s1, s2, tmp;
  float ***rho_, ***kappa_, ***buz_, ***bux_, ***buy_;
  float **d0, **_dres;
  char fname[sizeof("dsyn_0000")];
  FILE *fp;

  memset(g, 0, fwi->n*sizeof(float));

  rho_ = alloc3float(sim->n1, sim->n2, sim->n3);//m+dm
  kappa_ = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);//extended m+dm
  buz_ = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);//extended m+dm
  bux_ = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);//extended m+dm
  buy_ = alloc3float(sim->n1pad, sim->n2pad, sim->n3pad);//extended m+dm
  d0 = alloc2float(sim->nt, acq->nrec);//Ru0
  _dres = alloc2float(sim->nt, acq->nrec);//-dres=-R^H(delta_d-Rdu)

  for(ipar=0; ipar<fwi->npar; ipar++){
    if(fwi->idxpar[ipar]==1){//ln(vp)
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	    sim->vp[i3][i2][i1] = exp(x[j]);//velocity via m0=ln(vp), ip remains the same
	    sim->rho[i3][i2][i1] = sim->ip[i3][i2][i1]/sim->vp[i3][i2][i1];
	  }
	}
      }
    }//end if
  }//end for ipar
  for(ipar=0; ipar<fwi->npar; ipar++){
    if(fwi->idxpar[ipar]==2) {//dln(ip)
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	    tmp = exp(log(sim->ip[i3][i2][i1]) + x[j]);//exp(ln(ip)+dln(ip))=perturbed impedance
	    rho_[i3][i2][i1] = tmp/sim->vp[i3][i2][i1];//density under new ip, velocity remains the same
	  }
	}
      }
    }//end if
  }//end for ipar
  
  check_cfl(sim);
  cpml_init(sim);  
  extend_model_init(sim);
  fdtd_init(sim, 0);//flag=0, incident field
  fdtd_init(sim, 1);//flag=1, incident field
  fdtd_init(sim, 2);//flag=2, adjoint field
  fdtd_init(sim, 3);//flag=3, adjoint field
  fdtd_null(sim, 0);//flag=0, incident field
  fdtd_null(sim, 1);//flag=1, incident field
  fdtd_null(sim, 2);//flag=2, adjoint field
  fdtd_null(sim, 3);//flag=3, adjoint field
  decimate_interp_init(sim, 0);
  decimate_interp_init(sim, 1);
  computing_box_init(acq, sim, 0);
  computing_box_init(acq, sim, 1);
  extend_model(sim, sim->vp, rho_, kappa_, buz_, bux_, buy_);
  extend_model(sim, sim->vp, sim->rho, sim->kappa, sim->buz, sim->bux, sim->buy);
  
  /*--------------------------------------------------------------*/
  if(iproc==0) printf("----stage 1: forward modelling --------\n");
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    //A(m+dm)(u+du)=f
    decimate_interp_bndr(sim, 0, it, 0, sim->face1_, sim->face2_, sim->face3_);/* interp=0 */
    fdtd_update_v(sim, 0, it, 0, kappa_, buz_, bux_, buy_);
    fdtd_update_p(sim, 0, it, 0, kappa_, buz_, bux_, buy_);
    inject_source(sim, acq, sim->p0, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p0, sim->dcal, it);//dcal=R(u+du)

    //A(m)u=f
    decimate_interp_bndr(sim, 1, it, 0, sim->face1, sim->face2, sim->face3);/* interp=0 */
    fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p1, d0, it);//d=Ru
  }//end for it

  /*--------------------------------------------------------------*/
  fwi->fcost = 0;
  for(irec=0; irec<acq->nrec; irec++){
    for(it=0; it<sim->nt; it++){
      sim->dcal[irec][it] -= d0[irec][it];//Rdu=R(u+du) - Ru
      sim->dres[irec][it] = (sim->dobs[irec][it]-sim->dcal[irec][it])*acq->wdat[irec][it];
      fwi->fcost += sim->dres[irec][it]*sim->dres[irec][it];
      sim->dres[irec][it] *= acq->wdat[irec][it];//the weighting will mask diving waves
      _dres[irec][it] = -sim->dres[irec][it];//-R^H(\delta d - R\delta u)
    }//end for it
  }//end for irec
  MPI_Allreduce(&fwi->fcost, &tmp, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  fwi->fcost = 0.5*tmp*sim->dt;
  if(iproc==0) printf("fcost=%g\n", fwi->fcost);

  sprintf(fname, "dres_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr, "error opening file\n"); exit(1);}
  fwrite(&sim->dres[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);
  
  /*--------------------------------------------------------------*/
  if(iproc==0) printf("----stage 2: adjoint modelling -----------\n");
  sim->sign_dt = -1;
  for(it=sim->nt-1; it>=0; it--){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    //A^H(m+dm) lambda2=R^H(delta_d-R delta_u)
    inject_adjoint_source(sim, acq, sim->p3, sim->dres, it);
    fdtd_update_v(sim, 3, it, 1, kappa_, buz_, bux_, buy_);
    fdtd_update_p(sim, 3, it, 1, kappa_, buz_, bux_, buy_);

    inject_source(sim, acq, sim->p0, sim->stf[it]);
    fdtd_update_p(sim, 0, it, 0, kappa_, buz_, bux_, buy_);
    fdtd_update_v(sim, 0, it, 0, kappa_, buz_, bux_, buy_);
    decimate_interp_bndr(sim, 0, it, 1, sim->face1_, sim->face2_, sim->face3_);/* interp=1 */

    //A^H(m) lambda1=-R^H(delta_d-R delta_u)
    inject_adjoint_source(sim, acq, sim->p2, _dres, it);
    fdtd_update_v(sim, 2, it, 1, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_p(sim, 2, it, 1, sim->kappa, sim->buz, sim->bux, sim->buy);

    inject_source(sim, acq, sim->p1, sim->stf[it]);
    fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    decimate_interp_bndr(sim, 1, it, 1, sim->face1, sim->face2, sim->face3);/* interp=1 */

    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	i3_ = (sim->n3>1)?i3+sim->nb:0;
	in3 = (sim->n3>1)?(i3>=sim->order/2 && i3<sim->n3-sim->order):1;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2+sim->nb;
	  in2 =  (i2>=sim->order/2 && i2<sim->n2-sim->order);
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1+sim->nb;
	    in1 = (i1>=sim->order/2 && i1<sim->n1-sim->order);
	    in1 = in1 && (i1>fwi->ibathy[i3][i2]);

	    if(in1 && in2 && in3){//only compute below bathymetry
	      j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));

	      //G1=<lambda1|dA(m0)/dm0|u>
	      s1 = sim->p2[i3_][i2_][i1_]*sim->divv[i3_][i2_][i1_];//dJ/dln(kappa) under kappa-rho parametrization
	      s2 = (sim->vz2[i3_][i2_][i1_] + sim->vz2[i3_][i2_][i1_-1])*(sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1]);
	      s2 += (sim->vx2[i3_][i2_][i1_] + sim->vx2[i3_][i2_-1][i1_])*(sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_]);
	      if(sim->n3>1) s2 += (sim->vz2[i3_][i2_][i1_] + sim->vz2[i3_-1][i2_][i1_])*(sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_-1][i2_][i1_]);
	      s2 *= sim->rho[i3][i2][i1]*0.25;//dJ/dln(rho) under kappa-rho parametrization
	      if(fwi->idxpar[ipar]==1) g[j] += s1 - s2;//dJ/dln(vp)=dJ/dln(kappa)-dJ/dln(rho), add G1 only
	  
	      //G2=<lambda2|dA(m0+dm)/d[delta_m)]|u>
	      s1 = sim->p3[i3_][i2_][i1_]*sim->divv0[i3_][i2_][i1_];//dJ/dln(kappa) under kappa-rho parametrization
	      s2 = (sim->vz3[i3_][i2_][i1_] + sim->vz3[i3_][i2_][i1_-1])*(sim->dvzdt0[i3_][i2_][i1_] + sim->dvzdt0[i3_][i2_][i1_-1]);
	      s2 += (sim->vx3[i3_][i2_][i1_] + sim->vx3[i3_][i2_-1][i1_])*(sim->dvxdt0[i3_][i2_][i1_] + sim->dvxdt0[i3_][i2_-1][i1_]);
	      if(sim->n3>1) s2 += (sim->vz3[i3_][i2_][i1_] + sim->vz3[i3_-1][i2_][i1_])*(sim->dvzdt0[i3_][i2_][i1_] + sim->dvzdt0[i3_-1][i2_][i1_]);
	      s2 *= rho_[i3][i2][i1]*0.25;//dJ/dln(rho) under kappa-rho parametrization
	      if(fwi->idxpar[ipar]==1) g[j] += s1 - s2;//dJ/dln(vp)=dJ/dln(kappa)-dJ/dln(rho), add G2 also
	      if(fwi->idxpar[ipar]==2) g[j] += s1 + s2;//dJ/d[d(ln(ip))]=dJ/dln(ip)=dJ/dln(kappa)+dJ/dln(rho), add G2 only
	    }//end if
	  }//end for i1
	}//end for i2
      }//end for i3
    }//end for ipar

  }//end for it
  cpml_free(sim);
  extend_model_free(sim);
  fdtd_free(sim, 0);
  fdtd_free(sim, 1);
  fdtd_free(sim, 2);
  fdtd_free(sim, 3);
  decimate_interp_free(sim, 0);
  decimate_interp_free(sim, 1);
  computing_box_free(sim, 0);
  computing_box_free(sim, 1);

  tmp = sim->volume*sim->dt;
  for(j=0; j<fwi->n; j++) g[j] *= tmp;
  for(ipar=0; ipar<fwi->npar; ipar++){
    MPI_Allreduce(&g[ipar*sim->n123], &rho_[0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&g[ipar*sim->n123], &rho_[0][0][0], sim->n123*sizeof(float));
  }
  
  free3float(rho_);
  free3float(kappa_);
  free3float(buz_);
  free3float(bux_);
  free3float(buy_);
  free2float(d0);
  free2float(_dres);

  if(iproc==0){
    for(ipar=0; ipar<fwi->npar; ipar++){
      if(fwi->idxpar[ipar]==1){//grad_m0
	fp = fopen("gradient_fwi_m0", "wb");
	fwrite(&g[ipar*sim->n123], sim->n123*sizeof(float), 1, fp);
	fclose(fp);
      }
      if(fwi->idxpar[ipar]==2){//grad_dm
	fp = fopen("gradient_fwi_dm", "wb");
	fwrite(&g[ipar*sim->n123], sim->n123*sizeof(float), 1, fp);
	fclose(fp);
      }
    }
  }
  
  if(sim->mode==1 && fwi->firstgrad){
    /*
      s1 = fabs(x[0]);
      s2 = fabs(g[0]);
      for(j=0; j<fwi->n; j++){
      s1 = MAX(s1, fabs(x[j]));
      s2 = MAX(s2, fabs(g[j]));
      }
      fwi->alpha = 0.01*s1/s2;
    */
    s1 = 0;
    s2 = 0;
    for(j=0; j<fwi->n; j++){
      s1 += fabs(x[j]);
      s2 += fabs(g[j]);
    }
    if(s1==0) s1 = 1e5*fwi->n;//s1=0 happens if x=0 as input
    fwi->alpha = 5e-4*s1/s2;

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

