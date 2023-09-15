/* Linearized waveform inversion (Least-squares reverse time migration)
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
#include "fwi.h"
#include "mpi_info.h"
#include <mpi.h>

void read_data(sim_t *sim, acq_t *acq);
void write_data(sim_t *sim, acq_t *acq);
void setup_data_weight(acq_t *acq, sim_t *sim);

void check_cfl(sim_t *sim);

void fdtd_init(sim_t *sim, int flag);
void fdtd_null(sim_t *sim, int flag);
void fdtd_close(sim_t *sim, int flag);
void fdtd_update_v(sim_t *sim, int flag, int it, int adj);
void fdtd_update_p(sim_t *sim, int flag, int it, int adj);
void fdtd_update_v_scatter(sim_t *sim, int flag, int it, int adj);
void fdtd_update_p_scatter(sim_t *sim, int flag, int it, int adj);

void extend_model_init(sim_t *sim);
void extend_model(sim_t *sim);
void extend_model_close(sim_t *sim);

void computing_box_init(acq_t *acq, sim_t *sim, int adj);
void computing_box_close(sim_t *sim, int adj);

void cpml_init(sim_t *sim);
void cpml_close(sim_t *sim);

void decimate_interp_init(sim_t *sim);
void decimate_interp_close(sim_t *sim);
void decimate_interp_bndr(sim_t *sim, int interp, int it);

void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);
void inject_adjoint_source(sim_t *sim, acq_t *acq, float ***rp, float **dres, int it);

void precondition(sim_t *sim, fwi_t *fwi, float *x);

void do_lsrtm(sim_t *sim, acq_t *acq)
{
  char *bathyfile;
  int it, irec, ipar, j;
  int i1, i2, i3, i1_, i2_, i3_;
  float ****cg_b, ****cg_x, ****cg_p, ****cg_Ap, ****cg_r, ****cg_z;
  float ***g1, ***g2, ***h1, ***h2;
  float alpha, beta, pAp, rzold, rznew, rs;
  float s1, s2, tmp;
  fwi_t *fwi;
  FILE *fp;
  char fname[sizeof("dres_0000")];

  fwi = (fwi_t*)malloc(sizeof(fwi_t));
  fwi->bathy=alloc2float(sim->n2, sim->n3);
  fwi->ibathy=alloc2int(sim->n2, sim->n3);
  if(!getparstring("bathyfile",&bathyfile)){
    memset(fwi->bathy[0], 0, sim->n2*sim->n3*sizeof(float));
    memset(fwi->ibathy[0], 0, sim->n2*sim->n3*sizeof(int));
  }else{
    fp=fopen(bathyfile,"rb");
    if(fp==NULL) err("cannot open bathyfile=%s",bathyfile);
    if(fread(fwi->bathy[0],sizeof(float),sim->n2*sim->n3,fp)!=sim->n2*sim->n3) 
      err("error reading bathyfile=%s", bathyfile);
    fclose(fp);
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	fwi->ibathy[i3][i2]=NINT(fwi->bathy[i3][i2]/sim->d1);
      }
    }
  }
  if(!getparint("family", &fwi->family)) fwi->family = 1;//1=vp-rho; 2=vp-Ip
  if(!getparint("npar", &fwi->npar)) fwi->npar = 2;/* number of iterations */
  fwi->idxpar = alloc1int(fwi->npar);
  if(!getparint("idxpar", fwi->idxpar)) err("must provide idxpar=");//index of parameter in each family
  /*
    family 1: idxpar[0] =1, vp; idxpar[1]=2, rho
    family 2: idxpar[0] =1, vp; idxpar[1]=2, Ip
  */
  if(!getparint("niter", &fwi->niter)) fwi->niter = 30;/* number of iterations */
  if(!getparint("preco", &fwi->preco)) fwi->preco = 0;//1=precondition; 0=not

  fwi->n = fwi->npar*sim->n123;
  if(fwi->preco==2) fwi->hess = alloc1float(fwi->n);

  sim->dcal = alloc2float(sim->nt, acq->nrec);//synthetic data - direct wave
  sim->dobs = alloc2float(sim->nt, acq->nrec);
  sim->dres = alloc2float(sim->nt, acq->nrec);
  g1 = alloc3float(sim->n1, sim->n2, sim->n3);
  g2 = alloc3float(sim->n1, sim->n2, sim->n3);

  cg_b = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  cg_x = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  cg_r = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  cg_z = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  cg_p = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  cg_Ap = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);//A=L^H L

  read_data(sim, acq);
  setup_data_weight(acq, sim);//the muting will be used to remove direct waves
  
  check_cfl(sim);
  cpml_init(sim);  
  extend_model_init(sim);
  fdtd_init(sim, 0);//flag=0, scattering field
  fdtd_init(sim, 1);//flag=1, incident field
  fdtd_init(sim, 2);//flag=2, adjoint field
  decimate_interp_init(sim);
  extend_model(sim);
  computing_box_init(acq, sim, 0);
  computing_box_init(acq, sim, 1);

  h1 = g1;//reuse the same memory unit for H11
  h2 = g2;//reuse the same memory unit for H22
  memset(&h1[0][0][0], 0, sim->n123*sizeof(float));
  memset(&h2[0][0][0], 0, sim->n123*sizeof(float));
  memset(&cg_x[0][0][0][0], 0, fwi->npar*sim->n123*sizeof(float));
  if(iproc==0) printf("----Background field modelling --------\n");
  fdtd_null(sim, 1);//flag=1, incident field
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    decimate_interp_bndr(sim, 0, it);//interp=0
    fdtd_update_v(sim, 1, it, 0);
    fdtd_update_p(sim, 1, it, 0);
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);//create background field data

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
  tmp = 0;
  for(irec=0; irec<acq->nrec; irec++){
    for(it=0; it<sim->nt; it++){
      sim->dobs[irec][it] -= sim->dcal[irec][it];//only reflections are left in the observed data
      sim->dres[irec][it] = acq->wdat[irec][it]*sim->dobs[irec][it];//dref=Lm=0 when m=0
      tmp += sim->dres[irec][it]*sim->dres[irec][it];
    }
  }
  tmp *= 0.5*sim->dt;
  MPI_Allreduce(&tmp, &fwi->fcost_dat, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  if(iproc==0) printf("fcost=%g\n", fwi->fcost);

  sprintf(fname, "dres_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"error opening file\n"); exit(1);}
  fwrite(&sim->dres[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);
  
  if(iproc==0) printf("----migration: b=L^H (delta d)--------\n");
  memset(&g1[0][0][0], 0, sim->n123*sizeof(float));
  memset(&g2[0][0][0], 0, sim->n123*sizeof(float));
  fdtd_null(sim, 2);//flag=2, adjoint field
  sim->sign_dt = -1;
  for(it=sim->nt-1; it>=0; it--){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    //adjoint field modelling
    inject_adjoint_source(sim, acq, sim->p2, sim->dres, it);
    fdtd_update_v(sim, 2, it, 1);
    fdtd_update_p(sim, 2, it, 1);

    //reconstruction of the background field in reverse time order
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    fdtd_update_p(sim, 1, it, 0);
    fdtd_update_v(sim, 1, it, 0);
    decimate_interp_bndr(sim, 1, it);//interp=1

    for(i3=0; i3<sim->n3; i3++){
      i3_ = (sim->n3>1)?i3+sim->nb:0;
      for(i2=0; i2<sim->n2; i2++){
	i2_ = i2+sim->nb;
	for(i1=0; i1<sim->n1; i1++){
	  i1_ = i1+sim->nb;
	  g1[i3][i2][i1] += sim->divv[i3_][i2_][i1_]*sim->p2[i3_][i2_][i1_];
	  g2[i3][i2][i1] += sim->dvzdt[i3_][i2_][i1_]*sim->vz2[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1]*sim->vz2[i3_][i2_][i1_-1];
	  g2[i3][i2][i1] += sim->dvxdt[i3_][i2_][i1_]*sim->vx2[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_]*sim->vx2[i3_][i2_-1][i1_];
	  if(sim->n3>1) g2[i3][i2][i1] += sim->dvydt[i3_][i2_][i1_]*sim->vy2[i3_][i2_][i1_] + sim->dvydt[i3_-1][i2_][i1_]*sim->vy2[i3_-1][i2_][i1_];
	}
      }
    }
  }//end for it
  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      for(i1=0; i1<sim->n1; i1++){
	tmp = sim->volume*sim->dt;//kappa
	g1[i3][i2][i1] *= tmp;//dJ/dln(kappa)
	g2[i3][i2][i1] *= sim->rho[i3][i2][i1]*0.5*tmp;//dJ/dln(rho)
	if(i1<= fwi->ibathy[i3][i2]){//reset again to avoid leakage
	  g1[i3][i2][i1] = 0.;
	  g2[i3][i2][i1] = 0.;
	}
      }
    }
  }
  MPI_Allreduce(&g1[0][0][0], &cg_b[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(&g1[0][0][0], &cg_b[0][0][0][0], sim->n123*sizeof(float));
  MPI_Allreduce(&g2[0][0][0], &cg_b[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(&g2[0][0][0], &cg_b[0][0][0][0], sim->n123*sizeof(float));

  rs = 0;
  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(i1=0; i1<sim->n1; i1++){
	  if(fwi->family==1){//vp-rho
	    if(fwi->idxpar[ipar]==1) cg_b[ipar][i3][i2][i1] = 2.0*g1[i3][i2][i1];//dJ/dln(vp)=2*dJ/dln(kappa)
	    if(fwi->idxpar[ipar]==2) cg_b[ipar][i3][i2][i1] = g1[i3][i2][i1] + g2[i3][i2][i1];//dJ/dln(rho')=dJ/dln(rho) + dJ/dln(kappa)
	  }else if(fwi->family==2){//vp-Ip
	    if(fwi->idxpar[ipar]==1) cg_b[ipar][i3][i2][i1] = g1[i3][i2][i1] - g2[i3][i2][i1];//dJ/dln(vp)=dJ/dln(kappa) - dJ/dln(rho)
	    if(fwi->idxpar[ipar]==2) cg_b[ipar][i3][i2][i1] = g1[i3][i2][i1] + g2[i3][i2][i1];//dJ/dln(Ip)=dJ/dln(kappa) + dJ/dln(rho)
	  }//end if

	  cg_r[ipar][i3][i2][i1] = cg_b[ipar][i3][i2][i1];//r=b;
	  cg_p[ipar][i3][i2][i1] = cg_r[ipar][i3][i2][i1];//p=r;
	  rs += cg_p[ipar][i3][i2][i1]*cg_p[ipar][i3][i2][i1];
	}//end for i1
      }//end for i2
    }//end for i3
  }//end for ipar
  fwi->fcost = fwi->fcost_dat;
  rzold = rs;
  if(fwi->preco) {
    precondition(sim, fwi, &cg_p[0][0][0][0]);
    rzold = 0;
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    cg_z[ipar][i3][i2][i1] = cg_p[ipar][i3][i2][i1];
	    rzold += cg_r[ipar][i3][i2][i1]*cg_z[ipar][i3][i2][i1];
	  }//end for i1
	}//end for i2
      }//end for i3
    }//end for ipar
  }    

  if(iproc==0){
    fp = fopen("iterate.txt","w");
    fprintf(fp,"==========================================================\n");
    fprintf(fp,"Maximum number of CG iterations: %d\n", fwi->niter);
    fprintf(fp,"initial data residual: fwi->fcost_dat=0.5*|dobs|^2=%.4e\n", fwi->fcost_dat);
    fprintf(fp,"==========================================================\n");
    fprintf(fp,"iter    fcost_cg       fcost     ||r||\n");
    fclose(fp);
  }
  //start CG iterations
  for(fwi->iter=0; fwi->iter<fwi->niter; fwi->iter++) {
    if(iproc==0) {
      tmp = fwi->fcost-fwi->fcost_dat;//fcost_cg=0.5*x^H Ax-x^H b=fcost-0.5*|dobs|
      fp=fopen("iterate.txt","a");
      fprintf(fp,"%d    %.4e       %.4e     %.4e\n", fwi->iter, tmp, fwi->fcost, sqrt(rs));
      fclose(fp);
      printf("======= CG %d: fcost_cg=%.4e, ||r||=%.4e ========\n", fwi->iter, tmp, sqrt(rs));

      fp = fopen("param_final", "wb");
      fwrite(&cg_x[0][0][0][0], fwi->npar*sim->n123*sizeof(float), 1, fp);
      fclose(fp);
      if(fwi->iter==0) fp = fopen("param_iter", "wb");
      else             fp = fopen("param_iter", "ab");
      fwrite(&cg_x[0][0][0][0], fwi->npar*sim->n123*sizeof(float), 1, fp);
      fclose(fp);
    }
    
    if(iproc==0) printf("------ Born modelling/demigration, Lp--------\n");
    //use cg_p to construct [dmv,dmp] to excite Born modelling
    memset(&g1[0][0][0], 0, sim->n123*sizeof(float));
    memset(&g2[0][0][0], 0, sim->n123*sizeof(float));
    for(ipar=0; ipar<fwi->npar; ipar++){
      if(fwi->idxpar[ipar]==1) memcpy(&g1[0][0][0], &cg_p[ipar][0][0][0], sim->n123*sizeof(float));
      if(fwi->idxpar[ipar]==2) memcpy(&g2[0][0][0], &cg_p[ipar][0][0][0], sim->n123*sizeof(float));
    }//end for ipar
    fdtd_null(sim, 0);//flag=0, scattering field
    fdtd_null(sim, 1);//flag=1, incident field
    sim->sign_dt = 1;
    for(it=0; it<sim->nt; it++){
      if(iproc==0 && it%100==0) printf("it-----%d\n", it);

      //background field modelling, p1
      decimate_interp_bndr(sim, 0, it);//interp=0
      fdtd_update_v(sim, 1, it, 0);
      fdtd_update_p(sim, 1, it, 0);
      inject_source(sim, acq, sim->p1, sim->stf[it]);

      //scattering field modelling, p0 (using p1 to construct source)
      fdtd_update_v(sim, 0, it, 0);
      for(i3=0; i3<sim->n3; i3++){
	i3_ = (sim->n3>1)?i3+sim->nb:0;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2+sim->nb;
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1+sim->nb;

	    tmp = 0;
	    if(fwi->family==1) tmp = g2[i3][i2][i1];//vp-rho
	    if(fwi->family==2) tmp = g2[i3][i2][i1] - g1[i3][i2][i1];//vp-Ip
	    
	    sim->vz0[i3_][i2_][i1_] += -sim->dt*tmp*sim->dvzdt[i3_][i2_][i1_];
	    sim->vx0[i3_][i2_][i1_] += -sim->dt*tmp*sim->dvxdt[i3_][i2_][i1_];
	    if(sim->n3>1) sim->vy0[i3_][i2_][i1_] += -sim->dt*tmp*sim->dvydt[i3_][i2_][i1_];
	  }//end for i1
	}//end for i2
      }//end for i3
      fdtd_update_p(sim, 0, it, 0);
      for(i3=0; i3<sim->n3; i3++){
	i3_ = (sim->n3>1)?i3+sim->nb:0;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2+sim->nb;
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1+sim->nb;

	    tmp = 0;
	    if(fwi->family==1) tmp = 2.*g1[i3][i2][i1] + g2[i3][i2][i1];//vp-rho
	    if(fwi->family==2) tmp = g1[i3][i2][i1] + g2[i3][i2][i1];//vp-Ip

	    sim->p0[i3_][i2_][i1_] += -sim->dt*tmp*sim->kappa[i3_][i2_][i1_]*sim->divv[i3_][i2_][i1_];
	  }//end for i1
	}//end for i2
      }//end for i3
      extract_wavefield(sim, acq, sim->p0, sim->dcal, it);//reflection data, stored in dcal=L*dm
    }

    /* pAp = 0;//pAp=<Lp,Lp> */
    /* for(irec=0; irec<acq->nrec; irec++){ */
    /*   for(it=0; it<sim->nt; it++){ */
    /* 	sim->dres[irec][it] = acq->wdat[irec][it]*sim->dcal[irec][it];//dref=Lm=0 when m=0 */
    /* 	pAp += sim->dres[irec][it]*sim->dres[irec][it]; */
    /*   } */
    /* } */
    /* tmp = pAp *0.5*sim->dt; */
    /* MPI_Allreduce(&tmp, &pAp, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); */
    /* alpha = rzold/pAp; */
    
    if(iproc==0) printf("-------- Migration: Ap=L^H (Lp) ----------\n");
    memset(g1[0][0], 0, sim->n123*sizeof(float));
    memset(g2[0][0], 0, sim->n123*sizeof(float));
    fdtd_null(sim, 2);//flag=2, adjoint field
    sim->sign_dt = -1;
    for(it=sim->nt-1; it>=0; it--){
      if(iproc==0 && it%100==0) printf("it-----%d\n", it);

      //adjoint field modelling
      inject_adjoint_source(sim, acq, sim->p2, sim->dcal, it);
      fdtd_update_v(sim, 2, it, 1);
      fdtd_update_p(sim, 2, it, 1);

      //reconstruction of the background field in reverse time order
      inject_source(sim, acq, sim->p1, sim->stf[it]);
      fdtd_update_p(sim, 1, it, 0);
      fdtd_update_v(sim, 1, it, 0);
      decimate_interp_bndr(sim, 1, it);//interp=1

      for(i3=0; i3<sim->n3; i3++){
	i3_ = (sim->n3>1)?i3+sim->nb:0;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2+sim->nb;
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1+sim->nb;
	    g1[i3][i2][i1] += sim->divv[i3_][i2_][i1_]*sim->p2[i3_][i2_][i1_];
	    g2[i3][i2][i1] += sim->dvzdt[i3_][i2_][i1_]*sim->vz2[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1]*sim->vz2[i3_][i2_][i1_-1];
	    g2[i3][i2][i1] += sim->dvxdt[i3_][i2_][i1_]*sim->vx2[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_]*sim->vx2[i3_][i2_-1][i1_];
	    if(sim->n3>1) g2[i3][i2][i1] += sim->dvydt[i3_][i2_][i1_]*sim->vy2[i3_][i2_][i1_] + sim->dvydt[i3_-1][i2_][i1_]*sim->vy2[i3_-1][i2_][i1_];
	  }
	}
      }
    }//end for it
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(i1=0; i1<sim->n1; i1++){
	  tmp = -sim->volume*sim->dt;//note a missing minus sign in adjoint source - dcal
	  g1[i3][i2][i1] *= tmp;//dJ/dln(kappa)
	  g2[i3][i2][i1] *= sim->rho[i3][i2][i1]*0.5*tmp;//dJ/dln(rho)
	  if(i1<= fwi->ibathy[i3][i2]){//reset again to avoid leakage
	    g1[i3][i2][i1] = 0.;
	    g2[i3][i2][i1] = 0.;
	  }
	}
      }
    }
    MPI_Allreduce(&g1[0][0][0], &cg_Ap[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&g1[0][0][0], &cg_Ap[0][0][0][0], sim->n123*sizeof(float));
    MPI_Allreduce(&g2[0][0][0], &cg_Ap[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&g2[0][0][0], &cg_Ap[0][0][0][0], sim->n123*sizeof(float));

    pAp = 0;
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    if(fwi->family==1){//vp-rho
	      if(fwi->idxpar[ipar]==1) cg_Ap[ipar][i3][i2][i1] = g1[i3][i2][i1]*2.0;//dJ/dln(Vp)=2*dJ/dln(kappa)
	      if(fwi->idxpar[ipar]==2) cg_Ap[ipar][i3][i2][i1] = g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(rho')=dJ/dln(rho)+dJ/dln(kappa)
	    }else if(fwi->family==2){//vp-Ip
	      if(fwi->idxpar[ipar]==1) cg_Ap[ipar][i3][i2][i1] = -g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(vp)=-dJ/dln(rho)+dJ/dln(kappa)
	      if(fwi->idxpar[ipar]==2) cg_Ap[ipar][i3][i2][i1] = g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(Ip)=dJ/dln(rho)+dJ/dln(kappa)
	    }//end if
	    pAp += cg_p[ipar][i3][i2][i1]*cg_Ap[ipar][i3][i2][i1];
	  }//end for i1
	}//end for i2
      }//end for i3
    }//end for ipar
    alpha = rzold/pAp;
    //printf("alpha=%e |r|^2=%e pAp=%e \n", alpha, rzold, pAp);
    
    rs = 0;
    tmp = 0;
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    cg_x[ipar][i3][i2][i1] += alpha*cg_p[ipar][i3][i2][i1];
	    cg_r[ipar][i3][i2][i1] -= alpha*cg_Ap[ipar][i3][i2][i1];
	    rs += cg_r[ipar][i3][i2][i1]*cg_r[ipar][i3][i2][i1];
	    tmp += cg_x[ipar][i3][i2][i1]*(cg_r[ipar][i3][i2][i1] + cg_b[ipar][i3][i2][i1]);
	  }
	}
      }
    }
    tmp *= -0.5;
    fwi->fcost = tmp + fwi->fcost_dat;//0.5*||Lx-d||=0.5*x^H (L^H L)x-x^H d + 0.5*|d|^2=fcost_cg+0.5*|d|^2
    rznew = rs;
    //fcost += fwi->fcost_dat;

    memcpy(&cg_z[0][0][0][0], &cg_r[0][0][0][0], fwi->n*sizeof(float));
    if(fwi->preco) {
      precondition(sim, fwi, &cg_z[0][0][0][0]);
      rznew = 0;
      for(ipar=0; ipar<fwi->npar; ipar++){
	for(i3=0; i3<sim->n3; i3++){
	  for(i2=0; i2<sim->n2; i2++){
	    for(i1=0; i1<sim->n1; i1++){
	      rznew += cg_r[ipar][i3][i2][i1]*cg_z[ipar][i3][i2][i1];
	    }//end for i1
	  }//end for i2
	}//end for i3
      }//end for ipar
    }
    
    beta = rznew/rzold;    
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    cg_p[ipar][i3][i2][i1] = cg_z[ipar][i3][i2][i1] + beta*cg_p[ipar][i3][i2][i1];
	  }
	}
      }
    }
    rzold = rznew;
  }//end for iter

  cpml_close(sim);
  extend_model_close(sim);
  computing_box_close(sim, 0);
  computing_box_close(sim, 1);
  fdtd_close(sim, 0);
  fdtd_close(sim, 1);
  fdtd_close(sim, 2);
  decimate_interp_close(sim);

  free2float(sim->dcal);
  free2float(sim->dobs);
  free2float(sim->dres);
  free3float(g1);
  free3float(g2);

  free4float(cg_b);
  free4float(cg_x);
  free4float(cg_r);
  free4float(cg_z);
  free4float(cg_p);
  free4float(cg_Ap);

  if(fwi->preco==2) free1float(fwi->hess);
  free1int(fwi->idxpar);
  free(fwi);

}

