/* Data-domain Linearized waveform inversion (Least-squares RTM)
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
#include "acq.h"
#include "fwi.h"
 
#include <mpi.h>

void read_data(sim_t *sim, acq_t *acq);
void write_data(sim_t *sim, acq_t *acq);
void setup_data_weight(acq_t *acq, sim_t *sim);

void check_cfl(sim_t *sim);

void fdtd_init(sim_t *sim, int flag);
void fdtd_null(sim_t *sim, int flag);
void fdtd_free(sim_t *sim, int flag);
void fdtd_update_v(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy);
void fdtd_update_p(sim_t *sim, int flag, int it, int adj, float ***kappa, float ***buz, float ***bux, float ***buy);

void extend_model_init(sim_t *sim);
void extend_model(sim_t *sim, float ***vp, float ***rho, float ***kappa, float ***buz, float ***bux, float ***buy);
void extend_model_free(sim_t *sim);

void computing_box_init(acq_t *acq, sim_t *sim, int adj);
void computing_box_free(sim_t *sim, int adj);

void cpml_init(sim_t *sim);
void cpml_free(sim_t *sim);

void decimate_interp_init(sim_t *sim, int flag);
void decimate_interp_free(sim_t *sim, int flag);
void decimate_interp_bndr(sim_t *sim, int flag, int it, int interp, float **face1, float **face2, float **face3);

void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);
void inject_adjoint_source(sim_t *sim, acq_t *acq, float ***rp, float **dres, int it);


//PCGNR=preconditioned CGNR algorithm, see Algorithm 9.7 in Saad book
void do_lsrtm(sim_t *sim, acq_t *acq)
{
  char *bathyfile;
  int it, irec, ipar;
  int i1, i2, i3, i1_, i2_, i3_;
  float ****cg_x, ****cg_p, ****cg_z, ****cg_rt, ****hess;
  float **cg_b, **cg_r, **cg_Ap;
  float ***g1, ***g2, ***h1, ***h2;
  double alpha, beta, zrt_old, zrt_new, ws, rs;
  double s1, s2, tmp, maxval, rs0;
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
  if(!getparint("npar", &fwi->npar)) fwi->npar = 2;//number of parameters
  fwi->idxpar = alloc1int(fwi->npar);
  if(!getparint("idxpar", fwi->idxpar)) err("must provide idxpar=");//index of parameter in each family
  /*family 1: idxpar[0] =1, vp; idxpar[1]=2, rho
    family 2: idxpar[0] =1, vp; idxpar[1]=2, Ip  */
  if(!getparint("niter", &fwi->niter)) fwi->niter = 15;//number of iterations 
  if(!getparint("preco", &fwi->preco)) fwi->preco = 0;//1=precondition, 0=not

  fwi->n = fwi->npar*sim->n123;//number of unknowns
  sim->dcal = alloc2float(sim->nt, acq->nrec);//synthetic data - direct wave
  sim->dobs = alloc2float(sim->nt, acq->nrec);
  sim->dres = alloc2float(sim->nt, acq->nrec);
  g1 = alloc3float(sim->n1, sim->n2, sim->n3);
  g2 = alloc3float(sim->n1, sim->n2, sim->n3);

  cg_b = alloc2float(sim->nt, acq->nrec);
  cg_r = alloc2float(sim->nt, acq->nrec);
  cg_Ap = alloc2float(sim->nt, acq->nrec);
  cg_x = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  cg_z = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  cg_p = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  cg_rt = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  if(fwi->preco) hess = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  memset(&cg_x[0][0][0][0], 0, fwi->n*sizeof(float));

  read_data(sim, acq);
  setup_data_weight(acq, sim);//the muting will be used to remove direct waves
  
  check_cfl(sim);
  cpml_init(sim);  
  extend_model_init(sim);
  fdtd_init(sim, 0);//flag=0, scattering field
  fdtd_init(sim, 1);//flag=1, incident field
  fdtd_init(sim, 2);//flag=2, adjoint field
  decimate_interp_init(sim, 1);
  extend_model(sim, sim->vp, sim->rho, sim->kappa, sim->buz, sim->bux, sim->buy);
  computing_box_init(acq, sim, 0);
  computing_box_init(acq, sim, 1);

  //==============================================================
  if(iproc==0) printf("----Born modelling, p0=0, simulate background field p1 --------\n");
  if(fwi->preco){
    h1 = g1;
    h2 = g2;
    memset(&h1[0][0][0], 0, sim->n123*sizeof(float));
    memset(&h2[0][0][0], 0, sim->n123*sizeof(float));
  }
  fdtd_null(sim, 1);//flag=1, incident field
  fdtd_null(sim, 0);//flag=0, scattering field
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    //background field modelling, p1
    decimate_interp_bndr(sim, 1, it, 0, sim->face1, sim->face2, sim->face3);//interp=0
    fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);//direct wave + diving wave

    if(fwi->preco){
      for(i3=0; i3<sim->n3; i3++){
	i3_ = (sim->n3>1)?i3+sim->nb:0;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2+sim->nb;
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1+sim->nb;
	    h1[i3][i2][i1] += sim->divv[i3_][i2_][i1_]*sim->divv[i3_][i2_][i1_];
	    h2[i3][i2][i1] += (sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1])*(sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1]);
	    h2[i3][i2][i1] += (sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_])*(sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_]);
	    if(sim->n3>1) h2[i3][i2][i1] += (sim->dvydt[i3_][i2_][i1_] + sim->dvydt[i3_-1][i2_][i1_])*(sim->dvydt[i3_][i2_][i1_] + sim->dvydt[i3_-1][i2_][i1_]);
	  }//end for i1
	}//end for i2
      }//end for i3
    }//end if
  }
  if(fwi->preco){
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
    MPI_Allreduce(&h1[0][0][0], &hess[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&h1[0][0][0], &hess[0][0][0][0], sim->n123*sizeof(float));
    MPI_Allreduce(&h2[0][0][0], &hess[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&h2[0][0][0], &hess[0][0][0][0], sim->n123*sizeof(float));
    maxval = 0;
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    if(fwi->family==1){//vp-rho
	      if(fwi->idxpar[ipar]==1) hess[ipar][i3][i2][i1] = 4.0*h1[i3][i2][i1];//H_{ln(vp),ln(vp)}
	      if(fwi->idxpar[ipar]==2) hess[ipar][i3][i2][i1] = h1[i3][i2][i1] + h2[i3][i2][i1];//H_{ln(rho),ln(rho)}
	    }
	    if(fwi->family==2){//vp-ip
	      if(fwi->idxpar[ipar]==1) hess[ipar][i3][i2][i1] = h1[i3][i2][i1] + h2[i3][i2][i1];//H_{ln(vp),ln(vp)}
	      if(fwi->idxpar[ipar]==2) hess[ipar][i3][i2][i1] = h1[i3][i2][i1] + h2[i3][i2][i1];//H_{ln(ip),ln(ip)}
	    }
	    maxval = MAX(hess[ipar][i3][i2][i1], maxval);
	  }
	}
      }
    }
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    hess[ipar][i3][i2][i1] += 1e-3*maxval;//avoid division by zero
	  }
	}
      }
    } 
  }//end if
  tmp = 0;
  for(irec=0; irec<acq->nrec; irec++){
    for(it=0; it<sim->nt; it++){
      sim->dres[irec][it] = (double)(sim->dobs[irec][it]-sim->dcal[irec][it])*acq->wdat[irec][it];//muting direct wave should happen here
      cg_b[irec][it] = sim->dres[irec][it];//since dcal=Lx0=0 for x0=0
      cg_r[irec][it] = cg_b[irec][it];
      tmp += (double)cg_r[irec][it]*cg_r[irec][it];
      sim->dres[irec][it] *= (double)acq->wdat[irec][it];//prepare adjoint source
    }
  }
  tmp *= sim->dt;
  MPI_Allreduce(&tmp, &rs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  rs0 = rs;
  
  sprintf(fname, "b_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"error opening file\n"); exit(1);}
  fwrite(&cg_b[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);

  sprintf(fname, "d0_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"error opening file\n"); exit(1);}
  fwrite(&sim->dcal[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);
  
  if(iproc==0) printf("----migration: z0=L^H r0, r0=b-Lx0--------\n");
  memset(&g1[0][0][0], 0, sim->n123*sizeof(float));
  memset(&g2[0][0][0], 0, sim->n123*sizeof(float));
  fdtd_null(sim, 2);//flag=2, adjoint field
  sim->sign_dt = -1;
  for(it=sim->nt-1; it>=0; it--){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    //adjoint field modelling
    inject_adjoint_source(sim, acq, sim->p2, sim->dres, it);
    fdtd_update_v(sim, 2, it, 1, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_p(sim, 2, it, 1, sim->kappa, sim->buz, sim->bux, sim->buy);

    //reconstruction of the background field p1 in reverse time order
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    decimate_interp_bndr(sim, 1, it, 1, sim->face1, sim->face2, sim->face3);//interp=1

    for(i3=0; i3<sim->n3; i3++){
      i3_ = (sim->n3>1)?i3+sim->nb:0;
      for(i2=0; i2<sim->n2; i2++){
	i2_ = i2+sim->nb;
	for(i1=0; i1<sim->n1; i1++){
	  i1_ = i1+sim->nb;
	  g1[i3][i2][i1] += sim->divv[i3_][i2_][i1_]*sim->p2[i3_][i2_][i1_];
	  g2[i3][i2][i1] += (sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1])*(sim->vz2[i3_][i2_][i1_] + sim->vz2[i3_][i2_][i1_-1]);
	  g2[i3][i2][i1] += (sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_])*(sim->vx2[i3_][i2_][i1_] + sim->vx2[i3_][i2_-1][i1_]);
	  if(sim->n3>1) g2[i3][i2][i1] += (sim->dvydt[i3_][i2_][i1_] + sim->dvydt[i3_-1][i2_][i1_])*(sim->vy2[i3_][i2_][i1_] + sim->vy2[i3_-1][i2_][i1_]);
	}
      }
    }
  }//end for it
  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      for(i1=0; i1<sim->n1; i1++){
	tmp = -sim->volume*sim->dt;//kappa
	g1[i3][i2][i1] *= tmp;//dJ/dln(kappa)
	g2[i3][i2][i1] *= sim->rho[i3][i2][i1]*0.25*tmp;//dJ/dln(rho)
	if(i1<= fwi->ibathy[i3][i2]){//reset again to avoid leakage
	  g1[i3][i2][i1] = 0.;
	  g2[i3][i2][i1] = 0.;
	}
      }
    }
  }
  MPI_Allreduce(&g1[0][0][0], &cg_z[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(&g1[0][0][0], &cg_z[0][0][0][0], sim->n123*sizeof(float));
  MPI_Allreduce(&g2[0][0][0], &cg_z[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(&g2[0][0][0], &cg_z[0][0][0][0], sim->n123*sizeof(float));
  zrt_old = 0;
  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(i1=0; i1<sim->n1; i1++){
	  if(fwi->family==1){//vp-rho
	    if(fwi->idxpar[ipar]==1) cg_rt[ipar][i3][i2][i1] = 2.0*g1[i3][i2][i1];//dJ/dln(vp)=2*dJ/dln(kappa)
	    if(fwi->idxpar[ipar]==2) cg_rt[ipar][i3][i2][i1] = g1[i3][i2][i1] + g2[i3][i2][i1];//dJ/dln(rho')=dJ/dln(rho) + dJ/dln(kappa)
	  }else if(fwi->family==2){//vp-Ip
	    if(fwi->idxpar[ipar]==1) cg_rt[ipar][i3][i2][i1] = g1[i3][i2][i1] - g2[i3][i2][i1];//dJ/dln(vp)=dJ/dln(kappa) - dJ/dln(rho)
	    if(fwi->idxpar[ipar]==2) cg_rt[ipar][i3][i2][i1] = g1[i3][i2][i1] + g2[i3][i2][i1];//dJ/dln(Ip)=dJ/dln(kappa) + dJ/dln(rho)
	  }//end if
	  cg_z[ipar][i3][i2][i1] = cg_rt[ipar][i3][i2][i1];
	  if(fwi->preco) cg_z[ipar][i3][i2][i1] /= (double)hess[ipar][i3][i2][i1];
	  cg_p[ipar][i3][i2][i1] = cg_z[ipar][i3][i2][i1];//p=z
	  zrt_old += (double)cg_z[ipar][i3][i2][i1]*cg_rt[ipar][i3][i2][i1];
	}//end for i1
      }//end for i2
    }//end for i3
  }//end for ipar
  if(iproc==0){
    fp = fopen("param_final_rtm", "wb");
    fwrite(&cg_rt[0][0][0][0], fwi->n*sizeof(float), 1, fp);
    fclose(fp);
  }
  
  //start CG iterations
  for(fwi->iter=0; fwi->iter<fwi->niter; fwi->iter++){
    if(iproc==0){
      if(fwi->iter==0){
	fp = fopen("iterate.txt","w");
	fprintf(fp,"===========================================\n");
	fprintf(fp,"Number of PCGNR iterations: %d\n", fwi->niter);
	fprintf(fp,"fcost=0.5||Lx- delta_d||^2,  with x0=0\n");
	fprintf(fp,"preco=%d\n", fwi->preco);
	fprintf(fp,"===========================================\n");
	fprintf(fp,"iter     fcost\n");
	fclose(fp);
      }
      fp=fopen("iterate.txt","a");
      fprintf(fp,"%d    %.4e\n", fwi->iter, rs/rs0);
      fclose(fp);
      printf("======= iter=%d, fcost=%.4e ========\n", fwi->iter, rs/rs0);

      fp = fopen("param_final", "wb");
      fwrite(&cg_x[0][0][0][0], fwi->n*sizeof(float), 1, fp);
      fclose(fp);

      if(sim->n3==1){//we only store intermediate models for 2D case
	if(fwi->iter==0) fp = fopen("param_iter", "wb");
	else             fp = fopen("param_iter", "ab");
	fwrite(&cg_x[0][0][0][0], fwi->n*sizeof(float), 1, fp);
	fclose(fp);
      }
    }
    
    if(iproc==0) printf("------ Born modelling/demigration, w_k=L p_k--------\n");
    //use cg_p to construct [dmv,dmp] to excite Born modelling
    memset(&g1[0][0][0], 0, sim->n123*sizeof(float));
    memset(&g2[0][0][0], 0, sim->n123*sizeof(float));
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    if(i1>fwi->ibathy[i3][i2]){//reset again to avoid leakage
	      if(fwi->idxpar[ipar]==1) g1[i3][i2][i1] = cg_p[ipar][i3][i2][i1];
	      if(fwi->idxpar[ipar]==2) g2[i3][i2][i1] = cg_p[ipar][i3][i2][i1];
	    }
	  }
	}
      }
    }//end for ipar
    fdtd_null(sim, 0);//flag=0, scattering field
    fdtd_null(sim, 1);//flag=1, incident field
    memset(&sim->dcal[0][0], 0, sim->nt*acq->nrec*sizeof(float));
    sim->sign_dt = 1;
    for(it=0; it<sim->nt; it++){
      if(iproc==0 && it%100==0) printf("it-----%d\n", it);

      //background field modelling, p1
      decimate_interp_bndr(sim, 1, it, 0, sim->face1, sim->face2, sim->face3);//interp=0
      fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
      fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
      inject_source(sim, acq, sim->p1, sim->stf[it]);

      //scattering field modelling, p0 (using p1 to construct source)
      fdtd_update_v(sim, 0, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
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
      fdtd_update_p(sim, 0, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
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
      extract_wavefield(sim, acq, sim->p0, sim->dcal, it);//Born modelled data, stored in dcal=L*dm
    }//end for it

    tmp = 0;
    for(irec=0; irec<acq->nrec; irec++){
      for(it=0; it<sim->nt; it++){
	cg_Ap[irec][it] = (double)acq->wdat[irec][it]*sim->dcal[irec][it];//w=Ap
	tmp += (double)cg_Ap[irec][it]*cg_Ap[irec][it];
      }
    }
    tmp *= sim->dt;
    MPI_Allreduce(&tmp, &ws, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    alpha = zrt_old/ws;

    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    cg_x[ipar][i3][i2][i1] += alpha*cg_p[ipar][i3][i2][i1];
	  }//end for i1
	}//end for i2
      }//end for i3
    }//end for ipar
    tmp = 0;
    for(irec=0; irec<acq->nrec; irec++){
      for(it=0; it<sim->nt; it++){
	cg_r[irec][it] -= alpha*cg_Ap[irec][it];//dref=Lm=0 when m=0
	tmp += (double)cg_r[irec][it]*cg_r[irec][it];
	sim->dres[irec][it] = (double)cg_r[irec][it]*acq->wdat[irec][it];
      }
    }
    tmp *= sim->dt;
    MPI_Allreduce(&tmp, &rs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    sprintf(fname, "dres_%04d", acq->shot_idx[iproc]);
    fp=fopen(fname,"wb");
    if(fp==NULL) { fprintf(stderr,"error opening file\n"); exit(1);}
    fwrite(&cg_r[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
    fclose(fp);
    fflush(stdout);

    if(iproc==0) printf("-------- Migration: z_{k+1}=L^H r_{k+1} ----------\n");
    memset(g1[0][0], 0, sim->n123*sizeof(float));
    memset(g2[0][0], 0, sim->n123*sizeof(float));
    fdtd_null(sim, 2);//flag=2, adjoint field
    sim->sign_dt = -1;
    for(it=sim->nt-1; it>=0; it--){
      if(iproc==0 && it%100==0) printf("it-----%d\n", it);

      //adjoint field modelling
      inject_adjoint_source(sim, acq, sim->p2, sim->dres, it);
      fdtd_update_v(sim, 2, it, 1, sim->kappa, sim->buz, sim->bux, sim->buy);
      fdtd_update_p(sim, 2, it, 1, sim->kappa, sim->buz, sim->bux, sim->buy);

      //reconstruction of the background field in reverse time order
      inject_source(sim, acq, sim->p1, sim->stf[it]);
      fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
      fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
      decimate_interp_bndr(sim, 1, it, 1, sim->face1, sim->face2, sim->face3);//interp=1

      for(i3=0; i3<sim->n3; i3++){
	i3_ = (sim->n3>1)?i3+sim->nb:0;
	for(i2=0; i2<sim->n2; i2++){
	  i2_ = i2+sim->nb;
	  for(i1=0; i1<sim->n1; i1++){
	    i1_ = i1+sim->nb;
	    g1[i3][i2][i1] += sim->divv[i3_][i2_][i1_]*sim->p2[i3_][i2_][i1_];
	    g2[i3][i2][i1] += (sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1])*(sim->vz2[i3_][i2_][i1_] + sim->vz2[i3_][i2_][i1_-1]);
	    g2[i3][i2][i1] += (sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_])*(sim->vx2[i3_][i2_][i1_] + sim->vx2[i3_][i2_-1][i1_]);
	    if(sim->n3>1) g2[i3][i2][i1] += (sim->dvydt[i3_][i2_][i1_] + sim->dvydt[i3_-1][i2_][i1_])*(sim->vy2[i3_][i2_][i1_] + sim->vy2[i3_-1][i2_][i1_]);
	  }
	}
      }
    }//end for it
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(i1=0; i1<sim->n1; i1++){
	  tmp = -sim->volume*sim->dt;//note a missing minus sign in adjoint source - dcal
	  g1[i3][i2][i1] *= tmp;//dJ/dln(kappa)
	  g2[i3][i2][i1] *= sim->rho[i3][i2][i1]*0.25*tmp;//dJ/dln(rho)
	  if(i1<= fwi->ibathy[i3][i2]){//reset again to avoid leakage
	    g1[i3][i2][i1] = 0.;
	    g2[i3][i2][i1] = 0.;
	  }
	}
      }
    }
    MPI_Allreduce(&g1[0][0][0], &cg_z[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&g1[0][0][0], &cg_z[0][0][0][0], sim->n123*sizeof(float));
    MPI_Allreduce(&g2[0][0][0], &cg_z[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    memcpy(&g2[0][0][0], &cg_z[0][0][0][0], sim->n123*sizeof(float));

    zrt_new = 0;
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    if(fwi->family==1){//vp-rho
	      if(fwi->idxpar[ipar]==1) cg_rt[ipar][i3][i2][i1] = g1[i3][i2][i1]*2.0;//dJ/dln(Vp)=2*dJ/dln(kappa)
	      if(fwi->idxpar[ipar]==2) cg_rt[ipar][i3][i2][i1] = g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(rho')=dJ/dln(rho)+dJ/dln(kappa)
	    }else if(fwi->family==2){//vp-Ip
	      if(fwi->idxpar[ipar]==1) cg_rt[ipar][i3][i2][i1] = -g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(vp)=-dJ/dln(rho)+dJ/dln(kappa)
	      if(fwi->idxpar[ipar]==2) cg_rt[ipar][i3][i2][i1] = g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(Ip)=dJ/dln(rho)+dJ/dln(kappa)
	    }//end if
	    cg_z[ipar][i3][i2][i1] = cg_rt[ipar][i3][i2][i1];
	    if(fwi->preco) cg_z[ipar][i3][i2][i1] /= hess[ipar][i3][i2][i1];
	    zrt_new += cg_z[ipar][i3][i2][i1]*cg_rt[ipar][i3][i2][i1];
	  }//end for i1
	}//end for i2
      }//end for i3
    }//end for ipar
    beta = zrt_new/zrt_old;
    
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    cg_p[ipar][i3][i2][i1] = cg_z[ipar][i3][i2][i1] + beta*cg_p[ipar][i3][i2][i1];
	  }
	}
      }
    }
    zrt_old = zrt_new;
  }
  
  cpml_free(sim);
  extend_model_free(sim);
  computing_box_free(sim, 0);
  computing_box_free(sim, 1);
  fdtd_free(sim, 0);
  fdtd_free(sim, 1);
  fdtd_free(sim, 2);
  decimate_interp_free(sim, 1);

  free2float(sim->dcal);
  free2float(sim->dobs);
  free2float(sim->dres);
  free3float(g1);
  free3float(g2);
  
  free2float(cg_b);
  free2float(cg_r);
  free2float(cg_Ap);
  free4float(cg_x);
  free4float(cg_z);
  free4float(cg_p);
  free4float(cg_rt);
  if(fwi->preco) free4float(hess);

  free1int(fwi->idxpar);
  free(fwi);
}


