/* Image-domain least-squares RTM (ID-LSM, migration deconvolution-MD)
 *-----------------------------------------------------------------------
 * Note: Every cw1, cw2, cw3 in z, x and y directions are taken as the 
 *       sampling point!
 * 1. In PSF approach, we wish cw1,cw2,cw3 as small as possible to remove 
 *    checkbord effect;
 * 2. In iterative approach using PSF, we wish cw1>=nw1/2, cw2>=nw2/2 and 
 *    cw3>=nw3/2 to avoid PSF window takes energy from another point scatter.
 *------------------------------------------------------------------------
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
#include "opt.h"
#include <fftw3.h>

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

void cal_diagonal_hessian(sim_t *sim, fwi_t *fwi, float ****psf, float *hess);
void matmul_Hv(sim_t *sim, fwi_t *fwi, float ****psf, int adj, double *x, double *y);

//migration deconvolution
void do_psf_hessian(sim_t *sim, acq_t *acq)
{
  char *bathyfile;
  int it, irec, ipar;
  int i1, i2, i3, i1_, i2_, i3_;
  int isp1, isp2, isp3;//is scatter point or not (1=is, 0=not)
  float tmp;
  float ***g1, ***g2;
  float ****m1, ****m2, ****mr;
  char fname[sizeof("dres_0000")];
  fwi_t *fwi;
  FILE *fp;
  
  fwi = (fwi_t*)malloc(sizeof(fwi_t));
  fwi->bathy=alloc2float(sim->n2, sim->n3);
  fwi->ibathy=alloc2int(sim->n2, sim->n3);
  if(!getparstring("bathyfile",&bathyfile)){
    memset(fwi->bathy[0], 0, sim->n2*sim->n3*sizeof(float));
  }else{
    fp=fopen(bathyfile,"rb");
    if(fp==NULL) err("cannot open bathyfile=%s",bathyfile);
    if(fread(fwi->bathy[0],sizeof(float),sim->n2*sim->n3,fp)!=sim->n2*sim->n3) 
      err("error reading bathyfile=%s", bathyfile);
    fclose(fp);
  }
  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      fwi->ibathy[i3][i2]=NINT(fwi->bathy[i3][i2]/sim->d1);
    }
  }
  if(!getparint("family", &fwi->family)) fwi->family = 2;//1=vp-rho; 2=vp-Ip
  if(!getparint("npar", &fwi->npar)) fwi->npar = 1;//number of parameters
  if(!getparint("niter", &fwi->niter)) fwi->niter = 20;//number of iterations
  if(!getparint("preco", &fwi->preco)) fwi->preco = 1;//1=preconditioning; 0=not
  fwi->idxpar = alloc1int(fwi->npar);
  if(!getparint("idxpar", fwi->idxpar)) err("must provide idxpar=");
  /* //index of parameter in each family
     family 1: idxpar[0] =1, vp; idxpar[1]=2, rho
     family 2: idxpar[0] =1, vp; idxpar[1]=2, ip
  */
  if(!getparint("mdopt", &fwi->mdopt)) fwi->mdopt = 1;//1=iterative with PSF Hessian; 2=Weiner filtering
  if(!getparint("cw1", &sim->cw1)) sim->cw1 = 25;//select every cw1 points
  if(!getparint("cw2", &sim->cw2)) sim->cw2 = 25;//select every cw2 points
  sim->cw3 = 1;
  if(sim->n3>1){
    if(!getparint("cw3", &sim->cw3)) sim->cw3 = sim->nw3/2;//select every cw3 points
  }  
  if(iproc==0){
    printf("mdopt=%d (1=iterative with PSF; 2=Wiener filtering)\n", fwi->mdopt);
    printf("family=%d\n", fwi->family);
    printf("npar=%d\n", fwi->npar);
    if(fwi->mdopt==1) printf("niter=%d\n", fwi->niter);
    printf("preco=%d (1=preconditioning; 0=not)\n", fwi->preco);
    printf("[cw1, cw2, cw3]=[%d, %d, %d]\n", sim->cw1, sim->cw2, sim->cw3);
  }
  
  fwi->n = fwi->npar*sim->n123;//number of unknowns
  sim->dcal = alloc2float(sim->nt, acq->nrec);//synthetic data
  sim->dobs = alloc2float(sim->nt, acq->nrec);
  sim->dres = alloc2float(sim->nt, acq->nrec);
  g1 = alloc3float(sim->n1, sim->n2, sim->n3);
  g2 = alloc3float(sim->n1, sim->n2, sim->n3);
  mr = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  m1 = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  m2 = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  
  read_data(sim, acq);
  setup_data_weight(acq, sim);//the muting will be used to remove direct waves
  if(sim->muteopt==0) err("RTM must assign muteopt=1");
  
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

  if(iproc==0) printf("----0. Background field modelling --------\n");
  fdtd_null(sim, 1);//flag=1, incident field
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%sim->nt_verb==0) printf("it-----%d\n", it);

    decimate_interp_bndr(sim, 1, it, 0, sim->face1, sim->face2, sim->face3);//interp=0
    fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);//direct + diving wave
  }
  tmp = 0;
  for(irec=0; irec<acq->nrec; irec++){
    for(it=0; it<sim->nt; it++){
      sim->dres[irec][it] = sim->dobs[irec][it] - sim->dcal[irec][it];//remove direct + diving wave
      sim->dres[irec][it] *= acq->wdat[irec][it];//apply muting here
      tmp += sim->dres[irec][it]*sim->dres[irec][it];
      sim->dres[irec][it] *= acq->wdat[irec][it];//weighted residual of reflections
    }
  }
  tmp *= 0.5*sim->dt;
  MPI_Allreduce(&tmp, &fwi->fcost, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  if(iproc==0) printf("fcost=%g\n", fwi->fcost);

  sprintf(fname, "dres_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"error opening file\n"); exit(1);}
  fwrite(&sim->dres[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);

  sprintf(fname, "d0_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"error opening file\n"); exit(1);}
  fwrite(&sim->dcal[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);
  
  if(iproc==0) printf("----1. migration: b=L^H (delta d)--------\n");
  memset(&g1[0][0][0], 0, sim->n123*sizeof(float));
  memset(&g2[0][0][0], 0, sim->n123*sizeof(float));
  fdtd_null(sim, 2);//flag=2, adjoint field
  sim->sign_dt = -1;
  for(it=sim->nt-1; it>=0; it--){
    if(iproc==0 && it%sim->nt_verb==0) printf("it-----%d\n", it);

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
  	tmp = -sim->volume*sim->dt;//kappa
  	g1[i3][i2][i1] *= tmp;//dJ/dln(kappa)
  	g2[i3][i2][i1] *= sim->rho[i3][i2][i1]*0.25*tmp;//dJ/dln(rho)
      }
    }
  }
  MPI_Allreduce(&g1[0][0][0], &mr[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(&g1[0][0][0], &mr[0][0][0][0], sim->n123*sizeof(float));
  MPI_Allreduce(&g2[0][0][0], &mr[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(&g2[0][0][0], &mr[0][0][0][0], sim->n123*sizeof(float));
  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
  	for(i1=0; i1<sim->n1; i1++){
  	  if(fwi->family==1){//vp-rho
  	    if(fwi->idxpar[ipar]==1) mr[ipar][i3][i2][i1] = g1[i3][i2][i1]*2.0;//dJ/dln(Vp)=2*dJ/dln(kappa)
  	    if(fwi->idxpar[ipar]==2) mr[ipar][i3][i2][i1] = g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(rho')=dJ/dln(rho)+dJ/dln(kappa)
  	  }else if(fwi->family==2){//vp-Ip
  	    if(fwi->idxpar[ipar]==1) mr[ipar][i3][i2][i1] = -g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(vp)=-dJ/dln(rho)+dJ/dln(kappa)
  	    if(fwi->idxpar[ipar]==2) mr[ipar][i3][i2][i1] = g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(Ip)=dJ/dln(rho)+dJ/dln(kappa)
  	  }//end if
	  if(i1<=fwi->ibathy[i3][i2]) mr[ipar][i3][i2][i1] = 0;//mute image above bathymetry
  	}//end for i1
      }//end for i2
    }//end for i3
  }//end for ipar
  if(iproc==0){
    fp = fopen("param_final_rtm", "wb");
    fwrite(&mr[0][0][0][0], fwi->n*sizeof(float), 1, fp);
    fclose(fp);
  }
    
  if(iproc==0) printf("------2. Born modelling/demigration, dcal=L m1--------\n");
  memset(&m1[0][0][0][0], 0, fwi->n*sizeof(float));
  memset(&m2[0][0][0][0], 0, fwi->n*sizeof(float));
  memset(&g1[0][0][0], 0, sim->n123*sizeof(float));
  memset(&g2[0][0][0], 0, sim->n123*sizeof(float));
  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      isp3 = (sim->n3==1 || i3%sim->cw3==0)?1:0;//1=is a sampled point in axis 3
      for(i2=0; i2<sim->n2; i2++){
	isp2 = (i2%sim->cw2==0)?1:0;//1=is a sampled point in axis 2
	for(i1=0; i1<sim->n1; i1++){
	  isp1 = (i1%sim->cw1==0)?1:0;//1=is a sampled point in axis 1
	  //only consider scatters below bathymetry
	  if(i1>fwi->ibathy[i3][i2]){
	    if(fwi->mdopt==1){
	      if(fwi->idxpar[ipar]==1) g1[i3][i2][i1] = (isp1 && isp2 && isp3)?1:0;
	      if(fwi->idxpar[ipar]==2) g2[i3][i2][i1] = (isp1 && isp2 && isp3)?1:0;
	    }else if(fwi->mdopt==2){
	      if(fwi->idxpar[ipar]==1) g1[i3][i2][i1] = mr[ipar][i3][i2][i1];
	      if(fwi->idxpar[ipar]==2) g2[i3][i2][i1] = mr[ipar][i3][i2][i1];
	    }
	  }
	}
      }
    }
    if(fwi->idxpar[ipar]==1) memcpy(&m1[ipar][0][0][0], &g1[0][0][0], sim->n123*sizeof(float));
    if(fwi->idxpar[ipar]==2) memcpy(&m1[ipar][0][0][0], &g2[0][0][0], sim->n123*sizeof(float));
  }//end for ipar
  if(iproc==0){
    fp = fopen("param_final_m1", "wb");
    fwrite(&m1[0][0][0][0], fwi->n*sizeof(float), 1, fp);
    fclose(fp);
  }
  fdtd_null(sim, 0);//flag=0, scattering field
  fdtd_null(sim, 1);//flag=1, incident field
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%sim->nt_verb==0) printf("it-----%d\n", it);

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
    extract_wavefield(sim, acq, sim->p0, sim->dcal, it);//reflection data, stored in dcal=L*dm
  }

  //==================================================================
  if(iproc==0) printf("----3. migration: m2=L^H dcal=L^H L m1--------\n");
  memset(&g1[0][0][0], 0, sim->n123*sizeof(float));
  memset(&g2[0][0][0], 0, sim->n123*sizeof(float));
  fdtd_null(sim, 2);//flag=2, adjoint field
  sim->sign_dt = -1;
  for(it=sim->nt-1; it>=0; it--){
    if(iproc==0 && it%sim->nt_verb==0) printf("it-----%d\n", it);

    //adjoint field modelling
    inject_adjoint_source(sim, acq, sim->p2, sim->dcal, it);
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
  	tmp = -sim->volume*sim->dt;//kappa
  	g1[i3][i2][i1] *= tmp;//dJ/dln(kappa)
  	g2[i3][i2][i1] *= sim->rho[i3][i2][i1]*0.25*tmp;//dJ/dln(rho)
      }
    }
  }
  MPI_Allreduce(&g1[0][0][0], &m2[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(&g1[0][0][0], &m2[0][0][0][0], sim->n123*sizeof(float));
  MPI_Allreduce(&g2[0][0][0], &m2[0][0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(&g2[0][0][0], &m2[0][0][0][0], sim->n123*sizeof(float));
  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
  	for(i1=0; i1<sim->n1; i1++){
  	  if(fwi->family==1){//vp-rho
  	    if(fwi->idxpar[ipar]==1) m2[ipar][i3][i2][i1] = g1[i3][i2][i1]*2.0;//dJ/dln(Vp)=2*dJ/dln(kappa)
  	    if(fwi->idxpar[ipar]==2) m2[ipar][i3][i2][i1] = g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(rho')=dJ/dln(rho)+dJ/dln(kappa)
  	  }else if(fwi->family==2){//vp-Ip
  	    if(fwi->idxpar[ipar]==1) m2[ipar][i3][i2][i1] = -g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(vp)=-dJ/dln(rho)+dJ/dln(kappa)
  	    if(fwi->idxpar[ipar]==2) m2[ipar][i3][i2][i1] = g2[i3][i2][i1] + g1[i3][i2][i1];//dJ/dln(Ip)=dJ/dln(rho)+dJ/dln(kappa)
  	  }//end if
  	}//end for i1
      }//end for i2
    }//end for i3
  }//end for ipar
  if(iproc==0){
    fp = fopen("param_final_m2", "wb");
    fwrite(&m2[0][0][0][0], fwi->n*sizeof(float), 1, fp);
    fclose(fp);
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
  free4float(mr);
  free4float(m1);
  free4float(m2);

  free1int(fwi->idxpar);
  free(fwi);  
}

//FFT-based migration deconvolution
void do_mig_decon_fft(sim_t *sim, acq_t *acq)
{
  char *bathyfile;
  int ipar, i1, i2, i3, i1_, i2_, i3_;
  int isp1, isp2, isp3;//is scatter point or not (1=is, 0=not)
  float tmp, maxval;
  float ****m1, ****m2, ****mr, ****md, ****valence;
  fwi_t *fwi;
  FILE *fp;
  
  fwi = (fwi_t*)malloc(sizeof(fwi_t));
  fwi->bathy=alloc2float(sim->n2, sim->n3);
  fwi->ibathy=alloc2int(sim->n2, sim->n3);
  if(!getparstring("bathyfile",&bathyfile)){
    memset(fwi->bathy[0], 0, sim->n2*sim->n3*sizeof(float));
  }else{
    fp=fopen(bathyfile,"rb");
    if(fp==NULL) err("cannot open bathyfile=%s",bathyfile);
    if(fread(fwi->bathy[0],sizeof(float),sim->n2*sim->n3,fp)!=sim->n2*sim->n3) 
      err("error reading bathyfile=%s", bathyfile);
    fclose(fp);
  }
  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      fwi->ibathy[i3][i2]=NINT(fwi->bathy[i3][i2]/sim->d1);
    }
  }
  if(!getparint("family", &fwi->family)) fwi->family = 2;//1=vp-rho; 2=vp-Ip
  if(!getparint("npar", &fwi->npar)) fwi->npar = 1;//number of parameters
  fwi->idxpar = alloc1int(fwi->npar);
  if(!getparint("idxpar", fwi->idxpar)) err("must provide idxpar=");
  /* //index of parameter in each family
     family 1: idxpar[0] =1, vp; idxpar[1]=2, rho
     family 2: idxpar[0] =1, vp; idxpar[1]=2, ip
  */
  if(!getparint("nw1", &sim->nw1)) sim->nw1 = 25;//window length
  if(!getparint("nw2", &sim->nw2)) sim->nw2 = 25;//window length
  if(!getparint("cw1", &sim->cw1)) sim->cw1 = 1;//select every cw1 points
  if(!getparint("cw2", &sim->cw2)) sim->cw2 = 1;//select every cw2 points
  sim->nw3 = 1;
  sim->cw3 = 1;
  if(sim->n3>1){
    if(!getparint("nw3", &sim->nw3)) sim->nw3 = sim->nw2;//window length
    if(!getparint("cw3", &sim->cw3)) sim->cw3 = 1;//select every cw3 points
  }  
  if(iproc==0){
    printf("family=%d\n", fwi->family);
    printf("npar=%d\n", fwi->npar);
    printf("[nw1, nw2, nw3]=[%d, %d, %d]\n", sim->nw1, sim->nw2, sim->nw3);
    printf("[cw1, cw2, cw3]=[%d, %d, %d]\n", sim->cw1, sim->cw2, sim->cw3);
  }

  fwi->n = fwi->npar*sim->n123;//number of unknowns
  mr = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);//m_rtm
  m1 = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  m2 = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  md = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);//m_decon
  valence = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);//valence
  
  fp = fopen("param_final_rtm", "rb");
  fread(&mr[0][0][0][0], fwi->n*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("param_final_m1", "rb");
  fread(&m1[0][0][0][0], fwi->n*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("param_final_m2", "rb");
  fread(&m2[0][0][0][0], fwi->n*sizeof(float), 1, fp);
  fclose(fp);

  //======================================
  //FFT-based migration deconvolution
  //======================================
  int j1, j2, j3, j, jj;
  int n1fft, n2fft, n3fft;
  
  n1fft = 1;
  while(n1fft<=sim->nw1) n1fft *= 2;
  n2fft = 1;
  while(n2fft<=sim->nw2) n2fft *= 2;
  n3fft = 1;
  if(sim->n3>1){
    while(n3fft<=sim->nw3) n3fft *= 2;
  }
  if(iproc==0) printf("[n1fft,n2fft,n3fft]=[%d,%d,%d]\n", n1fft, n2fft, n3fft);

  int rank = 4;
  int n[4] = {n1fft, n2fft, n3fft, fwi->npar};
  int n123fft = n1fft*n2fft*n3fft*fwi->npar;
  fftw_complex *ft_tmp = fftw_malloc(sizeof(fftw_complex)*n123fft);
  fftw_complex *ft_m1 = fftw_malloc(sizeof(fftw_complex)*n123fft);
  fftw_complex *ft_m2 = fftw_malloc(sizeof(fftw_complex)*n123fft);
  fftw_plan fft = fftw_plan_dft(rank, n, ft_tmp, ft_tmp, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan ifft = fftw_plan_dft(rank, n, ft_tmp, ft_tmp, FFTW_BACKWARD, FFTW_MEASURE);

  float *win1, *win2, *win3, *win;
  win1 = alloc1float(sim->nw1);
  win2 = alloc1float(sim->nw2);
  win3 = alloc1float(sim->nw3);
  win = alloc1float(sim->nw1*sim->nw2*sim->nw3);

  for(j1=0; j1<sim->nw1; j1++){
    i1 = j1-sim->nw1/2;
    win1[j1] = cos(PI*i1/sim->nw1);
  }
  for(j2=0; j2<sim->nw2; j2++){
    i2 = j2-sim->nw2/2;
    win2[j2] = cos(PI*i2/sim->nw2);
  }
  win3[0] = 1.;
  if(sim->n3>1){
    for(j3=0; j3<sim->nw3; j3++){
      i3 = j3-sim->nw3/3;
      win3[j3] = cos(PI*i3/sim->nw3);
    }
  }
  for(j3=0; j3<sim->nw3; j3++){
    for(j2=0; j2<sim->nw2; j2++){
      for(j1=0; j1<sim->nw1; j1++){
	jj = j1 + sim->nw1*(j2 + sim->nw2*j3);
	win[jj] = win1[j1]*win2[j2]*win3[j3];
      }
    }
  }
  free1float(win1);
  free1float(win2);
  free1float(win3);
  
  memset(&md[0][0][0][0], 0, fwi->n*sizeof(float));
  memset(&valence[0][0][0][0], 0, fwi->n*sizeof(float));
  for(i3=0; i3<sim->n3; i3++){
    isp3 = (sim->n3==1 || (i3%sim->cw3==0||i3==sim->n3-1))?1:0;//1=is a sampled point in axis 3
    for(i2=0; i2<sim->n2; i2++){
      isp2 = (i2%sim->cw2==0||i2==sim->n2-1)?1:0;//1=is a sampled point in axis 2
      for(i1=0; i1<sim->n1; i1++){
	isp1 = (i1%sim->cw1==0||i1==sim->n1-1)?1:0;//1=is a sampled point in axis 1

	if(isp1 && isp2 && isp3){
	  //-------------------------------------------------	
	  memset(ft_tmp, 0, n123fft*sizeof(fftw_complex));
	  for(ipar=0; ipar<fwi->npar; ipar++){
	    for(j3=0; j3<sim->nw3; j3++){
	      i3_ = (sim->n3>1)?MIN(MAX(i3+j3-sim->nw3/2, 0), sim->n3-1):0;
	      for(j2=0; j2<sim->nw2; j2++){
		i2_ = MIN(MAX(i2+j2-sim->nw2/2, 0), sim->n2-1);
		for(j1=0; j1<sim->nw1; j1++){
		  i1_ = MIN(MAX(i1+j1-sim->nw1/2, 0), sim->n1-1);

		  j = j1 + n1fft*(j2 + n2fft*(j3 + n3fft*ipar));
		  ft_tmp[j] = m1[ipar][i3_][i2_][i1_];//take an image block
		}
	      }
	    }
	  }//end for ipar
	  fftw_execute(fft);
	  memcpy(ft_m1, ft_tmp, n123fft*sizeof(fftw_complex));

	  memset(ft_tmp, 0, n123fft*sizeof(fftw_complex));
	  for(ipar=0; ipar<fwi->npar; ipar++){
	    for(j3=0; j3<sim->nw3; j3++){
	      i3_ = (sim->n3>1)?MIN(MAX(i3+j3-sim->nw3/2, 0), sim->n3-1):0;
	      for(j2=0; j2<sim->nw2; j2++){
		i2_ = MIN(MAX(i2+j2-sim->nw2/2, 0), sim->n2-1);
		for(j1=0; j1<sim->nw1; j1++){
		  i1_ = MIN(MAX(i1+j1-sim->nw1/2, 0), sim->n1-1);

		  j = j1 + n1fft*(j2 + n2fft*(j3 + n3fft*ipar));
		  ft_tmp[j] = m2[ipar][i3_][i2_][i1_];//take an image block
		}
	      }
	    }
	  }
	  fftw_execute(fft);
	  memcpy(ft_m2, ft_tmp, n123fft*sizeof(fftw_complex));

	  maxval = 0;
	  for(ipar=0; ipar<fwi->npar; ipar++){
	    for(j3=0; j3<n3fft; j3++){
	      for(j2=0; j2<n2fft; j2++){
		for(j1=0; j1<n1fft; j1++){
		  j = j1 + n1fft*(j2 + n2fft*(j3 + n3fft*ipar));
		  tmp = creal(conj(ft_m2[j])*ft_m2[j]);
		  maxval = MAX(maxval, tmp);
		}
	      }
	    }
	  }
	  memset(ft_tmp, 0, n123fft*sizeof(fftw_complex));
	  for(ipar=0; ipar<fwi->npar; ipar++){
	    for(j3=0; j3<sim->nw3; j3++){
	      i3_ = (sim->n3>1)?MIN(MAX(i3+j3-sim->nw3/2, 0), sim->n3-1):0;
	      for(j2=0; j2<sim->nw2; j2++){
		i2_ = MIN(MAX(i2+j2-sim->nw2/2, 0), sim->n2-1);
		for(j1=0; j1<sim->nw1; j1++){
		  i1_ = MIN(MAX(i1+j1-sim->nw1/2, 0), sim->n1-1);

		  j = j1 + n1fft*(j2 + n2fft*(j3 + n3fft*ipar));
		  ft_tmp[j] = mr[ipar][i3_][i2_][i1_];//take an image block
		}
	      }
	    }
	  }
	  fftw_execute(fft);//F[m1]
	  for(ipar=0; ipar<fwi->npar; ipar++){
	    for(j3=0; j3<sim->nw3; j3++){
	      i3_ = (sim->n3>1)?MIN(MAX(i3+j3-sim->nw3/2, 0), sim->n3-1):0;
	      for(j2=0; j2<sim->nw2; j2++){
		i2_ = MIN(MAX(i2+j2-sim->nw2/2, 0), sim->n2-1);
		for(j1=0; j1<sim->nw1; j1++){
		  i1_ = MIN(MAX(i1+j1-sim->nw1/2, 0), sim->n1-1);

		  j = j1 + n1fft*(j2 + n2fft*(j3 + n3fft*ipar));
		  //filtering in frequency/wavenumber domain
		  ft_tmp[j] *= ft_m1[j]*conj(ft_m2[j])/(conj(ft_m2[j])*ft_m2[j] + 1e-3*maxval);
		}
	      }
	    }
	  }
	  fftw_execute(ifft);

	  for(ipar=0; ipar<fwi->npar; ipar++){
	    for(j3=0; j3<sim->nw3; j3++){
	      i3_ = (sim->n3>1)?MIN(MAX(i3+j3-sim->nw3/2, 0), sim->n3-1):0;
	      for(j2=0; j2<sim->nw2; j2++){
		i2_ = MIN(MAX(i2+j2-sim->nw2/2, 0), sim->n2-1);
		for(j1=0; j1<sim->nw1; j1++){
		  i1_ = MIN(MAX(i1+j1-sim->nw1/2, 0), sim->n1-1);

		  j = j1 + n1fft*(j2 + n2fft*(j3 + n3fft*ipar));
		  jj = j1 + sim->nw1*(j2 + sim->nw2*j3);
		  md[ipar][i3_][i2_][i1_] += win[jj]*creal(ft_tmp[j])/n123fft;
		  valence[ipar][i3_][i2_][i1_] += win[jj];//sum over weights
		}//end for j1
	      }//end for j2
	    }//end for j3
	  }//end for ipar
	}//end if
	//-------------------------------------------------	
      }//end for i1
    }//end for i2
  }//end for i3
  fftw_free(ft_tmp);
  fftw_free(ft_m1);
  fftw_free(ft_m2);
  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);

  if(iproc==0){
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    md[ipar][i3][i2][i1] /= valence[ipar][i3][i2][i1];//normalize to remove checkbord effect
	    if(i1<fwi->ibathy[i3][i2]) md[ipar][i3][i2][i1] = 0.;//mute reflectivity above bathymetry
	  }
	}
      }
    }
    fp = fopen("param_final_decon", "wb");
    fwrite(&md[0][0][0][0], fwi->n*sizeof(float), 1, fp);
    fclose(fp);
  }
  
  free4float(mr);
  free4float(m1);
  free4float(m2);
  free4float(md);
  free1float(win);
  
  free1int(fwi->idxpar);
  free(fwi);
}  

//migration deconvolution using PCGNR iterative method and PSFs
void do_mig_decon_pcgnr(sim_t *sim, acq_t *acq)
{
  char *bathyfile;
  int ipar, i, i1, i2, i3;
  float maxval;
  float ****m1, ****m2, ****mr, ****md;
  fwi_t *fwi;
  FILE *fp;
  
  fwi = (fwi_t*)malloc(sizeof(fwi_t));
  fwi->bathy=alloc2float(sim->n2, sim->n3);
  fwi->ibathy=alloc2int(sim->n2, sim->n3);
  if(!getparstring("bathyfile",&bathyfile)){
    memset(fwi->bathy[0], 0, sim->n2*sim->n3*sizeof(float));
  }else{
    fp=fopen(bathyfile,"rb");
    if(fp==NULL) err("cannot open bathyfile=%s",bathyfile);
    if(fread(fwi->bathy[0],sizeof(float),sim->n2*sim->n3,fp)!=sim->n2*sim->n3) 
      err("error reading bathyfile=%s", bathyfile);
    fclose(fp);
  }
  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      fwi->ibathy[i3][i2]=NINT(fwi->bathy[i3][i2]/sim->d1);
    }
  }
  if(!getparint("family", &fwi->family)) fwi->family = 2;//1=vp-rho; 2=vp-Ip
  if(!getparint("npar", &fwi->npar)) fwi->npar = 1;//number of parameters
  if(!getparint("niter", &fwi->niter)) fwi->niter = 20;//number of iterations
  if(!getparint("preco", &fwi->preco)) fwi->preco = 1;//1=preconditioning; 0=not
  fwi->idxpar = alloc1int(fwi->npar);
  if(!getparint("idxpar", fwi->idxpar)) err("must provide idxpar=");
  /* //index of parameter in each family
     family 1: idxpar[0] =1, vp; idxpar[1]=2, rho
     family 2: idxpar[0] =1, vp; idxpar[1]=2, ip
  */
  if(!getparint("mdopt", &fwi->mdopt)) fwi->mdopt = 1;//1=iterative with PSF Hessian; 2=Weiner filtering
  if(!getparint("cw1", &sim->cw1)) sim->cw1 = 25;//select every cw1 points
  if(!getparint("cw2", &sim->cw2)) sim->cw2 = 25;//select every cw2 points
  sim->cw3 = 1;
  if(sim->n3>1){
    if(!getparint("cw3", &sim->cw3)) sim->cw3 = sim->cw2;//select every cw3 points
  }  
  if(iproc==0){
    printf("mdopt=%d (1=iterative with PSF; 2=Wiener filtering)\n", fwi->mdopt);
    printf("family=%d\n", fwi->family);
    printf("npar=%d\n", fwi->npar);
    if(fwi->mdopt==1) printf("niter=%d\n", fwi->niter);
    printf("preco=%d (1=preconditioning; 0=not)\n", fwi->preco);
    printf("[cw1, cw2, cw3]=[%d, %d, %d]\n", sim->cw1, sim->cw2, sim->cw3);
  }
  
  fwi->n = fwi->npar*sim->n123;//number of unknowns
  mr = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  m1 = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  m2 = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  md = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);

  fp = fopen("param_final_rtm", "rb");
  fread(&mr[0][0][0][0], fwi->n*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("param_final_m1", "rb");
  fread(&m1[0][0][0][0], fwi->n*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("param_final_m2", "rb");
  fread(&m2[0][0][0][0], fwi->n*sizeof(float), 1, fp);
  fclose(fp);
    
  if(fwi->mdopt==1){
    printf("**** Migration-decon using PCGNR, preco=%d\n", fwi->preco);
    //========================================
    //solving min ||m_rtm-Hm||^2 using PCGNR
    //========================================
    double rs, rs0, ws, alpha, beta, zrt_new, zrt_old;
    float ****psf = m2;//PSF Hessian is stored in m2
    if(fwi->preco){
      fwi->hess = alloc1float(fwi->n); 
      cal_diagonal_hessian(sim, fwi, psf, fwi->hess);
      maxval = 0;
      for(i=0; i<fwi->n; i++) maxval = MAX(fwi->hess[i], maxval);
      for(i=0; i<fwi->n; i++) fwi->hess[i] += 1e-3*maxval;
    }

    double *x = alloc1double(fwi->n);
    double *r = alloc1double(fwi->n);
    double *rt = alloc1double(fwi->n);
    double *z = alloc1double(fwi->n);
    double *p = alloc1double(fwi->n);
    double *Ap = alloc1double(fwi->n);//A=H, Hessian matrix

    memset(x, 0, fwi->n*sizeof(double));
    rs = 0.;
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    i = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	    r[i] = mr[ipar][i3][i2][i1];//r=b-Ax
	    rs += r[i]*r[i];
	  }
	}
      }
    }

    matmul_Hv(sim, fwi, psf, 1, rt, r);//adj=1, z=A^t(b-Ax)
    memcpy(z, rt, fwi->n*sizeof(double));
    if(fwi->preco){
      for(i=0; i<fwi->n; i++) z[i] /= fwi->hess[i];
    }
    zrt_old = 0;
    for(i=0; i<fwi->n; i++) {
      p[i] = z[i];
      zrt_old += z[i]*rt[i];
    }
    rs0 = rs;
  
    for(fwi->iter=0; fwi->iter<fwi->niter; fwi->iter++){
      if(fwi->iter==0){
	fp = fopen("iterate.txt","w");
	fprintf(fp,"========================================\n");
	fprintf(fp,"Number of PCGNR iterations: %d\n", fwi->niter);
	fprintf(fp,"fcost=0.5*||H dm - m_rtm||^2\n");
	fprintf(fp,"preco=%d \n", fwi->preco);
	fprintf(fp,"========================================\n");
	fprintf(fp,"iter     fcost\n");
	fclose(fp);
      }
      fp=fopen("iterate.txt","a");
      fprintf(fp,"%d    %.4e\n", fwi->iter, rs/rs0);
      fclose(fp);
      printf("iter=%d, fcost=%.4e \n", fwi->iter, rs/rs0);
    
      matmul_Hv(sim, fwi, psf, 0, p, Ap);//adj=0, Ap
      ws = 0;
      for(i=0; i<fwi->n; i++) ws += Ap[i]*Ap[i];    
      alpha = zrt_old/ws;
      rs = 0;
      for(i=0; i<fwi->n; i++){
	x[i] += alpha*p[i];
	r[i] -= alpha*Ap[i];
	rs += r[i]*r[i];
      }
      matmul_Hv(sim, fwi, psf, 1, rt, r);//adj=1, z=A^t r
      memcpy(z, rt, fwi->n*sizeof(double));
      if(fwi->preco){
	for(i=0; i<fwi->n; i++) z[i] /= fwi->hess[i];
      }
      zrt_new = 0.;
      for(i=0; i<fwi->n; i++) zrt_new += z[i]*rt[i];
      beta = zrt_new/zrt_old;
      for(i=0; i<fwi->n; i++) p[i] = z[i] + beta*p[i];
      zrt_old = zrt_new;

      if(iproc==0){//only 1 processor is needed in the following
	for(ipar=0; ipar<fwi->npar; ipar++){
	  for(i3=0; i3<sim->n3; i3++){
	    for(i2=0; i2<sim->n2; i2++){
	      for(i1=0; i1<sim->n1; i1++){
		i = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
		md[ipar][i3][i2][i1] = (i1>fwi->ibathy[i3][i2])?x[i]:0;//x is double while mm is float
	      }//end for i1
	    }//end for i2
	  }//end for i3
	}//end for ipar
	fp = fopen("param_final_decon", "wb");
	fwrite(&md[0][0][0][0], fwi->n*sizeof(float), 1, fp);
	fclose(fp);
      }//end if

    }

    free1double(x);
    free1double(r);
    free1double(rt);
    free1double(z);
    free1double(p);
    free1double(Ap);
    if(fwi->preco) free1float(fwi->hess);
  }//end for
  free4float(mr);
  free4float(m1);
  free4float(m2);
  free4float(md);

  free1int(fwi->idxpar);
  free(fwi);
}

double quantile(int q, int n, double* a);

//migration deconvolution using steepest descent with L1 regularization
void do_mig_decon_l1reg(sim_t *sim, acq_t *acq)
{
  char *bathyfile;
  int ipar, i, i1, i2, i3;
  float perc;
  double threshold;
  float ****m1, ****m2, ****mr, ****md;
  fwi_t *fwi;
  FILE *fp;
  
  fwi = (fwi_t*)malloc(sizeof(fwi_t));
  fwi->bathy=alloc2float(sim->n2, sim->n3);
  fwi->ibathy=alloc2int(sim->n2, sim->n3);
  if(!getparstring("bathyfile",&bathyfile)){
    memset(fwi->bathy[0], 0, sim->n2*sim->n3*sizeof(float));
  }else{
    fp=fopen(bathyfile,"rb");
    if(fp==NULL) err("cannot open bathyfile=%s",bathyfile);
    if(fread(fwi->bathy[0],sizeof(float),sim->n2*sim->n3,fp)!=sim->n2*sim->n3) 
      err("error reading bathyfile=%s", bathyfile);
    fclose(fp);
  }
  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      fwi->ibathy[i3][i2]=NINT(fwi->bathy[i3][i2]/sim->d1);
    }
  }
  if(!getparint("family", &fwi->family)) fwi->family = 2;//1=vp-rho; 2=vp-Ip
  if(!getparint("npar", &fwi->npar)) fwi->npar = 1;//number of parameters
  if(!getparint("niter", &fwi->niter)) fwi->niter = 20;//number of iterations
  if(!getparint("preco", &fwi->preco)) fwi->preco = 1;//1=preconditioning; 0=not
  fwi->idxpar = alloc1int(fwi->npar);
  if(!getparint("idxpar", fwi->idxpar)) err("must provide idxpar=");
  if(!getparfloat("perc", &perc)) perc = 0.99;
  /* //index of parameter in each family
     family 1: idxpar[0] =1, vp; idxpar[1]=2, rho
     family 2: idxpar[0] =1, vp; idxpar[1]=2, ip
  */
  if(!getparint("mdopt", &fwi->mdopt)) fwi->mdopt = 1;//1=iterative with PSF Hessian; 2=Weiner filtering
  if(!getparint("cw1", &sim->cw1)) sim->cw1 = 25;//select every cw1 points
  if(!getparint("cw2", &sim->cw2)) sim->cw2 = 25;//select every cw2 points
  sim->cw3 = 1;
  if(sim->n3>1){
    if(!getparint("cw3", &sim->cw3)) sim->cw3 = sim->cw2;//select every cw3 points
  }  
  if(iproc==0){
    printf("mdopt=%d (1=iterative with PSF; 2=Wiener filtering)\n", fwi->mdopt);
    printf("family=%d\n", fwi->family);
    printf("npar=%d\n", fwi->npar);
    if(fwi->mdopt==1) printf("niter=%d\n", fwi->niter);
    printf("preco=%d (1=preconditioning; 0=not)\n", fwi->preco);
    printf("[cw1, cw2, cw3]=[%d, %d, %d]\n", sim->cw1, sim->cw2, sim->cw3);
  }
  
  fwi->n = fwi->npar*sim->n123;//number of unknowns
  mr = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  m1 = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  m2 = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);
  md = alloc4float(sim->n1, sim->n2, sim->n3, fwi->npar);

  fp = fopen("param_final_rtm", "rb");
  fread(&mr[0][0][0][0], fwi->n*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("param_final_m1", "rb");
  fread(&m1[0][0][0][0], fwi->n*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("param_final_m2", "rb");
  fread(&m2[0][0][0][0], fwi->n*sizeof(float), 1, fp);
  fclose(fp);
    
  if(fwi->mdopt==1){
    printf("**** Migration-decon using steepest descent, preco=%d\n", fwi->preco);
    //========================================
    //solving min ||m_rtm-Hm||^2 using PCGNR
    //========================================
    double rs, num, den, alpha, fcost0, freg;
    float ****psf = m2;//PSF Hessian is stored in m2

    double *b = alloc1double(fwi->n);
    double *x = alloc1double(fwi->n);
    double *r = alloc1double(fwi->n);
    double *g = alloc1double(fwi->n);
    double *Lg = alloc1double(fwi->n);

    memset(x, 0, fwi->n*sizeof(double));
    for(ipar=0; ipar<fwi->npar; ipar++){
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    i = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	    b[i] = mr[ipar][i3][i2][i1];
	  }
	}
      }
    }
  
    for(fwi->iter=0; fwi->iter<fwi->niter; fwi->iter++){    
      matmul_Hv(sim, fwi, psf, 0, x, r);//adj=0, Lg
      rs = 0;
      for(i=0; i<fwi->n; i++){
	r[i] -= b[i];
	rs += r[i]*r[i];
      }
      matmul_Hv(sim, fwi, psf, 1, g, r);//adj=1, g=L^t (Lx-b)=-L^t r
      matmul_Hv(sim, fwi, psf, 0, g, Lg);//adj=1, g=L^t (Lx-b)=-L^t r
      
      num = 0;
      den = 0;
      for(i=0; i<fwi->n; i++){
	num += g[i]*g[i];
	den += Lg[i]*Lg[i];
      }
      alpha = -num/den;
      freg = 0;
      for(i=0; i<fwi->n; i++){
	x[i] += alpha*g[i];
	Lg[i] = fabs(x[i]);//copy it into Lg[]
	freg += Lg[i];
      }

      //a soft thresholding
      threshold = quantile((int)(perc*fwi->n), fwi->n, Lg);
      for(i=0; i<fwi->n; i++){
	if(x[i]>threshold) x[i] -= threshold;
	else if(x[i]<-threshold) x[i] += threshold;
	else x[i] = 0;
      }
      fwi->fcost =  0.5*rs + threshold*freg;
      
      if(iproc==0){//only 1 processor is needed in the following
	for(ipar=0; ipar<fwi->npar; ipar++){
	  for(i3=0; i3<sim->n3; i3++){
	    for(i2=0; i2<sim->n2; i2++){
	      for(i1=0; i1<sim->n1; i1++){
		i = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
		md[ipar][i3][i2][i1] = (i1>fwi->ibathy[i3][i2])?x[i]:0;//x is double while mm is float
	      }//end for i1
	    }//end for i2
	  }//end for i3
	}//end for ipar
	fp = fopen("param_final_decon", "wb");
	fwrite(&md[0][0][0][0], fwi->n*sizeof(float), 1, fp);
	fclose(fp);
      }//end if

      if(fwi->iter==0){
	fp = fopen("iterate.txt","w");
	fprintf(fp,"========================================\n");
	fprintf(fp,"Number of PCGNR iterations: %d\n", fwi->niter);
	fprintf(fp,"fcost=0.5*||H dm - m_rtm||^2 + beta*||dm||_1\n");
	fprintf(fp,"preco=%d \n", fwi->preco);
	fprintf(fp,"========================================\n");
	fprintf(fp,"iter     fcost\n");
	fclose(fp);
	fcost0 = fwi->fcost;
      }
      fp = fopen("iterate.txt","a");
      fprintf(fp,"%d    %.4e\n", fwi->iter, fwi->fcost/fcost0);
      fclose(fp);
      printf("iter=%d, fcost=%.4e \n", fwi->iter, fwi->fcost/fcost0);
    }

    free1double(b);
    free1double(x);
    free1double(r);
    free1double(g);
    free1double(Lg);
  }//end for
  free4float(mr);
  free4float(m1);
  free4float(m2);
  free4float(md);

  free1int(fwi->idxpar);
  free(fwi);
}

//matrix vector product using PSF Hessian and its transpose (adj=0, y=Hx; adj=1; x=H^t y)
void matmul_Hv(sim_t *sim, fwi_t *fwi, float ****psf, int adj, double *x, double *y)
{
  int i1, i2, i3, ipar;
  int j1, j2, j3, jpar;
  int k1, k2, k3;
  int k1p1, k2p1, k3p1;
  int m1, m2, m3;
  int m1p1, m2p1, m3p1;
  int i1_, i2_, i3_;
  float w1, w2, w3, Hij;
  int ix, iy;
  
  if(adj) memset(x, 0, fwi->n*sizeof(double));
  else    memset(y, 0, fwi->n*sizeof(double));

  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      if(sim->n3>1){
	k3 = i3/sim->cw3;
	k3p1 = k3+1;
	w3 = (float)(i3-k3*sim->cw3)/sim->cw3;
      }else{
	k3 = 0;
	k3p1 = 0;
	w3 = 1.;
      }
      for(i2=0; i2<sim->n2; i2++){
	k2 = i2/sim->cw2;
	k2p1 = k2+1;
	w2 = (float)(i2-k2*sim->cw2)/sim->cw2;
	for(i1=0; i1<sim->n1; i1++){
	  k1 = i1/sim->cw1;
	  k1p1 = k1+1;
	  w1 = (float)(i1-k1*sim->cw1)/sim->cw1;

	  iy = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	  //-----------------------------------------------------
	  for(jpar=0; jpar<fwi->npar; jpar++){
	    for(j3=-sim->cw3/2; j3<(sim->cw3+1)/2; j3++){
	      if(sim->n3>1){
		i3_ = MIN(MAX(i3+j3, 0), sim->n3-1);
		m3 = MIN(MAX(k3*sim->cw3+j3, 0), sim->n3-1);
		m3p1 = MIN(MAX(k3p1*sim->cw3+j3, 0), sim->n3-1);
	      }else{
		i3_ = 0;
		m3 = 0;
		m3p1 = 0;
	      }
	      for(j2=-sim->cw2/2; j2<(sim->cw2+1)/2; j2++){
		i2_ = MIN(MAX(i2+j2, 0), sim->n2-1);
		m2 = MIN(MAX(k2*sim->cw2+j2, 0), sim->n2-1);
		m2p1 = MIN(MAX(k2p1*sim->cw2+j2,  0), sim->n2-1);
		for(j1=-sim->cw1/2; j1<(sim->cw1+1)/2; j1++){
		  i1_ = MIN(MAX(i1+j1, 0), sim->n1-1);
		  m1 = MIN(MAX(k1*sim->cw1+j1, 0), sim->n1-1);
		  m1p1 = MIN(MAX(k1p1*sim->cw1+j1, 0), sim->n1-1);

		  if(i1_>fwi->ibathy[i3_][i2_]){
		    ix = i1_ + sim->n1*(i2_ + sim->n2*(i3_ + sim->n3*jpar));
		    //-------------------------------------------
		    //interpolate Hessian using 8 corners
		    Hij = (1.-w1)*(1.-w2)*(1.-w3)*psf[jpar][m3][m2][m1]
		      + w1*(1.-w2)*(1.-w3)*psf[jpar][m3][m2][m1p1]
		      + (1.-w1)*w2*(1.-w3)*psf[jpar][m3][m2p1][m1]
		      + w1*w2*(1.-w3)*psf[jpar][m3][m2p1][m1p1]
		      + (1.-w1)*(1.-w2)*w3*psf[jpar][m3p1][m2][m1]
		      + w1*(1.-w2)*w3*psf[jpar][m3p1][m2][m1p1]
		      + (1.-w1)*w2*w3*psf[jpar][m3p1][m2p1][m1]
		      + w1*w2*w3*psf[jpar][m3p1][m2p1][m1p1];
		    if(ix==iy && Hij<0) Hij = 0;//in case the estimated diagonal is 0
		    /* if(adj) x[ix] += Hij*y[iy]; */
		    /* else    y[iy] += Hij*x[ix]; */
		    //I use a symmetrization of the PSF Hessian 0.5*(H+H^t)
		    if(adj){
		      x[ix] += Hij*y[iy];
		      x[iy] += Hij*y[ix];
		    }else{
		      y[iy] += Hij*x[ix];
		      y[ix] += Hij*x[iy];
		    }
		  }
		}//end for j1
	      }//end for j2
	    }//end for j3
	  }//end for jpar

	}//end for i1
      }//end for i2
    }//end for i3
  }//end for ipar

  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(i1=0; i1<sim->n1; i1++){
	  iy = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	  //symmetrization of the PSF Hessian 0.5*(H+H^t)
	  if(adj) x[iy] *= 0.5;
	  else    y[iy] *= 0.5;
	}//end for i1
      }//end for i2
    }//end for i3
  }//end for ipar

}

//calculate diagonal elements of Hessian matrix
void cal_diagonal_hessian(sim_t *sim, fwi_t *fwi, float ****psf, float *hess)
{
  int i1, i2, i3, ipar;
  int j1, j2, j3, jpar;
  int k1, k2, k3;
  int k1p1, k2p1, k3p1;
  int m1, m2, m3;
  int m1p1, m2p1, m3p1;
  int i1_, i2_, i3_;
  float w1, w2, w3, Hij;
  int ix, iy;
  
  for(ipar=0; ipar<fwi->npar; ipar++){
    for(i3=0; i3<sim->n3; i3++){
      if(sim->n3>1){
	k3 = i3/sim->cw3;
	k3p1 = k3+1;
	w3 = (float)(i3-k3*sim->cw3)/sim->cw3;
      }else{
	k3 = 0;
	k3p1 = 0;
	w3 = 1.;
      }
      for(i2=0; i2<sim->n2; i2++){
	k2 = i2/sim->cw2;
	k2p1 = k2+1;
	w2 = (float)(i2-k2*sim->cw2)/sim->cw2;
	for(i1=0; i1<sim->n1; i1++){
	  k1 = i1/sim->cw1;
	  k1p1 = k1+1;
	  w1 = (float)(i1-k1*sim->cw1)/sim->cw1;

	  iy = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	  //-----------------------------------------------------
	  for(jpar=0; jpar<fwi->npar; jpar++){
	    for(j3=-sim->cw3/2; j3<(sim->cw3+1)/2; j3++){
	      if(sim->n3>1){
		i3_ = MIN(MAX(i3+j3, 0), sim->n3-1);
		m3 = MIN(MAX(k3*sim->cw3+j3, 0), sim->n3-1);
		m3p1 = MIN(MAX(k3p1*sim->cw3+j3, 0), sim->n3-1);
	      }else{
		i3_ = 0;
		m3 = 0;
		m3p1 = 0;
	      }
	      for(j2=-sim->cw2/2; j2<(sim->cw2+1)/2; j2++){
		i2_ = MIN(MAX(i2+j2, 0), sim->n2-1);
		m2 = MIN(MAX(k2*sim->cw2+j2, 0), sim->n2-1);
		m2p1 = MIN(MAX(k2p1*sim->cw2+j2,  0), sim->n2-1);
		for(j1=-sim->cw1/2; j1<(sim->cw1+1)/2; j1++){
		  i1_ = MIN(MAX(i1+j1, 0), sim->n1-1);
		  m1 = MIN(MAX(k1*sim->cw1+j1, 0), sim->n1-1);
		  m1p1 = MIN(MAX(k1p1*sim->cw1+j1, 0), sim->n1-1);

		  if(i1_>fwi->ibathy[i3_][i2_]){
		    ix = i1_ + sim->n1*(i2_ + sim->n2*(i3_ + sim->n3*jpar));
		    //-------------------------------------------
		    if(ix==iy) {
		      //interpolate Hessian using 8 corners
		      Hij = (1.-w1)*(1.-w2)*(1.-w3)*psf[jpar][m3][m2][m1]
			+ w1*(1.-w2)*(1.-w3)*psf[jpar][m3][m2][m1p1]
			+ (1.-w1)*w2*(1.-w3)*psf[jpar][m3][m2p1][m1]
			+ w1*w2*(1.-w3)*psf[jpar][m3][m2p1][m1p1]
			+ (1.-w1)*(1.-w2)*w3*psf[jpar][m3p1][m2][m1]
			+ w1*(1.-w2)*w3*psf[jpar][m3p1][m2][m1p1]
			+ (1.-w1)*w2*w3*psf[jpar][m3p1][m2p1][m1]
			+ w1*w2*w3*psf[jpar][m3p1][m2p1][m1p1];
		      if(Hij<0) Hij = 0;
		      hess[ix] = Hij;
		    }
		  }//end if
		}//end for j1
	      }//end for j2
	    }//end for j3
	  }//end for jpar
	  
	}//end for i1
      }//end for i2
    }//end for i3
  }//end for ipar
    
}
