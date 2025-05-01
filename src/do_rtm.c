/* Reverse time migration (RTM)
 * (The generalized RTM image uses impendence kernel)
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *----------------------------------------------------------------------*/
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
void inject_adjoint_source(sim_t *sim, acq_t *acq, float ***rp, float **dres, int it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);

void do_rtm(sim_t *sim, acq_t *acq)
{
  char *bathyfile;
  int it, irec;
  int i1, i2, i3, i1_, i2_, i3_;
  float ***g1, ***g2, ***mm;
  float tmp;
  fwi_t *fwi;
  FILE *fp;
  char fname[sizeof("dres_0000")];

  fwi = (fwi_t*)malloc(sizeof(fwi_t));
  fwi->bathy = alloc2float(sim->n2, sim->n3);
  fwi->ibathy = alloc2int(sim->n2, sim->n3);
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
  fwi->family = 2;//1=vp-rho; 2=vp-Ip
  fwi->npar = 1;//number of parameters
  fwi->idxpar = alloc1int(fwi->npar);
  fwi->idxpar[0] = 2;//index of parameter in each family
  /*family 1: idxpar[0] =1, vp; idxpar[1]=2, rho
    family 2: idxpar[0] =1, vp; idxpar[1]=2, Ip  */
  
  fwi->n = fwi->npar*sim->n123;//number of unknowns
  sim->dcal = alloc2float(sim->nt, acq->nrec);//synthetic data - direct wave
  sim->dobs = alloc2float(sim->nt, acq->nrec);
  sim->dres = alloc2float(sim->nt, acq->nrec);
  g1 = alloc3float(sim->n1, sim->n2, sim->n3);
  g2 = alloc3float(sim->n1, sim->n2, sim->n3);
  mm = alloc3float(sim->n1, sim->n2, sim->n3);
  
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
  if(iproc==0) printf("----Simulate background field p1 --------\n");
  fdtd_null(sim, 1);//flag=1, incident field
  fdtd_null(sim, 0);//flag=0, scattering field
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%sim->nt_verb==0) printf("it-----%d\n", it);

    //background field modelling, p1
    decimate_interp_bndr(sim, 1, it, 0, sim->face1, sim->face2, sim->face3);//interp=0
    fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);//direct wave + diving wave
  }
  for(irec=0; irec<acq->nrec; irec++){
    for(it=0; it<sim->nt; it++){
      sim->dres[irec][it] = (double)(sim->dobs[irec][it]-sim->dcal[irec][it])*acq->wdat[irec][it];//muting direct wave should happen here
      sim->dres[irec][it] *= (double)acq->wdat[irec][it];//prepare adjoint source
    }
  }
  
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
  
  if(iproc==0) printf("----Migration for reflections--------\n");
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

    //reconstruction of the background field p1 in reverse time order
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    fdtd_update_p(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    fdtd_update_v(sim, 1, it, 0, sim->kappa, sim->buz, sim->bux, sim->buy);
    decimate_interp_bndr(sim, 1, it, 1, sim->face1, sim->face2, sim->face3);//interp=1

    //build density and bulk modulus kernels before applying chain rule to compute impedance kernel
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
	mm[i3][i2][i1] = g1[i3][i2][i1] + g2[i3][i2][i1];//dJ/dln(Ip)=dJ/dln(kappa) + dJ/dln(rho)
      }
    }
  }
  MPI_Allreduce(&mm[0][0][0], &g1[0][0][0], sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  if(iproc==0){
    fp = fopen("param_final_rtm", "wb");
    fwrite(&g1[0][0][0], fwi->n*sizeof(float), 1, fp);
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
  free3float(mm);

  free2float(fwi->bathy);
  free2int(fwi->ibathy);
  free1int(fwi->idxpar);
  free(fwi);
}
