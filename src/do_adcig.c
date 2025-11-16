/* Extracting angle domain common image gather (ADCIG) from RTM
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

void write_data(sim_t *sim, acq_t *acq);

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

void do_adcig(sim_t *sim, acq_t *acq)
{
  char *bathyfile;
  int it, irec;
  int i1, i2, i3, i1_, i2_, i3_;
  int ia, iw, ip, ix, iy;
  int nxa, nya;
  double tmp, a, b1, b2, s1, s2, g1, g2;
  float *xadcig, *yadcig, *win;
  float ****adcig, ****image;
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

  if(!(nxa = countparval("xadcig"))){
    nxa = sim->n2;
    xadcig = alloc1float(nxa);
    for(i2=0; i2<sim->n2; i2++) xadcig[i2] = acq->xmin + i2*sim->d2;
  }else{
    xadcig = alloc1float(nxa);
    getparfloat("xadcig", xadcig);
  }
  if(!(nya = countparval("yadcig"))){
    nya = sim->n3;
    yadcig = alloc1float(nya);  
    for(i3=0; i3<sim->n3; i3++) yadcig[i3] = acq->ymin + i3*sim->d3;
  }else{
    yadcig = alloc1float(nya);  
    getparfloat("yadcig", yadcig);
  }

  
  if(!getparint("na", &sim->na)) sim->na = 90;//imageber of discrete angles
  if(!getparint("awh", &sim->awh)) sim->awh = 1;//half window length for angle smoothing
  sim->da = 90./sim->na;//angle interval

  image = alloc4float(sim->n1, sim->n2, sim->n3, sim->na);
  adcig = alloc4float(sim->n1, sim->n2, sim->n3, sim->na);
  memset(&image[0][0][0][0], 0, sim->na*sim->n123*sizeof(float));
  memset(&adcig[0][0][0][0], 0, sim->na*sim->n123*sizeof(float));

  fwi->n = fwi->npar*sim->n123;//number of unknowns
    
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
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

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

    //build density and bulk modulus kernels before applying chain rule to compute impedance kernel
    for(iy=0; iy<nya; iy++){
      i3 = NINT((yadcig[iy]-acq->ymin)/sim->d3);
      i3_ = (sim->n3>1)?i3 + sim->nb:0;
      for(ix=0; ix<nxa; ix++){
	i2 = NINT((xadcig[ix]-acq->xmin)/sim->d2);
	i2_ = i2 + sim->nb;
	for(i1=0; i1<sim->n1; i1++){
	  i1_ = i1+sim->nb;
	  
	  g1 = sim->divv[i3_][i2_][i1_]*sim->p2[i3_][i2_][i1_];
	  g2 = (sim->dvzdt[i3_][i2_][i1_] + sim->dvzdt[i3_][i2_][i1_-1])*(sim->vz2[i3_][i2_][i1_] + sim->vz2[i3_][i2_][i1_-1]);
	  g2 += (sim->dvxdt[i3_][i2_][i1_] + sim->dvxdt[i3_][i2_-1][i1_])*(sim->vx2[i3_][i2_][i1_] + sim->vx2[i3_][i2_-1][i1_]);
	  if(sim->n3>1) g2 += (sim->dvydt[i3_][i2_][i1_] + sim->dvydt[i3_-1][i2_][i1_])*(sim->vy2[i3_][i2_][i1_] + sim->vy2[i3_-1][i2_][i1_]);
	  tmp =-sim->volume*sim->dt*(g1 + g2*sim->rho[i3][i2][i1]*0.25);//dJ/dln(Ip)=dJ/dln(kappa) + dJ/dln(rho)
	  
	  if(sim->n3==1 && i1>=fwi->ibathy[i3][i2]){
	    //source side Poynting vector Ss=(s1,s2) defined according to energy density flux
	    s1 = -sim->p1[i3_][i2_][i1_]*(sim->vz1[i3_][i2_][i1_] + sim->vz1[i3_][i2_][i1_-1]);//*0.5
	    s2 = -sim->p1[i3_][i2_][i1_]*(sim->vx1[i3_][i2_][i1_] + sim->vx1[i3_][i2_-1][i1_]);//*0.5	    
	    b1 = MAX(fabs(s1), fabs(s2));//stable division: , Ss<--Ss/|Ss|
	    if(b1>0){//Ss is not a zero vector
	      if(fabs(s1)<fabs(s2)){
		a = s1/s2;
		s2 = 1./sqrt(1. + a*a);
		s1 = a*s2;	      
	      }else{
		a = s2/s1;
		s1 = 1./sqrt(1. + a*a);
		s2 = a*s1;
	      }
	    }
	    //receiver side Poynting vector Sg=(g1,g2) defined according to energy density flux
	    g1 = -sim->p2[i3_][i2_][i1_]*(sim->vz2[i3_][i2_][i1_] + sim->vz2[i3_][i2_][i1_-1]);//*0.5
	    g2 = -sim->p2[i3_][i2_][i1_]*(sim->vx2[i3_][i2_][i1_] + sim->vx2[i3_][i2_-1][i1_]);//*0.5
	    b2 = MAX(fabs(g1), fabs(g2));
	    if(b2>0){//Sg is not a zero vector, stable division: Sg<--Sg/|Sg|
	      if(fabs(g1)<fabs(g2)){
		a = g1/g2;
		g2 = 1./sqrt(1. + a*a);
		g1 = a*g2;	      
	      }else{
		a = g2/g1;
		g1 = 1./sqrt(1. + a*a);
		g2 = a*g1;
	      }
	    }
	    if(b1>0 && b2>0){//both Ss and Sg are non-zero vectors
	      a = s1*g1 + s2*g2; //<Ss,Sg>
	      a = 0.5*acos(a);//theta
	      a = a*180./PI;//convert rad to degree
	      ia = NINT(a/sim->da);//index of angle
	      ia = MAX(ia, 0);
	      ia = MIN(ia, sim->na-1);
	      adcig[ia][i3][i2][i1] += tmp;
	    }
	  }//end if

	}
      }
    }
    
  }//end for it
  MPI_Allreduce(&adcig[0][0][0][0], &image[0][0][0][0], sim->na*sim->n123, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  memcpy(&adcig[0][0][0][0], &image[0][0][0][0], sim->na*sim->n123*sizeof(float));
  
  win = alloc1float(2*sim->awh + 1);
  s1 = 0;
  for(iw=-sim->awh; iw<=sim->awh; iw++){
    win[iw + sim->awh] = exp(-iw*iw);
    s1 += win[iw + sim->awh];
  }
  for(iw=-sim->awh; iw<=sim->awh; iw++) win[iw + sim->awh] /= s1;

  memset(&image[0][0][0][0], 0, sim->na*sim->n123*sizeof(float));
  for(ia=0; ia<sim->na; ia++){
    for(iw=-sim->awh; iw<=sim->awh; iw++){
      //Gauss-windowed smoothing in first Fresnel zone
      ip = ia + iw;
      ip = MAX(ip, 0);
      ip = MIN(ip, sim->na-1);
      for(i3=0; i3<sim->n3; i3++){
	for(i2=0; i2<sim->n2; i2++){
	  for(i1=0; i1<sim->n1; i1++){
	    image[ip][i3][i2][i1] += adcig[ia][i3][i2][i1]*win[iw];
	  }//end for i1
	}//end for i2
      }//end for i3
    }//end for iw
  }//end for ia
  if(iproc==0){
    fp = fopen("adcig", "wb");
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	for(ia=0; ia<sim->na; ia++){
	  for(i1=0; i1<sim->n1; i1++){
	    fwrite(&image[ia][i3][i2][i1], sizeof(float), 1, fp);
	  }
	}
      }
    }
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

  free1float(win);
  free4float(image);
  free4float(adcig);

  free2float(fwi->bathy);
  free2int(fwi->ibathy);
  free1int(fwi->idxpar);
  free(fwi);
}
