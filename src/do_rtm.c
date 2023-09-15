/* 2D/3D seismic modelling, RTM and FWI code
 *-----------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------
 */
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
void inject_adjoint_source(sim_t *sim, acq_t *acq, float ***rp, float **dres, int it);

void laplace_filter(sim_t *sim, float ***in, float ***out);

void do_rtm(sim_t *sim, acq_t *acq)
{
  int it, irec, i1, i2, i3, i1_, i2_, i3_;
  float maxval, ***image, ***den, ***num;
  fwi_t *fwi;
  FILE *fp;
  char *bathyfile, fname[sizeof("dref_0000")];

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
  
  sim->freesurf = 0;//do not use free surface for RTM, only reflections should be imaged
  
  image = alloc3float(sim->n1, sim->n2, sim->n3);
  den = alloc3float(sim->n1, sim->n2, sim->n3);/* denominator of normalized X-correlation */
  num = alloc3float(sim->n1, sim->n2, sim->n3);/* numerator of normalized X-correlation */
  memset(&image[0][0][0], 0, sim->n123*sizeof(float));
  memset(&den[0][0][0], 0, sim->n123*sizeof(float));
  memset(&num[0][0][0], 0, sim->n123*sizeof(float));
  
  sim->dcal = alloc2float(sim->nt, acq->nrec);
  sim->dobs = alloc2float(sim->nt, acq->nrec);
  sim->dres = alloc2float(sim->nt, acq->nrec);
  read_data(sim, acq);
  setup_data_weight(acq, sim);
  
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
  if(iproc==0) printf("----stage 1: modelling in true model!\n");
  sim->sign_dt = 1;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    decimate_interp_bndr(sim, 0, it);/* interp=0 */
    fdtd_update_v(sim, 1, it, 0);
    fdtd_update_p(sim, 1, it, 0);
    inject_source(sim, acq, sim->p1, sim->stf[it]);
  }
  
  /*--------------------------------------------------------------*/
  for(irec=0; irec<acq->nrec; irec++){
    for(it=0; it<sim->nt; it++){
      sim->dres[irec][it] = sim->dobs[irec][it]*acq->wdat[irec][it];
    }
  }
  sprintf(fname, "dref_%04d", acq->shot_idx[iproc]);
  fp = fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"error opening file\n"); exit(1);}
  fwrite(&sim->dres[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fclose(fp);
  fflush(stdout);


  /*--------------------------------------------------------------*/
  if(iproc==0) printf("----stage 2: backpropagate observed data!\n");
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
    
    for(i3=0; i3<sim->n3; i3++){
      i3_ = (sim->n3>1)?i3+sim->nb:0;
      for(i2=0; i2<sim->n2; i2++){
	i2_ = i2+sim->nb;
	for(i1=0; i1<sim->n1; i1++){
	  i1_ = i1+sim->nb;
	  num[i3][i2][i1] += sim->p1[i3_][i2_][i1_]*sim->p2[i3_][i2_][i1_];
	  den[i3][i2][i1] += sim->p1[i3_][i2_][i1_]*sim->p1[i3_][i2_][i1_];
	}
      }
    }
  }
  
  maxval = 0;
  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      for(i1=0; i1<sim->n1; i1++){
	if(i1<= fwi->ibathy[i3][i2]) num[i3][i2][i1] = 0.;//mute image above seafloor
	maxval = MAX(maxval, den[i3][i2][i1]);
      }
    }
  }
  for(i3=0; i3<sim->n3; i3++){
    for(i2=0; i2<sim->n2; i2++){
      for(i1=0; i1<sim->n1; i1++){
	image[i3][i2][i1] = num[i3][i2][i1]/(den[i3][i2][i1]+1e-5*maxval);
      }
    }
  }

  MPI_Allreduce(image[0][0], sim->vp[0][0], sim->n123,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  laplace_filter(sim, sim->vp, image);
  if(iproc==0){
    fp = fopen("image_normalized_xcorr", "wb");
    fwrite(&image[0][0][0], sim->n123*sizeof(float), 1, fp);
    fclose(fp);
  }

  MPI_Allreduce(num[0][0], sim->vp[0][0], sim->n123,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  laplace_filter(sim, sim->vp, num);
  if(iproc==0){
    fp = fopen("image_xcorr", "wb");
    fwrite(&num[0][0][0], sim->n123*sizeof(float), 1, fp);
    fclose(fp);
  }
  

  
  cpml_close(sim);
  extend_model_close(sim);
  computing_box_close(sim, 0);
  computing_box_close(sim, 1);
  fdtd_close(sim, 1);
  fdtd_close(sim, 2);
  decimate_interp_close(sim);

  free2float(sim->dcal);
  free2float(sim->dobs);
  free2float(sim->dres);
  free3float(image);

  free2float(fwi->bathy);
  free2int(fwi->ibathy);
  free(fwi);
}
