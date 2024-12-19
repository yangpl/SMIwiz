/* 2D/3D seismic modelling, RTM and FWI code
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *---------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
#include "acq.h"
#include "mpi_info.h"
#include <mpi.h>

int iproc, nproc, ierr;

void acq_init(sim_t *sim, acq_t *acq);
void acq_close(sim_t *sim, acq_t *acq);

void do_updown(sim_t *sim, acq_t *acq);
void do_modelling(sim_t *sim, acq_t *acq);
void do_fwi(sim_t *sim, acq_t *acq);
void do_rtm(sim_t *sim, acq_t *acq);
void do_lsrtm(sim_t *sim, acq_t *acq);
void do_invert_source(sim_t *sim, acq_t *acq);
void do_psf_hessian(sim_t *sim, acq_t *acq);
void do_mig_decon_pcgnr(sim_t *sim, acq_t *acq);
void do_mig_decon_fft(sim_t *sim, acq_t *acq);
void do_mig_decon(sim_t *sim, acq_t *acq);
void do_adcig(sim_t *sim, acq_t *acq);


int main(int argc, char* argv[])
{
  char current_time[128];
  time_t      t;
  struct tm*  ptm;
  char *stffile, *vpfile, *rhofile;
  FILE *fp;
  acq_t *acq;
  sim_t *sim;
    
  MPI_Init(&argc,&argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
  initargs(argc, argv);  
  acq = (acq_t *)malloc(sizeof(acq_t));
  sim = (sim_t *)malloc(sizeof(sim_t));
  
  if(!getparint("mode", &sim->mode)) sim->mode=0;
  if(iproc==0){
    t = time(NULL);
    ptm = localtime(&t);
    strftime(current_time, 128, "%d-%b-%Y %H:%M:%S", ptm);
    printf("  Current date and time: %s\n", current_time);
    printf("=====================================================\n");
    printf("   Welcome to SMIwiz software: an integrated toolbox \n");
    printf("   for seismic modeling, RTM imaging, linearized and \n");
    printf("         nonlinear full waveform inversion           \n");
    printf("             Copyright (c) Pengliang Yang            \n");
    printf("       2021, Harbin Institute of Technology, China   \n");
    printf("=====================================================\n");
    if(sim->mode==0) printf(" Forward modeling \n");
    else if(sim->mode==1) printf(" FWI in the time domain \n");
    else if(sim->mode==2) printf(" RTM for reflectivity \n");
    else if(sim->mode==3) printf(" Data-domain LSRTM\n");
    else if(sim->mode==4) printf(" FWI gradient building \n");
    else if(sim->mode==5) printf(" Source inversion \n");
    else if(sim->mode==6) printf(" Extracting angle-domain common-image-gather (ADCIG)\n");
    else if(sim->mode==7) printf(" Compute PSF Hessian\n");
    else if(sim->mode==8) printf(" Iterative migration deconvolution via PSF\n");
    else if(sim->mode==9) printf(" Migration deconvolution via FFT-Wiener filter\n");
    printf("=====================================================\n");
  }
  if(!getparint("order",&sim->order)) sim->order = 4;//only accepts 4 or 8-th order FD
  if(!getparint("freesurf", &sim->freesurf)) sim->freesurf = 1;// 1=free surface; 0=no freesurf
  if(!getparint("nb", &sim->nb)) sim->nb = 20;   //number of layers for PML absorbing boundary
  if(!getparfloat("freq",&sim->freq)) sim->freq = 15;//reference frequency for PML
  if(!getparfloat("dt",&sim->dt)) err("must give dt= "); //temporal sampling
  if(!getparint("nt", &sim->nt)) err("must give nt= "); //total number of time steps
  if(!getparint("n1",&sim->n1)) err("must give n1= for FD grid");
  if(!getparint("n2",&sim->n2)) err("must give n2= for FD grid");
  if(!getparfloat("d1",&sim->d1)) err("must give d1= for FD grid"); 
  if(!getparfloat("d2",&sim->d2)) err("must give d2= for FD grid");
  if(!getparint("nt_verb", &sim->nt_verb)) sim->nt_verb = 100;//verbose display every nt_verb timesteps
  sim->n1pad = sim->n1+2*sim->nb;
  sim->n2pad = sim->n2+2*sim->nb;
  if(!getparint("n3",&sim->n3)) { //default, 2D
    sim->n3=1;
    sim->d3=1;
    sim->n3pad=1;
  }
  if(sim->n3>1) {//ny>1, 3D
    if(!getparfloat("d3",&sim->d3)) err("must give d3= for FD grid");
    sim->n3pad = sim->n3+2*sim->nb;
    sim->volume = sim->d1*sim->d2*sim->d3;
  }else  sim->volume = sim->d1*sim->d2;
  sim->n123 = sim->n1*sim->n2*sim->n3;
  sim->n123pad = sim->n1pad*sim->n2pad*sim->n3pad;
  sim->ibox = 1;//by default, computing box should be applied
  sim->ri = sim->order/2;//interpolation radius of Bessel I0 function for sinc
  if(!getparint("ps", &sim->ps)) sim->ps = 0; //1=PS decomposition; 0=no PS decomposition
  
  if(sim->mode!=0){
    if(!getparint("dr", &sim->dr)) sim->dr = 1;/* decimation ratio */
    if(sim->n3>1) sim->dr = 10;//r=5*vmax/vmin, assume vmax/vmin>=2
    sim->mt = sim->nt/sim->dr;/* decimation ratio */
    if(sim->mt*sim->dr!=sim->nt) err("nt must be multiple of dr!");
  }


  if(!getparfloat("zmin", &acq->zmin)) acq->zmin = 0;
  if(!getparfloat("zmax", &acq->zmax)) acq->zmax = acq->zmin+(sim->n1-1)*sim->d1;
  if(!getparfloat("xmin", &acq->xmin)) acq->xmin = 0;
  if(!getparfloat("xmax", &acq->xmax)) acq->xmax = acq->xmin+(sim->n2-1)*sim->d2;
  if(!getparfloat("ymin", &acq->ymin)) acq->ymin = 0;
  if(!getparfloat("ymax", &acq->ymax)) acq->ymax = acq->ymin+(sim->n3-1)*sim->d3;
  if(iproc==0){
    printf("freesurf=%d (1=with free surface; 0=no free surface)\n", sim->freesurf);
    printf("ri=%d (interpolation radius)\n", sim->ri);
    printf("dt=%g (time step)\n", sim->dt);
    printf("nt=%d (number of time steps)\n", sim->nt);
    printf("nb=%d (number of boundary layers)\n", sim->nb);
    printf("[zmin, zmax]=[%g, %g]\n", acq->zmin, acq->zmax);
    printf("[xmin, xmax]=[%g, %g]\n", acq->xmin, acq->xmax);
    printf("[ymin, ymax]=[%g, %g]\n", acq->ymin, acq->ymax);
    printf("[d1, d2, d3]=[%g, %g, %g]\n", sim->d1, sim->d2, sim->d3);
    printf("[n1, n2, n3]=[%d, %d, %d]\n", sim->n1, sim->n2, sim->n3);
    printf("[n1pad, n2pad, n3pad]=[%d, %d, %d]\n", sim->n1pad, sim->n2pad, sim->n3pad);
  }
  
  if(!getparstring("stffile",&stffile)) err("mute give stffile= ");
  if(!getparstring("vpfile",&vpfile)) err("must give vpfile= ");
  if(!getparstring("rhofile",&rhofile)) err("must give rhofile= ");
  
  sim->stf = alloc1float(sim->nt); //source wavelet
  sim->vp = alloc3float(sim->n1, sim->n2, sim->n3);
  sim->rho = alloc3float(sim->n1, sim->n2, sim->n3);

  if(sim->mode!=5){
    fp=fopen(stffile,"rb");
    if(fp==NULL) err("cannot open stffile=%s", stffile);
    if(fread(sim->stf,sizeof(float),sim->nt,fp)!=sim->nt) 
      err("error reading stffile=%s, size unmatched",stffile);
    fclose(fp);
  }
  
  fp=fopen(rhofile,"rb");
  if(fp==NULL) err("cannot open rhofile=%s", rhofile);
  if(fread(sim->rho[0][0],sizeof(float),sim->n123,fp)!=sim->n123)
    err("error reading rhofile=%s, size unmatched",rhofile);
  fclose(fp);
  
  fp=fopen(vpfile,"rb");
  if(fp==NULL) err("cannot open vpfile=%s", vpfile);
  if(fread(sim->vp[0][0],sizeof(float),sim->n123,fp)!=sim->n123)
    err("error reading vpfile=%s, size unmatched",vpfile);
  fclose(fp);

  ierr = MPI_Barrier(MPI_COMM_WORLD);
  acq_init(sim, acq);
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  if(iproc==0) printf("------------- acquisition init done ---------------\n");

  //====================do the job here========================
  if(sim->mode==0)      do_modelling(sim, acq); 
  else if(sim->mode==1) do_fwi(sim, acq);
  else if(sim->mode==2) do_rtm(sim, acq);
  else if(sim->mode==3) do_lsrtm(sim, acq);
  else if(sim->mode==4) do_fwi(sim, acq); 
  else if(sim->mode==5) do_invert_source(sim, acq);
  //else if(sim->mode==6) do_adcig(sim, acq);
  else if(sim->mode==7) do_psf_hessian(sim, acq);
  else if(sim->mode==8) do_mig_decon_pcgnr(sim, acq);
  else if(sim->mode==9) do_mig_decon_fft(sim, acq);
  else if(sim->mode==10) do_updown(sim, acq);
  //===========================================================

  ierr = MPI_Barrier(MPI_COMM_WORLD);

  free(sim->stf);
  free3float(sim->vp);
  free3float(sim->rho);

  acq_close(sim, acq);
  
  free(sim);
  free(acq);
  
  MPI_Finalize();

  return 0;
}

