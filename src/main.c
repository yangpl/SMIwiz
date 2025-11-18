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
#include <mpi.h>

int iproc, nproc, ierr;

void acq_init(sim_t *sim, acq_t *acq);
void acq_free(sim_t *sim, acq_t *acq);

void read_data(sim_t *sim, acq_t *acq);
void setup_data_mask(acq_t *acq, sim_t *sim);

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
    else if(sim->mode==10) printf(" Up-down wavefield separation\n");
    printf("=====================================================\n");
  }
  
  //=========== specify parameters for simulation ============
  if(!getparint("nt", &sim->nt)) err("must give nt= "); //total number of time steps
  if(!getparfloat("dt",&sim->dt)) err("must give dt= "); //temporal sampling
  if(!getparint("nb", &sim->nb)) sim->nb = 20;   //number of layers for PML absorbing boundary
  if(!getparint("n1",&sim->n1)) err("must give n1= for FD grid");
  if(!getparint("n2",&sim->n2)) err("must give n2= for FD grid");
  if(!getparfloat("d1",&sim->d1)) err("must give d1= for FD grid"); 
  if(!getparfloat("d2",&sim->d2)) err("must give d2= for FD grid");
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
  if(!getparint("order",&sim->order)) sim->order = 4;//only accepts 4 or 8-th order FD
  sim->ri = sim->order/2;//interpolation radius of Bessel I0 function for sinc
  sim->ibox = 1;//by default, computing box should be applied
  sim->n123 = sim->n1*sim->n2*sim->n3;
  sim->n123pad = sim->n1pad*sim->n2pad*sim->n3pad;
  if(sim->mode!=0){
    if(!getparint("dr", &sim->dr)) sim->dr = (sim->n3>1)?10:1;//decimation ratio dr=5*vmax/vmin, assume vmax/vmin>=2
    sim->mt = sim->nt/sim->dr;/* decimation ratio */
    if(sim->mt*sim->dr!=sim->nt) err("nt must be multiple of dr!");
  }
  
  if(!getparint("eachopt", &sim->eachopt)) sim->eachopt = 0;//1=each shot use different source wavelet
  if(!getparint("freesurf", &sim->freesurf)) sim->freesurf = 1;// 1=free surface; 0=no freesurf
  if(!getparfloat("freq",&sim->freq)) sim->freq = 15;//reference frequency for PML
  if(iproc==0){
    printf("nt=%d (number of time steps)\n", sim->nt);
    printf("dt=%g (time step)\n", sim->dt);
    printf("nb=%d (number of boundary layers)\n", sim->nb);
    printf("[d1, d2, d3]=[%g, %g, %g]\n", sim->d1, sim->d2, sim->d3);
    printf("[n1, n2, n3]=[%d, %d, %d]\n", sim->n1, sim->n2, sim->n3);
    printf("[n1pad, n2pad, n3pad]=[%d, %d, %d]\n", sim->n1pad, sim->n2pad, sim->n3pad);
    printf("order=%d (FD order=4 or 8)\n", sim->order);
    printf("ri=%d (interpolation radius ri=order/4)\n", sim->ri);
    printf("ibox=%d (1=apply computing box; 0=no computing box)\n", sim->ibox);
    if(sim->mode!=0) printf("dr=%d (ratio for boundary decimation + interpolation)\n", sim->dr);
    printf("eachopt=%d (1=one source per shot; 0=one source for all shots)\n", sim->eachopt);
    printf("freesurf=%d (1=with free surface; 0=no free surface)\n", sim->freesurf);
  }

  //=================== specify acquisition ================
  if(!getparint("suopt", &acq->suopt)) acq->suopt = 0;//0=default, 1=for real data precessing using RTM,FWI,LSRTM
  if(!getparfloat("zmin", &acq->zmin)) acq->zmin = 0;
  if(!getparfloat("zmax", &acq->zmax)) acq->zmax = acq->zmin+(sim->n1-1)*sim->d1;
  if(!getparfloat("xmin", &acq->xmin)) acq->xmin = 0;
  if(!getparfloat("xmax", &acq->xmax)) acq->xmax = acq->xmin+(sim->n2-1)*sim->d2;
  if(!getparfloat("ymin", &acq->ymin)) acq->ymin = 0;
  if(!getparfloat("ymax", &acq->ymax)) acq->ymax = acq->ymin+(sim->n3-1)*sim->d3;
  if(iproc==0){
    printf("-------- input model range -----------\n");
    printf("[zmin, zmax]=[%g, %g]\n", acq->zmin, acq->zmax);
    printf("[xmin, xmax]=[%g, %g]\n", acq->xmin, acq->xmax);
    if(sim->n3>1) printf("[ymin, ymax]=[%g, %g]\n", acq->ymin, acq->ymax);
  }
  acq->shot_idx = alloc1int(nproc);
  int nsrc = countparval("shots");
  if(nsrc>0){
    if(nsrc<nproc) err("nproc > number of shot indices! ");
    getparint("shots", acq->shot_idx);/* a list of source index separated by comma */
  }
  if(nsrc==0){
    for(int j=0; j<nproc; j++) acq->shot_idx[j] = j+1;//index starts from 1
  }
  
  if(!acq->suopt) acq_init(sim, acq);//read acquisition file if suopt==0
  if(sim->mode>=1 && sim->mode<=7){//mode=1,2,3,4,5,6,7,9 requires reading data
    read_data(sim, acq);//read data in binary or SU format
    setup_data_mask(acq, sim);//the muting will be used to remove direct waves
  }
  ierr = MPI_Barrier(MPI_COMM_WORLD);

  //=================== read vp,rho,stf files =================
  if(!getparstring("vpfile",&vpfile)) err("must give vpfile= ");
  sim->vp = alloc3float(sim->n1, sim->n2, sim->n3);
  fp=fopen(vpfile, "rb");
  if(fp==NULL) err("cannot open vpfile=%s", vpfile);
  if(fread(&sim->vp[0][0][0], sizeof(float), sim->n123, fp)!=sim->n123)
    err("error reading vpfile=%s,  size unmatched", vpfile);
  fclose(fp);
  
  if(!getparstring("rhofile",&rhofile)) err("must give rhofile= ");
  sim->rho = alloc3float(sim->n1, sim->n2, sim->n3);
  fp = fopen(rhofile, "rb");
  if(fp==NULL) err("cannot open rhofile=%s", rhofile);
  if(fread(&sim->rho[0][0][0], sizeof(float), sim->n123, fp)!=sim->n123)
    err("error reading rhofile=%s,  size unmatched", rhofile);
  fclose(fp);

  sim->stf = alloc1float(sim->nt); //source wavelet
  if(sim->mode!=5){
    if(!getparstring("stffile",&stffile)) err("must give stffile= ");
    if(sim->eachopt){//read each source wavelet for every shot
      char number[sizeof("0000")];
      char fname[10];
      sprintf(number, "%04d", acq->shot_idx[iproc]);
      //sources will be named stf_0001, stf_0002, ..., where stffile='stf'
      snprintf(fname, sizeof(fname), "%s_%s", stffile, number);
    
      fp=fopen(fname,"rb");
      if(fp==NULL) err("cannot open stffile=%s", fname);
      if(fread(sim->stf, sizeof(float), sim->nt, fp)!=sim->nt) 
	err("error reading stffile=%s,  size unmatched", fname);
      fclose(fp);
    }else{//read the same wavelet for all shots
      fp=fopen(stffile, "rb");
      if(fp==NULL) err("cannot open stffile=%s", stffile);
      if(fread(sim->stf, sizeof(float), sim->nt, fp)!=sim->nt) 
	err("error reading stffile=%s,  size unmatched", stffile);
      fclose(fp);
    }
  }

  //====================do the job here========================
  if(sim->mode==0) do_modelling(sim, acq); 
  else if(sim->mode==1) do_fwi(sim, acq);
  else if(sim->mode==2) do_rtm(sim, acq);
  else if(sim->mode==3) do_lsrtm(sim, acq);
  else if(sim->mode==4) do_fwi(sim, acq); 
  else if(sim->mode==5) do_invert_source(sim, acq);
  else if(sim->mode==6) do_adcig(sim, acq);
  else if(sim->mode==7) do_psf_hessian(sim, acq);
  else if(sim->mode==8) do_mig_decon_pcgnr(sim, acq);
  else if(sim->mode==9) do_mig_decon_fft(sim, acq);
  else if(sim->mode==10) do_updown(sim, acq);
  
  //===========================================================
  free(sim->stf);
  free3float(sim->vp);
  free3float(sim->rho);

  if(sim->mode!=8 && sim->mode!=9) acq_free(sim, acq);  
  free(sim);
  free(acq);
  
  MPI_Finalize();

  return 0;
}

