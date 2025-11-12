/* 2D/3D seismic modelling, RTM and FWI code
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
#include "acq.h"

double kaiser_windowed_sinc(double x, double dx, int r);

void acq_init(sim_t *sim, acq_t *acq)
/*< read acquisition file to initialize acquisition geometry >*/
{
  char *acquifile;
  int nrec_max;
  float *rx1, *rx2, *rx3;
  float zs, xs, ys;
  float zz, xx, yy, dip, azimuth,  tmp,  frac;
  float zmin, zmax, xmin, xmax, ymin, ymax;
  int isreceiver, isrc, irec, iseof, j;
  FILE *fp;
  
  if(!getparfloat("zmin", &acq->zmin)) acq->zmin = 0;
  if(!getparfloat("zmax", &acq->zmax)) acq->zmax = acq->zmin+(sim->n1-1)*sim->d1;
  if(!getparfloat("xmin", &acq->xmin)) acq->xmin = 0;
  if(!getparfloat("xmax", &acq->xmax)) acq->xmax = acq->xmin+(sim->n2-1)*sim->d2;
  if(!getparfloat("ymin", &acq->ymin)) acq->ymin = 0;
  if(!getparfloat("ymax", &acq->ymax)) acq->ymax = acq->ymin+(sim->n3-1)*sim->d3;
  if(!getparint("nrec_max", &nrec_max)) nrec_max = 100000;//maximum dimensions/receivers per shot
  if(iproc==0){
    printf("[zmin, zmax]=[%g, %g]\n", acq->zmin, acq->zmax);
    printf("[xmin, xmax]=[%g, %g]\n", acq->xmin, acq->xmax);
    if(sim->n3>1) printf("[ymin, ymax]=[%g, %g]\n", acq->ymin, acq->ymax);
    printf("nrec_max=%d \n", nrec_max);
  }
  if(!getparstring("acquifile", &acquifile)) err("must give acquifile= ");

  rx1 = alloc1float(nrec_max);
  rx2 = alloc1float(nrec_max);
  rx3 = alloc1float(nrec_max);
  
  fp = fopen(acquifile,"r");
  if(fp==NULL) err("file %s does not exist!", acquifile); 
  iseof = fscanf(fp, "%*[^\n]\n");//skip a line at the beginning of the file
  isrc = 0;
  while(1){
    iseof = fscanf(fp,"%f %f %f %f %f %d",&zz,&xx,&yy,&dip,&azimuth,&isreceiver);
    if(iseof==EOF)
      break;
    else{
      if(isreceiver==0){// a source line, origin of axes stripped outsh
	isrc++;
	if(acq->shot_idx[iproc]==isrc){
	  zs = zz;
	  xs = xx;
	  ys = yy;
	  acq->nrec = 0;//start to count number of receivers from 0
	}
      }else{
	if(acq->shot_idx[iproc]==isrc){//a receiver line associated with source-isrc
	  rx1[acq->nrec] = zz;
	  rx2[acq->nrec] = xx;
	  rx3[acq->nrec] = yy;
	  acq->nrec++;
	}
      }
    }
  }
  fclose(fp);
  acq->nsrc = 1;//by default, each process handles one shot
  printf("isrc=%d, nrec=%d, (zs,xs,ys)=(%.2f,%.2f,%.2f)\n",
	 acq->shot_idx[iproc], acq->nrec, zs, xs, ys);
  
  zmin = acq->zmin;
  zmax = acq->zmax;
  xmin = acq->xmin;
  xmax = acq->xmax;
  ymin = acq->ymin;
  ymax = acq->ymax;

  acq->src_x1 = alloc1float(acq->nsrc);
  acq->src_x2 = alloc1float(acq->nsrc);
  acq->src_x3 = alloc1float(acq->nsrc);
  acq->src_w1 = alloc2float(2*sim->ri+1, acq->nsrc);//positive src depth below free surface
  acq->src_w1m = alloc2float(2*sim->ri+1, acq->nsrc);//negative src depth above free surface
  acq->src_w2 = alloc2float(2*sim->ri+1, acq->nsrc);
  acq->src_w3 = alloc2float(2*sim->ri+1, acq->nsrc);
  acq->src_i1 = alloc1int(acq->nsrc);
  acq->src_i1m = alloc1int(acq->nsrc);
  acq->src_i2 = alloc1int(acq->nsrc);
  acq->src_i3 = alloc1int(acq->nsrc);
  acq->src_nm = alloc1int(acq->nrec);
  for(isrc=0; isrc< acq->nsrc; isrc++){
    acq->src_x1[isrc] = zs;
    acq->src_x2[isrc] = xs;
    acq->src_x3[isrc] = ys;
    //reset the acquisition limits
    zmin = MIN(zmin, zs);
    zmax = MAX(zmax, zs);
    xmin = MIN(xmin, xs);
    xmax = MAX(xmax, xs);
    if(sim->n3>1){
      ymin = MIN(ymin, ys);
      ymax = MAX(ymax, ys);
    }
    
    //point at true depth
    tmp = (acq->src_x1[isrc]-acq->zmin)/sim->d1;
    frac = tmp-NINT(tmp);
    acq->src_i1[isrc] = NINT(tmp) + sim->nb;
    acq->src_nm[isrc] = 0;
    for(j=-sim->ri; j<=sim->ri; j++){
      acq->src_w1[isrc][j+sim->ri] = kaiser_windowed_sinc(-frac+j, 1.0, sim->ri);
      if(tmp-sim->ri+j<0.) acq->src_nm[isrc]++;//count number of points above free surface
    }
    //point at mirror location around free surface
    tmp = -tmp;
    frac = tmp-NINT(tmp);
    acq->src_i1m[isrc] = NINT(tmp) + sim->nb;
    for(j=-sim->ri; j<=sim->ri; j++)
      acq->src_w1m[isrc][j+sim->ri] = kaiser_windowed_sinc(-frac+j, 1.0, sim->ri);
    
    tmp = (acq->src_x2[isrc]-acq->xmin)/sim->d2;
    frac = tmp-NINT(tmp);
    acq->src_i2[isrc] = NINT(tmp) + sim->nb;
    for(j=-sim->ri; j<=sim->ri; j++)
      acq->src_w2[isrc][j+sim->ri] = kaiser_windowed_sinc(-frac+j, 1.0, sim->ri);
    
    tmp = (acq->src_x3[isrc]-acq->ymin)/sim->d3;
    frac = tmp-NINT(tmp);
    acq->src_i3[isrc] = NINT(tmp) + sim->nb;
    for(j=-sim->ri; j<=sim->ri; j++)
      acq->src_w3[isrc][j+sim->ri] = kaiser_windowed_sinc(-frac+j, 1.0, sim->ri);

  }
  if(zmin<acq->zmin) err("source location: z<zmin");
  if(zmax>acq->zmax) err("source location: z>zmax");
  if(xmin<acq->xmin) err("source location: x<xmin");
  if(xmax>acq->xmax) err("source location: x>xmax");
  if(sim->n3>1){
    if(ymin<acq->ymin) err("source location: y<ymin");
    if(ymax>acq->ymax) err("source location: y>ymax");
  }

  acq->rec_x1 = alloc1float(acq->nrec);
  acq->rec_x2 = alloc1float(acq->nrec);
  acq->rec_x3 = alloc1float(acq->nrec);
  acq->rec_w1 = alloc2float(2*sim->ri+1, acq->nrec);
  acq->rec_w1m = alloc2float(2*sim->ri+1, acq->nrec);
  acq->rec_w2 = alloc2float(2*sim->ri+1, acq->nrec);
  acq->rec_w3 = alloc2float(2*sim->ri+1, acq->nrec);
  acq->rec_i1 = alloc1int(acq->nrec);
  acq->rec_i1m = alloc1int(acq->nrec);
  acq->rec_i2 = alloc1int(acq->nrec);
  acq->rec_i3 = alloc1int(acq->nrec);
  acq->rec_nm = alloc1int(acq->nrec);
  for(irec=0; irec< acq->nrec; irec++){
    acq->rec_x1[irec] = rx1[irec];
    acq->rec_x2[irec] = rx2[irec];
    acq->rec_x3[irec] = rx3[irec];

    //reset the acquisition limits
    zmin = MIN(zmin, rx1[irec]);
    zmax = MAX(zmax, rx1[irec]);
    xmin = MIN(xmin, rx2[irec]);
    xmax = MAX(xmax, rx2[irec]);
    if(sim->n3>1){
      ymin = MIN(ymin, rx3[irec]);
      ymax = MAX(ymax, rx3[irec]);
    }

    //point at true depth
    tmp = (acq->rec_x1[irec]-acq->zmin)/sim->d1;
    frac = tmp-NINT(tmp);
    acq->rec_i1[irec] = NINT(tmp) + sim->nb;
    acq->rec_nm[irec] = 0;
    for(j=-sim->ri; j<=sim->ri; j++){
      acq->rec_w1[irec][j+sim->ri] = kaiser_windowed_sinc(-frac+j, 1.0, sim->ri);
      if(tmp-sim->ri+j<0.) acq->rec_nm[irec]++;//count number of points above free surface
    }
    //point at mirror location around free surface
    tmp = -tmp;
    frac = tmp-NINT(tmp);
    acq->rec_i1m[irec] = NINT(tmp) + sim->nb;
    for(j=-sim->ri; j<=sim->ri; j++)
      acq->rec_w1m[irec][j+sim->ri] = kaiser_windowed_sinc(-frac+j, 1.0, sim->ri);

    tmp = (acq->rec_x2[irec]-acq->xmin)/sim->d2;
    frac = tmp-NINT(tmp);
    acq->rec_i2[irec] = NINT(tmp) + sim->nb;
    for(j=-sim->ri; j<=sim->ri; j++)
      acq->rec_w2[irec][j+sim->ri] = kaiser_windowed_sinc(-frac+j, 1.0, sim->ri);

    tmp = (acq->rec_x3[irec]-acq->ymin)/sim->d3;
    frac = tmp-NINT(tmp);
    acq->rec_i3[irec] = NINT(tmp) + sim->nb;
    for(j=-sim->ri; j<=sim->ri; j++)
      acq->rec_w3[irec][j+sim->ri] = kaiser_windowed_sinc(-frac+j, 1.0, sim->ri);
  }

  if(zmin<acq->zmin) err("receiver location: z<zmin");
  if(zmax>acq->zmax) err("receiver location: z>zmax");
  if(xmin<acq->xmin) err("receiver location: x<xmin");
  if(xmax>acq->xmax) err("receiver location: x>xmax");
  if(sim->n3>1){
    if(ymin<acq->ymin) err("receiver location: y<ymin");
    if(ymax>acq->ymax) err("receiver location: y>ymax");
  }
  if(iproc==0) printf("------------- acquisition init done ---------------\n");

  free1float(rx1);
  free1float(rx2);
  free1float(rx3);

  //allocate memory for observed and synthetic data, and data residual
  sim->dobs = alloc2float(sim->nt,acq->nrec);
  sim->dcal = alloc2float(sim->nt,acq->nrec);
  sim->dres = alloc2float(sim->nt,acq->nrec);
  memset(&sim->dobs[0][0], 0, sim->nt*acq->nrec*sizeof(float));
  memset(&sim->dcal[0][0], 0, sim->nt*acq->nrec*sizeof(float));
  memset(&sim->dres[0][0], 0, sim->nt*acq->nrec*sizeof(float));
}


void acq_free(sim_t *sim, acq_t *acq)
/*< free the allocated variables for acquisition >*/
{
  free1float(acq->src_x1);
  free1float(acq->src_x2);
  free1float(acq->src_x3);
  free1float(acq->rec_x1);
  free1float(acq->rec_x2);
  free1float(acq->rec_x3);

  free2float(acq->src_w1);
  free2float(acq->src_w1m);
  free2float(acq->src_w2);
  free2float(acq->src_w3);
  free2float(acq->rec_w1);
  free2float(acq->rec_w1m);
  free2float(acq->rec_w2);
  free2float(acq->rec_w3);

  free1int(acq->src_i1);
  free1int(acq->src_i1m);
  free1int(acq->src_i2);
  free1int(acq->src_i3);
  free1int(acq->rec_i1);
  free1int(acq->rec_i1m);
  free1int(acq->rec_i2);
  free1int(acq->rec_i3);

  free1int(acq->src_nm);
  free1int(acq->rec_nm);
  free1int(acq->shot_idx);

  free2float(sim->dobs);
  free2float(sim->dcal);
  free2float(sim->dres);
}
