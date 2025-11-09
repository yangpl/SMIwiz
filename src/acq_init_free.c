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
#include "segy.h"

double kaiser_windowed_sinc(double x, double dx, int r);

void acq_file_init(sim_t *sim, acq_t *acq);
void su_file_init(char *fname, sim_t *sim, acq_t *acq);

void acq_init(sim_t *sim, acq_t *acq)
{
  char fname[sizeof "dat_0000"];
  if(!getparint("suopt", &acq->suopt)) acq->suopt = 0;//1=for RTM,FWI,LSRTM
  if(acq->suopt){
    sprintf(fname, "dat_%04d", acq->shot_idx[iproc]);
    su_file_init(fname, sim, acq);
  }else{
    if(!getparfloat("zmin", &acq->zmin)) acq->zmin = 0;
    if(!getparfloat("zmax", &acq->zmax)) acq->zmax = acq->zmin+(sim->n1-1)*sim->d1;
    if(!getparfloat("xmin", &acq->xmin)) acq->xmin = 0;
    if(!getparfloat("xmax", &acq->xmax)) acq->xmax = acq->xmin+(sim->n2-1)*sim->d2;
    if(!getparfloat("ymin", &acq->ymin)) acq->ymin = 0;
    if(!getparfloat("ymax", &acq->ymax)) acq->ymax = acq->ymin+(sim->n3-1)*sim->d3;
    if(iproc==0){
      printf("[zmin, zmax]=[%g, %g]\n", acq->zmin, acq->zmax);
      printf("[xmin, xmax]=[%g, %g]\n", acq->xmin, acq->xmax);
      printf("[ymin, ymax]=[%g, %g]\n", acq->ymin, acq->ymax);
    }
    acq_file_init(sim, acq);
  }

}

void acq_file_init(sim_t *sim, acq_t *acq)
/*< read acquisition file to initialize acquisition geometry >*/
{
  char *acquifile;
  int nrec_max, nsrc;
  float *rx1, *rx2, *rx3;
  float zs, xs, ys;
  float zz, xx, yy, dip, azimuth,  tmp,  frac;
  float zmin, zmax, xmin, xmax, ymin, ymax;
  int isreceiver, isrc, irec, iseof, j;
  FILE *fp;
  
  acq->shot_idx = alloc1int(nproc);
  nsrc = countparval("shots");
  if(nsrc>0){
    if(nsrc<nproc) err("nproc > number of shot indices! ");
    getparint("shots", acq->shot_idx);/* a list of source index separated by comma */
  }
  if(nsrc==0){
    for(j=0; j<nproc; j++) acq->shot_idx[j] = j+1;//index starts from 1
  }
  
  if(!getparint("nrec_max", &nrec_max)) nrec_max = 100000;//maximum dimensions/receivers per shot
  if(iproc==0) printf("nrec_max=%d \n", nrec_max);
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

  free1float(rx1);
  free1float(rx2);
  free1float(rx3);
}

void su_file_init(char *fname, sim_t *sim, acq_t *acq)
/*< read SU file to initialize acquisition geometry >*/
{
  unsigned short ns;
  unsigned long int filesize;
  trace_header *trhdr;
  float **trace, tmp, frac;
  float zmin, zmax, xmin, xmax, ymin, ymax;
  int j, isrc,irec;
  FILE *fp;

  acq->nsrc=1;//by default 1 shot per process
  fp=fopen(fname,"rb");
  if(fp==NULL) err("File does not exist!");
  fseek(fp, 114, SEEK_SET );//pointer skips first 114 bytes
  fread(&ns, 2, 1,fp); // ns 115-116 byte in the trace header

  fseek(fp, 0, SEEK_END);//pointer goes to the end of the file
  filesize=ftell(fp);//ftell tells you the number of bytes the file has
  //each trace starts with 240 bytes header,
  //then proceeds with ns samples in float point.
  acq->nrec=filesize/(240+ns*4);

  rewind(fp);//pointer goes to the beginning of the file
  trhdr=(trace_header*)malloc(acq->nrec*sizeof(trace_header));
  trace=alloc2float(ns,acq->nrec);
  for(irec=0; irec<acq->nrec; irec++){
    fread(trhdr,240,1,fp);
    fread(trace[irec],ns*sizeof(float),1,fp);
  }
  fclose(fp);
  //neeed to check whether nt==ns: if not, interpolate the traces

  //allocate memory for source and receivers
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
  for(isrc=0; isrc< acq->nsrc; isrc++){
    //source position from SU header
    if(trhdr[isrc].scalel==0) tmp=1.;
    else if(trhdr[isrc].scalel>0) tmp=trhdr[isrc].scalel;
    else tmp=1./fabs(trhdr[isrc].scalel);
    acq->src_x1[isrc] = (trhdr[isrc].selev + trhdr[isrc].sdepth)*tmp;

    if(trhdr[isrc].scalco==0) tmp=1.;
    else if(trhdr[isrc].scalco>0) tmp=trhdr[isrc].scalco;
    else tmp=1./fabs(trhdr[isrc].scalco);
    acq->src_x2[isrc] = trhdr[isrc].sx * tmp;
    acq->src_x3[isrc] = trhdr[isrc].sy * tmp;

    //coordinate of sources, origin stripped out
    //z0=zs-acq->zmin;
    //x0=xs-acq->xmin;

    //reset the acquisition limits
    zmin = (isrc==0)? acq->src_x1[isrc] : MIN(zmin, acq->src_x1[isrc]);
    zmax = (isrc==0)? acq->src_x1[isrc] : MAX(zmax, acq->src_x1[isrc]);
    xmin = (isrc==0)? acq->src_x2[isrc] : MIN(xmin, acq->src_x2[isrc]);
    xmax = (isrc==0)? acq->src_x2[isrc] : MAX(xmax, acq->src_x2[isrc]);
    if(sim->n3>1){
      ymin = (isrc==0)? acq->src_x3[isrc] : MIN(ymin, acq->src_x3[isrc]);
      ymax = (isrc==0)? acq->src_x3[isrc] : MAX(ymax, acq->src_x3[isrc]);
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
  /* if(zmin<acq->zmin) err("source location: z<zmin"); */
  /* if(zmax>acq->zmax) err("source location: z>zmax"); */
  /* if(xmin<acq->xmin) err("source location: x<xmin"); */
  /* if(xmax>acq->xmax) err("source location: x>xmax"); */
  /* if(sim->n3>1){ */
  /*   if(ymin<acq->ymin) err("source location: y<ymin"); */
  /*   if(ymax>acq->ymax) err("source location: y>ymax"); */
  /* } */


  //receiver position from SU header
  for(irec=0; irec<acq->nrec; irec++){
    if(trhdr[irec].scalel==0) tmp=1.;
    else if(trhdr[irec].scalel>0) tmp=trhdr[irec].scalel;
    else tmp=1./fabs(trhdr[irec].scalel);
    acq->rec_x1[irec] = -trhdr[irec].gelev*tmp;

    if(trhdr[irec].scalco==0) tmp=1.;
    else if(trhdr[irec].scalco>0) tmp=trhdr[irec].scalco;
    else tmp=1./fabs(trhdr[irec].scalco);
    acq->rec_x2[irec] = trhdr[irec].gx * tmp;
    acq->rec_x3[irec] = trhdr[irec].gy * tmp;

    //coorxinate of sources, origin stripped out
    //z0=zs-acq->zmin;
    //x0=xs-acq->xmin;

    //reset the acquisition limits
    zmin = (irec==0)? acq->rec_x1[irec] : MIN(zmin, acq->rec_x1[irec]);
    zmax = (irec==0)? acq->rec_x1[irec] : MAX(zmax, acq->rec_x1[irec]);
    xmin = (irec==0)? acq->rec_x2[irec] : MIN(xmin, acq->rec_x2[irec]);
    xmax = (irec==0)? acq->rec_x2[irec] : MAX(xmax, acq->rec_x2[irec]);
    if(sim->n3>1){
      ymin = (irec==0)? acq->rec_x3[irec] : MIN(ymin, acq->rec_x3[irec]);
      ymax = (irec==0)? acq->rec_x3[irec] : MAX(ymax, acq->rec_x3[irec]);
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
  /* if(zmin<acq->zmin) err("source location: z<zmin"); */
  /* if(zmax>acq->zmax) err("source location: z>zmax"); */
  /* if(xmin<acq->xmin) err("source location: x<xmin"); */
  /* if(xmax>acq->xmax) err("source location: x>xmax"); */
  /* if(sim->n3>1){ */
  /*   if(ymin<acq->ymin) err("source location: y<ymin"); */
  /*   if(ymax>acq->ymax) err("source location: y>ymax"); */
  /* } */

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
}
