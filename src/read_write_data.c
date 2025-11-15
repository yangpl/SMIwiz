/* read/write seismic data 
 *---------------------------------------------------------------------------
 *  Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
 *  Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
 *  Homepage: https://yangpl.wordpress.com
 *  E-mail: ypl.2100@gmail.com
 *-------------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"
#include "acq.h"
#include "segy.h"

double kaiser_windowed_sinc(double x, double dx, int r);

//read SU/binary data file to initialize sim->dobs
void read_data(sim_t *sim, acq_t *acq)
{
  unsigned short ns;
  unsigned long int filesize;
  float dt_trace, tmp, frac;
  float zmin, zmax, xmin, xmax, ymin, ymax;
  int it, j, isrc, irec;
  FILE *fp;

  if(acq->suopt){//deal with observed in SU format
    //======================================================
    char fname[sizeof("dat_0000.su")];
    sprintf(fname, "dat_%04d.su", acq->shot_idx[iproc]);
    
    fp = fopen(fname,"rb");
    if(fp==NULL) { fprintf(stderr,"error in read_data in SU format"); exit(1); }
    fseek(fp, 114, SEEK_SET );//pointer skips first 114 bytes
    fread(&ns, 2, 1,fp); // ns 115-116 byte in the trace header
    if(iproc==0) printf("ns=%d\n", ns);
    
    fseek(fp, 0, SEEK_END);//pointer goes to the end of the file
    filesize = ftell(fp);//ftell tells you the number of bytes the file has
    //each trace starts with 240 bytes header, proceeds with ns samples in float point.
    acq->nrec = filesize/(240+ns*4);
    acq->nsrc = 1;//by default 1 shot per process
    if(iproc==0) printf("nrec=%d\n", acq->nrec);

    sim->dcal = alloc2float(sim->nt, acq->nrec);
    sim->dobs = alloc2float(sim->nt, acq->nrec);
    sim->dres = alloc2float(sim->nt, acq->nrec);

    trace_header *trhdr = malloc(acq->nrec*sizeof(trace_header));
    float *trace = malloc(ns*sizeof(float));

    rewind(fp);//pointer goes to the beginning of the file
    fread(&trhdr[0], 240, 1, fp);//240 Byte for trace header
    dt_trace = trhdr[0].dt*1e-6;//convert micro-seconds to seconds
    if(iproc==0){
      printf("dt_trace=%g\n", dt_trace);
      if(ns*dt_trace<sim->nt*sim->dt) printf("traces will be extrapolated longer\n");
      if(ns*dt_trace>sim->nt*sim->dt) printf("traces will be cut shorter!\n");
    }

    rewind(fp);//pointer goes to the beginning of the file
    for(irec=0; irec<acq->nrec; irec++){
      fread(&trhdr[irec], 240, 1, fp);//sizeof(trace_header)=240 Byte
      fread(trace, ns*sizeof(float), 1, fp);//then ns float numbers as one trace

      //linearly interpolate trace obtained from SU file
      for(it=0; it<sim->nt; it++){
	tmp = it*sim->dt;
	j = (int)(tmp/dt_trace);
	frac = (tmp-j*dt_trace)/dt_trace;
	sim->dobs[irec][it] = (1.-frac)*trace[MIN(j,ns-1)] + frac*trace[MIN(j+1,ns-1)];      
      }
    }
    fclose(fp);
    fflush(stdout);
    
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

      //reset the acquisition limits
      zmin = (irec==0)? acq->rec_x1[irec] : MIN(zmin, acq->rec_x1[irec]);
      zmax = (irec==0)? acq->rec_x1[irec] : MAX(zmax, acq->rec_x1[irec]);
      xmin = (irec==0)? acq->rec_x2[irec] : MIN(xmin, acq->rec_x2[irec]);
      xmax = (irec==0)? acq->rec_x2[irec] : MAX(xmax, acq->rec_x2[irec]);
      ymin = (irec==0)? acq->rec_x3[irec] : MIN(ymin, acq->rec_x3[irec]);
      ymax = (irec==0)? acq->rec_x3[irec] : MAX(ymax, acq->rec_x3[irec]);

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
    if(iproc==0){
      printf("receiver range:\n");
      printf("[zmin, zmax]=[%g, %g]\n", zmin, zmax);
      printf("[xmin, xmax]=[%g, %g]\n", xmin, xmax);
      printf("[ymin, ymax]=[%g, %g]\n", ymin, ymax);
    }

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
    for(isrc=0; isrc<acq->nsrc; isrc++){
      //source position from SU header
      if(trhdr[0].scalel==0) tmp=1.;
      else if(trhdr[0].scalel>0) tmp=trhdr[0].scalel;
      else tmp=1./fabs(trhdr[0].scalel);
      acq->src_x1[isrc] = -(trhdr[0].selev + trhdr[0].sdepth)*tmp;

      if(trhdr[0].scalco==0) tmp=1.;
      else if(trhdr[0].scalco>0) tmp=trhdr[0].scalco;
      else tmp=1./fabs(trhdr[0].scalco);
      acq->src_x2[isrc] = trhdr[0].sx * tmp;
      acq->src_x3[isrc] = trhdr[0].sy * tmp;

      //reset the acquisition limits
      zmin = (isrc==0)? acq->src_x1[isrc] : MIN(zmin, acq->src_x1[isrc]);
      zmax = (isrc==0)? acq->src_x1[isrc] : MAX(zmax, acq->src_x1[isrc]);
      xmin = (isrc==0)? acq->src_x2[isrc] : MIN(xmin, acq->src_x2[isrc]);
      xmax = (isrc==0)? acq->src_x2[isrc] : MAX(xmax, acq->src_x2[isrc]);
      ymin = (isrc==0)? acq->src_x3[isrc] : MIN(ymin, acq->src_x3[isrc]);
      ymax = (isrc==0)? acq->src_x3[isrc] : MAX(ymax, acq->src_x3[isrc]);
    
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
    if(iproc==0){
      printf("source range:\n");
      printf("[zmin, zmax]=[%g, %g]\n", zmin, zmax);
      printf("[xmin, xmax]=[%g, %g]\n", xmin, xmax);
      printf("[ymin, ymax]=[%g, %g]\n", ymin, ymax);
    }

    free(trace);
    free(trhdr);

  }else{//deal with observed data in binary after modelling
    //======================================================
    char fname[sizeof("dat_0000")];
    sprintf(fname, "dat_%04d", acq->shot_idx[iproc]);
    
    fp = fopen(fname,"rb");
    if(fp==NULL) { fprintf(stderr,"error in read_data in binary format\n"); exit(1);}
    fread(&sim->dobs[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
    fclose(fp);
  }

}


void write_data(sim_t *sim, acq_t *acq)
/*< write synthetic data for each shot/process >*/
{
  char fname[sizeof "dat_0000"];
  FILE *fp;

  sprintf(fname, "dat_%04d", acq->shot_idx[iproc]);
  fp=fopen(fname,"wb");
  if(fp==NULL) { fprintf(stderr,"write_data, error opening file\n"); exit(1);}
  fwrite(&sim->dcal[0][0], sim->nt*acq->nrec*sizeof(float), 1, fp);
  fflush(fp);
  fclose(fp);

}
