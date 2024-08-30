/* set data weights to mimic sumute in Seismic Unix */
/*
  Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
  Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
  Homepage: https://yangpl.wordpress.com
  E-mail: ypl.2100@gmail.com
*/
#include "cstd.h"
#include "sim.h"
#include "acq.h"
#include "mpi_info.h"

void setup_data_weight(acq_t *acq, sim_t *sim)
/*< set up data weighting used in RTM and FWI >*/
{
  const float PI=3.141592653589793238462643;
  int ntaper;//number of points for cosine taper
  int nxwdat;//number of points in x-axis for data weights
  float dxwdat;//spatial interval in x-axis for data weights
  int nxmute1,ntmute1,nxmute2,ntmute2,nxtmute;
  float *xwdat,*xmute1,*tmute1,*xmute2,*tmute2;

  int it,is,ir,id,k,itmute,itstart,itend,ntaper_true;
  float d1,d2,d3,dist,s,tmute;
  FILE *fp;

  int nt = sim->nt;
  float dt = sim->dt;
  //data weighting parameters
  if(!getparfloat("dxwdat",&dxwdat)) dxwdat = 100;
  if(!(nxwdat = countparval("xwdat"))) nxwdat = 6;
  xwdat = alloc1float(nxwdat);
  if(!getparfloat("xwdat",xwdat)){
    for(k=0; k<nxwdat; k++) xwdat[k] = 1.;
  }

  acq->xweight=alloc1float(acq->nrec);
  acq->wdat=alloc2float(nt,acq->nrec);
  is=0;
  for(ir=0; ir<acq->nrec; ir++){
    //compute distance between source and receiver
    d1=acq->src_x1[is]-acq->rec_x1[ir];
    d2=acq->src_x2[is]-acq->rec_x2[ir];
    d3=acq->src_x3[is]-acq->rec_x3[ir];
    if(sim->n3>1) dist=sqrt(d1*d1+d2*d2+d3*d3);//offset
    else dist=sqrt(d1*d1+d2*d2);//offset
    
    s=dist/dxwdat;
    id = (int)s;
    s=s-id;//reduced coordinate between 0 and 1
    //linear interpolation +extrapolation of the last value if the table is not big enough
    acq->xweight[ir] = xwdat[MIN(id,nxwdat-1)]*(1.-s) + xwdat[MIN(id+1,nxwdat-1)]*s;
    for(it=0; it<nt; it++) acq->wdat[ir][it] = acq->xweight[ir];
  }

  //data muting parameters
  if(!getparint("muteopt", &sim->muteopt)) sim->muteopt = 0;
  //0, do not mute; 1,mute above; 2, mute below; 3, mute above and below
  if(!getparint("ntaper", &ntaper)) ntaper = 10;
  if(sim->muteopt==1||sim->muteopt==3){
    if (!(nxmute1 = countparval("xmute1"))) err("must give xmute1= vector");
    if (!(ntmute1 = countparval("tmute1"))) err("must give tmute1= vector");
    if (nxmute1 != ntmute1) err("length of xmute1, tmute1 must be the same");
    nxtmute=nxmute1;
    xmute1=alloc1float(nxtmute);
    tmute1=alloc1float(nxtmute);
    getparfloat("xmute1", xmute1);
    getparfloat("tmute1", tmute1);
    for(ir=0; ir<acq->nrec; ir++){
      d1=acq->src_x1[is]-acq->rec_x1[ir];
      d2=acq->src_x2[is]-acq->rec_x2[ir];
      d3=acq->src_x3[is]-acq->rec_x3[ir];
      if(sim->n3>1) dist=sqrt(d1*d1+d2*d2+d3*d3);//offset
      else dist=sqrt(d1*d1+d2*d2);//offset
      id=nxtmute-1;
      for(k=0; k<nxtmute; k++){
	if(dist<xmute1[k]) { 
	  id=k-1; break; 
	}
      }
      if(id==-1){
	tmute = tmute1[0];
      }else if(id==nxtmute-1){
	tmute=tmute1[nxtmute-1];
      }else{
	s=(dist-xmute1[id])/(xmute1[id+1]-xmute1[id]);
	tmute=tmute1[id]*(1.-s)+tmute1[id+1]*s;//mute time
      }

      itmute=MIN(MAX(0,NINT(tmute/dt)),nt-1);//index of mute time
      itstart=MAX(0,itmute-ntaper);
      if(itmute==nt-1) itstart=nt-1;
      ntaper_true=itmute-itstart+1;//true ntaper
      for(it=0; it<itstart; it++) acq->wdat[ir][it]=0.;//apply front mute
      k=ntaper_true;
      for(it=itstart; it<=itmute; it++){
	s=0.5*(1.+cos(k*PI/ntaper_true));
	acq->wdat[ir][it] *= s;//apply the taper
	k--;
      }
    }
    free1float(xmute1);
    free1float(tmute1);
  }
  if(sim->muteopt==2||sim->muteopt==3){
    if (!(nxmute2 = countparval("xmute2"))) err("must give xmute2= vector");
    if (!(ntmute2 = countparval("tmute2"))) err("must give tmute2= vector");
    if (nxmute2 != ntmute2) err("length of xmute2, tmute2 must be the same");
    nxtmute = nxmute2;
    xmute2 = alloc1float(nxtmute);
    tmute2 = alloc1float(nxtmute);
    getparfloat("xmute2", xmute2);
    getparfloat("tmute2", tmute2);
    for(ir=0; ir<acq->nrec; ir++){
      d1=acq->src_x1[is]-acq->rec_x1[ir];
      d2=acq->src_x2[is]-acq->rec_x2[ir];
      d3=acq->src_x3[is]-acq->rec_x3[ir];
      if(sim->n3>1) dist=sqrt(d1*d1+d2*d2+d3*d3);//offset
      else dist=sqrt(d1*d1+d2*d2);//offset
      id=nxtmute-1;
      for(k=0; k<nxtmute; k++){
	if(dist<xmute2[k]) { 
	  id=k-1; break; 
	}
      }
      if(id==nxtmute-1){
	tmute=tmute2[nxtmute-1];
      }else{
	s=(dist-xmute2[id])/(xmute2[id+1]-xmute2[id]);
	tmute=tmute2[id]*(1.-s)+tmute2[id+1]*s;//mute time
      }

      itmute=MIN(MAX(0,NINT(tmute/dt)),nt-1);//index of mute time
      itend=MIN(nt-1,itmute+ntaper);
      ntaper_true=itend-itmute+1;//true ntaper
      k=1;
      for(it=itmute; it<itend; it++){
	s=0.5*(1.+cos(k*PI/ntaper_true));
	acq->wdat[ir][it]*=s;//apply the taper
	k++;
      }
      for(it=itend; it<nt; it++) acq->wdat[ir][it]=0.;//tail mute
    }
    free1float(xmute2);
    free1float(tmute2);
  }

  if(iproc==0){
    fp = fopen("data_weight", "wb");
    fwrite(&acq->wdat[0][0], acq->nrec*sim->nt*sizeof(float), 1, fp);
    fclose(fp);
  }
}

