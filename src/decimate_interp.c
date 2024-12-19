/* 2D/3D seismic modelling, RTM and FWI code
 *-----------------------------------------------------------------------
 *  Reference:
 *  Yang, P., Brossier, R., and Virieux, J. (2016). Wavefield reconstruction
 *  from significantly decimated boundaries. Geophysics, 81 (5) T197–T209
 *
 * Remark: it is possible to reduce one additional layer in storing boundaries
 *   over staggered grid, as I show in 2014 Computers & Geosciences. It is 
 *   not implemented for the ease of grid manipulation.
 *-----------------------------------------------------------------------
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com/
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "sim.h"

double sinc(double x);
double bessi0(double x);

void decimate_interp_init(sim_t *sim, int flag)
{
  if(flag==1){
    if(sim->order==8){
      //allocate memory to store decimated boundaries, 8 layers on each side
      sim->face1 = alloc2float(sim->mt, 16*sim->n2*sim->n3);
      sim->face2 = alloc2float(sim->mt, 16*sim->n1*sim->n3);
      memset(sim->face1[0], 0, sim->mt*16*sim->n2*sim->n3*sizeof(float));
      memset(sim->face2[0], 0, sim->mt*16*sim->n1*sim->n3*sizeof(float));
      if(sim->n3>1){
	sim->face3 = alloc2float(sim->mt, 16*sim->n1*sim->n2);
	memset(sim->face3[0], 0, sim->mt*16*sim->n1*sim->n2*sizeof(float));
      }
    }else if(sim->order==4){
      sim->face1 = alloc2float(sim->mt, 8*sim->n2*sim->n3);
      sim->face2 = alloc2float(sim->mt, 8*sim->n1*sim->n3);
      memset(sim->face1[0], 0, sim->mt*8*sim->n2*sim->n3*sizeof(float));
      memset(sim->face2[0], 0, sim->mt*8*sim->n1*sim->n3*sizeof(float));
      if(sim->n3>1){
	sim->face3 = alloc2float(sim->mt, 8*sim->n1*sim->n2);
	memset(sim->face3[0], 0, sim->mt*8*sim->n1*sim->n2*sizeof(float));
      }
    }
  }else if(flag==0){
    if(sim->order==8){
      //allocate memory to store decimated boundaries, 8 layers on each side
      sim->face1_ = alloc2float(sim->mt, 16*sim->n2*sim->n3);
      sim->face2_ = alloc2float(sim->mt, 16*sim->n1*sim->n3);
      memset(sim->face1_[0], 0, sim->mt*16*sim->n2*sim->n3*sizeof(float));
      memset(sim->face2_[0], 0, sim->mt*16*sim->n1*sim->n3*sizeof(float));
      if(sim->n3>1){
	sim->face3_ = alloc2float(sim->mt, 16*sim->n1*sim->n2);
	memset(sim->face3_[0], 0, sim->mt*16*sim->n1*sim->n2*sizeof(float));
      }
    }else if(sim->order==4){
      sim->face1_ = alloc2float(sim->mt, 8*sim->n2*sim->n3);
      sim->face2_ = alloc2float(sim->mt, 8*sim->n1*sim->n3);
      memset(sim->face1_[0], 0, sim->mt*8*sim->n2*sim->n3*sizeof(float));
      memset(sim->face2_[0], 0, sim->mt*8*sim->n1*sim->n3*sizeof(float));
      if(sim->n3>1){
	sim->face3_ = alloc2float(sim->mt, 8*sim->n1*sim->n2);
	memset(sim->face3_[0], 0, sim->mt*8*sim->n1*sim->n2*sizeof(float));
      }
    }
  }
  
}

void decimate_interp_close(sim_t *sim, int flag)
{
  if(flag==1){
    free2float(sim->face1);
    free2float(sim->face2);
    if(sim->n3>1) free2float(sim->face3);
  }else if(flag==0){
    free2float(sim->face1_);
    free2float(sim->face2_);
    if(sim->n3>1) free2float(sim->face3_);
  }
}

double kwsinc(double x, int l, int r)
{
  static float b = 6.31;
  static float PI = 3.141592653589793238462643;
  double a = x/(l*r);
  if(fabs(a)>1) return 0;

  double s = bessi0(b*sqrt(1.0-a*a))/bessi0(b);
  a = PI*x/r;
  if(fabs(a)>1.e-15) s *= sin(a)/a;
    
  return s;
}

void decimate_interp_bndr(sim_t *sim, int flag, int it, int interp, float **face1, float **face2, float **face3)
{
  int id[8];//8-point Kaiser windowed sinc interpolation, half lenght l=4
  float kw[8];
  int i1,i2,i3,i,j,k,m,m_;
  int jb = sim->order/2;//minimum and maximum bounds for j is [-jb,jb-1]
  float ***vz, ***vx, ***vy;

  if(flag==0){
    vz = sim->vz0;
    vx = sim->vx0;
    vy = sim->vy0;
  }else if(flag==1){
    vz = sim->vz1;
    vx = sim->vx1;
    vy = sim->vy1;
  }

  m = (it+1)/sim->dr;//range from 0 to sim->mt
  m_ = MAX(1, m)-1;//index range from 0 to sim->mt-1
  if(interp){
    //printf("it=%d\t", it);
    //interpolation from stored boundaries, do it every step
    if(m_<3) {
      for(i=0; i<8; i++) {
	id[i] = i;	
	kw[i] = kwsinc((double)(it+1-(i+1)*sim->dr), 4, sim->dr);
      }
    }else if(m_+4>sim->mt-1){
      for(i=sim->mt-8; i<=sim->mt-1; i++){
	k = i-(sim->mt-8);
	id[k] = i;
	kw[k] = kwsinc((double)(it+1-(i+1)*sim->dr), 4, sim->dr);
      }
    }else{//m_-3>=0; m_+4<=mt-1
      for(i=m_-3; i<=m_+4; i++){
	k = i-(m_-3);
	id[k] = i;
	kw[k] = kwsinc((double)(it+1-(i+1)*sim->dr), 4, sim->dr);
      }
    }

    //========================================================
    if(sim->n3>1){/* 3D */
      //face1
      k=0;
      for(i3=sim->nb; i3<sim->nb+sim->n3; i3++){
	for(i2=sim->nb; i2<sim->nb+sim->n2; i2++){
	  for(j=-jb; j<jb; j++){
	    //top
	    vz[i3][i2][sim->nb+j]
	      = kw[0]*face1[k][id[0]] + kw[1]*face1[k][id[1]]
	      + kw[2]*face1[k][id[2]] + kw[3]*face1[k][id[3]]
	      + kw[4]*face1[k][id[4]] + kw[5]*face1[k][id[5]]
	      + kw[6]*face1[k][id[6]] + kw[7]*face1[k][id[7]];
	    k++;
	    //bottom
	    vz[i3][i2][sim->nb+sim->n1+j]
	      = kw[0]*face1[k][id[0]] + kw[1]*face1[k][id[1]]
	      + kw[2]*face1[k][id[2]] + kw[3]*face1[k][id[3]]
	      + kw[4]*face1[k][id[4]] + kw[5]*face1[k][id[5]]
	      + kw[6]*face1[k][id[6]] + kw[7]*face1[k][id[7]];
	    k++;
	  }
	}
      }
      //face2
      k=0;
      for(i3=sim->nb; i3<sim->nb+sim->n3; i3++){
	for(i1=sim->nb; i1<sim->nb+sim->n1; i1++){
	  for(j=-jb; j<jb; j++){
	    //left
	    vx[i3][sim->nb+j][i1]
	      = kw[0]*face2[k][id[0]] + kw[1]*face2[k][id[1]]
	      + kw[2]*face2[k][id[2]] + kw[3]*face2[k][id[3]]
	      + kw[4]*face2[k][id[4]] + kw[5]*face2[k][id[5]]
	      + kw[6]*face2[k][id[6]] + kw[7]*face2[k][id[7]];
	    k++;
	    //right
	    vx[i3][sim->nb+sim->n2+j][i1]
	      = kw[0]*face2[k][id[0]] + kw[1]*face2[k][id[1]]
	      + kw[2]*face2[k][id[2]] + kw[3]*face2[k][id[3]]
	      + kw[4]*face2[k][id[4]] + kw[5]*face2[k][id[5]]
	      + kw[6]*face2[k][id[6]] + kw[7]*face2[k][id[7]];
	    k++;
	  }
	}
      }
      //face3
      k = 0;
      for(i2=sim->nb; i2<sim->nb+sim->n2; i2++){
	for(i1=sim->nb; i1<sim->nb+sim->n1; i1++){
	  for(j=-jb; j<jb; j++){
	    //left
	    vy[sim->nb+j][i2][i1]
	      = kw[0]*face3[k][id[0]] + kw[1]*face3[k][id[1]]
	      + kw[2]*face3[k][id[2]] + kw[3]*face3[k][id[3]]
	      + kw[4]*face3[k][id[4]] + kw[5]*face3[k][id[5]]
	      + kw[6]*face3[k][id[6]] + kw[7]*face3[k][id[7]];
	    k++;
	    //right
	    vy[sim->nb+sim->n3+j][i2][i1]
	      = kw[0]*face3[k][id[0]] + kw[1]*face3[k][id[1]]
	      + kw[2]*face3[k][id[2]] + kw[3]*face3[k][id[3]]
	      + kw[4]*face3[k][id[4]] + kw[5]*face3[k][id[5]]
	      + kw[6]*face3[k][id[6]] + kw[7]*face3[k][id[7]];
	    k++;
	  }

	}
      }
      
    }else{/* 2D */
      i3=0;
      //face1
      k=0;
      for(i2=sim->nb; i2<sim->nb+sim->n2; i2++){
	for(j=-jb; j<jb; j++){
	  //top
	  vz[i3][i2][sim->nb+j]
	    = kw[0]*face1[k][id[0]] + kw[1]*face1[k][id[1]]
	    + kw[2]*face1[k][id[2]] + kw[3]*face1[k][id[3]]
	    + kw[4]*face1[k][id[4]] + kw[5]*face1[k][id[5]]
	    + kw[6]*face1[k][id[6]] + kw[7]*face1[k][id[7]];
	  k++;
	  //bottom
	  vz[i3][i2][sim->nb+sim->n1+j]
	    = kw[0]*face1[k][id[0]] + kw[1]*face1[k][id[1]]
	    + kw[2]*face1[k][id[2]] + kw[3]*face1[k][id[3]]
	    + kw[4]*face1[k][id[4]] + kw[5]*face1[k][id[5]]
	    + kw[6]*face1[k][id[6]] + kw[7]*face1[k][id[7]];
	  k++;
	}
      }
      //face2
      k=0;
      for(i1=sim->nb; i1<sim->nb+sim->n1; i1++){
	for(j=-jb; j<jb; j++){
	  //left
	  vx[i3][sim->nb+j][i1]
	    = kw[0]*face2[k][id[0]] + kw[1]*face2[k][id[1]]
	    + kw[2]*face2[k][id[2]] + kw[3]*face2[k][id[3]]
	    + kw[4]*face2[k][id[4]] + kw[5]*face2[k][id[5]]
	    + kw[6]*face2[k][id[6]] + kw[7]*face2[k][id[7]];
	  k++;
	  //right
	  vx[i3][sim->nb+sim->n2+j][i1]
	    = kw[0]*face2[k][id[0]] + kw[1]*face2[k][id[1]]
	    + kw[2]*face2[k][id[2]] + kw[3]*face2[k][id[3]]
	    + kw[4]*face2[k][id[4]] + kw[5]*face2[k][id[5]]
	    + kw[6]*face2[k][id[6]] + kw[7]*face2[k][id[7]];
	  k++;
	}
      }
    }

  }else if((it+1)%sim->dr==0){
    //decimate the boundary, do it only at Nyquist step
    if(sim->n3>1){
      //face1
      k=0;
      for(i3=sim->nb; i3<sim->nb+sim->n3; i3++){
	for(i2=sim->nb; i2<sim->nb+sim->n2; i2++){
	  for(j=-jb; j<jb; j++){
	    //top
	    face1[k][m_] = vz[i3][i2][sim->nb+j];
	    k++;
	    //bottom
	    face1[k][m_] = vz[i3][i2][sim->nb+sim->n1+j];
	    k++;
	  }
	}
      }
      //face2
      k=0;
      for(i3=sim->nb; i3<sim->nb+sim->n3; i3++){
	for(i1=sim->nb; i1<sim->nb+sim->n1; i1++){
	  for(j=-jb; j<jb; j++){
	    //left
	    face2[k][m_] = vx[i3][sim->nb+j][i1];
	    k++;
	    //right
	    face2[k][m_] = vx[i3][sim->nb+sim->n2+j][i1];
	    k++;
	  }
	}
      }
      //face3
      k=0;
      for(i2=sim->nb; i2<sim->nb+sim->n2; i2++){
	for(i1=sim->nb; i1<sim->nb+sim->n1; i1++){
	  for(j=-jb; j<jb; j++){
	    //left
	    face3[k][m_] = vy[sim->nb+j][i2][i1];
	    k++;
	    //right
	    face3[k][m_] = vy[sim->nb+sim->n3+j][i2][i1];
	    k++;
	  }
	}
      }
      
    }else{
      i3=0;
      //face1
      k=0;
      for(i2=sim->nb; i2<sim->nb+sim->n2; i2++){
	for(j=-jb; j<jb; j++){
	  //top
	  face1[k][m_] = vz[i3][i2][sim->nb+j];
	  k++;
	  //bottom
	  face1[k][m_] = vz[i3][i2][sim->nb+sim->n1+j];
	  k++;
	}
      }
      //face2
      k=0;
      for(i1=sim->nb; i1<sim->nb+sim->n1; i1++){
	for(j=-jb; j<jb; j++){
	  //left
	  face2[k][m_] = vx[i3][sim->nb+j][i1];
	  k++;
	  //right
	  face2[k][m_] = vx[i3][sim->nb+sim->n2+j][i1];
	  k++;
	}
      }
    }//end if

  }//end if
  
}
