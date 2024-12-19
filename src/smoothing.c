/*  Claerbout's box-triangle smoothing adapted for 2D/3D!
 *---------------------------------------------------------------------------
 *  Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
 *  Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
 *  Homepage: https://yangpl.wordpress.com
 *  E-mail: ypl.2100@gmail.com
 *-------------------------------------------------------------------------*/
#include "cstd.h"

float *pp, *qq, *bb1, *bb2, *yy;
int np, nq;

void triangle_init(int nbox, int nd)
{
  np = nbox + nd -1;
  nq = nbox + np -1;  
  pp = alloc1float(np);
  qq = alloc1float(nq);
  bb1 = alloc1float(np);
  bb2 = alloc1float(nq);
  yy = alloc1float(nd);
  
}

void triangle_close()
{
  free1float(pp);
  free1float(qq);
  free1float(bb1);
  free1float(bb2);
  free1float(yy);
}

void boxconv_lop(int nbox, int nx, float *xx, float *yy, float *bb)
{
  int ny;
  int i;

  ny  = nx + nbox-1;
  for(i=0; i<ny; i++) bb[i] = 0.;
  bb[0] = xx[0];
  for(i=1; i<nx; i++) bb[i] = bb[i-1] + xx[i]; //make B(z) = X(z)/(1-z)
  for(i=nx; i<ny; i++) bb[i] = bb[i-1];
  for(i=0; i<nbox; i++) yy[i] = bb[i];
  for(i=nbox; i<ny; i++) yy[i] = bb[i]-bb[i-nbox];//make Y(z) = B(z)*(1-z**nbox)
  for(i=0; i<ny; i++) yy[i] /= nbox;
}

void triangle_lop(int nbox, int nd, float *xx)
{
  int i;
  
  boxconv_lop(nbox, nd, xx, pp, bb1);
  boxconv_lop(nbox, np, pp, qq, bb2);
  for(i=0; i<nd; i++) yy[i] = qq[i+nbox-1];
  for(i=0; i<nbox-1; i++) yy[i] += qq[nbox-2-i]; //fold back near end
  for(i=0; i<nbox-1; i++) yy[nd-i-1] += qq[nd +(nbox-1)+i];//fold back far end
  memcpy(xx, yy, nd*sizeof(float));
  
}


void triangle_smoothing(float ***mod, int n1, int n2, int n3, int r1, int r2, int r3, int repeat)
{
  int i, i1, i2, i3;
  float *tmp, *tmp2;
  
  tmp = alloc1float(n2);
  if(n3>1) tmp2 = alloc1float(n3);

  for(i=0; i<repeat; i++){
    /*-----------------------------------------*/
    triangle_init(r1, n1);
    for(i3=0; i3<n3; i3++)
      for(i2=0; i2<n2; i2++)
	triangle_lop(r1, n1, mod[i3][i2]);
    triangle_close();

    /*-----------------------------------------*/
    triangle_init(r2, n2);
    for(i3=0; i3<n3; i3++){
      for(i1=0; i1<n1; i1++){
	for(i2=0; i2<n2; i2++) tmp[i2] = mod[i3][i2][i1];
	triangle_lop(r2, n2, tmp);
	for(i2=0; i2<n2; i2++) mod[i3][i2][i1] = tmp[i2];
      }
    }
    triangle_close();

    /*-----------------------------------------*/
    if(n3>1){
      triangle_init(r3, n3);
      for(i2=0; i2<n2; i2++){
	for(i1=0; i1<n1; i1++){
	  for(i3=0; i3<n3; i3++) tmp2[i3] = mod[i3][i2][i1];
	  triangle_lop(r3, n3, tmp2);
	  for(i3=0; i3<n3; i3++) mod[i3][i2][i1] = tmp2[i3];
	}
      }
      triangle_close();
    }
  }

  free1float(tmp);
  if(n3>1) free1float(tmp2);

}
  
