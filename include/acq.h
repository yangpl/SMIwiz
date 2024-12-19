#ifndef _acq_h_
#define _acq_h_

typedef struct {
  int nsrc, nrec;//number of sources and receivers on each processor
  float zmin, zmax, xmin, xmax, ymin, ymax;//domain [zmin,zmax]*[xmin,xmax]*[ymin,ymax]

  int *shot_idx;//shot index for each processor
  
  float *src_x1, *src_x2, *src_x3;//physical coordinates for sources
  float *rec_x1, *rec_x2, *rec_x3;//physical coordinates for receivers

  float **src_w1, **src_w1m, **src_w2, **src_w3;//interpolation weights
  float **rec_w1, **rec_w1m, **rec_w2, **rec_w3;//interpolation weights

  int *src_i1, *src_i1m, *src_i2, *src_i3;//incerpolacion center
  int *rec_i1, *rec_i1m, *rec_i2, *rec_i3;//incerpolacion center

  int *src_nm, *rec_nm;//number of points to mirror

  float *xweight;//1D data weights only in x-axis
  float **wdat;//data weights of size [nt,nrec]
} acq_t;

#endif
