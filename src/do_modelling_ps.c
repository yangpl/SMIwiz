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
#include "mpi_info.h"

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

void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);

void write_data(sim_t *sim, acq_t *acq);

void do_modelling_ps(sim_t *sim, acq_t *acq)
{
  int it, itt, nrt;
  int i1, i2, i3;
  int i1_, i2_, i3_;
  int i1_nb, i2_nb, i3_nb;
  int i1_int, i2_int, i3_int;
  int i1p1, i2p1, i3p1;  
  float w1, w2, w3;
  char *subvpfile= "";//model perturbation
  char *subrhofile= "";//model perturbation
  FILE *fp;
  float ***dm1, ***dm2;
  float dvzdt, dvxdt, dvydt, divv, tmp, vpmax1, vpmax2;
  sim_t *sim2;

  sim2 = (sim_t *)malloc(sizeof(sim_t));
  if(!getparint("nr", &sim->nr)) sim->nr = 2; //ratio between fine and coarse grid interval
  if(!getparint("i1start", &sim->i1start)) err("must give i1start= ");
  if(!getparint("i2start", &sim->i2start)) err("must give i2start= ");
  if(!getparint("i1end", &sim->i1end)) err("must give i1end= ");
  if(!getparint("i2end", &sim->i2end)) err("must give i2end= ");

  sim2->ibox = 0;//no computing box for secondary field simulation
  sim2->fm = sim->fm;
  sim2->order = sim->order;
  sim2->nb = sim->nb;
  sim2->n1 = (sim->i1end - sim->i1start)*sim->nr + 1;
  sim2->n2 = (sim->i2end - sim->i2start)*sim->nr + 1;
  sim2->d1 = sim->d1/sim->nr;
  sim2->d2 = sim->d2/sim->nr;
  sim2->n1pad = sim2->n1 + 2*sim2->nb;
  sim2->n2pad = sim2->n2 + 2*sim2->nb;
  sim2->n3 = 1;
  sim2->d3 = 1;
  sim2->n3pad = 1;
  sim2->i3start = 0;
  sim2->i3end = 0;
  if(sim->n3>1){//ny>1, 3D
    if(!getparint("i3start", &sim->i3start)) err("must give i3start= ");
    if(!getparint("i3end", &sim->i3end)) err("must give i3end= ");
    sim2->n3 = (sim->i3end - sim->i3start)*sim->nr + 1;
    sim2->n3pad = sim2->n3 + 2*sim2->nb;
    sim2->d3 = sim->d3/sim->nr;
  }
  sim2->n123 = sim2->n1*sim2->n2*sim2->n3;
  sim2->n123pad = sim2->n1pad*sim2->n2pad*sim2->n3pad;
  if(iproc==0){
    printf("[i1start, i1end]=[%d, %d]\n", sim->i1start, sim->i1end);
    printf("[i2start, i2end]=[%d, %d]\n", sim->i2start, sim->i2end);
    printf("[i3start, i3end]=[%d, %d]\n", sim->i3start, sim->i3end);
    printf("---------------submodel size-------------\n");
    printf("[n1sub, n2sub, n3sub]=[%d, %d, %d]\n", sim2->n1, sim2->n2, sim2->n3);
    printf("[n1padsub, n2padsub, n3padsub]=[%d, %d, %d]\n", sim2->n1pad, sim2->n2pad, sim2->n3pad);
  }
  
  if(!getparstring("subvpfile",&subvpfile)) err("must give subvpfile= ");
  if(!getparstring("subrhofile",&subrhofile)) err("must give subrhofile= ");
  sim2->vp = alloc3float(sim2->n1, sim2->n2, sim2->n3);
  sim2->rho = alloc3float(sim2->n1, sim2->n2, sim2->n3);
  
  fp=fopen(subvpfile,"rb");
  if(fp==NULL) err("cannot open subvpfile=%s", subvpfile);
  if(fread(sim2->vp[0][0],sizeof(float),sim2->n123,fp)!=sim2->n123)
    err("error reading vpfile=%s, size unmatched", subvpfile);
  fclose(fp);

  fp=fopen(subrhofile,"rb");
  if(fp==NULL) err("cannot open subrhofile=%s", subrhofile);
  if(fread(sim2->rho[0][0],sizeof(float),sim2->n123,fp)!=sim2->n123)
    err("error reading rhofile=%s, size unmatched", subrhofile);
  fclose(fp);
  if(iproc==0) printf("-------reading subgrid vp and rho done!-----\n");

  vpmax1 = 0;
  vpmax2 = 0;
  dm1 = alloc3float(sim2->n1, sim2->n2, sim2->n3);//dln(vp)
  dm2 = alloc3float(sim2->n1, sim2->n2, sim2->n3);//dln(rho)
  for(i3_=0; i3_<sim2->n3; i3_++){
    if(sim2->n3==1){
      i3_int = 0;
      i3 = 0;
      i3p1 = 0;
      w3 = 1.;
    }else{
      i3_int = i3_/sim->nr;//integer part
      i3 = sim->i3start + i3_int;
      i3p1 = i3+1;
      w3 = 1. - ((float)i3_/sim->nr-i3_int);
    }
    for(i2_=0; i2_<sim2->n2; i2_++){
      i2_int = i2_/sim->nr;//integer part
      i2 = sim->i2start + i2_int;
      i2p1 = i2+1;
      w2 = 1. - ((float)i2_/sim->nr-i2_int);
      for(i1_=0; i1_<sim2->n1; i1_++){
	i1_int = i1_/sim->nr;//integer part
	i1 = sim->i1start + i1_int;
	i1p1 = i1+1;
	w1 = 1. - ((float)i1_/sim->nr-i1_int);
	//printf("[i3,i2,i1]=[%d,%d,%d]\n", i3, i2, i1);
	
	//trilinear interpolation of the velocity model at coarse grid
	tmp = w3*( w2*(w1*sim->vp[i3][i2][i1] + (1.-w1)*sim->vp[i3][i2][i1p1])
		   + (1.-w2)*( w1*sim->vp[i3][i2p1][i1] + (1.-w1)*sim->vp[i3][i2p1][i1p1]) )
	  +(1.-w3)*( w2*(w1*sim->vp[i3p1][i2][i1] + (1.-w1)*sim->vp[i3p1][i2][i1p1])
		     + (1.-w2)*( w1*sim->vp[i3p1][i2p1][i1] + (1.-w1)*sim->vp[i3p1][i2p1][i1p1]) );
	dm1[i3_][i2_][i1_] = log(sim2->vp[i3_][i2_][i1_]/tmp);//velocity perturbation
	dm2[i3_][i2_][i1_] = 0.;//density perturbation
	//dm1[i3_][i2_][i1_] = (sim2->vp[i3_][i2_][i1_]-tmp)/tmp;//velocity perturbation
	//dm2[i3_][i2_][i1_] = 0.;//density perturbation
	vpmax1 = MAX(vpmax1, sim2->vp[i3_][i2_][i1_]);
	vpmax2 = MAX(vpmax2, tmp);
      }
    }
  }
  nrt = ceil(sim->nr*vpmax1/vpmax2);
  sim2->nt = sim->nt*nrt;
  sim2->dt = sim->dt/nrt;

  
  sim->dcal = alloc2float(sim->nt, acq->nrec);  

  sim->sign_dt = 1;
  if(iproc==0) printf("-------coarse grid info-------\n");
  check_cfl(sim);
  cpml_init(sim);  
  extend_model_init(sim);
  extend_model(sim);
  computing_box_init(acq, sim, 0);
  fdtd_init(sim, 1);//flag=1, incident field
  fdtd_null(sim, 1);//flag=1, incident field

  sim2->sign_dt = 1;
  if(iproc==0) printf("-------fine grid info-------\n");
  check_cfl(sim2);
  cpml_init(sim2);  
  extend_model_init(sim2);
  extend_model(sim2);
  fdtd_init(sim2, 1);//flag=1, incident field

  fp = fopen("secondary.bin", "wb");
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    //1. primary field modelling:
    inject_source(sim, acq, sim->p1, sim->stf[it]);
    fdtd_update_v(sim, 1, it, 0);//flag=1
    fdtd_update_p(sim, 1, it, 0);//flag=1
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);

    //2. secondary field modelling:
    fdtd_null(sim2, 1);//flag=1, incident field
    //2.1 multi-step forward modelling, each small dt' we inject the same source
    // this means the first temporal and spatial derivatives of primary field
    // are coarsely sampled in time (dt=r'*dt'), however their Fourier integration
    //should be equivalent. It is similar to compute/integration gradient every r' time steps
    for(itt=0; itt<nrt; itt++){
      fdtd_update_v(sim2, 1, itt, 0);//flag=1
      //2.1.1 inject source excited by primary field
      /* for(i3_=0; i3_<sim2->n3; i3_++){ */
      /* 	if(sim2->n3>1){ */
      /* 	  i3_int = i3_/sim->nr;//integer part */
      /* 	  i3 = sim->i3start + i3_int + sim->nb; */
      /* 	  i3p1 = i3+1; */
      /* 	  w3 = 1. - ((float)i3_/sim->nr-i3_int); */
      /* 	  i3_nb = i3_ + sim2->nb; */
      /* 	}else{ */
      /* 	  i3_int = 0; */
      /* 	  i3 = 0; */
      /* 	  i3p1 = 0; */
      /* 	  w3 = 1.; */
      /* 	  i3_nb = 0; */
      /* 	} */
      /* 	for(i2_=0; i2_<sim2->n2; i2_++){ */
      /* 	  i2_int = i2_/sim->nr;//integer part */
      /* 	  i2 = sim->i2start + i2_int + sim->nb; */
      /* 	  i2p1 = i2+1; */
      /* 	  w2 = 1. - ((float)i2_/sim->nr-i2_int); */
      /* 	  i2_nb = i2_ + sim2->nb; */
      /* 	  for(i1_=0; i1_<sim2->n1; i1_++){ */
      /* 	    i1_int = i1_/sim->nr;//integer part */
      /* 	    i1 = sim->i1start + i1_int + sim->nb; */
      /* 	    i1p1 = i1+1; */
      /* 	    w1 = 1. - ((float)i1_/sim->nr-i1_int); */
      /* 	    i1_nb = i1_ + sim2->nb; */

      /* 	    //every coarse-grid time step, a constant interpolation at every fine-grid time step */
      /* 	    tmp = dm2[i3_][i2_][i1_];//rho-vp */
      /* 	    dvzdt = w3*( w2*(w1*sim->dvzdt[i3][i2][i1] + (1.-w1)*sim->dvzdt[i3][i2][i1p1]) */
      /* 			 + (1.-w2)*(w1*sim->dvzdt[i3][i2p1][i1] + (1.-w1)*sim->dvzdt[i3][i2p1][i1p1]) ) */
      /* 	      + (1.- w3)*( w2*(w1*sim->dvzdt[i3p1][i2][i1] + (1.-w1)*sim->dvzdt[i3p1][i2][i1p1]) */
      /* 			   + (1.-w2)*(w1*sim->dvzdt[i3p1][i2p1][i1] + (1.-w1)*sim->dvzdt[i3p1][i2p1][i1p1]) ); */
      /* 	    sim2->vz1[i3_nb][i2_nb][i1_nb] += -sim2->dt*tmp*dvzdt; */
      /* 	    dvxdt = w3*( w2*(w1*sim->dvxdt[i3][i2][i1] + (1.-w1)*sim->dvxdt[i3][i2][i1p1]) */
      /* 			 + (1.-w2)*(w1*sim->dvxdt[i3][i2p1][i1] + (1.-w1)*sim->dvxdt[i3][i2p1][i1p1]) ) */
      /* 	      + (1.- w3)*( w2*(w1*sim->dvxdt[i3p1][i2][i1] + (1.-w1)*sim->dvxdt[i3p1][i2][i1p1]) */
      /* 			   + (1.-w2)*(w1*sim->dvxdt[i3p1][i2p1][i1] + (1.-w1)*sim->dvxdt[i3p1][i2p1][i1p1]) ); */
      /* 	    sim2->vx1[i3_nb][i2_nb][i1_nb] += -sim2->dt*tmp*dvxdt; */
      /* 	    if(sim->n3>1) { */
      /* 	      dvydt = w3*( w2*(w1*sim->dvydt[i3][i2][i1] + (1.-w1)*sim->dvydt[i3][i2][i1p1]) */
      /* 			   + (1.-w2)*(w1*sim->dvydt[i3][i2p1][i1] + (1.-w1)*sim->dvydt[i3][i2p1][i1p1]) ) */
      /* 		+ (1.- w3)*( w2*(w1*sim->dvydt[i3p1][i2][i1] + (1.-w1)*sim->dvydt[i3p1][i2][i1p1]) */
      /* 			     + (1.-w2)*(w1*sim->dvydt[i3p1][i2p1][i1] + (1.-w1)*sim->dvydt[i3p1][i2p1][i1p1]) ); */
      /* 	      sim2->vy1[i3_nb][i2_nb][i1_nb] += -sim2->dt*tmp*dvydt; */
      /* 	    } */
      /* 	  }//end for i1 */
      /* 	}//end for i2 */
      /* }//end for i3 */
      fdtd_update_p(sim2, 1, itt, 0);//flag=1
      //2.1.2 inject source excited by primary field
      for(i3_=0; i3_<sim2->n3; i3_++){
      	if(sim2->n3>1){
      	  i3_int = i3_/sim->nr;//integer part
      	  i3 = sim->i3start + i3_int + sim->nb;
      	  i3p1 = i3+1;
      	  w3 = 1. - ((float)i3_/sim->nr-i3_int);
      	  i3_nb = i3_ + sim2->nb;
      	}else{
      	  i3 = 0;
      	  i3p1 = 0;
      	  w3 = 1.;
      	  i3_nb = 0;
      	}
      	for(i2_=0; i2_<sim2->n2; i2_++){
      	  i2_int = i2_/sim->nr;//integer part
      	  i2 = sim->i2start + i2_int + sim->nb;
      	  i2p1 = i2+1;
      	  w2 = 1. - ((float)i2_/sim->nr-i2_int);
      	  i2_nb = i2_ + sim2->nb;
      	  for(i1_=0; i1_<sim2->n1; i1_++){
      	    i1_int = i1_/sim->nr;//integer part
      	    i1 = sim->i1start + i1_int + sim->nb;
      	    i1p1 = i1+1;
      	    w1 = 1. - ((float)i1_/sim->nr-i1_int);
      	    i1_nb = i1_ + sim2->nb;

      	    //every coarse-grid time step, a constant interpolation at every fine-grid time step
      	    tmp = 2.*dm1[i3_][i2_][i1_] + dm2[i3_][i2_][i1_];
	    tmp *= sim->rho[i3_][i2_][i1_]*sim->vp[i3_][i2_][i1_]*sim->vp[i3_][i2_][i1_];
      	    divv = w3*( w2*(w1*sim->divv[i3][i2][i1] + (1.-w1)*sim->divv[i3][i2][i1p1])
      			 + (1.-w2)*(w1*sim->divv[i3][i2p1][i1] + (1.-w1)*sim->divv[i3][i2p1][i1p1]) )
      	      + (1.- w3)*( w2*(w1*sim->divv[i3p1][i2][i1] + (1.-w1)*sim->divv[i3p1][i2][i1p1])
      			   + (1.-w2)*(w1*sim->divv[i3p1][i2p1][i1] + (1.-w1)*sim->divv[i3p1][i2p1][i1p1]) );
      	    sim2->p1[i3_nb][i2_nb][i1_nb] += -sim->dt*tmp*divv;
      	  }//end for i1
      	}//end for i2
      }//end for i3

    }//end for itt
    fwrite(&sim2->p1[0][0][0], sim2->n123*sizeof(float), 1, fp);

    //2.2. adding them together
    for(i3=sim->i3start; i3<=sim->i3end; i3++){
      i3_nb = (sim->n3>1)?i3+sim->nb:0;
      i3_ = (sim->n3>1)?(i3 - sim->i3start)*sim->nr + sim2->nb:0;
      for(i2=sim->i2start; i2<=sim->i2end; i2++){
    	i2_nb = i2+sim->nb;
    	i2_ = (i2 - sim->i2start)*sim->nr + sim2->nb;
    	for(i1=sim->i1start; i1<=sim->i1end; i1++){
    	  i1_nb = i1 + sim->nb;
    	  i1_ = (i1 - sim->i1start)*sim->nr + sim2->nb;
	  
    	  sim->p1[i3_nb][i2_nb][i1_nb] += sim2->p1[i3_][i2_][i1_];
    	  sim->vz1[i3_nb][i2_nb][i1_nb] += sim2->vz1[i3_][i2_][i1_];
    	  sim->vx1[i3_nb][i2_nb][i1_nb] += sim2->vx1[i3_][i2_][i1_];
    	  if(sim->n3>1) sim->vy1[i3_nb][i2_nb][i1_nb] += sim2->vy1[i3_][i2_][i1_];
    	}
      }
    }
    
  }//end for it
  fclose(fp);
  write_data(sim, acq);
  
  extend_model_close(sim);
  fdtd_close(sim, 1);//flag=1, incident field
  cpml_close(sim);
  computing_box_close(sim, 0);

  extend_model_close(sim2);
  fdtd_close(sim2, 1);//flag=1, incident field
  cpml_close(sim2);
  free3float(sim2->vp);
  free3float(sim2->rho);
  free(sim2);

  free3float(dm1);
  free3float(dm2);
  free2float(sim->dcal);
}
