/* 2D/3D isotropic acoustic forward modelling
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
#include <mpi.h>

void check_cfl(sim_t *sim);

void fdtd_init(sim_t *sim, int flag);
void fdtd_null(sim_t *sim, int flag);
void fdtd_free(sim_t *sim, int flag);
void fdtd_update_v(sim_t *sim, int flag, int it, int adj);
void fdtd_update_p(sim_t *sim, int flag, int it, int adj);

void extend_model_init(sim_t *sim);
void extend_model_free(sim_t *sim);

void computing_box_init(acq_t *acq, sim_t *sim, int adj);
void computing_box_free(sim_t *sim, int adj);

void cpml_init(sim_t *sim);
void cpml_free(sim_t *sim);

void inject_source(sim_t *sim, acq_t *acq, float ***sp, float stf_it);
void extract_wavefield(sim_t *sim, acq_t *acq, float ***sp, float **dat, int it);

void write_data(sim_t *sim, acq_t *acq);

void do_modelling(sim_t *sim, acq_t *acq)
{
  double t0, t_update_v, t_update_p, t_inject_src, t_extract_field;
  float tmp;
  int i1, i2, i3, i1_, i2_, i3_, it;
  FILE *fp;

  if(!getparint("itcheck", &sim->itcheck)) sim->itcheck = sim->nt/2;

  sim->sign_dt = 1;
  check_cfl(sim);
  cpml_init(sim);
  extend_model_init(sim);
  computing_box_init(acq, sim, 0);
  fdtd_init(sim, 1);//flag=1, incident field
  fdtd_null(sim, 1);//flag=1, incident field

  t_update_v = 0.;
  t_update_p = 0.;
  t_inject_src = 0.;
  t_extract_field = 0.;
  for(it=0; it<sim->nt; it++){
    if(iproc==0 && it%100==0) printf("it-----%d\n", it);

    if(iproc==0) t0 = MPI_Wtime();
    fdtd_update_v(sim, 1, it, 0);
    if(iproc==0) t_update_v += MPI_Wtime()-t0;
    
    if(iproc==0) t0 = MPI_Wtime();
    fdtd_update_p(sim, 1, it, 0);
    if(iproc==0) t_update_p += MPI_Wtime()-t0;

    if(iproc==0) t0 = MPI_Wtime();
    if(sim->aniso){
      tmp = sim->stf[it]/3.;
      inject_source(sim, acq, sim->pv1, tmp);
      tmp = 2.*sim->stf[it]/3.;
      inject_source(sim, acq, sim->ph1, tmp);
    }else
      inject_source(sim, acq, sim->p1, sim->stf[it]);
    if(iproc==0) t_inject_src += MPI_Wtime()-t0;

    if(iproc==0) t0 = MPI_Wtime();
    extract_wavefield(sim, acq, sim->p1, sim->dcal, it);
    if(iproc==0) t_extract_field += MPI_Wtime()-t0;
    
    if(iproc==0 && it==sim->itcheck){
      fp = fopen("snapshot.bin", "wb");
      for(i3=0; i3<sim->n3; i3++){
      	i3_ = (sim->n3>1)?i3 + sim->nb:0;
      	for(i2=0; i2<sim->n2; i2++){
      	  i2_ = i2 + sim->nb;
      	  for(i1=0; i1<sim->n1; i1++){
      	    i1_ = i1 + sim->nb;
      	    fwrite(&sim->p1[i3_][i2_][i1_], sizeof(float), 1, fp);
      	  }
      	}
      }
      fclose(fp);
    }//end if
  }
  write_data(sim, acq);

  if(iproc==0) {
    t0 = t_update_v + t_update_p + t_inject_src + t_extract_field;
    FILE *fp = fopen("time_info.txt", "w");
    fprintf(fp, "update_v      \t %e\n", t_update_v);
    fprintf(fp, "update_p      \t %e\n", t_update_p);
    fprintf(fp, "inject_src    \t %e\n", t_inject_src);
    fprintf(fp, "extract_field \t %e\n", t_extract_field);
    fclose(fp);
    
    printf("-------------- elapsed time --------------------\n");
    printf("update_v      \t %e\n", t_update_v);
    printf("update_p      \t %e\n", t_update_p);
    printf("inject_src    \t %e\n", t_inject_src);
    printf("extract_field \t %e\n", t_extract_field);
    printf("total time    \t %e\n", t0);
    printf("------------------------------------------------\n");
  }
  
  extend_model_free(sim);
  fdtd_free(sim, 1);//flag=1, incident field
  cpml_free(sim);
  computing_box_free(sim, 0);

}
