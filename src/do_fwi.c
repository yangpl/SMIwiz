/* nonlinear optimization for full waveform inversion (FWI)
 *---------------------------------------------------------------------------
 *  Copyright (c) Pengliang Yang, 2020, Harbin Institute of Technology, China
 *  Copyright (c) Pengliang Yang, 2018, University Grenoble Alpes, France
 *  Homepage: https://yangpl.wordpress.com
 *  E-mail: ypl.2100@gmail.com
 *--------------------------------------------------------------------------*/
#include <mpi.h>
#include "cstd.h"
#include "acq.h"
#include "opt.h"
#include "sim.h"
#include "fwi.h"

float l2norm(int n, float *a);
void flipsign(int n, float *a, float *b);
void lbfgs_save(int n, float *x, float *g, float **sk, float **yk, opt_t *opt);
void lbfgs_update(int n, float *x, float *g, float **sk, float **yk, opt_t *opt);
void lbfgs_descent(int n, float *g, float *d, float **sk, float **yk, float *q, float *rho, float *alp, opt_t *opt);
bool lbfgs_descent1(int n, float *g, float *q, float *rho, float *alp, float **sk, float **yk, opt_t *opt);
void lbfgs_descent2(int n, float *g, float *q, float *rho, float *alp, float **sk, float **yk, opt_t *opt);
void boundx(float *x, int n, float *xmin, float *xmax);
void line_search(int n, //dimension of x
		 float *x, //input vector x
		 float *g, //gradient of misfit function
		 float *d, //descent direction
		 opt_fg fg, //subroutine to evaluation function and gradient
		 opt_t *opt); //pointer of l-BFGS optimization parameters
void cg_solve(int n, //dimension of x
	      float *x, //input vector x
	      float *g, //gradient of misfit function
	      float *d, //descent direction
	      opt_Hv Hv, //subroutine to evaluation function and gradient
	      opt_t *opt); //pointer of l-BFGS optimization parameters


void fg_fwi_init(sim_t *sim_, acq_t *acq_, fwi_t *fwi_);
void fg_fwi_free(sim_t *sim);
void fg_mod_reg(sim_t *sim, fwi_t *fwi, float ***grad, float *x, float *g);
float fg_fwi(float *x, float *g);
void precondition(sim_t *sim, fwi_t *fwi, float *x);


void do_fwi(sim_t *sim, acq_t *acq)
/*< perform FWI using l-BFGS optimization >*/
{
  int i1, i2, i3, ipar, j;
  float tmp, fcost;
  opt_t *opt; //pointer for opt_t parameters
  fwi_t *fwi;
  char *bathyfile;
  FILE *fp;

  opt = malloc(sizeof(opt_t));
  if(!getparint("niter", &opt->niter)) opt->niter=50;//maximum number of iterations
  if(!getparint("nls", &opt->nls)) opt->nls=20;//maximum number of line searches
  if(!getparfloat("tol", &opt->tol)) opt->tol=1e-8;//convergence tolerance 
  if(!getparint("npair", &opt->npair)) opt->npair=5; //l-BFGS memory length
  if(!getparfloat("c1", &opt->c1)) opt->c1=1e-4; //Nocedal value for Wolfe condition
  if(!getparfloat("c2", &opt->c2)) opt->c2=0.9;  //Nocedal value for Wolfe condition
  if(!getparint("bound", &opt->bound)) opt->bound=1;//use bounds or not
  if(!getparint("preco", &opt->preco)) opt->preco=0;//1=precondition; 0=not
  opt->verb = (iproc==0)?1:0; //other process are silent.

  fwi = malloc(sizeof(fwi_t));
  fwi->bathy=alloc2float(sim->n2, sim->n3);
  fwi->ibathy=alloc2int(sim->n2, sim->n3);
  if(!getparstring("bathyfile",&bathyfile)){
    memset(fwi->bathy[0], 0, sim->n2*sim->n3*sizeof(float));
    memset(fwi->ibathy[0], 0, sim->n2*sim->n3*sizeof(int));
  }else{
    fp=fopen(bathyfile,"rb");
    if(fp==NULL) err("cannot open bathyfile=%s",bathyfile);
    if(fread(fwi->bathy[0],sizeof(float),sim->n2*sim->n3,fp)!=sim->n2*sim->n3) 
      err("error reading bathyfile=%s", bathyfile);
    fclose(fp);
    for(i3=0; i3<sim->n3; i3++){
      for(i2=0; i2<sim->n2; i2++){
	fwi->ibathy[i3][i2]=NINT(fwi->bathy[i3][i2]/sim->d1);
      }
    }
  }
  if(!getparfloat("transition", &fwi->transition)) fwi->transition=0;//thinkness of transition zone
  fwi->itransition = NINT(fwi->transition/sim->d1);
  
  if(!getparint("family", &fwi->family)) fwi->family = 1;//1=vp-rho'; 2=vp-Ip
  if(!getparint("npar", &fwi->npar)) fwi->npar = 1;//number of unknown parameters
  if(!(j=countparval("idxpar"))) err("must give idxpar= vector");  
  if(j!= fwi->npar) err("must have length[idxpar]=%d", fwi->npar);
  fwi->idxpar = alloc1int(fwi->npar);
  getparint("idxpar", fwi->idxpar);
  //family=1: idxpar[0] =1, vp; idxpar[1]=2, rho; family=2: idxpar[0] =1, vp; idxpar[1]=2, Ip
  if(!getparfloat("rhomin", &fwi->rhomin)) fwi->rhomin = 1000;
  if(!getparfloat("rhomax", &fwi->rhomax)) fwi->rhomax = 3000;
  if(!getparfloat("vpmin", &fwi->vpmin)) fwi->vpmin = 1000;
  if(!getparfloat("vpmax", &fwi->vpmax)) fwi->vpmax = 6000;
  
  if(!getparint("rwi", &fwi->rwi)) fwi->rwi = 0;//0=FWI, 1=RWI
  if(fwi->rwi==1){
    fwi->family = 2;//vp-ip parametrization
    opt->preco = 0;
  }
  fwi->preco = opt->preco;//copy option of precondtioner
  if(iproc==0) {
    printf("rwi=%d (1=RWI; 0=FWI)\n", fwi->rwi);
    printf("family=%d (1=vp-rho; 2=vp-ip)\n", fwi->family);
    printf("npar=%d\n", fwi->npar);
    printf("parameter bound: [rhomin, rhomax]=[%g,%g]\n", fwi->rhomin, fwi->rhomax);
    printf("parameter bound: [vpmin, vpmax]=[%g,%g]\n", fwi->vpmin, fwi->vpmax);
  }
  
  fwi->n = sim->n123*fwi->npar;
  opt->x = alloc1float(fwi->n);
  opt->g = alloc1float(fwi->n);
  opt->d = alloc1float(fwi->n);
  if(opt->preco) opt->pg = alloc1float(fwi->n);
  opt->sk= alloc2float(fwi->n, opt->npair);
  opt->yk= alloc2float(fwi->n, opt->npair);
  if(opt->bound){
    opt->xmin = alloc1float(fwi->n);
    opt->xmax = alloc1float(fwi->n);
    if(fwi->rwi){
      for(ipar=0; ipar<fwi->npar; ipar++){
	for(j=0; j<sim->n123; j++){
	  if(fwi->idxpar[ipar]==1){//min and max for ln(vp)
	    opt->xmin[j+ipar*sim->n123] = log(fwi->vpmin);
	    opt->xmax[j+ipar*sim->n123] = log(fwi->vpmax);
	  }else if(fwi->idxpar[ipar]==2){//min and max for d[ln(vp)]
	    opt->xmin[j+ipar*sim->n123] = -0.4;//min value for dlnvp=dm/vp
	    opt->xmax[j+ipar*sim->n123] = 0.4;//max value for dlnvp=dm/vp
	  }//end if
	}//end for j
      }//end for ipar
    }else{
      for(ipar=0; ipar<fwi->npar; ipar++){
	for(j=0; j<sim->n123; j++){
	  if(fwi->family==1){//1=vp; 2=rho
	    if(fwi->idxpar[ipar]==1) {
	      opt->xmin[j+ipar*sim->n123] = log(fwi->vpmin);
	      opt->xmax[j+ipar*sim->n123] = log(fwi->vpmax);
	    }
	    if(fwi->idxpar[ipar]==2) {
	      opt->xmin[j+ipar*sim->n123] = log(fwi->rhomin);
	      opt->xmax[j+ipar*sim->n123] = log(fwi->rhomax);
	    }
	  }else if(fwi->family==2){//1=vp; 2=ip
	    if(fwi->idxpar[ipar]==1) {
	      opt->xmin[j+ipar*sim->n123] = log(fwi->vpmin);
	      opt->xmax[j+ipar*sim->n123] = log(fwi->vpmax);
	    }
	    if(fwi->idxpar[ipar]==2) {
	      opt->xmin[j+ipar*sim->n123] = log(fwi->rhomin*fwi->vpmin);
	      opt->xmax[j+ipar*sim->n123] = log(fwi->rhomax*fwi->vpmax);
	    }
	  }//end if
	}//end for j
      }//end for ipar
    }//end if
  }//end if
  fg_fwi_init(sim, acq, fwi);
	
  if(sim->mode==1||sim->mode==4) {//FWI or simply output FWI gradient
    //initialize opt->x[]
    if(fwi->rwi){
      for(ipar=0; ipar<fwi->npar; ipar++){
	for(i3=0; i3<sim->n3; i3++){
	  for(i2=0; i2<sim->n2; i2++){
	    for(i1=0; i1<sim->n1; i1++){
	      j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	      if(fwi->idxpar[ipar]==1) opt->x[j] = log(sim->vp[i3][i2][i1]);
	      if(fwi->idxpar[ipar]==2) opt->x[j] = sim->dm[i3][i2][i1];//dln(ip)
	    }//end for i1
	  }//end for i2
	}//end for i3
      }//end for ipar
    }else{
      for(ipar=0; ipar<fwi->npar; ipar++){
	for(i3=0; i3<sim->n3; i3++){
	  for(i2=0; i2<sim->n2; i2++){
	    for(i1=0; i1<sim->n1; i1++){
	      j = i1 + sim->n1*(i2 + sim->n2*(i3 + sim->n3*ipar));
	      if(fwi->family==1){//vp-rho
		if(fwi->idxpar[ipar]==1) opt->x[j] = log(sim->vp[i3][i2][i1]);
		if(fwi->idxpar[ipar]==2) opt->x[j] = log(sim->rho[i3][i2][i1]);
	      }
	      if(fwi->family==2){//vp-ip
		if(fwi->idxpar[ipar]==1) opt->x[j] = log(sim->vp[i3][i2][i1]);
		if(fwi->idxpar[ipar]==2) opt->x[j] = log(sim->rho[i3][i2][i1]*sim->vp[i3][i2][i1]);
	      }
	    }//end for i1
	  }//end for i2
	}//end for i3
      }//end for ipar
    }
    fcost = fg_fwi(opt->x, opt->g);
  }

  if(sim->mode==1){//FWI
    //initialize all counters
    opt->f0 = fcost;
    opt->fk = fcost;
    opt->igrad = 0;
    opt->kpair = 0;
    opt->ils = 0;
    if(opt->verb){
      opt->gk_norm = l2norm(fwi->n, opt->g);
      fp=fopen("iterate.txt","w");
      fprintf(fp,"==========================================================\n");
      fprintf(fp,"l-BFGS memory length: %d\n",opt->npair);
      fprintf(fp,"Maximum number of iterations: %d\n",opt->niter);
      fprintf(fp,"Convergence tolerance: %3.2e\n", opt->tol);
      fprintf(fp,"maximum number of line search: %d\n",opt->nls);
      fprintf(fp,"initial step length: alpha=%g\n",opt->alpha);
      fprintf(fp,"==========================================================\n");
      fprintf(fp,"iter    fk       fk/f0      ||gk||    alpha    nls   ngrad\n");
      fclose(fp);
    }
    //l-BFGS optimization 
    for(fwi->iter=0; fwi->iter<opt->niter; fwi->iter++){
      if(opt->verb){
	printf("==========================================================\n");
	printf("# iter=%d  fk/f0=%g\n", fwi->iter,opt->fk/opt->f0);
	opt->gk_norm = l2norm(fwi->n, opt->g);
	fp=fopen("iterate.txt","a");
	fprintf(fp,"%3d   %3.2e  %3.2e   %3.2e  %3.2e  %3d  %4d\n",
		fwi->iter,opt->fk,opt->fk/opt->f0,opt->gk_norm,opt->alpha,opt->ils,opt->igrad);
	fclose(fp);
      }
      if(fwi->iter==0){//first iteration, no stored gradient
	if(opt->preco){
	  memcpy(opt->pg, opt->g, fwi->n*sizeof(float));
	  precondition(sim, fwi, opt->pg);
	  flipsign(fwi->n, opt->pg, opt->d);
	}else
	  flipsign(fwi->n, opt->g, opt->d);//descent direction=-gradient
      }else{
	// allocate vector q and  rho, alpha in lbfgs_descent1()
	// they will be freed in lbfgs_descent2()
	opt->q = alloc1float(fwi->n);
	opt->rho = alloc1float(opt->kpair);
	opt->alp = alloc1float(opt->kpair);

	lbfgs_update(fwi->n, opt->x, opt->g, opt->sk, opt->yk, opt);
	//lbfgs_descent(n, opt->g, opt->d, opt->sk, opt->yk, opt);
	//---------------------------------------------------option 1-----------
	//1st loop of two-loop recursion
	opt->loop1=lbfgs_descent1(fwi->n, opt->g, opt->q, opt->rho, opt->alp, opt->sk, opt->yk, opt);
	if(opt->preco) precondition(sim, fwi, opt->q);
	//2nd loop of two-loop recursion if 1st loop was done
	if(opt->loop1) lbfgs_descent2(fwi->n, opt->g, opt->q, opt->rho, opt->alp, opt->sk, opt->yk, opt);
	flipsign(fwi->n, opt->q, opt->d); //descent direction d=-q where q=H^{-1}g

	//free allocated vector q and rho, alpha in lbfgs_descent2()
	free1float(opt->q);
	free1float(opt->alp);
	free1float(opt->rho);
      } 
      lbfgs_save(fwi->n, opt->x, opt->g, opt->sk, opt->yk, opt);
      line_search(fwi->n, opt->x, opt->g, opt->d, fg_fwi, opt);
      
      if(opt->ls_fail){
	if(opt->verb) {
	  fp=fopen("iterate.txt","a");
	  fprintf(fp, "==>Line search failed!\n");
	  fclose(fp);
	}
	break;
      }
      //not break, then line search succeeds or descent direction accepted
      if(opt->verb) {
	if(fwi->rwi){
	  for(ipar=0; ipar<fwi->npar; ipar++){
	    if(fwi->idxpar[ipar]==1){
	      fp = fopen("param_final_m0","wb");
	      fwrite(&opt->x[ipar*sim->n123], sim->n123*sizeof(float), 1, fp);
	      fclose(fp);
	    }
	    if(fwi->idxpar[ipar]==2){
	      fp = fopen("param_final_dm","wb");
	      fwrite(&opt->x[ipar*sim->n123], sim->n123*sizeof(float), 1, fp);
	      fclose(fp);
	    }
	  }

	}else{
	  fp = fopen("param_final","wb");
	  for(j=0; j<fwi->n; j++){
	    tmp = exp(opt->x[j]);
	    fwrite(&tmp, sizeof(float), 1, fp);
	  }
	  fclose(fp);
	}

	fp=fopen("gradient_final","wb");
	fwrite(opt->g, fwi->n*sizeof(float),1,fp);
	fclose(fp);
	
	if(sim->n3==1){//we only store intermediate models and gradients in 2D 
	  if(fwi->iter==0) fp = fopen("param_iter","wb");
	  else             fp = fopen("param_iter","ab");
	  for(j=0; j<fwi->n; j++){	    
	    tmp = exp(opt->x[j]);
	    if(fwi->rwi) tmp = opt->x[j];
	    fwrite(&tmp, sizeof(float), 1, fp);
	  }
	  fclose(fp);

	  if(fwi->iter==0) fp = fopen("gradient_iter","wb");
	  else             fp = fopen("gradient_iter","ab");
	  fwrite(opt->g, fwi->n*sizeof(float),1,fp);
	  fclose(fp);
	}//end if n3>1
      }

      if(opt->fk < opt->tol * opt->f0){//here we assume misfit function is always positive
	if(opt->verb){
	  fp=fopen("iterate.txt","a");
	  fprintf(fp, "==>Convergence reached!\n");
	  fclose(fp);
	}
	break;
      }
      fflush(stdout);
    } 
    if(opt->verb && fwi->iter==opt->niter) {
      fp=fopen("iterate.txt","a");
      fprintf(fp, "==>Maximum iteration number reached!\n");
      fclose(fp);
    }
  }

  fg_fwi_free(sim);
  free1int(fwi->idxpar);
  free2float(fwi->bathy);
  free2int(fwi->ibathy);
  free(fwi);
  
  free1float(opt->x);
  free1float(opt->g);
  free1float(opt->d);
  if(opt->preco) free1float(opt->pg);
  free2float(opt->sk);
  free2float(opt->yk);
  if(opt->bound){
    free1float(opt->xmin);
    free1float(opt->xmax);
  }
  free(opt);
}

