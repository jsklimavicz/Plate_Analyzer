#include "multimin.h"


void f(const size_t n, const double *x,void *fparams,double *fval){
  
  double *p = (double *) fparams;
  
  *fval=p[2]*(x[0]-p[0])*(x[0]-p[0]) + p[3]*(x[1]-p[1])*(x[1]-p[1]) + p[4];

}



void df(const size_t n, const double *x,void *fparams,double *grad)
{
  
    double *p = (double *) fparams;
  
    grad[0]=2*p[2]*(x[0]-p[0]);

    grad[1]=2*p[3]*(x[1]-p[1]);

}


void fdf(const size_t n, const double *x,void *fparams,double *fval,double *grad)
{
  
  double *p = (double *) fparams;
  
  *fval=p[2]*(x[0]-p[0])*(x[0]-p[0]) + p[3]*(x[1]-p[1])*(x[1]-p[1]) + p[4];
  
  grad[0]=2*p[2]*(x[0]-p[0]);
  
  grad[1]=2*p[3]*(x[1]-p[1]);
  
}


int main(){
  
  double par[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 };
  
  double x[2]={0,0};
  double minimum;
  
  double xmin[2], xmax[2];
  unsigned type[2];
  
  struct multimin_params optim_par = {.1,1e-2,100,1e-3,1e-5,2,0};
  
  /* unconstrained minimization */
  multimin(2, x, &minimum, NULL,NULL,NULL,&f,&df,&fdf, (void *) par, optim_par);
  
  printf("unconstrained minimum %e at x=%e y=%e\n",minimum,x[0],x[1]);
  
  /* minimum constrained in [-1,2]x(-5,5] */
  type[0]=3;
  xmin[0]=-1;
  xmax[0]=2;
  x[0]=0;

  type[1]=2;
  xmin[1]=-5;
  xmax[1]=5;
  x[1]=0;

  /* increase verbosity */
  optim_par.verbosity=1;

  multimin(2,x,&minimum,type,xmin,xmax,&f,&df,&fdf,par,optim_par);
  
  printf("constrained minimum %e at x=%e y=%e\n",minimum,x[0],x[1]);
  
  
  return 0;
}