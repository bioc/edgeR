#include <Rmath.h>
#include "ql.h"

/* compute the adjusted deviance and degree of freedom for QL method */

/*
  inputs:
  y raw count matrix
  mu fitted value matrix
  nt number of genes
  nl number of samples
  design design matrix
  nv number of columns of design matrix
  disp dispersion vector

  outputs:
  df adjusted degree of freedom 
  dev adjusted deviance or total.dev
  s2 quasi dispersion estimate by adjusted dev and df
  hatmat hatvalue matrix
  devmat unit deviance matrix
  dfmat individual df matrix
*/

const double thresholdzero=1e-4;

void c_compute_adjust_s2 (double *y, double *mu, int *nt, int *nl, double *design, int *nv, double *disp, double *work, double *df, double *dev, double *s2, double *hatmat, double *devmat, double *dfmat)
{
  int ntag=(*nt), nlib=(*nl), nvar=(*nv);
  double *yptr, *uptr, *hptr, *dptr, *fptr; yptr=y, uptr=mu, hptr=hatmat, dptr=devmat, fptr=dfmat; // set up points and prepare the loop
  double xdesign[nlib*nvar], hatvalues[nlib], weight[nlib], prior=(*work); // prepare the working matrix for the hat values
  double wpt[2], hdp, udp; // store the weights, complementary hat
  for(int tag=0; tag<ntag; ++tag){
    /*
      The preparation of working weight matrix is separated in two steps:
      I. compute the weight for each sample
      II. compute the working weight matrix by multiply the weight looped by columns (trick [i%nlib])
    */
    for(int lib=0; lib<nlib; ++lib) weight[lib]=sqrt((*(uptr+lib))/(1+(*(uptr+lib)*(disp[tag]/prior)))), hatvalues[lib]=0; // compute the weights and reset hatvalues
    for(int i=0; i< (nlib*nvar); ++i) xdesign[i]=design[i] * weight[i % nlib];
    QR_hat(xdesign, nlib, nvar, hatvalues); // compute hat values
    /*
      The loop is based on that a matrix is stored in C as a vector with column based.
      So the tranpose of y and mu matrix are prepared in R level
      in this sense it is easy to visit one gene by simply ++yptr and ++uptr
    */
    dev[tag]=0, df[tag]=0;
    for(int lib=0; lib<nlib; ++lib, ++yptr,++uptr, ++dptr,++fptr, ++hptr){
      compute_weight(*(uptr), disp[tag], prior, wpt);
      udp = compute_unit_nb_deviance(*(yptr), *(uptr), disp[tag]/prior);
      hdp = 1.0 - hatvalues[lib];
      if(hdp < thresholdzero) udp=0.0, hdp=0.0;
      (*dptr) = udp*wpt[0]; // dev matrix
      (*fptr) = hdp*wpt[1]; //df matrix
      (*hptr) = hatvalues[lib]; // hat matrix
      dev[tag] += (*dptr), df[tag] += (*fptr);
    }
    if(df[tag] < thresholdzero){
      s2[tag] = 0, df[tag] = 0;
    }
    else{
      s2[tag] = dev[tag]/df[tag];
    } 
  }
  return;
}

void c_compute_adjust_s2_non (double *y, double *mu, int *nt, int *nl, double *design, int *nv, double *disp, double *work, double *df, double *dev, double *s2)
{
  int ntag=(*nt), nlib=(*nl), nvar=(*nv);
  double *uptr, *yptr; yptr=y, uptr=mu;
  double xdesign[nlib*nvar], hatvalues[nlib], weight[nlib], prior=(*work); 
  double wpt[2], hdp, udp;
  for(int tag=0; tag<ntag; ++tag){
    for(int lib=0; lib<nlib; ++lib) weight[lib]=sqrt((*(uptr+lib))/(1+(*(uptr+lib)*(disp[tag]/prior)))), hatvalues[lib]=0;  
    for(int i=0; i< (nlib*nvar); ++i) xdesign[i]=design[i] * weight[i % nlib];
    QR_hat(xdesign, nlib, nvar, hatvalues); 
    dev[tag]=0, df[tag]=0;
    for(int lib=0; lib<nlib; ++lib, ++yptr, ++uptr){
      compute_weight(*(uptr), disp[tag], prior, wpt);
      udp = compute_unit_nb_deviance(*(yptr), *(uptr), disp[tag]/prior);
      hdp = 1.0 - hatvalues[lib];
      if(hdp < thresholdzero) udp=0.0, hdp=0.0;
      dev[tag] += udp*wpt[0], df[tag] += hdp*wpt[1];
    }
    if(df[tag] < thresholdzero){
      s2[tag] = 0, df[tag] = 0;
    }
    else{
      s2[tag] = dev[tag]/df[tag];
    } 
  }
  return;
}
