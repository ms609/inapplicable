#define USE_RINTERNALS


#include <Rmath.h>
#include <math.h>
#include <R.h> 

#include <Rinternals.h>

void fitchinapp(int *dat1, int *dat2, int *nr, int *pars, double *weight, int *inapp, double *w){
  int k, tmp;
  for(k = 0; k < (*nr); k++){
    tmp = dat1[k] & dat2[k];
    if(!tmp){
      tmp = dat1[k] | dat2[k];
      if ((dat1[k] == *inapp) || (dat2[k] == *inapp)) { // Do not increment pscore
      } else {
        if ((*inapp) & tmp) { // Only delete inapplicable token if it is present
          tmp = tmp - (*inapp);
        }
        (pars[k])++;
        (*w) += weight[k];
      }
    }
    dat1[k] = tmp;
  } 
}

void fitchbridge(int *dat, int *nr, int *pars, int *node, int *edge, int *nl, double *weight, int *inapp, double *pvec, double *pscore)
{
  int i, ni, k;
  ni = 0;
  for (i=0; i< *nl; i++) {
    if (ni == node[i]){
      pvec[ni-1L] += pvec[edge[i]-1L];
      fitchinapp(&dat[(ni-1L) * (*nr)], &dat[(edge[i]-1L) * (*nr)], nr, pars, weight, inapp, &pvec[(ni-1L)]);
    } else {
      ni = node[i];
      pvec[(ni-1L)] += pvec[(edge[i]-1L)];         
      for(k = 0; k < (*nr); k++) dat[(ni-1L)*(*nr) + k] = dat[(edge[i]-1L)*(*nr) + k];                     
    }
  }
  pscore[0]=pvec[ni-1];
}

SEXP FITCHI(SEXP dat, SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP weight, SEXP mx, SEXP q, SEXP inapp){   
  int *data, *nr=INTEGER(nrx), m=INTEGER(mx)[0], i, n=INTEGER(q)[0];   
  double *pvtmp;
  SEXP DAT, pars, pvec, pscore, RESULT;
  PROTECT(RESULT = allocVector(VECSXP, 4L));
  PROTECT(pars = allocVector(INTSXP, *nr));
  PROTECT(pscore = allocVector(REALSXP, 1L));
  PROTECT(DAT = allocMatrix(INTSXP, nr[0], m));
  PROTECT(pvec = allocVector(REALSXP, m));
  pvtmp = REAL(pvec);
  data = INTEGER(DAT);
  for(i=0; i<m; i++) pvtmp[i] = 0.0;
  for(i=0; i<*nr; i++) INTEGER(pars)[i] = 0L;
  REAL(pscore)[0]=0.0;
  for(i=0; i<(*nr * n); i++) data[i] = INTEGER(dat)[i];
  
  fitchbridge(data, nr, INTEGER(pars), INTEGER(node), INTEGER(edge), INTEGER(l), REAL(weight), INTEGER(inapp), pvtmp, REAL(pscore));
  
  SET_VECTOR_ELT(RESULT, 0, pscore);
  SET_VECTOR_ELT(RESULT, 1, pars);
  SET_VECTOR_ELT(RESULT, 2, DAT);
  SET_VECTOR_ELT(RESULT, 3, pvec);
  UNPROTECT(5);
  return(RESULT); 
}