/* History: Feb 20 2019 Initial coding
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define DEBUG 0
#define EQUAL_EPS 1e-6

#define MAX( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define MIN( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}

/* Function to allocate memory for a double vector */
static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */

static int * iVec_alloc(n, initFlag, initVal)
int n, initFlag, initVal;
{
  int i, *ret, *p;

  ret = (int *) malloc(n*sizeof(int));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: iVec_alloc */

static void copy_iVec(v1, v2, n)
int *v1, *v2, n;
{
  int i, *p1, *p2;
  
  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) *p1 = *p2;

} /* END: copy_iVec */

static void permute_iVec(vec, n, ret)
int *vec, *ret, n;
{
  /*
  To shuffle an array a of n elements (indices 0..n-1):
  for i from n-1 downto 1 do
     j = random integer such that 0 <= j <= i
     exchange a[j] and a[i]
  */

  int i, j, tmp;

  copy_iVec(ret, vec, n);
  
  for (i=n-1; i>0; i--) {
    j = floor(runif(0.0, i+1.0));
    if (j > i) j = i;

    tmp    = ret[i];
    ret[i] = ret[j];
    ret[j] = tmp;
  }

} /* END: permute_iVec */

static void M_fn(t, t_all, N, Nfail, retVec)
double *t, *t_all, *retVec;
int N, Nfail;
{
  int i, k, ii=0; 
  double tmin, tk, tkm1, ti, test;

  /* Return by row */
  for (k=1; k<Nfail+1; k++) {
    tk   = t_all[k];
    tkm1 = t_all[k-1];
    test = tkm1 - EQUAL_EPS;

    for (i=0; i<N; i++) {
      ti = t[i];
      if (ti >= test) {
        tmin = MIN(tk, ti);
        retVec[ii] = tmin - tkm1;
      } else {
        retVec[ii] = 0.0;
      }
      ii++;
    }
  }

} /* END: M_fn */

void C_M_fn(t, t_all, pN, pNfail, retCode, retVec)
double *t, *t_all, *retVec;
int *pN, *pNfail, *retCode;
{
  *retCode = -1;

  /* retVec is by row */
  M_fn(t, t_all, *pN, *pNfail, retVec);

  *retCode = 0;

  return;

} /* END: C_M_fn */

static double loglik_p(n, X, Zt, exst1, nevents, xt1, xt12, trt1LogEffect, 
             trt1Log1mEffect, lambda, h_nr)
int n, nevents, *exst1;
double *X, *Zt, *xt1, *xt12, *trt1LogEffect, *trt1Log1mEffect, lambda, h_nr;
{
  int i;
  double ret, r1=0.0, r2, r3=0.0, r4=0.0, r5=0.0, r6=0.0, r7=0.0, Zti, oneMinusZti;

  r2 = (double) nevents;
  
  for (i=0; i<n; i++) {
    Zti         = Zt[i];
    oneMinusZti = 1.0 - Zti;
    r1 += Zti*exst1[i];
    r3 += Zti*xt1[i];
    r4 += Zti*xt12[i];
    r5 += oneMinusZti*X[i];
    r6 += Zti*trt1LogEffect[i];
    r7 += oneMinusZti*trt1Log1mEffect[i];
  }

  ret = r1*log(lambda) + r2*log(h_nr) - h_nr*(r3 + lambda*r4 + r5) + r6 + r7;

  return(ret);

} /* END: loglik_p */

static void EM_setupObjP(n, X, eventStatus, t1, xt1, xt12, exst1, Xminust1, nevents)
int n, *eventStatus, *exst1, *nevents;
double *X, t1, *xt1, *xt12, *Xminust1;
{
  int i, tmp, evi, sum=0;
  double Xi;

  for (i=0; i<n; i++) { 
     Xi                  = X[i];
     tmp                 = 0;
     evi                 = eventStatus[i];
     if (Xi > t1) tmp    = 1;
     Xminust1[i]         = Xi - t1;
     exst1[i]            = tmp*evi;
     sum                += evi;
     if (Xi < t1) {
       xt1[i] = Xi;
     } else {
       xt1[i] = t1;
     }
     xt12[i]             = Xminust1[i]*tmp;
  }
  *nevents = sum;

  return;

} /* END: EM_setupObjP */

static void EM_setupTrtP(n, trt, logEff, log1mEff, trt1LogEffect, trt1Log1mEffect)
int n, *trt;
double logEff, log1mEff, *trt1LogEffect, *trt1Log1mEffect;
{
  int i, trti;

  for (i=0; i<n; i++) { 
    trti                = trt[i];
    trt1LogEffect[i]    = trti*logEff;
    trt1Log1mEffect[i]  = trti*log1mEff;
  }

  return;

} /* END: EM_setupTrtP */

static void setupZt(n, trt, ret)
int n, *trt;
double *ret;
{
  int i;

  for (i=0; i<n; i++) { 
    if (trt[i]) {
      ret[i] = 0.5;
    } else {
      ret[i] = 0.0;
    }
  }

  return;

} /* END: EM_setupTrtP */


static int EM_loopP(n, nevents, maxiter, trt, eventStatus, exst1, X, xt1, xt12, 
                    effect_p, trt1LogEffect, trt1Log1mEffect, stopTol, t1, print, 
                    Xminust1, ret_h_nr, ret_lambda, Zt)
int n, nevents, maxiter, *trt, *eventStatus, *exst1, print;
double *X, *xt1, *xt12, effect_p, *trt1LogEffect, *trt1Log1mEffect, stopTol, *Zt,
       t1, *Xminust1, *ret_h_nr, *ret_lambda;
{
  int i, iter=0, conv=0, evi;
  double r1, r2, r3, r4, r5, lambda, h_nr, Zti, Xi, dtmp, vec, vec2, pdfr, pdfnr,
         oneMinusEffect, loglik, loglik0=-1.0, diff;

  r2             = (double) nevents;
  oneMinusEffect = 1.0 - effect_p;

  while(iter < maxiter) {
    iter++;
    r1 = 0.0; 
    r3 = 0.0;
    r4 = 0.0;
    r5 = 0.0;
    for (i=0; i<n; i++) {
      Zti = Zt[i];
      r1 += Zti*exst1[i];
      r3 += Zti*xt1[i];
      r4 += Zti*xt12[i];
      r5 += (1.0 - Zti)*X[i];   
    }
    h_nr   = (r2 - r1)/(r3 + r5);
    lambda = r1/(r4*h_nr);

    for (i=0; i<n; i++) {
      if (trt[i]) {
        Xi   = X[i];
        evi  = eventStatus[i];
        dtmp = lambda*h_nr;
        if (!evi) {
          vec = 1.0;
        } else {
          vec  = dtmp;
        }
        if (Xi < t1) {
          vec  = h_nr;
          vec2 = -h_nr*Xi;
        } else {
          vec2 = -h_nr*t1 - dtmp*Xminust1[i];
        }  
        pdfr  = vec*exp(vec2);
        vec   = exp(-h_nr*Xi);
        if (evi) vec = vec*h_nr;
        pdfnr = vec;
        vec   = effect_p*pdfr;
        Zt[i] = vec/(vec + oneMinusEffect*pdfnr);
      }
    }

    loglik = loglik_p(n, X, Zt, exst1, nevents, xt1, xt12, trt1LogEffect, 
                      trt1Log1mEffect, lambda, h_nr);

    if (iter > 1) {
      diff = fabs(loglik - loglik0);
      if (diff <= stopTol) {  
        conv = 1;
        break;
      }
    }
    loglik0 = loglik;
    if (print > 1) {
      if (iter > 1) {
        Rprintf("Iter=%d, loglike=%11.4f, diff=%g\n", iter, loglik, diff);
      } else {
        Rprintf("Iter=%d, loglike=%11.4f\n", iter, loglik);
      }     
    }    
  
  } /* END: while */

  if (print) {
    if (conv) {
      Rprintf("EM algoritm converged in %d iteration(s)\n", iter);
    } else {
      Rprintf("EM algoritm did not converge\n");
    }
  }

  *ret_lambda = lambda;
  *ret_h_nr   = h_nr;

  return(conv);

} /* END: EM_loopP */

/* Zt is input and output */
static int EM_mainP(n, X, trt, eventStatus, effect_p, t1, stopTol, maxiter, print, 
                    ret_lambda, ret_h_nr, Zt)
int n, print, maxiter, *trt, *eventStatus;
double *X, *Zt, effect_p, t1, stopTol, *ret_lambda, *ret_h_nr;
{
  double *Xminust1, *xt1, *xt12, *trt1LogEffect, *trt1Log1mEffect, logEff, log1mEff;
  int nevents=0, conv=0, *exst1;

  Xminust1        = dVec_alloc(n, 0, 0.0);
  exst1           = iVec_alloc(n, 0, 0.0);
  xt1             = dVec_alloc(n, 0, t1);
  xt12            = dVec_alloc(n, 0, 0.0);
  trt1LogEffect   = dVec_alloc(n, 0, 0.0);
  trt1Log1mEffect = dVec_alloc(n, 0, 0.0);
  logEff          = log(effect_p);
  log1mEff        = log(1.0 - effect_p);
 
  EM_setupObjP(n, X, eventStatus, t1, xt1, xt12, exst1, Xminust1, &nevents);
  EM_setupTrtP(n, trt, logEff, log1mEff, trt1LogEffect, trt1Log1mEffect);

  conv = EM_loopP(n, nevents, maxiter, trt, eventStatus, exst1, X, xt1, xt12, 
                    effect_p, trt1LogEffect, trt1Log1mEffect, stopTol, t1, print, 
                    Xminust1, ret_h_nr, ret_lambda, Zt);

  free(Xminust1);
  free(exst1);
  free(xt1);
  free(xt12);
  free(trt1LogEffect);
  free(trt1Log1mEffect);

  return(conv);

} /* END: EM_mainP */

void C_EM_mainP(pn, X, trt, eventStatus, peffect_p, pt1, pstopTol, pmaxiter, pprint, 
                    ret_conv, ret_lambda, ret_h_nr, Zt)
int *pn, *pprint, *pmaxiter, *trt, *eventStatus, *ret_conv;
double *X, *Zt, *peffect_p, *pt1, *pstopTol, *ret_lambda, *ret_h_nr;
{

  *ret_conv = EM_mainP(*pn, X, trt, eventStatus, *peffect_p, *pt1, *pstopTol, 
                *pmaxiter, *pprint, ret_lambda, ret_h_nr, Zt);

  return;

} /* END: C_EM_mainP */


static int ReRandP(num_rand, n, X, trt, eventStatus, effect_p, t1, stopTol, maxiter, print, 
                   lambda_obs, ret_p)
int num_rand, n, *trt, *eventStatus, maxiter, print;
double *X, effect_p, t1, stopTol, *ret_p, lambda_obs;
{
  double *Xminust1, *xt1, *xt12, *trt1LogEffect, *trt1Log1mEffect, logEff, log1mEff;
  double *Zt, lambda, h_nr;
  int iter, conv, nevents=0, *exst1, *trtPermute, sumNrand=0, sumGTobs=0;

  *ret_p          = -1.0;
  Xminust1        = dVec_alloc(n, 0, 0.0);
  exst1           = iVec_alloc(n, 0, 0.0);
  xt1             = dVec_alloc(n, 0, t1);
  xt12            = dVec_alloc(n, 0, 0.0);
  trt1LogEffect   = dVec_alloc(n, 0, 0.0);
  trt1Log1mEffect = dVec_alloc(n, 0, 0.0);
  trtPermute      = iVec_alloc(n, 0, 0.0);
  Zt              = dVec_alloc(n, 0, 0.0); 
  logEff          = log(effect_p);
  log1mEff        = log(1.0 - effect_p);
 
  EM_setupObjP(n, X, eventStatus, t1, xt1, xt12, exst1, Xminust1, &nevents);

  for (iter=0; iter<num_rand; iter++) {
    permute_iVec(trt, n, trtPermute);
    EM_setupTrtP(n, trtPermute, logEff, log1mEff, trt1LogEffect, trt1Log1mEffect);
    setupZt(n, trtPermute, Zt);
    conv = EM_mainP(n, X, trtPermute, eventStatus, effect_p, t1, stopTol, maxiter, print, 
                    &lambda, &h_nr, Zt);
    if (conv) {
      sumNrand++;
      if (lambda > lambda_obs) sumGTobs++;
    }
  }

  if (sumNrand) *ret_p = 1.0 - ((double) sumGTobs)/((double) sumNrand);

  free(Xminust1);
  free(exst1);
  free(xt1);
  free(xt12);
  free(trt1LogEffect);
  free(trt1Log1mEffect);
  free(trtPermute);
  free(Zt);

  return(sumNrand);

} /* END: ReRandP */

void C_ReRandP(pnum_rand, pn, X, trt, eventStatus, peffect_p, pt1, pstopTol, pmaxiter,
               pprint, plambda_obs, ret_nrand, ret_p)
int *pnum_rand, *pn, *trt, *eventStatus, *pmaxiter, *pprint, *ret_nrand;
double *X, *peffect_p, *pt1, *pstopTol, *plambda_obs, *ret_p;
{

  /* For random number generation */
  GetRNGstate();

  *ret_nrand = ReRandP(*pnum_rand, *pn, X, trt, eventStatus, *peffect_p, *pt1, 
                       *pstopTol, *pmaxiter, *pprint, *plambda_obs, ret_p);

  PutRNGstate();  

  return;

} /* END: C_ReRandP */
