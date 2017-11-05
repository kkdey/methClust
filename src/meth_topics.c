#include <stdlib.h>
#include <assert.h>
#include <Rmath.h>
#include "latools.h"
#include "rhelp.h"
# ifdef _OPENMP
#include <omp.h>
# endif


double wllhd(int nwrd, double *m, double *u, double *q){
  double l = 0.0;
  int j;
  for(j=0; j<nwrd; j++) l += m[j]*log(q[j]) + u[j]*log(1 - q[j]);
  return l;
}

void wrdprob(double *q, int nwrd, int K, int *wrd, double **freq, double *w){
  int j, k;
  zero_dvec(q,nwrd);
  for(j=0; j<nwrd; j++)
    for(k=0; k<K; k++)
      q[j] += w[k]*freq[k][wrd[j]];
}

void wgrad(double *grad, int nwrd, int K, int *wrd, double *m,
           double *u, double *q, double **freq, double *w, int nef){
  int j, k;
  zero_dvec(grad,K);
  for(k=0; k<K; k++)
  { for(j=0; j<nwrd; j++) {
    grad[k] += (freq[k][wrd[j]]*m[j]/q[j] - freq[k][wrd[j]]*u[j]/(1 - q[j]));
  }
    if(nef) grad[k] += 1.0/(w[k]*((double) K)); }
}

void wneghess_lowtri(double *nH, int nwrd, int K, int *wrd, double *m,
                     double *u, double *q, double **freq, double *w, int nef){
  int j, h, k;
  zero_dvec(nH, K*K);
  for(k=0; k<K; k++){
    if(nef) nH[k*(K+1)] += 1.0/(w[k]*w[k]*((double) K));
    for(h=k; h<K; h++)
      for(j=0; j<nwrd; j++)
        nH[k*K+h] += m[j]*freq[k][wrd[j]]*freq[h][wrd[j]]/(q[j]*q[j]) - u[j]*freq[k][wrd[j]]*freq[h][wrd[j]]/(q[j]*q[j]);
  }
}

/* find a new active set  */
int wactivate(int K, double *w, int *active){
  int h, k, nact;
  nact = sum_ivec(active, K);

  for(h=0; h<K; h++)
    if(!active[h]){
      if(w[h] <= 0.0)
      { w[h]=0.0; active[h] = 1; nact++; }
      if(w[h] >= 1.0)
      { zero_dvec(w,K); w[h] = 1.0;
      for(k=0; k<K; k++) active[k] = 1;
      return 0; } }

    return K - nact;
}

/* update function for a given solution */
int wupdate(int K, double *w, double *B, int *active){
  double delmin = 1.0;
  double delta;
  int h, d;

  d = 0;
  for(h=0; h<K; h++)
    if(!active[h])
    { assert(w[h]>0);
      delta = 1.0;
      if(B[d] < -w[h]) delta = -w[h]/B[d];
      if(B[d] > 1.0-w[h]) delta = (1.0-w[h])/B[d];
      if(delta < delmin) delmin = delta;
      d++; }

    d = 0;
    for(h=0; h<K; h++)
      if(!active[h])
      { w[h] += delmin*B[d]; d++; }

      return  wactivate(K, w, active);
}

/* main program: sequential quadratic programming for topic weights */
// NB: assumed prior concentration of 1/K if in NEF parametrization
int sqpw(int p, int nwrd, int K,  int *wrd, double *m, double *u,
         double **freq, double *w, int nef, double tol, int tmax, int verb){

  int h,k,c,d;

  // single-word solution
  if(nwrd==1){
    zero_dvec(w,K);
    k = 0;
    for(h=1;h<K;h++)
      if(freq[h][*wrd] > freq[k][*wrd]) k = h;
      w[k] = 1.0;
      return 1; }

  double diff = tol+1.0;

  double *q = new_dvec(nwrd);

  double *A = new_dzero( (K+1)*(K+1) );
  double *B = new_dzero(K+1);
  double *nhess = new_dvec(K*K);
  double *grad = new_dvec(K);
  double *wold = new_dup_dvec(w, K);

  int *active = new_izero(K);
  int nw = K;
  int info;
  int dosysv = 1;
  int iter = 1;

  while(diff > tol && iter < tmax && nw > 0){

    /* update gradient, and negative hessian */
    wrdprob(q, nwrd, K, wrd, freq, w);
    wgrad(grad, nwrd, K, wrd, m, u, q, freq, w, nef);
    wneghess_lowtri(nhess, nwrd, K, wrd, m, u, q, freq, w, nef);

    /* build linear system for current active set */
    d = 0;
    for(k=0; k<K; k++)
      if(!active[k]){
        B[d] = grad[k];
        c = d;
        for(h=k; h<K; h++)
          if(!active[h]){ A[d*(nw+1) + c] = nhess[k*K + h]; c++; }
          d++;
      }
      assert(d==nw);
      for(c=0; c<nw; c++) A[c*(nw+1) + nw] = 1.0;
      A[(nw+1)*(nw+1)-1] = 0.0;
      B[nw] = 0.0;


      /* solve it  (sysv can be unstable for small samples; switch to gesv if necessary) */
      if(dosysv){ info = la_dsysv(nw+1, 1, A, B);
        if(info!=0){ dosysv = 0; continue; } }
      else{
        for(d=1; d<=nw; d++) for(c=0; c<d; c++) A[d*(nw+1) + c] = A[c*(nw+1) +d];
        info = la_dgesv(nw+1, 1, A, B);
      }

      /* check and update */
      if(fabs(sum_dvec(B, nw)) > .001 || (info !=0 && dosysv == 0)){
        myprintf(mystdout, "trouble in wsolve\n");
        nw = wactivate(K, w, active);
        return 0;
      }
      else nw = wupdate(K, w, B, active);


      /* iterate */
      diff = 0.0;
      for(h=0; h<K; h++) diff += fabs(wold[h]-w[h]);
      copy_dvec(wold, w, K);

      if(verb==1){
        myprintf(mystdout,
                 "Omega Fit: iter = %d, unique words = %d, active=%d, diff = %g \n",
                 iter, nwrd, K-nw,  diff);
        myprintf(mystdout, "          W = "); print_dvec(w, K, mystdout); }

      iter++;
  }

  free(wold);
  free(active);
  free(q);
  free(B);
  free(grad);
  free(A);
  free(nhess);

  return 1;
}

void Romega( int *n, int *p, int *K,
             int *doc, int *wrd, double *m, double *u,
             double *freq_vec, double *W, int *nef,
             double *tol, int *tmax, int *verb){

  int speakup = *verb;

  double **freq = new_mat_fromv(*p, *K, freq_vec);
  int i;

# ifdef _OPENMP
  speakup=0; // otherwise it'll crash R-GUI trying to print all at once
#pragma omp parallel for private(i)
# endif
  for(i=0; i<(*n); i++)
    if(sqpw(*p, doc[i+1]-doc[i], *K, &wrd[doc[i]], &m[doc[i]], &u[doc[i]],  freq,  &W[i*(*K)], *nef, *tol, *tmax, speakup) == 0)
      myprintf(mystdout, "Failed to converge for omega at i = %d\n", i+1);

    delete_mat(freq);
}

void RmixQ( int *n_in, int *p_in, int *K_in, int *N_in, int *B_in,
            double *m, double *u,  int *doc, int *wrd, int *grp,
            double *omega_vec, double *freq_vec, double *Q){

  int l, i, k, n, p, K, N, B;

  n = *n_in;
  p = *p_in;
  K = *K_in;
  N = *N_in;
  B = *B_in;

  for(l=0; l<N; l++)
    for(k=0; k<K; k++)
      Q[k*n + doc[l]] += m[l]*log(freq_vec[k*p + wrd[l]]) + u[l]*log(1 - freq_vec[k*p + wrd[l]]);

  for(i=0; i<n; i++)
    for(k=0; k<K; k++)
      Q[k*n + i] += log(omega_vec[k*B + grp[i]]);

}



void RcalcQ(int *n_in, int *p_in, int *K_in, int *doc, int *wrd,
            int *N, double *omega, double *freq, double *q ){

  int n, p, K, l, h;

  n = *n_in;
  p = *p_in;
  K = *K_in;

  for(l=0; l<(*N); l++){
    q[l] = 0.0;
    for(h=0; h<K; h++) q[l] += omega[h*n + doc[l]]*freq[h*p + wrd[l]];
  }
}


