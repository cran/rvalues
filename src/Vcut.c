#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP Vcut(SEXP Vmat, SEXP lamfun, SEXP nunits, SEXP ngrid) {
  /* Vmat - nunits x ngrid matrix */
  /* lam_fun - vector of length ngrid */
  SEXP ans;
  int nrow, ncol, ii, jj;
  double tst;
  
  nrow = INTEGER(nunits)[0];
  ncol = INTEGER(ngrid)[0];
  
  PROTECT(ans = allocVector(INTSXP, nrow));
  
  for (ii=0; ii < nrow; ii++) {
     for(jj=0; jj < ncol; jj++)  {
         tst = REAL(Vmat)[ii + nrow*jj] - REAL(lamfun)[jj];
         if(tst > 0) {
            INTEGER(ans)[ii] = jj + 1;
            break;
         }
         if(jj == ncol - 1) {
            INTEGER(ans)[ii] = jj + 1;
            // If no root exists across columns return max 
         }
     }
  }
  
  UNPROTECT(1); 
  return(ans);
} 
