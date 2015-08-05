#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#define R_NO_REMAP 1
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

static SEXP grouped_updateCovMatrix(SEXP covExpr, SEXP xtExpr, SEXP x_tExpr, SEXP parsExpr, SEXP sig_f_sqExpr)
{
  int* dims = INTEGER(Rf_getAttrib(xtExpr, R_DimSymbol));
  size_t p = (size_t) dims[0];
  size_t n = (size_t) dims[1];
  size_t m = (size_t) INTEGER(Rf_getAttrib(x_tExpr, R_DimSymbol))[1];
  
  
  double* cov  = REAL(covExpr);
  double* xt   = REAL(xtExpr);
  double* x_t  = REAL(x_tExpr);
  double* pars = REAL(parsExpr);
  
  double sig_f_sq = REAL(sig_f_sqExpr)[0];
  
  if (xt == x_t) {
    for (size_t row = 0; row < (n - 1); ++row) {
      double* x_Row = x_t + row * p;
      for (size_t col = (row + 1); col < n; ++col) {
        double* xRow = xt + col * p;
        
        double temp1 = 0.0;
        for (size_t i = 0; i < p; ++i) {
          double temp2 = (xRow[i] - x_Row[i]) / pars[i];
          temp1 += temp2 * temp2;
        }
        temp1 = sig_f_sq * exp(-temp1);
        
        cov[row + col * n] = temp1;
        cov[col + row * n] = temp1;
      }
      cov[row + row * n] = sig_f_sq;
    }
    cov[n * n - 1] = sig_f_sq;
  } else {
    for (size_t row = 0; row < n; ++row) {
      double* xRow = xt + row * p;
      for (size_t col = 0; col < m; ++col) {
        double* x_Row = x_t + col * p;
        
        double temp1 = 0.0;
        for (size_t i = 0; i < p; ++i) {
          double temp2 = (xRow[i] - x_Row[i]) / pars[i];
          temp1 += temp2 * temp2;
        }
        cov[row + col * n] = sig_f_sq * exp(-temp1);
      }
    }
  }
  
  return R_NilValue;
}

static SEXP grouped_updateLeftFactor(SEXP LExpr, SEXP covExpr)
{
  size_t dim = (size_t) INTEGER(Rf_getAttrib(LExpr, R_DimSymbol))[1];
  
  double* L         = REAL(LExpr);
  const double* cov = REAL(covExpr);
  
  memcpy(L, cov, dim * dim * sizeof(double));
  
  for (size_t i = 0; i < dim; ++i) L[i + i * dim] += 1.0;
  
  char triangleType = 'L';
  
  int i_dim = (int) dim;
  
  int lapackResult;
  
  F77_CALL(dpotrf)(&triangleType, &i_dim, L, &i_dim, &lapackResult);
  
  if (lapackResult < 0) return Rf_ScalarInteger(EINVAL);
  if (lapackResult > 0) return Rf_ScalarInteger(EDOM);
  
  for (size_t row = 0; row < dim - 1; ++row) {
    for (size_t col = row + 1; col < dim; ++col) {
      L[row + col * dim] = 0.0;
    }
  }
  
  return Rf_ScalarInteger(0.0);
}

#define defineFunction(_N_, _F_, _A_) { _N_, (DL_FUNC) (&_F_), _A_ }

static R_CallMethodDef R_callMethods[] = {
  defineFunction("npci_grouped_updateCovMatrix", grouped_updateCovMatrix, 5),
  defineFunction("npci_grouped_updateLeftFactor", grouped_updateLeftFactor, 2),
  { NULL, NULL, 0 }
};

#undef defineFunction

typedef struct {
  const char* name;
  DL_FUNC function;
} C_CallMethodDef;

void R_init_npci(DllInfo* info)
{
  R_registerRoutines(info, NULL, R_callMethods, NULL, NULL);
  R_useDynamicSymbols(info, (Rboolean) false);
}

