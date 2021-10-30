#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#ifdef HAVE_ALLOCA_H
#  include <alloca.h>
#elif !defined alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# elif defined(__DECC)
#  define alloca __ALLOCA
# elif defined(_MSC_VER)
#  include <malloc.h>
#  define alloca _alloca
# elif defined(__sun)
#  include <alloca.h>
# else
#  include <stdlib.h>
# endif
#endif

#include <Rversion.h>

#if R_VERSION >= R_Version(3, 6, 2)
#  define USE_FC_LEN_T
#endif

#if R_VERSION <= R_Version(3, 3, 1)
#  define NO_C_HEADERS
#endif

#define R_NO_REMAP 1
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

#undef NO_C_HEADERS
#undef R_NO_REMAP
#undef USE_FC_LEN_T

#ifndef FCONE
# define FCONE
#endif

static SEXP squaredExponential_updateCovMatrix(SEXP covExpr, SEXP xtExpr, SEXP x_tExpr, SEXP parsExpr, SEXP sig_f_sqExpr)
{
  int* dims = INTEGER(Rf_getAttrib(xtExpr, R_DimSymbol));
  size_t p = (size_t) dims[0];
  size_t n = (size_t) dims[1];
  size_t m = (size_t) INTEGER(Rf_getAttrib(x_tExpr, R_DimSymbol))[1];
  
  
  double* cov  = REAL(covExpr);
  double* xt   = REAL(xtExpr);
  double* x_t  = REAL(x_tExpr);
  double* pars = REAL(parsExpr);
  
  double* scales = (double*) alloca(p * sizeof(double));
  for (size_t i = 0; i < p; ++i) scales[i] = exp(-0.5 * pars[i]);
  
  double sig_f_sq = REAL(sig_f_sqExpr)[0];
  
  if (xt == x_t) {
    for (size_t row = 0; row < (n - 1); ++row) {
      double* x_Row = x_t + row * p;
      for (size_t col = (row + 1); col < n; ++col) {
        double* xRow = xt + col * p;
        
        double temp1 = 0.0;
        for (size_t i = 0; i < p; ++i) {
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
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
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
          temp1 += temp2 * temp2;
        }
        cov[row + col * n] = sig_f_sq * exp(-temp1);
      }
    }
  }
  
  return R_NilValue;
}

static SEXP matern_updateCovMatrix(SEXP covExpr, SEXP xtExpr, SEXP x_tExpr, SEXP parsExpr, SEXP sig_f_sqExpr)
{
  int* dims = INTEGER(Rf_getAttrib(xtExpr, R_DimSymbol));
  size_t p = (size_t) dims[0];
  size_t n = (size_t) dims[1];
  size_t m = (size_t) INTEGER(Rf_getAttrib(x_tExpr, R_DimSymbol))[1];
  
  
  double* cov  = REAL(covExpr);
  double* xt   = REAL(xtExpr);
  double* x_t  = REAL(x_tExpr);
  double* pars = REAL(parsExpr);
  
  double nu = exp(pars[0]);
  double* scales = (double*) alloca(p * sizeof(double));
  for (size_t i = 0; i < p; ++i) scales[i] = exp(-0.5 * pars[i + 1]);

  double sig_f_sq = REAL(sig_f_sqExpr)[0];
  
  double constTerm = (1.0 - nu) * M_LN2 - Rf_lgammafn(nu);
  if (xt == x_t) {
    for (size_t row = 0; row < (n - 1); ++row) {
      double* x_Row = x_t + row * p;
      for (size_t col = (row + 1); col < n; ++col) {
        double* xRow = xt + col * p;
        
        double temp1 = 0.0;
        for (size_t i = 0; i < p; ++i) {
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
          temp1 += temp2 * temp2;
        }
        
        temp1 = sqrt(2.0 * nu * temp1);
        temp1 = log(Rf_bessel_k(temp1, nu, 1.0)) - temp1 + nu * log(temp1);
        
        temp1 = sig_f_sq * exp(constTerm + temp1);
        
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
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
          temp1 += temp2 * temp2;
        }
        temp1 = sqrt(2.0 * nu * temp1);
        temp1 = log(Rf_bessel_k(temp1, nu, 1.0)) - temp1 + nu * log(temp1);
        
        cov[row + col * n] = sig_f_sq * exp(constTerm + temp1);
      }
    }
  }
  
  return R_NilValue;
}

static SEXP exponential_updateCovMatrix(SEXP covExpr, SEXP xtExpr, SEXP x_tExpr, SEXP parsExpr, SEXP sig_f_sqExpr)
{
  int* dims = INTEGER(Rf_getAttrib(xtExpr, R_DimSymbol));
  size_t p = (size_t) dims[0];
  size_t n = (size_t) dims[1];
  size_t m = (size_t) INTEGER(Rf_getAttrib(x_tExpr, R_DimSymbol))[1];
  
  
  double* cov  = REAL(covExpr);
  double* xt   = REAL(xtExpr);
  double* x_t  = REAL(x_tExpr);
  double* pars = REAL(parsExpr);
  
  double* scales = (double*) alloca(p * sizeof(double));
  for (size_t i = 0; i < p; ++i) scales[i] = exp(-0.5 * pars[i]);
  
  double sig_f_sq = REAL(sig_f_sqExpr)[0];
  
  if (xt == x_t) {
    for (size_t row = 0; row < (n - 1); ++row) {
      double* x_Row = x_t + row * p;
      for (size_t col = (row + 1); col < n; ++col) {
        double* xRow = xt + col * p;
        
        double temp1 = 0.0;
        for (size_t i = 0; i < p; ++i) {
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
          temp1 += temp2 * temp2;
        }
        temp1 = sig_f_sq * exp(-sqrt(temp1));
        
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
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
          temp1 += temp2 * temp2;
        }
        cov[row + col * n] = sig_f_sq * exp(-sqrt(temp1));
      }
    }
  }
  
  return R_NilValue;
}

static SEXP gammaExponential_updateCovMatrix(SEXP covExpr, SEXP xtExpr, SEXP x_tExpr, SEXP parsExpr, SEXP sig_f_sqExpr)
{
  int* dims = INTEGER(Rf_getAttrib(xtExpr, R_DimSymbol));
  size_t p = (size_t) dims[0];
  size_t n = (size_t) dims[1];
  size_t m = (size_t) INTEGER(Rf_getAttrib(x_tExpr, R_DimSymbol))[1];
  
  
  double* cov  = REAL(covExpr);
  double* xt   = REAL(xtExpr);
  double* x_t  = REAL(x_tExpr);
  double* pars = REAL(parsExpr);
  
  double gamma = exp(pars[0]);
  gamma = 2.0 * gamma / (1.0 + gamma);
  double* scales = (double*) alloca(p * sizeof(double));
  for (size_t i = 0; i < p; ++i) scales[i] = exp(-0.5 * pars[i + 1]);
  
  double sig_f_sq = REAL(sig_f_sqExpr)[0];
  
  if (xt == x_t) {
    for (size_t row = 0; row < (n - 1); ++row) {
      double* x_Row = x_t + row * p;
      for (size_t col = (row + 1); col < n; ++col) {
        double* xRow = xt + col * p;
        
        double temp1 = 0.0;
        for (size_t i = 0; i < p; ++i) {
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
          temp1 += temp2 * temp2;
        }
        temp1 = sig_f_sq * exp(-pow(sqrt(temp1), gamma));
        
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
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
          temp1 += temp2 * temp2;
        }
        cov[row + col * n] = sig_f_sq * exp(-pow(sqrt(temp1), gamma));
      }
    }
  }
  
  return R_NilValue;
}

static SEXP rationalQuadratic_updateCovMatrix(SEXP covExpr, SEXP xtExpr, SEXP x_tExpr, SEXP parsExpr, SEXP sig_f_sqExpr)
{
  int* dims = INTEGER(Rf_getAttrib(xtExpr, R_DimSymbol));
  size_t p = (size_t) dims[0];
  size_t n = (size_t) dims[1];
  size_t m = (size_t) INTEGER(Rf_getAttrib(x_tExpr, R_DimSymbol))[1];
  
  
  double* cov  = REAL(covExpr);
  double* xt   = REAL(xtExpr);
  double* x_t  = REAL(x_tExpr);
  double* pars = REAL(parsExpr);
  
  double alpha = exp(pars[0]);
  double* scales = (double*) alloca(p * sizeof(double));
  for (size_t i = 0; i < p; ++i) scales[i] = exp(-0.5 * pars[i + 1]);
  
  double sig_f_sq = REAL(sig_f_sqExpr)[0];
  
  if (xt == x_t) {
    for (size_t row = 0; row < (n - 1); ++row) {
      double* x_Row = x_t + row * p;
      for (size_t col = (row + 1); col < n; ++col) {
        double* xRow = xt + col * p;
        
        double temp1 = 0.0;
        for (size_t i = 0; i < p; ++i) {
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
          temp1 += temp2 * temp2;
        }
        temp1 = sig_f_sq * pow(1.0 + temp1 / (2.0 * alpha), -alpha);
        
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
          double temp2 = (xRow[i] - x_Row[i]) / scales[i];
          temp1 += temp2 * temp2;
        }
        cov[row + col * n] = sig_f_sq * pow(1.0 + temp1 / (2.0 * alpha), -alpha);
      }
    }
  }
  
  return R_NilValue;
}

static SEXP neuralNetwork_updateCovMatrix(SEXP covExpr, SEXP xtExpr, SEXP x_tExpr, SEXP parsExpr, SEXP sig_f_sqExpr)
{
  int* dims = INTEGER(Rf_getAttrib(xtExpr, R_DimSymbol));
  size_t p = (size_t) dims[0];
  size_t n = (size_t) dims[1];
  size_t m = (size_t) INTEGER(Rf_getAttrib(x_tExpr, R_DimSymbol))[1];
  
  
  double* cov  = REAL(covExpr);
  double* xt   = REAL(xtExpr);
  double* x_t  = REAL(x_tExpr);
  double* pars = REAL(parsExpr);
  
  double offset = exp(pars[0]);
  double* scales = (double*) alloca(p * sizeof(double));
  for (size_t i = 0; i < p; ++i) scales[i] = exp(-0.5 * pars[i + 1]);
  
  double sig_f_sq = REAL(sig_f_sqExpr)[0];
  
  if (xt == x_t) {
    for (size_t row = 0; row < (n - 1); ++row) {
      double* x_Row = x_t + row * p;
      for (size_t col = (row + 1); col < n; ++col) {
        double* xRow = xt + col * p;
        
        double crossSum = offset;
        double xSum = offset;
        double x_Sum = offset;
        for (size_t i = 0; i < p; ++i) {
          crossSum += ( xRow[i] * x_Row[i]) / scales[i];
          xSum     += ( xRow[i] *  xRow[i]) / scales[i];
          x_Sum    += (x_Row[i] * x_Row[i]) / scales[i];
        }
        
        crossSum = sig_f_sq * M_2_PI * asin(crossSum / sqrt((1.0 + xSum) * (1.0 + x_Sum)));
        cov[row + col * n] = crossSum;
        cov[col + row * n] = crossSum;
      }
      cov[row + row * n] = sig_f_sq;
    }
    cov[n * n - 1] = sig_f_sq;
  } else {
    for (size_t row = 0; row < n; ++row) {
      double* xRow = xt + row * p;
      for (size_t col = 0; col < m; ++col) {
        double* x_Row = x_t + col * p;
        
        double crossSum = offset;
        double xSum = offset;
        double x_Sum = offset;
        for (size_t i = 0; i < p; ++i) {
          crossSum += ( xRow[i] * x_Row[i]) / scales[i];
          xSum     += ( xRow[i] *  xRow[i]) / scales[i];
          x_Sum    += (x_Row[i] * x_Row[i]) / scales[i];
        }
        
        crossSum = sig_f_sq * M_2_PI * asin(crossSum / sqrt((1.0 + xSum) * (1.0 + x_Sum)));
        cov[row + col * n] = crossSum;
      }
    }
  }
  
  return R_NilValue;
}


static SEXP updateLeftFactor(SEXP LExpr, SEXP covExpr)
{
  size_t dim = (size_t) INTEGER(Rf_getAttrib(LExpr, R_DimSymbol))[1];
  
  double* L         = REAL(LExpr);
  const double* cov = REAL(covExpr);
  
  memcpy(L, cov, dim * dim * sizeof(double));
  
  for (size_t i = 0; i < dim; ++i) L[i + i * dim] += 1.0;
  
  char triangleType = 'L';
  
  int i_dim = (int) dim;
  
  int lapackResult;
  
  F77_CALL(dpotrf)(&triangleType, &i_dim, L, &i_dim, &lapackResult FCONE);
  
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
  defineFunction("npci_squaredExponential_updateCovMatrix", squaredExponential_updateCovMatrix, 5),
  defineFunction("npci_matern_updateCovMatrix", matern_updateCovMatrix, 5),
  defineFunction("npci_exponential_updateCovMatrix", exponential_updateCovMatrix, 5),
  defineFunction("npci_gammaExponential_updateCovMatrix", gammaExponential_updateCovMatrix, 5),
  defineFunction("npci_rationalQuadratic_updateCovMatrix", rationalQuadratic_updateCovMatrix, 5),
  defineFunction("npci_neuralNetwork_updateCovMatrix", neuralNetwork_updateCovMatrix, 5),
  defineFunction("npci_updateLeftFactor", updateLeftFactor, 2),
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

