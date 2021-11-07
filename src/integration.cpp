/*
 * =====================================================================================
 *
 *       Filename:  integration.cpp
 *
 *    Description:  Integration routine based on Cuba.
 *
 *        Version:  1.0
 *        Created:  17/02/2021 23:40:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rademananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#include "../include/integration.h"

struct IntData {
  std::function<double(double, void *)> fun;
  void *param;
};

struct IntData1 {
  std::function<double(double, double, void *)> fun;
  void *param;
};

struct IntData2 {
  std::function<double(double, double, double, void *)> fun;
  void *param;
};

int integration::CubaIntegrator(
    int method, int(Func)(int *, double *, int *, double *, void *), int ndim,
    int ncomp, double prec, double **res, int *fail, double **error,
    double **prob, void *par, int verbos) {
  static const int last = 4;
  static const int seed = 0;
  static const int gridno = 0;
  // static const int nmin = 2;
  static const int key1 = 13;
  static const int key2 = 13;
  static const int key3 = 1;
  static const int ngiven = 0;
  static const int nextra = 0;
  static const int key = 13;
  static const int nvec = 1;
  // static const int spin = -1;
  static const int maxpass = 2;
  static const int mineval = 0;
  static const int nnew = 1e5;
  static const int nmin = 1e5;
  static const int nstart = 1e5;
  static const int nbatch = 1e4;
  static const int nincrease = 5e4;
  static const int ldxgiven = ndim;
  static const int maxeval = 1e9;
  static const int verbose = verbos;
  static const double border = 0.1;
  static const double maxchisq = 10.;
  static const double epsabs = 1e-100;
  static const double flatness = 1;
  static const double mindeviation = 0.25;
  const char *statefile = NULL;

  int neval, nregions;
  double *ris = NULL;
  double *err = NULL;
  double *proba = NULL;
  ris = new double[ncomp];
  err = new double[ncomp];
  proba = new double[ncomp];
  integrand_t func = (integrand_t)Func;

  switch (method) {
    case 0:
      Vegas(ndim, ncomp, func, par, nvec, prec, epsabs, verbose, seed, mineval,
            maxeval, nstart, nincrease, nbatch, gridno, statefile, NULL, &neval,
            fail, ris, err, proba);
      break;
    case 1:
      Suave(ndim, ncomp, func, par, nvec, prec, epsabs, verbose, seed, mineval,
            maxeval, nnew, nmin, flatness, statefile, NULL, &nregions, &neval,
            fail, ris, err, proba);
      break;
    case 2:
      Divonne(ndim, ncomp, func, par, nvec, prec, epsabs, verbose, seed,
              mineval, maxeval, key1, key2, key3, maxpass, border, maxchisq,
              mindeviation, ngiven, ldxgiven, NULL, nextra, NULL, statefile,
              NULL, &nregions, &neval, fail, ris, err, proba);
      break;
    case 3:
      Cuhre(ndim, ncomp, func, par, nvec, prec, epsabs, verbose | last, mineval,
            maxeval, key, statefile, NULL, &nregions, &neval, fail, ris, err,
            proba);
      break;
  }
  *error = err;
  *prob = proba;
  *res = ris;

  return 0;
}

//==============================================================================================//
//                  Routines that perform the integration for partoic xsec. //
//----------------------------------------------------------------------------------------------//
int DeltaPartonicIntegrand(int *ndim, double *x, int *ncomp, double *y,
                           void *p) {
  IntData inparam = *reinterpret_cast<IntData *>(p);
  double zz = x[0];
  y[0] = inparam.fun(zz, inparam.param);
  return 0;
}

double integration::IntegrateDeltaPartonic(int method,
                                           double(Func)(double, void *),
                                           void *pp, double *error) {
  int fail;
  double prec = 1e-3;
  double *res = NULL;
  double *err = NULL;
  double *prb = NULL;

  res = new double[1];
  err = new double[1];
  prb = new double[1];

  int ndim = 1;
  int ncomp = 1;
  int verbose = 0;

  IntData UserData;
  UserData.fun = Func;
  UserData.param = pp;

  CubaIntegrator(method, DeltaPartonicIntegrand, ndim, ncomp, prec, &res, &fail,
                 &err, &prb, &UserData, verbose);

  *error = err[0];
  return res[0];
}

int DistrPartonicIntegrand(int *ndim, double *x, int *ncomp, double *y,
                           void *p) {
  IntData1 inparam = *reinterpret_cast<IntData1 *>(p);
  double zz1 = x[0];
  double zz2 = x[1];
  y[0] = inparam.fun(zz1, zz2, inparam.param);
  return 0;
}

double integration::IntegrateDistrPartonic(int method,
                                           double(Func)(double, double, void *),
                                           void *pp, double *error) {
  int fail;
  double prec = 1e-3;
  double *res = NULL;
  double *err = NULL;
  double *prb = NULL;

  res = new double[1];
  err = new double[1];
  prb = new double[1];

  int ndim = 2;
  int ncomp = 1;
  int verbose = 0;

  IntData1 UserData;
  UserData.fun = Func;
  UserData.param = pp;

  CubaIntegrator(method, DistrPartonicIntegrand, ndim, ncomp, prec, &res, &fail,
                 &err, &prb, &UserData, verbose);

  *error = err[0];
  return res[0];
}

//==============================================================================================//

//==============================================================================================//
//          Routines that perform the full hadronic integration for Delta terms.
//          //
//----------------------------------------------------------------------------------------------//
int DeltaIntegrand(int *ndim, double *x, int *ncomp, double *y, void *p) {
  IntData1 inparam = *reinterpret_cast<IntData1 *>(p);
  double t = x[0];
  double yc = x[1];
  y[0] = inparam.fun(t, yc, inparam.param);
  return 0;
}

double integration::IntegrateDeltaHadronic(int method,
                                           double(Func)(double, double, void *),
                                           void *pp, double *error) {
  int fail;
  double prec = 1e-3;
  double *res = NULL;
  double *err = NULL;
  double *prb = NULL;

  res = new double[1];
  err = new double[1];
  prb = new double[1];

  int ndim = 2;
  int ncomp = 1;
  int verbose = 0;

  IntData1 UserData;
  UserData.fun = Func;
  UserData.param = pp;

  CubaIntegrator(method, DeltaIntegrand, ndim, ncomp, prec, &res, &fail, &err,
                 &prb, &UserData, verbose);

  *error = err[0];
  return res[0];
}
//==============================================================================================//

//==============================================================================================//
//            Routines that perform the full hadronic integration for regular
//            terms.            //
//----------------------------------------------------------------------------------------------//
int DistrIntegrand(int *ndim, double *x, int *ncomp, double *y, void *p) {
  IntData2 inparam = *reinterpret_cast<IntData2 *>(p);
  double zz1 = x[0];
  double zz2 = x[1];
  double zz3 = x[2];
  y[0] = inparam.fun(zz1, zz2, zz3, inparam.param);
  return 0;
}

double integration::IntegrateDistrHadronic(int method,
                                           double(Func)(double, double, double,
                                                        void *),
                                           void *pp, double *error) {
  int fail;
  double prec = 1e-3;
  double *res = NULL;
  double *err = NULL;
  double *prb = NULL;

  res = new double[1];
  err = new double[1];
  prb = new double[1];

  int ndim = 3;
  int ncomp = 1;
  int verbose = 0;

  IntData2 UserData;
  UserData.fun = Func;
  UserData.param = pp;

  CubaIntegrator(method, DistrIntegrand, ndim, ncomp, prec, &res, &fail, &err,
                 &prb, &UserData, verbose);

  *error = err[0];
  return res[0];
}
//==============================================================================================//
