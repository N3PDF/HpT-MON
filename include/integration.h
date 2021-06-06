/*
 * =====================================================================================
 *
 *       Filename:  integration.h
 *
 *    Description:  Header file for the integration.cpp file which implements
 * the integration routines using cuba.
 *
 *        Version:  1.0
 *        Created:  17/02/2021 23:40:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#pragma once

#include <cuba.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <vector>

namespace integration {
extern int CubaIntegrator(int method,
                          int(Func)(int *, double *, int *, double *, void *),
                          int ndim, int ncomp, double prec, double **res,
                          int *fail, double **error, double **prob, void *par,
                          int verbos = 0);

double IntegrateDeltaPartonic(int method, double(Func)(double, void *),
                              void *pp, double *error);

double IntegrateDistrPartonic(int method, double(Func)(double, double, void *),
                              void *pp, double *error);

double IntegrateDeltaHadronic(int method, double(Func)(double, double, void *),
                              void *pp, double *error);

double IntegrateDistrHadronic(int method,
                              double(Func)(double, double, double, void *),
                              void *pp, double *error);
}  // namespace integration