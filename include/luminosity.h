/*
 * =====================================================================================
 *
 *       Filename:  luminosity.h
 *
 *    Description:  Header file for the luminosity.cpp file.
 *
 *        Version:  1.0
 *        Created:  25/02/2021 20:48:01
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#pragma once

#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_zeta.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class Luminosity {
 public:
  typedef LHAPDF::PDF* pdfptr;

  Luminosity(std::string pdfname, double MUF, int NF);
  virtual ~Luminosity();

  // flavour based luminosity
  double Lumgg(double x1, double x2);
  double Lumgq(double x1, double x2);
  double Lumqg(double x1, double x2);
  double Lumqq(double x1, double x2);
  double LumqQ(double x1, double x2);
  double Lumqqb(double x1, double x2);
  double LumqQb(double x1, double x2);

  // compute pdf
  double xfxQ(int pidflv, double x);
  double fxQ(double x, double pidflv);

  // Extract alphas
  double get_alphas(double mur);

 private:
  LHAPDF::PDF* _pdf;
  double _nf, _muf;
};
