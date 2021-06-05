/*
 * =====================================================================================
 *
 *       Filename:  hadronic.h
 *
 *    Description:  Header file for the hadronic.cpp file which computes the
 * full hadronic cross section in x space.
 *
 *        Version:  1.0
 *        Created:  24/02/2021 08:46:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <vector>

#include "./higgspt.h"
#include "./integration.h"

class HadronicHiggs {
 public:
  HadronicHiggs(int order, int channel, std::string pdfname, void *params);
  virtual ~HadronicHiggs();

  // Main function that performs the
  // integration over rapidities
  std::vector<double> higgsdpt(double pt, double y1, double y2);

  double ExtractDeltaHadronic(double pt, double tm, double t, double y1,
                              double y2, double yc);

  double ExtractDistrHadronic(double pt, double tm, double y1, double y2,
                              double zz1, double zz2, double zz3);

 private:
  double NF, SROOT;
  double MH2, q2, shad;
  double aass, SIGMA0, xtau;
  double ORD, CHANNEL;

  // Instantiate HiggsDpT class
  HiggsDpT higgs;
};
