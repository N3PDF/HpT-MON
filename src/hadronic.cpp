/*
 * =====================================================================================
 *
 *       Filename:  hadronic.cpp
 *
 *    Description:  Integration over the momentum fractions in order to get the
 * final hadronic results.
 *
 *        Version:  1.0
 *        Created:  24/02/2021 09:31:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#include "../include/hadronic.h"

HadronicHiggs::HadronicHiggs(int order, int channel, std::string pdfname,
                             void *params)
    : higgs(order, channel, pdfname, params) {
  PhysParams param = *reinterpret_cast<PhysParams *>(params);

  NF = param.nf;
  ORD = order;
  SIGMA0 = param.sigma0;
  aass = param.alphas / M_PI;
  SROOT = param.sroot;
  MH2 = std::pow(param.mh, 2);

  q2 = MH2;
  shad = std::pow(SROOT, 2);
  xtau = q2 / shad;
}

HadronicHiggs::~HadronicHiggs() {}

struct HadronicStruct {
  HadronicHiggs *hadrhiggs;
  double pt;
  double tm;
  double y1;
  double y2;
};

//==============================================================================================//
//                  Routines that perform the integration over rapidity. //
//----------------------------------------------------------------------------------------------//
double HadronicHiggs::ExtractDeltaHadronic(double pt, double tm, double t,
                                           double y1, double y2, double yc) {
  double xsection = higgs.delta(pt, tm, t, y1, y2, yc) * SIGMA0 * (aass / 2.);
  return xsection;
}

double deltaIntegrand(double t, double yc, void *p) {
  HadronicStruct par = *reinterpret_cast<HadronicStruct *>(p);
  return par.hadrhiggs->ExtractDeltaHadronic(par.pt, par.tm, t, par.y1, par.y2,
                                             yc);
}
//==============================================================================================//

//==============================================================================================//
//                  Routines that perform the full hadronic integration. //
//----------------------------------------------------------------------------------------------//
double HadronicHiggs::ExtractDistrHadronic(double pt, double tm, double y1,
                                           double y2, double zz1, double zz2,
                                           double zz3) {
  double crosres;
  double distres = higgs.distr(pt, tm, y1, y2, zz1, zz2, zz3);

  if (y2 == -y1) {
    crosres = distres;
  } else {
    crosres = higgs.distrcross(pt, tm, y1, y2, zz1, zz2, zz3);
  }

  double xsection = (distres + crosres) * SIGMA0 * (aass / 2.);

  return xsection;
}

double distrIntegrand(double zz1, double zz2, double zz3, void *p) {
  HadronicStruct par = *reinterpret_cast<HadronicStruct *>(p);
  return par.hadrhiggs->ExtractDistrHadronic(par.pt, par.tm, par.y1, par.y2,
                                             zz1, zz2, zz3);
}
//==============================================================================================//

std::vector<double> HadronicHiggs::higgsdpt(double pt, double y1, double y2) {
  //////////////////////////////////////////////////
  // Descpription:                                //
  // pt:  value of the pt                         //
  // y1:  lower boundary of the rapidity yh       //
  // y2:  upper boundary of the rapidity yh       //
  //////////////////////////////////////////////////

  int method = 0;
  double deltres = 0, distres = 0, delterr = 0, disterr = 0;
  double yy1, yy2;
  std::vector<double> results;

  double q = sqrt(q2);
  double tm = sqrt(q2 + std::pow(pt, 2));  // Transverse mass
  double tmpx = (q2 + shad) / SROOT / tm;
  double xr = std::pow(1 - xtau, 2) -
              4 * xtau * std::pow(pt / q, 2);  // pt kinematic limit
  double ymax =
      log((tmpx + sqrt(std::pow(tmpx, 2) - 4)) / 2.);  // y kinematic limit

  if (xr < 0) {
    return {0., 0.};
  }

  // Define rapidity boundaries
  if (y1 > ymax || y2 < -ymax) {
    return {0., 0.};
  } else if (y1 < -ymax && y2 > ymax) {
    yy1 = -ymax;
    yy2 = ymax;
  } else if (y1 < -ymax && y2 > -ymax && y2 < ymax) {
    yy1 = -ymax;
    yy2 = y2;
  } else if (y1 < -ymax && y2 > -ymax && y2 < ymax) {
    yy1 = -ymax;
    yy2 = y2;
  } else if (y1 > -ymax && y1 < ymax && y2 > ymax) {
    yy1 = y1;
    yy2 = ymax;
  } else {
    yy1 = y1;
    yy2 = y2;
  }

  // Pass parameters
  HadronicStruct finalparams;
  finalparams.pt = pt;
  finalparams.tm = tm;
  finalparams.y1 = yy1;
  finalparams.y2 = yy2;
  finalparams.hadrhiggs = this;

  // LO
  deltres = integration::IntegrateDeltaHadronic(method, deltaIntegrand,
                                                &finalparams, &delterr);

  // NLO
  if (ORD >= 1) {
    distres = integration::IntegrateDistrHadronic(method, distrIntegrand,
                                                  &finalparams, &disterr);

    distres *= aass / 2.;
    disterr *= aass / 2.;
  }

  results.push_back(deltres + distres);
  results.push_back(delterr + disterr);

  return results;
}
