/*
 * =====================================================================================
 *
 *       Filename:  partonic.cpp
 *
 *    Description:  This file provide the partonic results in N space for a
 * given channel.
 *
 *        Version:  1.0
 *        Created:  20/02/2021 22:21:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#include "../include/partonic.h"

CrossHiggs::CrossHiggs(int order, int channel, std::string pdfname,
                       void *params)
    : higgs(order, channel, pdfname, params) {
  PhysParams param = *reinterpret_cast<PhysParams *>(params);

  NF = param.nf;
  ORD = order;
  SIGMA0 = param.sigma0;
  aass = param.alphas;
  MH2 = std::pow(param.mh, 2);
}

CrossHiggs::~CrossHiggs() {}

struct PartStruct {
  CrossHiggs *croshiggs;
  double pt;
  double nn;
  double zz;
};

//==============================================================================================//
//                  Routines that perform the integration over rapidity. //
//----------------------------------------------------------------------------------------------//

double CrossHiggs::ExtractDeltaPartonic(double pt, double nn, double zz) {
  double xsection =
      higgs.deltapartonic(pt, nn, zz) * SIGMA0 * (aass / (2. * M_PI));
  return xsection;
}

double PartonicDeltaIntegrand(double zz, void *p) {
  PartStruct par = *reinterpret_cast<PartStruct *>(p);
  return par.croshiggs->ExtractDeltaPartonic(par.pt, par.nn, zz);
}

//==============================================================================================//
//                  Routines that perform the full partonic integration. //
//----------------------------------------------------------------------------------------------//
double CrossHiggs::ExtractDistrPartonic(double pt, double nn, double zz1,
                                        double zz2) {
  // A factor two is here to account for (th,a)<->(uh,b)
  // TODO: figure out why we also need a factor 2 for the  jacobian, to get good
  // agreement with
  //       threshold expanded and full LO.
  // TODO: also add a `cross' distribution to accomodate for the qg-channels
  double distres = higgs.distrpartonic(pt, nn, zz1, zz2);
  distres += higgs.distrpartoniccross(pt, nn, zz1, zz2);
  double xsection = SIGMA0 * (aass / (2. * M_PI)) * distres;
  return xsection;
}

double PartonicDistrIntegrand(double zz1, double zz2, void *p) {
  PartStruct par = *reinterpret_cast<PartStruct *>(p);
  return par.croshiggs->ExtractDistrPartonic(par.pt, par.nn, zz1, zz2);
}

//==============================================================================================//
//           Main function that computes the partonic part for a given channel
//           //
//----------------------------------------------------------------------------------------------//
std::vector<double> CrossHiggs::partonichiggsdpt(double pt, double nn) {
  //////////////////////////////////////////////////
  // Descpription:                                //
  // pt:  value of the pt                         //
  // y1:  lower boundary of the rapidity yh       //
  // y2:  upper boundary of the rapidity yh       //
  //////////////////////////////////////////////////

  int method = 1;
  double deltres, delterr;
  double distres = 0, disterr = 0;
  std::vector<double> results;

  // Pass parameters
  PartStruct finalparams;
  finalparams.pt = pt;
  finalparams.nn = nn;
  finalparams.croshiggs = this;

  // LO
  deltres = integration::IntegrateDeltaPartonic(method, PartonicDeltaIntegrand,
                                                &finalparams, &delterr);

  // NLO
  if (ORD >= 1) {
    distres = integration::IntegrateDistrPartonic(
        method, PartonicDistrIntegrand, &finalparams, &disterr);

    distres *= aass / (2. * M_PI);
    disterr *= aass / (2. * M_PI);
  }

  results.push_back(deltres + distres);
  results.push_back(delterr + disterr);

  return results;
}
//==============================================================================================//
