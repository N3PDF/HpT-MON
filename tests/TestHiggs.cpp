/*
 * =====================================================================================
 *
 *       Filename:  TestHiggs.cpp
 *
 *    Description:  This file performs the benchmark of the full hadronic expression
 *                  and all the expressions it depends on with the HqT code
 *                  http://theory.fi.infn.it/grazzini/codes.html.
 *
 *        Version:  1.0
 *        Created:  18/02/2021 00:15:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#define CATCH_CONFIG_MAIN

#include <LHAPDF/LHAPDF.h>
#include "catch.hpp"
#include "../include/params.h"
#include "../include/higgspt.h"


// Settings in HqT (Match the physparams)
double X1 = 0.1;
double X2 = 0.2;
double Q = 125;
double MH = 125;
double MUR = 165;
double MUF = 165;
double SROOT = 13e3;
const char* PDF = "MSTW2008nnlo68cl";


TEST_CASE("Test distrcross class") {
    LHAPDF::setVerbosity(0);

    PhysParams physparam;
    physparam.alphas = 0.107536;    // alphas(Q=165GeV)
    physparam.nf     = 5.;
    physparam.mh     = 125.;
    physparam.mur    = 165.;
    physparam.muf    = 165.;
    physparam.sroot  = 13000.;
    physparam.sigma0 = 0.;

    // Init. higgspt class
    int order = 0, channel = 0;
    HiggsDpT higgs(order, channel, PDF, &physparam);

    // Define Phys. Variables
    double pt = 2.;
    double q2 = std::pow(physparam.mh,2);
    double tm = sqrt(q2+std::pow(pt,2));
    double y1 = -1.0;
    double y2 = 0.8;
    double zz1 = 0.10703379542565672;
    double zz2 = 0.24052913077693314;
    double zz3 = 0.77898332085577449;

    // Compute the splittings
    double distrcross = higgs.distrcross(pt,tm,y1,y2,zz1,zz2,zz3);

    // Result from HqT
    double Ddistrcross = -152114.75789147534;

    REQUIRE(Ddistrcross == Approx(distrcross).epsilon(1e-3));
}

TEST_CASE("Test distr class") {
    LHAPDF::setVerbosity(0);

    PhysParams physparam;
    physparam.alphas = 0.107536;    // alphas(Q=165GeV)
    physparam.nf     = 5.;
    physparam.mh     = 125.;
    physparam.mur    = 165.;
    physparam.muf    = 165.;
    physparam.sroot  = 13000.;
    physparam.sigma0 = 0.;

    // Init. higgspt class
    int order = 0, channel = 0;
    HiggsDpT higgs(order, channel, PDF, &physparam);

    // Define Phys. Variables
    double pt = 2.;
    double q2 = std::pow(physparam.mh,2);
    double tm = sqrt(q2+std::pow(pt,2));
    double y1 = -1.0;
    double y2 = 1.0;
    double zz1 = 0.10703379542565672;
    double zz2 = 0.24052913077693314;
    double zz3 = 0.77898332085577449;

    // Compute the splittings
    double distr = higgs.distr(pt,tm,y1,y2,zz1,zz2,zz3);

    // Results from HqT
    double Ddistr = -203250.61489061269;

    REQUIRE(Ddistr == Approx(distr).epsilon(1.5e-3));
}


TEST_CASE("Test delta function") {
    LHAPDF::setVerbosity(0);

    PhysParams physparam;
    physparam.alphas = 0.107536;    // alphas(Q=165GeV)
    physparam.nf     = 5.;
    physparam.mh     = 125.;
    physparam.mur    = 165.;
    physparam.muf    = 165.;
    physparam.sroot  = 13000.;
    physparam.sigma0 = 0.;

    // Init. higgspt class
    int order = 1, channel = 3;
    HiggsDpT higgs(order, channel, PDF, &physparam);

    // Define Phys. Variables
    double pt = 2.;
    double q2 = std::pow(physparam.mh,2);
    double tm = sqrt(q2+std::pow(pt,2));
    double t = 0.99489933627730553;
    double y1 = -4.6442628918517483;
    double y2 = 4.6442628918517483;
    double yc = 0.51081027482208352;

    // Compute the splittings
    double delta = higgs.delta(pt,tm,t,y1,y2,yc);

    // Results from HqT
    double Ddelta = 150.11378589543421;

    REQUIRE(Ddelta == Approx(delta).epsilon(1e-3));
}


TEST_CASE("Test channel-specific non-singular Functions") {
    LHAPDF::setVerbosity(0);

    PhysParams physparam;
    physparam.alphas = 0.107536;    // alphas(Q=165GeV)
    physparam.nf     = 5;
    physparam.mh     = 125;
    physparam.mur    = 165;
    physparam.muf    = 165;
    physparam.sroot  = 13e3;
    physparam.sigma0 = 0.;

    // Init. higgspt class
    int order = 0, channel = 0;
    HiggsDpT higgs(order, channel, PDF, &physparam);

    // Define Phys. Variables
    double pt = 2.;
    double s = 189510.29071559146;
    double t = -33859.927765699482;
    double u = -44228.706220675893;

    // Compute the splittings
    double Rgq = higgs.Rgq(pt,u,t,s);
    double Rqqb = higgs.Rqqb(pt,u,t,s);

    // Results from HqT
    double RRqqb = 66.556233071183797;

    REQUIRE(Rgq == Approx(Rgq).epsilon(1e-4));
    REQUIRE(RRqqb == Approx(Rqqb).epsilon(1e-4));
}


TEST_CASE("Test B Functions") {
    LHAPDF::setVerbosity(0);

    PhysParams physparam;
    physparam.alphas = 0.107536;    // alphas(Q=165GeV)
    physparam.nf     = 5;
    physparam.mh     = 125;
    physparam.mur    = 165;
    physparam.muf    = 165;
    physparam.sroot  = 13e3;
    physparam.sigma0 = 0.;

    // Init. higgspt class
    int order = 0, channel = 0;
    HiggsDpT higgs(order, channel, PDF, &physparam);

    // Define Phys. Variables
    double pt = 2.;
    double s = 189510.29071559146;
    double t = -33859.927765699482;
    double u = -44228.706220675893;

    // Compute the splittings
    double B1pmB1pp = higgs.B1pmB1pp(pt,u,t,s);
    double B2pmB2pp = higgs.B2pmB2pp(pt,u,t,s);

    // Results from HqT
    double BB1pmB1pp = 1.0471600598868515E+017;
    double BB2pmB2pp = -1.2753537934175438E+017;

    REQUIRE(BB1pmB1pp == Approx(B1pmB1pp).epsilon(1e-4));
    REQUIRE(BB2pmB2pp == Approx(B2pmB2pp).epsilon(1e-4));
}


TEST_CASE("Test A Functions") {
    LHAPDF::setVerbosity(0);

    PhysParams physparam;
    physparam.alphas = 0.107536;    // alphas(Q=165GeV)
    physparam.nf     = 5;
    physparam.mh     = 125;
    physparam.mur    = 165;
    physparam.muf    = 165;
    physparam.sroot  = 13e3;
    physparam.sigma0 = 0.;

    // Init. higgspt class
    int order = 0, channel = 0;
    HiggsDpT higgs(order, channel, PDF, &physparam);

    // Define Phys. Variables
    double pt = 2.;
    double s = 189510.29071559146;
    double t = -33859.927765699482;
    double u = -44228.706220675893;

    // Compute the splittings
    double A1234 = higgs.A1234(pt,u,t,s);
    double A3412 = higgs.A3412(pt,u,t,s);
    double A1324 = higgs.A1324(pt,u,t,s);
    double A3241 = higgs.A3241(pt,u,t,s);

    // Results from HqT
    double AA1234 = 1934269973561736.0;
    double AA3412 = 5.9838541894288346E+017;
    double AA1324 = 1944656843036286.5;
    double AA3241 = -16002762776884736.;

    REQUIRE(AA1234 == Approx(A1234).epsilon(1e-4));
    REQUIRE(AA3412 == Approx(A3412).epsilon(1e-4));
    REQUIRE(AA1324 == Approx(A1324).epsilon(1e-4));
    REQUIRE(AA3241 == Approx(A3241).epsilon(1e-4));
}


TEST_CASE("Test Splitting Functions") {
    LHAPDF::setVerbosity(0);

    PhysParams physparam;
    physparam.alphas = 0.107536;    // alphas(Q=165GeV)
    physparam.nf     = 5;
    physparam.mh     = 125;
    physparam.mur    = 165;
    physparam.muf    = 165;
    physparam.sroot  = 13e3;
    physparam.sigma0 = 0.;

    double YH = -3.0882480144500732;
    double TM = 125.01599884033203;

    // Init. higgspt class
    int order = 0, channel = 0;
    HiggsDpT higgs(order, channel, PDF, &physparam);

    // Define Phys. Variables
    double pt = 2.;
    double q2 = std::pow(physparam.mh,2);
    double tm = sqrt(q2+std::pow(pt,2));

    // Compute the splittings
    double gg0 = higgs.gg0(X1,X2,tm,YH);
    double qg0 = higgs.qg0(X1,X2,tm,YH);
    double gq0 = higgs.gq0(X1,X2,tm,YH);
    double qqb0 = higgs.qqb0(X1,X2,tm,YH);

    // Results from HqT
    double GG0 = -89411809094.238815;
    double QG0 = -39595825124.684952;
    double GQ0 = 4291030.0243432121;
    double QQB0 = 13255986.327664498;

    // Perform checks with tolerance
    REQUIRE(TM == Approx(tm).epsilon(1e-4));
    REQUIRE(GG0 == Approx(gg0).epsilon(1e-4));
    REQUIRE(QG0 == Approx(qg0).epsilon(1e-4));
    REQUIRE(GQ0 == Approx(gq0).epsilon(1e-4));
    REQUIRE(QQB0 == Approx(qqb0).epsilon(1e-4));
}
