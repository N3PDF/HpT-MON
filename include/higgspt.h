/*
 * =====================================================================================
 *
 *       Filename:  higgspt.h
 *
 *    Description:  Header file for the higgspt.cpp file.
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

#include <tuple>
#include <string>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_result.h>
#include <iostream>
#include <complex>
#include <functional>

#include "./utils.h"
#include "./params.h"
#include "./luminosity.h"
#include "./splittings.h"


class HiggsDpT {
 public:
    HiggsDpT(
            int order,
            int channel,
            std::string pdfname,
            void *params);

    virtual ~HiggsDpT();

    // Partonic functions
    double gg0(double x1, double x2, double tm, double yh);
    double gq0(double x1, double x2, double tm, double yh);
    double qg0(double x1, double x2, double tm, double yh);
    double qqb0(double x1, double x2, double tm, double yh);

    // A functions
    double A1234(double pt, double u, double t, double s);
    double A3412(double pt, double u, double t, double s);
    double A1324(double pt, double u, double t, double s);
    double A3241(double pt, double u, double t, double s);

    // B functions
    double B1pmB1pp(double pt, double u, double t, double s);
    double B2pmB2pp(double pt, double u, double t, double s);

    // Regular part of the channel specific coefficient
    // functions
    double Rgq(double pt, double u, double t, double s);
    double Rqqb(double pt, double u, double t, double s);

    // Terms that enter in the
    // hadronic cross-sections
    double delta(
            double pt,
            double tm,
            double t,
            double y1,
            double y2,
            double yc);

    double distr(
            double pt,
            double tm,
            double y1,
            double y2,
            double zz1,
            double zz2,
            double zz3);

    double distrcross(
            double pt,
            double tm,
            double y1,
            double y2,
            double zz1,
            double zz2,
            double zz3);

 private:
    Luminosity *Lumi;

    int NF, ORD, CHANNEL;
    double MH2, MUF2, MUR2, SROOT;
    double beta0, aass, SIGMA0;
    double xtau, shad, q2;

    void coeff(double pt, double uh, double th, double sh);
    double big1, big2, big3, big4, big5;

    void REG(double pt, double uh, double th, double sh);
    double REGgg, REGgq, REGqg, REGqqb, REGqqpb, REGqq;
};
