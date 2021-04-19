/*
 * =====================================================================================
 *
 *       Filename:  higgsptpartonic.h
 *
 *    Description:  Header file for the higgsptpartonic.cpp file.
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_result.h>

#include <iostream>
#include <complex>
#include <functional>
#include <tuple>
#include <string>

#include "./utils.h"
#include "./params.h"
#include "./splittings.h"


class HiggsDpTpartonic{
 public:
    HiggsDpTpartonic(
        int order,
        int channel,
        std::string pdfname,
        void *params);

    virtual ~HiggsDpTpartonic();

    // Partonic functions
    double gg0(double sh, double th, double uh, double tauh);
    double gq0(double sh, double th, double uh);
    double qg0(double sh, double th, double uh);
    double qqb0(double sh, double th, double uh);

    // A functions
    double A1234(double pt, double u, double t, double s, double q2);
    double A3412(double pt, double u, double t, double s, double q2);
    double A1324(double pt, double u, double t, double s, double q2);
    double A3241(double pt, double u, double t, double s, double q2);

    // B functions
    double B1pmB1pp(double pt, double u, double t, double s, double q2);
    double B2pmB2pp(double pt, double u, double t, double s, double q2);

    // Regular part of the channel specific coefficient
    // functions
    double Rgq(double pt, double u, double t, double s, double q2);
    double Rqqb(double pt, double u, double t, double s, double q2);

    // Terms that enter in the
    // partonic cross-section
    double deltapartonic(double pt, double nn, double zz);
    double distrpartonic(double pt, double nn, double zz1, double zz2);

 private:
    int NF, ORD, CHANNEL;
    double MH2, MUF2, MUR2;
    double beta0, aass, SIGMA0;

    void coeff(double pt, double uh, double th, double sh, double q2);
    double big1, big2, big3, big4, big5;

    void REG(double pt, double uh, double th, double sh, double q2);
    double REGgg, REGgq, REGqg, REGqqb, REGqqpb, REGqq;
};
