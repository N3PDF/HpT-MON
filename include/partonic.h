/*
 * =====================================================================================
 *
 *       Filename:  partonic.h
 *
 *    Description:  Header file for the partonic.cpp file.
 *
 *        Version:  1.0
 *        Created:  20/02/2021 22:21:31
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
#include <iostream>
#include <string>
#include <vector>

#include "./higgsptpartonic.h"
#include "./integration.h"


class CrossHiggs {
 public:
    CrossHiggs(
        int order,
        int channel,
        std::string pdfname,
        void *params);
    virtual ~CrossHiggs();

    // Main function that performs the
    // integration over rapidities
    std::vector<double> partonichiggsdpt(double pt, double nn);
    double ExtractDeltaPartonic(double pt, double nn, double zz);
    double ExtractDistrPartonic(double pt, double nn, double zz1, double zz2);

 private:
    double NF, MH2, aass, SIGMA0, ORD, CHANNEL;

    // Instantiate HiggsDpT class
    HiggsDpTpartonic higgs;
};
