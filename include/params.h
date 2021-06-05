/*
 * =====================================================================================
 *
 *       Filename:  params.h
 *
 *    Description:  Definition of physical parameters.
 *
 *        Version:  1.0
 *        Created:  17/02/2021 23:16:58
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#pragma once

#include <string>

struct PhysParams{
    int nc;
    int nf;
    double mh;
    double mur;
    double muf;
    double alphas;
    double sroot;
    double sigma0;
};

extern PhysParams globalStruct;
