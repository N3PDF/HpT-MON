/*
 * =====================================================================================
 *
 *       Filename:  splittings.cpp
 *
 *    Description:  This file contains the expression of the splitting
 * functions.
 *
 *        Version:  1.0
 *        Created:  18/02/2021 12:34:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#include "../include/splittings.h"

double pgg(double z) {
  return 3. * (1 + std::pow(z, 4) + std::pow(1 - z, 4)) / z;
}

double pqq(double z) { return 4. / 3. * (1 + std::pow(z, 2)); }

double Pgq(double z) { return 4. / 3. * (1 + std::pow(1. - z, 2)) / z; }

double Pqg(double z) { return 0.5 * (std::pow(z, 2) + std::pow(1 - z, 2)); }
