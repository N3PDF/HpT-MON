/*
 * =====================================================================================
 *
 *       Filename:  utils.cpp
 *
 *    Description:  This file contains various routines.
 *
 *        Version:  1.0
 *        Created:  20/02/2021 22:07:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#include "../include/utils.h"

//==============================================================================================//
//                                    Some math functions //
//----------------------------------------------------------------------------------------------//
double Li2(double val) { return gsl_sf_dilog(val); }
