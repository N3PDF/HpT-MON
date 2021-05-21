/*
 * =====================================================================================
 *
 *       Filename:  luminosity.cpp
 *
 *    Description:  This file computes the product of two PDFs, for convenience, we call
 *                  the product luminosity.
 *
 *        Version:  1.0
 *        Created:  25/02/2021 21:42:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tanjona R. Rabemananjara, Roy Stegeman 
 *   Organization:  N3PDF
 *
 * =====================================================================================
 */

#include "../include/luminosity.h"


Luminosity::Luminosity(
        std::string pdfname,
        double MUF, int NF
) {
    _nf  = NF;
    _muf = MUF;
    _pdf = LHAPDF::mkPDF(pdfname);
}


Luminosity::~Luminosity() {}


// compute pdf //
double Luminosity::fxQ(double x, double pidflv) {
    return (_pdf->xfxQ(pidflv,x,_muf))/x;
}

// Extract alphas //
double Luminosity::get_alphas(double mur) {
    return _pdf->alphasQ(mur);
}


//// compute luminosity ////
double Luminosity::Lumgg(double x1, double x2) {
    return fxQ(x1,0)*fxQ(x2,0);
}


double Luminosity::Lumgq(double x1, double x2) {
    double lgq = 0;

    for (int i = 0; i < _nf; i++) {
        lgq += (fxQ(x2,(i+1))+fxQ(x2,-(i+1))) \
               *fxQ(x1,0);
    }
    return lgq;
}


double Luminosity::Lumqg(double x1, double x2) {
    return Lumgq(x2,x1);
}


double Luminosity::Lumqq(double x1, double x2) {
    double lqq = 0;

    for (int i = 0; i < _nf; i++) {
        lqq += fxQ(x1,(i+1))*fxQ(x2,(i+1)) \
                +fxQ(x1,-(i+1))*fxQ(x2,-(i+1));
    }
    return lqq;
}


double Luminosity::Lumqqb(double x1, double x2) {
    double lqqb = 0;

    for (int i = 0; i < _nf; i++) {
        lqqb += fxQ(x1,(i+1))*fxQ(x2,-(i+1)) \
                + fxQ(x2,(i+1))*fxQ(x1,-(i+1));
    }
    return lqqb;
}


double Luminosity::LumqQ(double x1, double x2) {
    double lqQ = 0;

    for (int i = 0; i < _nf; i++) {
        for (int j = 0; j < _nf; j++) {
            if (i != j) \
            lqQ += fxQ(x1,(i+1))*fxQ(x2,(j+1)) \
                    +fxQ(x1,-(i+1))*fxQ(x2,-(j+1));
        }
    }
    return lqQ;
}


double Luminosity::LumqQb(double x1, double x2) {
    double lqQb = 0;

    for (int i = 0; i < _nf; i++) {
        for (int j = 0; j < _nf; j++) {
            if (i != j) \
            lqQb += fxQ(x1,-(i+1))*fxQ(x2,(j+1)) \
                    +fxQ(x1,(i+1))*fxQ(x2,-(j+1));
        }
    }
    return lqQb;
}


//// compute TRUE luminosity ////
// gg channel
double f_gg(double x1, void *p)
{
	Luminosity::par_struct* p_struct = (Luminosity::par_struct *) p;
	double x2 = p_struct -> t /x1;
	double Muf = p_struct -> Muf;
	Luminosity::PDF_ptr PDF = p_struct -> PDF;
	return (PDF-> xfxQ(0,x1,Muf)/x1*PDF->xfxQ(0,x2,Muf)/x2)/x1;
}


double Luminosity::LumNgg(double x) {
	par_struct param;
	param.Muf = _muf;
	param.PDF = _pdf;
	param.t = x;
	param.Nf = _nf;

	double reslt = 0.;
	double error = 0.;
	double precision = 1e-4;

	gsl_function Integrand;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	Integrand.function = f_gg;
	Integrand.params = param;

	gsl_integration_qags(&Integrand,x,1.0,0,precision,1000,w,&reslt,&error);
	gsl_integration_workspace_free(w);
	return reslt;
}
