/*
 * =====================================================================================
 *
 *       Filename:  higgspt.cpp
 *
 *    Description:  This file contains all the terms contribuing to the Higgs fixed
 *                  order results.
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

#include "../include/higgspt.h"


HiggsDpT::HiggsDpT(
        int order,
        int channel,
        std::string pdfname,
        void *params
) {
    PhysParams param = *reinterpret_cast<PhysParams*>(params);

    NF = param.nf;
    MH2  = std::pow(param.mh, 2);
    MUR2 = std::pow(param.mur, 2);
    MUF2 = std::pow(param.muf, 2);
    aass = param.alphas/M_PI;
    SROOT = param.sroot;
    SIGMA0 = param.sigma0;

    ORD = order;
    CHANNEL = channel;

    q2 = MH2;
    beta0 = (33-2*NF)/6.;
    shad = std::pow(SROOT,2);
    xtau = q2/shad;

    Lumi = new Luminosity(
            pdfname,
            param.mur,
            param.nf);
}


HiggsDpT::~HiggsDpT() {}


//==============================================================================================//
//                          Partonic Functions (Eqs. [2.9])                                     //
//----------------------------------------------------------------------------------------------//
double HiggsDpT::gg0(double x1, double x2, double tm, double yh) {
    double sh = x1*x2*shad;
    double th = q2-SROOT*x2*tm*exp(yh);
    double uh = q2-SROOT*x1*tm*exp(-yh);
    return 3*(std::pow(q2,4)+std::pow(sh,4)+std::pow(th,4)+std::pow(uh,4))/(uh*th*sh);
}


double HiggsDpT::qg0(double x1, double x2, double tm, double yh) {
    double sh = x1*x2*shad;
    double th = q2-SROOT*x2*tm*exp(yh);
    double uh = q2-SROOT*x1*tm*exp(-yh);
    return 4./3.*(std::pow(uh,2)+std::pow(sh,2))/(-th);
}


double HiggsDpT::gq0(double x1, double x2, double tm, double yh) {
    double sh = x1*x2*shad;
    double th = q2-SROOT*x2*tm*exp(yh);
    double uh = q2-SROOT*x1*tm*exp(-yh);
    return 4./3.*(std::pow(th,2)+std::pow(sh,2))/(-uh);
}


double HiggsDpT::qqb0(double x1, double x2, double tm, double yh) {
    double sh = x1*x2*shad;
    double th = q2-SROOT*x2*tm*exp(yh);
    double uh = q2-SROOT*x1*tm*exp(-yh);
    return 2.*16./9.*(std::pow(th,2)+std::pow(uh,2))/sh;
}
//==============================================================================================//



//==============================================================================================//
//                        Regular Part of Coefficient Functions                                 //
//----------------------------------------------------------------------------------------------//
void HiggsDpT::REG(double pt, double uh, double th, double sh) {
    ///////////////////////////////////////////////////////////
    // Below are terms that contribute to the NonSingular    //
    // Real Contribution in Appendix A of G & S. Secifically //
    // these are given by Eqs. (A.1),(A.17),(A.27),(A.35)    //
    ///////////////////////////////////////////////////////////

    double pt2 = std::pow(pt,2);
    double QQ2 = uh+th+sh-q2;
    double xa = th/(th-QQ2);
    double xb = uh/(uh-QQ2);
    double xQt2 = QQ2+pt2;

    double E2 = sh*pt2*QQ2*(QQ2-pt2)/std::pow(pt2+QQ2,2)*(std::pow(QQ2-th,2)+std::pow(QQ2-uh,2)) \
                +2.*std::pow(sh,2)*pt2*QQ2;
    double E4 = 2.*sh*QQ2*pt2*(std::pow(sh,2)+std::pow(QQ2,2))/(pt2+QQ2)*log(pt2/QQ2) \
                +4.*std::pow(sh,2)*pt2*QQ2*log(pt2/(pt2+QQ2));

    // eq. (A.4)
    double A0 = (std::pow(th/xa,4)+std::pow(uh/xb,4))*pt2*QQ2/std::pow(xQt2,2) \
                *(5.-7.*QQ2/xQt2+20./3.*std::pow(QQ2/xQt2,2))+std::pow(sh,2)*QQ2 \
                *pt2*(17./3.+4.*log(pt2/xQt2));
    double Aep = 4.*std::pow(sh*pt2*q2*QQ2,2)*(1./std::pow(th,4)+1./std::pow(uh,4));


    // Regular part for specific channel
    REGgg = 1./std::pow(sh,2)/pt2/QQ2*(9.*(A0+Aep+A1234(pt,uh,th,sh)+A1234(pt,th,uh,sh) \
            +A3412(pt,uh,th,sh)+A3412(pt,th,uh,sh)+A1324(pt,uh,th,sh)+A1324(pt,th,uh,sh) \
            +A3241(pt,uh,th,sh)+A3241(pt,th,uh,sh))+4./3.*NF*(B1pmB1pp(pt,uh,th,sh) \
            +B1pmB1pp(pt,th,uh,sh))+3.*NF*(B2pmB2pp(pt,uh,th,sh)+B2pmB2pp(pt,th,uh,sh)));
    REGgq = Rgq(pt,uh,th,sh);
    REGqg = Rgq(pt,th,uh,sh);
    REGqqb = Rqqb(pt,uh,th,sh)+Rqqb(pt,th,uh,sh);
    REGqqpb = 16./9./std::pow(sh,2)/pt2/QQ2*E2;
    REGqq = 16./9./std::pow(sh,2)/pt2/QQ2*(E2+E4/3.);
}
//==============================================================================================//



//==============================================================================================//
//                                      Coefficients                                            //
//----------------------------------------------------------------------------------------------//
void HiggsDpT::coeff(double pt, double uh, double th, double sh) {
    double tiny = 1e-4;
    double QQ2 = uh+sh+th-q2;
    double xa = th/(th-QQ2);
    double xb = uh/(uh-QQ2);
    double tmp = tiny*std::pow(pt,2);

    ////// big1 & big2 & big3 are the terms that contribute to the gg///////
    // Term that gets multiplied by NC^2/2 in Eq. (3,17) (inside the curly brackets)
    big1 = 0.5*9.*((std::pow(q2,4)+std::pow(sh,4)+std::pow(QQ2,4)+std::pow(uh,4)+std::pow(th,4)) \
            +xa*xb*(std::pow(q2,4)+std::pow(sh,4)+std::pow(QQ2,4)+std::pow(uh/xb,4) \
            +std::pow(th/xa,4)))/sh/uh/th;

    // Term that gets multiplied by Beat0*NC in Eq. (3,17) (inside the curly brackets)
    big2 = beta0*3./2.*\
            (std::pow(q2,4)+std::pow(sh,4)+xa*xb*(std::pow(uh/xb,4)+std::pow(th/xa,4)))/sh/uh/th;

    // Terms that get multiplied by NC^2 in the square brackets
    if (QQ2 > tmp) {
    big3 = 9.*((std::pow(q2,4)+std::pow(sh,4)+std::pow(QQ2,4)+std::pow(uh/xb,4)+std::pow(th/xa,4)) \
            *(2*QQ2+std::pow(pt,2))/std::pow(sh,2)/QQ2/(QQ2+std::pow(pt,2))+2*std::pow(q2,2) \
            *(std::pow(q2-th,4)+std::pow(q2-uh,4)+std::pow(uh,4)+std::pow(th,4))/sh/uh/th/(q2-uh) \
            /(q2-th))/std::pow(pt,2)*log(std::pow(pt,2)/(std::pow(pt,2)+QQ2));
    } else {
        big3 = -9.*((std::pow(q2,4)+std::pow(sh,4)+std::pow(uh,4)+std::pow(th,4))/std::pow(sh,2) \
                /std::pow(pt,4));
    }

    ////// big4 & big5 are the terms that contribute to the qg///////
    // Term that gets multiplied by the first NC*CF in Eq. (3.20) (inside square bracket)
    big4 = 4.*((-std::pow(sh,3)*th-sh*std::pow(th,3)+std::pow(QQ2,3)*th+QQ2*std::pow(th,3)) \
            /sh/uh/th+xa*xb*(-std::pow(sh,3)*th/xa-sh*std::pow(th/xa,3)-std::pow(QQ2,3)*uh/xb \
            -QQ2*std::pow(uh/xb,3))/(sh*uh*th));

    // Term that gets multiplied by the seond NC*CF in Eq. (3.20) (inside square bracket)
    if (QQ2 > tmp) {
        big5 = 4.*((-std::pow(sh,3)*th/xa-sh*std::pow(th/xa,3)-std::pow(QQ2,3)*uh/xb \
                -QQ2*std::pow(uh/xb,3))*(2.*QQ2+std::pow(pt,2))/std::pow(sh,2) \
                /(QQ2+std::pow(pt,2))-2.*QQ2*std::pow(q2,2)*(std::pow(q2-th,2) \
                +std::pow(th,2))/sh/uh/(q2-uh))*log(std::pow(pt,2)/(std::pow(pt,2)+QQ2)) \
                /(QQ2*std::pow(pt,2));
    } else {
        big5 = -4.*(-std::pow(sh,3)*th/xa-sh*std::pow(th/xa,3))/std::pow(sh,2)/std::pow(pt,4);
    }
}



//==============================================================================================//
//                                      A Functions                                             //
//----------------------------------------------------------------------------------------------//
double HiggsDpT::A1234(double pt, double u, double t, double s) {
    double pt2 = std::pow(pt,2);
    double q4 = std::pow(q2,2);
    double QQ2 = u+t+s-q2;
    double x1 = t/(t-QQ2);
    double x2 = u/(u-QQ2);
    double t2 = std::pow(t,2);
    double t4 = std::pow(t,4);
    double u2 = std::pow(u,2);
    double u3 = std::pow(u,3);
    double u4 = std::pow(u,4);

    // Eqs. (A.2)-(A.3) G & S
    double A = s+q2-QQ2;
    double B = sqrt(std::pow(A,2)-4.*q2*s);
    double L1a = log(q2/s/std::pow(x1,2));
    double L1b = log(q2/s/std::pow(x2,2));
    double L2a = log(q2*s/std::pow(A-s*x1,2));
    double L2b = log(q2*s/std::pow(A-s*x2,2));
    double L3 = log((A+B)/(A-B));

    double result = -0.5*((std::pow(s*pt2/t,4)+std::pow(q2*QQ2/t,4))*L1a+std::pow(u,4)*L1b) \
                    +0.5*s*q4*QQ2*u3/t/(u-q2)/(t-q2)*(L1a+L1b)\
                    +0.5*s*pt2*QQ2*u3/A/(t-q2)*(L2b-L1b)\
                    +0.5*s*pt2*QQ2*(std::pow(q4,2)+std::pow(u-q2,4))/A/t/(u-q2)*(L2a-L1a) \
                    +s*pt2*q2*QQ2*(u4/std::pow(B*t,2)+u2/std::pow(B,2))/2.\
                    +std::pow(s*pt2*q2*QQ2,2)*(-6./std::pow(B,4)-4./t4+8./std::pow(B*t,2))\
                    +L3*(s*pt2*u3*(u+t)/B/t\
                    +std::pow(s*pt2*q2*QQ2,2)*(3.*A/std::pow(B,5)-1./A/std::pow(B,3))\
                    -s*pt2*q2*QQ2*((t2+t*u+4*std::pow(u,2)-2.*q2*QQ2)/B/t\
                    +A*(t2+3.*t*u+3.*std::pow(u,2)-6.*q2*QQ2)/2./std::pow(B,3)\
                    +(t2+t*u+7.*u2-2.*q2*QQ2)/(2*A*B)));

    return result;
}


double HiggsDpT::A3412(double pt, double u, double t, double s) {
    double pt2 = std::pow(pt,2);
    double q4 = std::pow(q2,2);
    double q6 = std::pow(q2,3);
    double q8 = std::pow(q2,4);
    double QQ2 = u+t+s-q2;
    double x1 = t/(t-QQ2);
    double x2 = u/(u-QQ2);
    double xQt2 = QQ2+pt2;
    double t2 = std::pow(t,2);
    double t3 = std::pow(t,3);
    double u2 = std::pow(u,2);
    double u3 = std::pow(u,3);
    double u4 = std::pow(u,4);
    double s2 = std::pow(s,2);
    double s3 = std::pow(s,3);
    double s4 = std::pow(s,4);

    // Eqs. (A.2)-(A.3) G & S
    double A = s+q2-QQ2;
    double B = sqrt(std::pow(A,2)-4.*q2*s);
    double L1b = log(q2/s/std::pow(x2,2));
    double L2a = log(q2*s/std::pow(A-s*x1,2));
    double L3 = log((A+B)/(A-B));

    double result = s*pt2*QQ2*std::pow(A,3)/2./t/(u-q2)*(L2a+L1b)+s*pt2*(u+t)/16/u/t/B*(std::pow(A,4) \
                    +6.*std::pow(A,2)*std::pow(B,2)+std::pow(B,4))*L3+(-s*pt2/2/u/t*(std::pow(s-QQ2,4) \
                    +q8+2*QQ2*A*std::pow(s-QQ2,2)-2.*QQ2*q6)-std::pow(s*pt2*q2*QQ2,2)/u4+2.*s*pt2*q2 \
                    *QQ2*(std::pow(A,2)-s*q2)/u2)*L1b+s*pt2/(8*u*t)*(std::pow(QQ2-u,3)/(QQ2-t) \
                    *((QQ2-t)*s+QQ2*u)+std::pow(QQ2-t,3)/(QQ2-u)*((QQ2-u)*s+QQ2*t))*(4./3+2*pt2 \
                    /xQt2+4.*std::pow(pt2/xQt2,2)-44./3.*std::pow(pt2/xQt2,3))+s*pt2*std::pow(QQ2-u,2) \
                    /4./u/t/(QQ2-t)*(-3.*(t-q2)*((QQ2-t)*s+QQ2*u)-QQ2*(q2*(t-q2)+QQ2*(u-q2))) \
                    *(1.+2.*pt2/xQt2-6.*std::pow(pt2/xQt2,2))+s*pt2*std::pow(QQ2-t,2)/4./u/t/(QQ2-u) \
                    *(-3.*(u-q2)*((QQ2-u)*s+QQ2*t)+3*QQ2*(q2*(u-q2)+QQ2*(t-q2))+4.*u*s2)*(1.+2. \
                    *pt2/xQt2-6.*std::pow(pt2/xQt2,2))+s*pt2*(QQ2-u)/2./u/t/(QQ2-t)*(3.*std::pow(t-q2,2) \
                    *((QQ2-t)*s+QQ2*u)+3.*(t-q2)*QQ2*(q2*(t-q2)+QQ2*(u-q2))+QQ2*u*(q2*(t-QQ2)+QQ2 \
                    *(u-q2)))*(1.-2.*pt2/xQt2)+s*pt2*(QQ2-t)/2/u/t/(QQ2-u)*(3*std::pow(u-q2,2) \
                    *((QQ2-u)*s+QQ2*t)+8.*u*t*s2+2.*u*s3-2.*QQ2*u*std::pow(u-QQ2,2)-3.*q2*QQ2 \
                    *std::pow(t-q2,2)-3.*QQ2*(q2-QQ2)*t2-QQ2*u*(4.*u*t-u*q2-QQ2*t+2.*t2-4.*q4)+3. \
                    *q2*std::pow(QQ2,2)*(t-q2)+q2*QQ2*u*(t-QQ2))*(1.-2.*pt2/xQt2)-4.*std::pow(s* \
                    pt2*q2*QQ2,2)/u4+s*pt2*q2*QQ2*(std::pow(B,2))/(2.*u2)+s2*pt2*q4/6.*((s+QQ2)/u \
                    /t+QQ2/u2+QQ2/t2)+2.*s2*pt2*QQ2*q4/u2+s2*pt2*q4/u-s2*pt2/12./u/t*(30.*q6+54. \
                    *(std::pow(QQ2,2))*q2+8.*std::pow(QQ2,3))+s*pt2/12./u/t*(11.*s4+17.*q8+QQ2* \
                    (61.*u2*t+17.*u3+73.*u*t2+29.*t3)+q2*(24.*u2*t+6.*u3+36.*u*t2+18.*t3)+std::pow(QQ2,2) \
                    *(-21.*u2-33.*t2-52.*u*t)+q2*QQ2*(-73.*u2-109.*t2-170.*u*t)+q4*(-23.*u2-35.*t2-52.*u \
                    *t)+q4*QQ2*(134.*t+110.*u)+4.*std::pow(QQ2,4)+52.*q2*std::pow(QQ2,3)+20.*q4* \
                    std::pow(QQ2,2)-22.*q6*QQ2);

    return result;
}


double HiggsDpT::A1324(double pt, double u, double t, double s) {
    double pt2 = std::pow(pt,2);
    double q4 = std::pow(q2,2);
    double q8 = std::pow(q2,4);
    double QQ2 = u+t+s-q2;
    double x1 = t/(t-QQ2);
    double x2 = u/(u-QQ2);
    double t2 = std::pow(t,2);
    double t4 = std::pow(t,4);
    double u3 = std::pow(u,3);
    double u4 = std::pow(u,4);
    double s2 = std::pow(s,2);

    // Eqs. (A.2)-(A.3) G & S
    double A = s+q2-QQ2;
    double B = sqrt(std::pow(A,2)-4.*q2*s);
    double L1a = log(q2/s/std::pow(x1,2));
    double L1b = log(q2/s/std::pow(x2,2));
    double L2a = log(q2*s/std::pow(A-s*x1,2));
    double L2b = log(q2*s/std::pow(A-s*x2,2));
    double L3 = log((A+B)/(A-B));

    double result = -0.5*((std::pow(s*pt2/t,4)+std::pow(q2*QQ2/t,4))*L1a+u4*L1b)+s2*pt2*q4*QQ2 \
                    /t2*L1a+(s*q4*QQ2*u3/t/(u-q2)/(t-q2)+s*pt2*u3/t)/2*(L1a+L1b)+s2*pt2*(1-x2) \
                    *u3/A/2./(t-q2)*(L2b-L1b)+s2*pt2*(1-x1)*(q8+std::pow(u-q2,4))/2./A/t/(u-q2) \
                    *(L2a-L1a)+s2*pt2*q4*QQ2/A/B*L3+s*pt2*q2*QQ2/2./t4*(std::pow(s*pt2,2)-6.*s \
                    *pt2*q2*QQ2+q4*std::pow(QQ2,2));

    return result;
}


double HiggsDpT::A3241(double pt, double u, double t, double s) {
    double pt2 = std::pow(pt,2);
    double q4 = std::pow(q2,2);
    double QQ2 = u+t+s-q2;
    double x1 = t/(t-QQ2);
    double xQt2 = QQ2+pt2;
    double t2 = std::pow(t,2);
    double t4 = std::pow(t,4);
    double u2 = std::pow(u,2);
    double s2 = std::pow(s,2);

    // Eqs. (A.2)-(A.3) G & S
    double A = s+q2-QQ2;
    double B = sqrt(std::pow(A,2)-4.*q2*s);
    double L1a = log(q2/s/std::pow(x1,2));
    double L2a = log(q2*s/std::pow(A-s*x1,2));

    double result = s2*pt2*std::pow(A,3)*(1-x1)/2./t/(u-q2)*(L2a-L1a)\
                    +(-std::pow(s*pt2*q2*QQ2,2)/t4+s*pt2*q4*std::pow(QQ2,2)/u/t\
                    -s*pt2*q2*QQ2*std::pow(A,4)/2./u/t/(u-q2)/(t-q2)\
                    +s*pt2*QQ2*q2*(u+t)*(2.*std::pow(A,2)-s*q2)/u/t2)*L1a\
                    +s2*pt2*QQ2*std::pow(QQ2-u,2)/2./u/t/std::pow(QQ2-t,2)\
                    *(-u*t-std::pow(QQ2-t,2))*(-3.+10.*QQ2/xQt2-6.*std::pow(QQ2/xQt2,2))\
                    +s2*pt2*QQ2*(QQ2-u)/u/t/std::pow(QQ2-t,2)*(u*t*(QQ2-u)-std::pow(QQ2-t,3)\
                    -q2*std::pow(QQ2-t,2)-q2*(QQ2-t)*(QQ2-u))*(-1.+2.*QQ2/xQt2) \
                    +s*pt2*q2*QQ2*(std::pow(B,2)/2./t2-2.*q2*QQ2/t2+std::pow(u+t,2)/2./u/t)\
                    -4.*std::pow(s*pt2*q2*QQ2,2)/t4+s2*pt2*QQ2/4./u/t\
                    *(std::pow(t+u,2)-(t+u)*(6.*QQ2+4.*q2)+6.*std::pow(QQ2,2)+8.*q2*QQ2)\
                    +s2*pt2*QQ2*q4*std::pow(t+u,2)/4./u2/t2;

    return result;
}
//==============================================================================================//



//==============================================================================================//
//                                      B Functions                                             //
//----------------------------------------------------------------------------------------------//
double HiggsDpT::B1pmB1pp(double pt, double u, double t, double s) {
    double pt2 = std::pow(pt,2);
    double QQ2 = u+t+s-q2;
    double x1 = t/(t-QQ2);
    double xQt2 = QQ2+pt2;
    double t2 = std::pow(t,2);
    double s2 = std::pow(s,2);

    double B1pm = std::pow(s2,2)*pt2*x1*std::pow(1-x1,3)/t+std::pow(s2,2) \
                *std::pow(pt2,3)*std::pow(x1,3)*(1-x1)/std::pow(t,3)+4 \
                *std::pow(s2,2)*std::pow(pt2,2)*std::pow(x1,2)*std::pow(1-x1,2)/t2 \
                -s2*pt2*QQ2*(1+log(pt2/xQt2));
    double B1pp = std::pow(s2,2)*pt2*std::pow(x1,3)*(1-x1)/t+std::pow(s2,2) \
                *std::pow(pt2,3)*x1*std::pow(1-x1,3)/std::pow(t,3)+4*std::pow(s2,2) \
                *std::pow(pt2,2)*std::pow(x1,2)*std::pow(1-x1,2)/t2-s2*pt2 \
                *QQ2/std::pow(xQt2,2)/u/t*(std::pow(u*t+pt2*QQ2,2)+2*s*pt2*QQ2*xQt2) \
                +s2*pt2*QQ2/u/t*(s2+std::pow(QQ2,2))*log(xQt2/QQ2);

    return B1pm+B1pp;
}

double HiggsDpT::B2pmB2pp(double pt, double u, double t, double s) {
    double pt2 = std::pow(pt,2);
    double q4 = std::pow(q2,2);
    double q6 = std::pow(q2,3);
    double QQ2 = u+t+s-q2;
    double x1 = t/(t-QQ2);
    double xQt2 = QQ2+pt2;
    double t2 = std::pow(t,2);
    double u2 = std::pow(u,2);
    double s2 = std::pow(s,2);
    double s3 = std::pow(s,3);
    double s4 = std::pow(s,4);

    double B2pm = 1./3.*std::pow(t/x1,4)*pt2/xQt2*(std::pow(pt2,3)-std::pow(QQ2,3)-std::pow(xQt2,3)) \
                    /std::pow(xQt2,3)-s2*pt2*QQ2/3.;

    double B2pp = -s2*pt2*QQ2/(2.*u*t)*(s2+std::pow(QQ2,2))*log(xQt2/QQ2)\
                +s*pt2*std::pow(QQ2-u,3)/(2.*u*t)/(QQ2-t)*((QQ2-t)*s+QQ2*u)\
                *(2./3.+QQ2/xQt2-10./3.*std::pow(QQ2/xQt2,3))\
                -s*pt2*std::pow(QQ2-u,2)/2./u/t/std::pow(QQ2-t,2)*(3.*std::pow(QQ2-t,3)*QQ2 \
                +(QQ2-t)*QQ2*(2.*u*t+q4)+std::pow(QQ2-t,2)\
                *(s2+4.*q2*QQ2-u*(QQ2+q2))-u2*QQ2*QQ2+u*t2*q2)\
                *(1.-2.*std::pow(QQ2/xQt2,2))\
                +s*pt2*(QQ2-u)/2./u/t/(QQ2-t)\
                *(3.*QQ2*s*(QQ2+s)*(QQ2-t)-t*s3+q2*QQ2*s2+QQ2*u*std::pow(q2-QQ2,2))\
                *(1.-2.*QQ2/xQt2) \
                +s*pt2/12./u/t*(-2.*s4+6.*s*q2*t*(t-q2)\
                +2.*s*q6+8.*QQ2*s*std::pow(s-QQ2,2)-2.*u*t*s*QQ2+7.*s2*q2*QQ2\
                -2.*s*std::pow(QQ2,2)*q2-q6*QQ2+3.*q2*std::pow(QQ2,3)-4.*u*t*q2*QQ2)\
                +11./6.*s3*pt2*std::pow(QQ2,2)/u/t-s2*pt2*q4*QQ2/3./t2;

    return B2pm+B2pp;
}
//==============================================================================================//



//==============================================================================================//
//                          Regular channel-specific coefficients                               //
//----------------------------------------------------------------------------------------------//
double HiggsDpT::Rgq(double pt, double u, double t, double s) {
    //////////////////////////////////////////////////////////////
    // gq->H+X terms in Eq. (A.17)                              //
    //////////////////////////////////////////////////////////////

    double pt2 = std::pow(pt,2);
    double q4 = std::pow(q2,2);
    double q6 = std::pow(q2,3);
    double QQ2 = u+t+s-q2;
    double x1 = t/(t-QQ2);
    double x2 = u/(u-QQ2);
    double xQt2 = QQ2+pt2;
    double t2 = std::pow(t,2);
    double u2 = std::pow(u,2);
    double u4 = std::pow(u,4);
    double s2 = std::pow(s,2);
    double s4 = std::pow(s,4);

    // Eqs. (A.2)-(A.3) G & S
    double A = s+q2-QQ2;
    double B = sqrt(std::pow(A,2)-4.*q2*s);
    double L1a = log(q2/s/std::pow(x1,2));
    double L1b = log(q2/s/std::pow(x2,2));
    double L2a = log(q2*s/std::pow(A-s*x1,2));
    double L2b = log(q2*s/std::pow(A-s*x2,2));
    double L3 = log((A+B)/(A-B));

    // Eqs. (A.18)-(A.26) G & S
    double C1pm = -2.*s2*pt2*QQ2*log(pt2/xQt2);
    double C2pm = 0.;
    double C1mp = s2*pt2*QQ2-1.5*s2*pt2*t2/u+pt2/2./xQt2*(s*std::pow(QQ2-t,3)\
                +QQ2*std::pow(QQ2-u,3))*(-3.+10.*QQ2/xQt2-6*std::pow(QQ2/xQt2,2));
    double C2mp = pt2*QQ2/std::pow(xQt2,2)*(s*std::pow(QQ2-t,3)+QQ2*std::pow(QQ2-u,3))*(-2.+3.*QQ2 \
                /xQt2)+2.*s2*pt2*QQ2+4.*s2*pt2*QQ2*log(pt2/xQt2);
    double C1pp = -1.5*s4*pt2/u-s*pt2*QQ2*A*A/u*L2b+s*pt2/u/t2*L1a\
                *(std::pow(A-q2,2)*s*t2-2.*QQ2*q2*u*t*A-QQ2*q4*(QQ2-t)*u)\
                +s*pt2/u/B*L3*((s+QQ2-q2)*(s*std::pow(A-q2,2)+QQ2*A*A)\
                -4.*s*QQ2*A*(A-q2))+0.5*s*pt2*(std::pow(QQ2-t,2)*(s/(QQ2-u)-QQ2/u)\
                +std::pow(QQ2-u,2)*(QQ2/(QQ2-t)+s/u))*(-3.+10.*QQ2/xQt2-6.*std::pow(QQ2/xQt2,2))\
                +s*pt2*(QQ2-t)/u/(QQ2-u)*(2*s*u*(s+t)\
                -QQ2*(4.*q2*QQ2-QQ2*t-q2*u-2*u*t))*(-1+2*QQ2/xQt2)\
                +s*pt2*(QQ2-u)/u/(QQ2-t)*(t*s2-2.*u*t*s+2.*QQ2*u*(QQ2-t))\
                *(-1.+2.*QQ2/xQt2)+s*pt2*QQ2*q2*(u+t)/t-2.*s*pt2*QQ2*QQ2*q4/t2\
                +s2*pt2*QQ2*q4/2./u2+s*pt2/2/u*(-2.*(QQ2+q2)*u*s+2.*q2*s2\
                +q4*s+q2*QQ2*(2.*(s-QQ2)+3.*q2-u)+5.*QQ2*s*(s-QQ2));
    double C2pp = 0.5*s*pt2*A*A*(1-x1)*(L2a-L1a)+s*pt2*(q2-t)*A*A\
                *(1.-x2)/2./u*(L1b-L2b)+0.5*s*pt2/u*(L1b-L1a)*std::pow(s-QQ2,3)\
                +s*pt2*QQ2*A*A/(q2-u)*(L1b-L2a)+s*pt2*QQ2/u2*L1b\
                *(2.*u*std::pow(s-QQ2,2)+4.*q2*(q2-t)*A-2.*q4*(QQ2-t)-q6)\
                -0.5*s*pt2*(std::pow(QQ2-t,2)*(s/(QQ2-u)-QQ2/u)\
                +std::pow(QQ2-u,2)*(QQ2/(QQ2-t)+s/u))*(-3.+10.*QQ2/xQt2-6.*std::pow(QQ2/xQt2,2))\
                +0.5*s*pt2*(QQ2-t)/u/(QQ2-u)*((-3.*s*pt2+u2-QQ2*QQ2)*(s+QQ2)\
                -q2*s*u+2.*QQ2*QQ2*(QQ2-u))*(-1.+2.*QQ2/xQt2)\
                +0.5*s*pt2*(QQ2-u)/u/(QQ2-t)*(3.*s*pt2*(s+QQ2)-u*QQ2*(QQ2-t)\
                +3.*QQ2*q2*(QQ2-u)+s*t*(u+s))*(-1.+2.*QQ2/xQt2)\
                +0.5*s*pt2*q2*QQ2/u2*(2.*std::pow(s-QQ2,2)-2.*q2*(s-q2)\
                -u*(QQ2-u)-4.*q2*QQ2)+s2*pt2*(u-t)*(q2+QQ2)/(2.*u)\
                -2.*std::pow(s*pt2*q2*QQ2,2)/u4*(4.+L1b);
    double C1mm = s2*pt2*t2/u*L1a-s*pt2*QQ2*std::pow(q2-t,2)/u*L2b+s*pt2*q2*QQ2/std::pow(B,2) \
                *(t*(u+t)-2.*q2*QQ2)+s*pt2/u/B*(t2*std::pow(B,2)-q2*t2*(u+t)+2.*std::pow(QQ2,2) \
                *q4+QQ2*q4*(3.*t-u)+QQ2*q4*u/std::pow(B,2)*(-t*(u+t)+2.*q2*QQ2+QQ2*(t-u)))*L3;
    double C2mm = s*pt2*t2*QQ2/(2.*u)*(L1a+3.*L1b)\
                +0.5*s*pt2*t2*(1-x1)*(L2a-L1a)+s*pt2*QQ2*std::pow(q2-t,3)*x2/(2.*u2)\
                *(L2b-L1b)+s*pt2*t2*QQ2/(q2-u)*(L1b+L2a)\
                +s2*pt2*t2/2./u*(L1b-L1a)+s*pt2*q2*QQ2/u2\
                *(4.*t*(t-q2)+q4)*L1b-2.*std::pow(s*pt2*q2*QQ2,2)/u4*(L1b+4.)\
                +s*pt2*t2*q2*QQ2/u2;
    double C2ep = 4./u4*std::pow(s*pt2*q2*QQ2,2);

    double result = 1./s2/pt2/QQ2*(16./9.*(C1pm+C1mp+C1pp+C1mm)+4.*(C2pm+C2mp+C2pp+C2mm+C2ep));

    return result;
}


double HiggsDpT::Rqqb(double pt, double u, double t, double s) {
    //////////////////////////////////////////////////////////////
    // qqbar->H+X terms in Eq. (A.27)                           //
    //////////////////////////////////////////////////////////////
    double x1, x2;
    double D1pm, D1pp, D2pm, D2pp, E1, E2, E3;
    double u2, t2, s2, s3, QQ2, q4;
    double xQt2;
    double A, B, L1a, L1b, L2a, L2b, L3;
    double pt2;
    double result;

    pt2 = std::pow(pt,2);
    q4 = std::pow(q2,2);
    QQ2 = u+t+s-q2;
    x1 = t/(t-QQ2);
    x2 = u/(u-QQ2);
    xQt2 = QQ2+pt2;
    t2 = std::pow(t,2);
    s2 = std::pow(s,2);
    s3 = std::pow(s,3);
    u2 = std::pow(u,2);

    // Eqs. (A.2)-(A.3) G & S
    A = s+q2-QQ2;
    B = sqrt(std::pow(A,2)-4.*q2*s);
    L1a = log(q2/s/std::pow(x1,2));
    L1b = log(q2/s/std::pow(x2,2));
    L2a = log(q2*s/std::pow(A-s*x1,2));
    L2b = log(q2*s/std::pow(A-s*x2,2));
    L3 = log((A+B)/(A-B));

    // Eqs. (A.28)-(A.34) G & S
    D1pm = -s2*pt2*QQ2*(1.+log(pt2/xQt2))-s3*std::pow(pt,2)*x1*(1.-x1)/t-s3*pt2*std::pow(1-x1,2);
    D2pm = -s2*pt2*QQ2/3.-s*pt2*t2/6./std::pow(x1,2)*(11.-12.*QQ2/xQt2+3.*std::pow(QQ2/xQt2,2)) \
            +11.*s*pt2*t2/6.;
    D1pp = s2*pt2*u2*(1-x2)*(L2b-L1b)/A+s2*pt2*std::pow(q2-u,2)*(1-x1)/A*(L2a-L1a)+s*pt2*q2*QQ2 \
           *(s*pt2+u*t)/t2*L1a-2.*s2*pt2*q4*QQ2/A/B*L3+s*pt2*q2*QQ2*(2.*s*pt2-u*t)/t2;

    D2pp = s*pt2*u2*(q2-t)*(1.-x2)/(2.*A)*(L1b-L2b)+s*pt2*std::pow(q2-u,3)\
            *(1-x1)/(2*A)*(L1a-L2a)-0.5*s*pt2*u2*(L1a+L1b)\
            +6.*std::pow(s*pt2*q2*QQ2,2)/std::pow(B,4)-s*pt2*q2*QQ2*u2/std::pow(B,2)\
            +L3*(s*pt2*u2*(u+t)/B+std::pow(s*pt2*q2*QQ2,2)*(1./A/std::pow(B,3)-3.*A/std::pow(B,5))\
            +s*pt2*q2*QQ2*((t-3.*u)/2./B+A*(std::pow(B,2)+2*u2)/4./std::pow(B,3) \
            +(t2-6.*u*t+7.*u2)/4./A/B));

    E1 = 4./3.*(2.*s2*pt2*QQ2-s*pt2*q2*QQ2);
    E2 = s*pt2*QQ2*(QQ2-pt2)/std::pow(xQt2,2)*(std::pow(QQ2-t,2)+(std::pow(QQ2-u,2)))+2.*s2*pt2*QQ2;
    E3 = -2.*s*pt2*(std::pow(u+t-2.*QQ2,2)-2.*s*pt2)*log(pt2/xQt2)-s*pt2*QQ2*(2.*xQt2+QQ2) \
        /std::pow(xQt2,2)*(std::pow(QQ2-t,2)+(std::pow(QQ2-u,2)))-6.*s2*pt2*QQ2;

    result  = (128./27.*(D1pm+D1pp)+32./3.*(D2pm+D2pp)+0.5*(16./9.*NF*E1+16./9.*E2+16./27.*E3)) \
            /s2/pt2/QQ2;

    return result;
}
//==============================================================================================//



//**********************************************************************************************//
// The below compute the delta and regular terms that enter into the hadronic cross section     //
// i.e. the partonic part convoluted with the PDFs and integrated over the momentunm fractions. //
// The integration w.r.t. the momentum fractions follow Section B of G & S.                     //
//----------------------------------------------------------------------------------------------//

//==============================================================================================//
//                                        Delta Terms                                           //
//----------------------------------------------------------------------------------------------//
double HiggsDpT::delta(
        double pt,
        double tm,
        double t,
        double y1,
        double y2,
        double yc
) {
    ///////////////////////////////////////////////////////////////
    // This function computes the terms proportional to delat(Q) //
    // and delta(Q^2) in the SINGULAR part of G&S. These are    //
    // given in Eqs. (3.17), (3.20), (3.24), (3.28).             //
    ///////////////////////////////////////////////////////////////
    // yc is the integration variable that gets passed to CUBA

    double yh = 0;
    double result = 0;
    double jac = 1;

    if (y1 != y2) {
        jac *= (y2-y1);
        yh = y1+(y2-y1)*yc;
    } else {
        yh = y1;
    }

    double x2min = (q2-SROOT*tm*exp(-yh))/(SROOT*tm*exp(yh)-shad);
    if (x2min < 0. || x2min > 1.) {
        std::cout << "ERROR IN x2min !" << std::endl;
        exit(EXIT_FAILURE);
    }

    double esp = 8;
    double x2 = exp((1.-std::pow(t,esp))*(log(x2min)));
    double x1 = (SROOT*tm*exp(yh)*x2-q2)/(x2*shad-SROOT*tm*exp(-yh));
    jac *= x2*log(x2min)*esp*std::pow(t,esp-1.);

    // Make sure x1&x2 are within the physical region
    double tiny = 1e-8;
    if (1.-tiny < x1 || 1.-tiny < x2) {
        return 0.;
    }

    double sh = x1*x2*shad;
    double th = q2-SROOT*x2*tm*exp(yh);
    double uh = q2-SROOT*x1*tm*exp(-yh);
    double jj = abs(x2*shad-SROOT*tm*exp(-yh));

    // LO
    if (ORD >= 0) {
        result += gg0(x1,x2,tm,yh)*Lumi->Lumgg(x1,x2) \
                +qg0(x1,x2,tm,yh)*Lumi->Lumqg(x1,x2) \
                +gq0(x1,x2,tm,yh)*Lumi->Lumgq(x1,x2) \
                +qqb0(x1,x2,tm,yh)*Lumi->Lumqqb(x1,x2);
    }

    // NLO
    if (ORD >= 1) {
        // gg-channel //
        double uu = 0.5*std::pow(log(uh/th),2)+std::pow(M_PI,2)/3. \
                    -log(sh/q2)*log(-th/q2)-log(sh/q2)*log(-uh/q2) \
                    -log(-uh/q2)*log(-th/q2)\
                    +std::pow(log(q2/sh),2)+std::pow(log(q2/(q2-th)),2)\
                    +std::pow(log(q2/(q2-uh)),2)+2.*Li2(1.-q2/sh)\
                    +2.*Li2(q2/(q2-th))+2.*Li2(q2/(q2-uh));
        double de = 1.5*beta0*(log(-MUR2/th)+log(-MUR2/uh))+67./6. \
                    -5./9.*NF;

        double ddgg = ((11.+de+3.*uu)*gg0(x1,x2,tm,yh)\
                    +(3.-NF)*(std::pow(q2,2)/sh+std::pow(q2,2)/th \
                    +std::pow(q2,2)/uh+q2))*Lumi->Lumgg(x1,x2);

        // gq-channel //
        // V functions Eqs. (3.22)-(3.24)
        double V1 = 0.5*(std::pow(log(uh/th),2)+std::pow(log(-sh/uh),2) \
                    -std::pow(log(-sh/th),2))+log(sh/q2)*log(-th/q2) \
                    -log(sh/q2)*log(-uh/q2)-log(-th/q2)*log(-uh/q2)+2. \
                    *Li2(q2/(q2-uh))+std::pow(log(q2/(q2-uh)),2) \
                    +std::pow(M_PI,2);
        double V2 = std::pow(log(q2/sh),2)+std::pow(log(q2/(q2-th)),2) \
                    -2.*log(sh/q2)*log(-th/q2)+2.*Li2(1.-q2/sh)+2. \
                    *Li2(q2/(q2-th))-3.5-2.*std::pow(M_PI,2)/3.;
        double V3 = beta0*(2*log(-MUR2/uh)+log(-MUR2/th))+67./3.-10.*NF/9.;

        double ddgq = ((11.+3*V1+4.*V2/3.+V3)*gq0(x1,x2,tm,yh) \
                    +20./9.*(std::pow(sh,2)+std::pow(th,2) \
                    +std::pow(uh,2)-uh*q2)/(-uh))*Lumi->Lumgq(x1,x2);

        // qg-channel //
        // V functions Eqs. (3.22)-(3.24)
        double V1c = 0.5*(std::pow(log(th/uh),2)+std::pow(log(-sh/th),2) \
                    -std::pow(log(-sh/uh),2))+log(sh/q2)*log(-uh/q2) \
                    -log(sh/q2)*log(-th/q2)-log(-uh/q2)*log(-th/q2)+2. \
                    *Li2(q2/(q2-th))+std::pow(log(q2/(q2-th)),2) \
                    +std::pow(M_PI,2);
        double V2c = std::pow(log(q2/sh),2)+std::pow(log(q2/(q2-uh)),2) \
                    -2*log(sh/q2)*log(-uh/q2)+2.*Li2(1.-q2/sh)+2. \
                    *Li2(q2/(q2-uh))-3.5-2.*std::pow(M_PI,2)/3.;
        double V3c = beta0*(2*log(-MUR2/th)+log(-MUR2/uh))+67./3.-10.*NF/9.;

        double ddqg = ((11.+3.*V1c+4.*V2c/3.+V3c)*qg0(x1,x2,tm,yh) \
                    +20./9.*(std::pow(sh,2)+std::pow(th,2)+std::pow(uh,2) \
                    -th*q2)/(-th))*Lumi->Lumqg(x1,x2);

        // qqbar-channel //
        // W functions Eqs. (3.26)-(3.28)
        double W1 = log(-uh/q2)*log(-th/q2)-log(sh/q2)*log(-uh/q2)-log(sh/q2) \
                    *log(-th/q2)+2.*Li2(1.-q2/sh)+std::pow(log(q2/sh),2)-0.5 \
                    *std::pow(log(uh/th),2)-5.*std::pow(M_PI,2)/3.;
        double W2 = 1.5*log(std::pow(sh,2)/th/uh)+std::pow(log(uh/th),2) \
                    -2.*log(-uh/q2)*log(-th/q2)+std::pow(log(q2/(q2-uh)),2) \
                    +std::pow(log(q2/(q2-th)),2)+2.*Li2(q2/(q2-uh))+2. \
                    *Li2(q2/(q2-th))-7.+2.*std::pow(M_PI,2);
        double W3 = beta0/2*(4*log(MUR2/sh)+log(-MUR2/uh)+log(-MUR2/th)) \
                    +(67./2.-5.*NF/3.);
            // NOTE: are we missing the lumi here?
        double ddqqb = ((11.+3.*W1+4.*W2/3.+W3)*qqb0(x1,x2,tm,yh)+160./27. \
                *(std::pow(th,2)+std::pow(uh,2)+std::pow(sh,2)-sh*q2)/sh);

        // Full NLO results
        result += aass/2.*(ddgg+ddgq+ddqg+ddqqb);
    }

    result *= -2.*pt*jac/sh/jj;

    return result;
}



//==============================================================================================//
//                                    Non-Delta Terms                                           //
//----------------------------------------------------------------------------------------------//
double HiggsDpT::distr(
        double pt,
        double tm,
        double y1,
        double y2,
        double zz1,
        double zz2,
        double zz3
) {
    ////////////////////////////////////////////////////////////////
    // This function computes the remaining terms that are not    //
    // multiplied by either delta(Q) or delta(Q^2) (both singular //
    // and non-Singular)                                          //
    // zz1 == z1                                                  //
    // zz2 == z2                                                  //
    // zz3 == yh                                                  //
    ////////////////////////////////////////////////////////////////

    double yh;
    double jacy = 1;

    if (y1 != y2) {
        jacy *= (y2-y1);
        yh = y1+(y2-y1)*zz3;
    } else {
        yh = y1;
    }

    double x10 = tm/SROOT*exp(yh);
    double x20 = tm/SROOT*exp(-yh);
    double dcut = x10*std::pow(pt/tm,2)/(1.-x10*(1.-std::pow(pt/tm,2)));

    double z1 = x20+(1.-dcut-x20)*zz1;
    double jac = 1.-dcut-x20;
    double lb = std::pow(pt,2)/std::pow(tm,2)*z1/(1.-z1);
    double z2 = x10*(1.+lb)+(1.-x10*(1.+lb))*zz2;
    jac *= 1.-x10*(1.+lb);
    double x1 = x10/z2*(1.+lb);
    double x2 = x20/z1;

    double sh = std::pow(tm,2)/z1/z2*(1.+lb);
    double th = -std::pow(tm,2)/z1*(1.-z1)*(1.+lb);
    double uh = q2-std::pow(tm,2)/z2*(1.+lb);

    double QQ2 = std::pow(tm,2)/z1/z2*(1.-z2)*(1.-z1)*(1.+lb);
    double pt2 = std::pow(pt,2);
    double xQt2 = QQ2+pt2;

    double pre1 = std::pow(tm,2)*(1.+lb)/std::pow(z1*z2,2);

    // a1:: (log(1-z1)/(1-z1))+ terms
    // b1:: 1/(1-z1))+ terms
    // c1:: regular terms

    // gg-channel
    coeff(pt,uh,th,sh);
    double a1 = (1./(-th)*pgg(z2)*gg0(z2*x1,x2,tm,yh)+(-z2/th)*big1) \
                *Lumi->Lumgg(x1,x2)/sh;
    double b1 = (1./th*pgg(z2)*log(-MUF2*z2/th)*gg0(z2*x1,x2,tm,yh) \
                +z2/th*big1*log((QQ2+pt*pt)*z2/(-th))+z2/th*big2) \
                *Lumi->Lumgg(x1,x2)/sh;
    double c1 = (1/(-th)*(-2*NF*pqg(z2)*log(MUF2/QQ2)+2*NF*z2*(1-z2)) \
                *qg0(z2*x1,x2,tm,yh)+0.5*big3)*Lumi->Lumgg(x1,x2)/sh;

    // gq-channel
    a1 += (1./(-th)*pgg(z2)*gq0(z2*x1,x2,tm,yh)+(-z2/th)*big4) \
            *Lumi->Lumgq(x1,x2)/sh;
    b1 += (1./th*pgg(z2)*log(-MUF2*z2/th)*gq0(z2*x1,x2,tm,yh) \
            +z2/th*big4*log((QQ2+pt*pt)*z2/(-th)))*Lumi->Lumgq(x1,x2)/sh;
    c1 += (1/(-th)*(-pqg(z2)*log(MUF2/QQ2)+z2*(1-z2))\
            *qqb0(z2*x1,x2,tm,yh)+0.5*big5)*Lumi->Lumgq(x1,x2)/sh;

    // qg-channel
    coeff(pt,th,uh,sh);
    a1 += -1./th*pqq(z2)*qg0(z2*x1,x2,tm,yh)*Lumi->Lumqg(x1,x2)/sh;
    b1 += (1./th*pqq(z2)*log(-MUF2*z2/th)*qg0(z2*x1,x2,tm,yh) \
        +z2/th*8./3.*(std::pow(uh,2)+std::pow(sh,2))/(-th)) \
        *Lumi->Lumqg(x1,x2)/sh;
    c1 += (-1/th*(4./3.*(1-z2)*qg0(z2*x1,x2,tm,yh)\
        +(-pgq(z2)*log(MUF2/QQ2)+4/3*z2)*gg0(z2*x1,x2,tm,yh))\
        +0.5*big5)*Lumi->Lumqg(x1,x2)/sh;

    // qqb-channel
    a1 += 1./(-th)*(pqq(z2)*qqb0(z2*x1,x2,tm,yh)-z2*16./27.\
        *(std::pow(th,2)+std::pow(uh,2)+std::pow(QQ2-th,2) \
        +std::pow(QQ2-uh,2))/sh)*Lumi->Lumqqb(x1,x2)/sh;
    b1 += (1./th*pqq(z2)*log(-MUF2*z2/th)*qqb0(z2*x1,x2,tm,yh) \
        -z2/th*log((QQ2+std::pow(pt,2))*z2/(-th))*16./27. \
        *(std::pow(th,2)+std::pow(uh,2)+std::pow(QQ2-th,2) \
        +std::pow(QQ2-uh,2))/sh+z2/th*16./9. \
        *beta0*(std::pow(uh,2)+std::pow(th,2))/sh) \
        *Lumi->Lumqqb(x1,x2)/sh;
    c1 += (-1/th*(4./3*(1-z2)*qqb0(z2*x1,x2,tm,yh)+(-pgq(z2) \
        *log(MUF2/QQ2)+4./3*z2)*gq0(z2*x1,x2,tm,yh))+16./9 \
        *(std::pow(sh-QQ2,2)+std::pow(uh+th-2*QQ2,2))/sh \
        *log(std::pow(pt,2)/(std::pow(pt,2)+QQ2))/std::pow(pt,2)) \
        *Lumi->Lumqqb(x1,x2)/sh;

    // qQ, qQb and qq - channel (same flavours)
    c1 += (-1./th*(-pgq(z2)*log(MUF2/QQ2)+4./3*z2)*gq0(z2*x1,x2,tm,yh)\
        +16./9*(std::pow(sh-QQ2,2)+std::pow(uh+th-2*QQ2,2))/sh \
        *log(std::pow(pt,2)/(std::pow(pt,2)+QQ2))/std::pow(pt,2)) \
        *(Lumi->Lumqq(x1,x2)+Lumi->LumqQb(x1,x2)+Lumi->LumqQ(x1,x2))/sh;


    // Non-Singular terms //
    REG(pt,uh,th,sh);
    double nosingular = 0.5*(0.+REGgg*Lumi->Lumgg(x1,x2) \
                +REGgq*Lumi->Lumgq(x1,x2) \
                +REGqg*Lumi->Lumqg(x1,x2) \
                +REGqqb*Lumi->Lumqqb(x1,x2) \
                +REGqqpb*(Lumi->LumqQ(x1,x2) \
                +Lumi->LumqQb(x1,x2)) \
                +REGqq*Lumi->Lumqq(x1,x2))/sh;

    x1 = x10*(1.+lb);
    x2 = x20/z1;

    sh = std::pow(tm,2)/z1*(1.+lb);
    th = -std::pow(tm,2)/z1*(1.-z1)*(1.+lb);
    uh = q2-std::pow(tm,2)*(1.+lb);
    QQ2 = 0;
    xQt2 = QQ2+pt2;

    double pre10 = std::pow(tm,2)*(1.+lb)/std::pow(z1,2);

    // gg-channel
    coeff(pt,uh,th,sh);
    double a10 = (1/(-th)*pgg(1.)*gg0(x1,x2,tm,yh)+(-1./th)*big1) \
                *Lumi->Lumgg(x1,x2)/sh;
    double b10 = (6./th*log(-MUF2/th)*gg0(x1,x2,tm,yh)+big1/th*log((xQt2) \
                /(-th))+big2/th)*Lumi->Lumgg(x1,x2)/sh;
    double d10 = 1/th*beta0*log(MUF2/(-th))*gg0(x1,x2,tm,yh) \
                *Lumi->Lumgg(x1,x2)/sh;

    // gq-channel
    a10 += (1/(-th)*pgg(1.)*gq0(x1,x2,tm,yh)+(-1./th)*big4) \
           *Lumi->Lumgq(x1,x2)/sh;
    b10 += (1./th*pgg(1.)*log(-MUF2/th)*gq0(x1,x2,tm,yh)+1./th*big4 \
           *log((xQt2)/(-th)))*Lumi->Lumgq(x1,x2)/sh;
    d10 += 1./th*beta0*log(MUF2/(-th))*gq0(x1,x2,tm,yh) \
           *Lumi->Lumgq(x1,x2)/sh;

    // qg-channel
    coeff(pt,th,uh,sh);
    a10 += -1./th*pqq(1.)*qg0(x1,x2,tm,yh)*Lumi->Lumqg(x1,x2)/sh;
    b10 += (1./th*pqq(1.)*log(-MUF2/th)*qg0(x1,x2,tm,yh)+1./th*8./3. \
           *(std::pow(uh,2)+std::pow(sh,2))/(-th))*Lumi->Lumqg(x1,x2)/sh;
    d10 += 1./th*2.*log(MUF2/(-th))*qg0(x1,x2,tm,yh)*Lumi->Lumqg(x1,x2)/sh;

    // qqb-channel
    a10 += 1./(-th)*(pqq(1)*qqb0(x1,x2,tm,yh)-16./27.*(std::pow(th,2) \
        +std::pow(uh,2)+std::pow(QQ2-th,2)+std::pow(QQ2-uh,2))/sh) \
        *Lumi->Lumqqb(x1,x2)/sh;
    b10 += (1./th*pqq(1.)*log(-MUF2/th)*qqb0(x1,x2,tm,yh)-1./th*log((xQt2) \
        /(-th))*16./27.*(std::pow(th,2)+std::pow(uh,2)+std::pow(QQ2-th,2) \
        +std::pow(QQ2-uh,2))/sh+1./th*16./9.*beta0*(std::pow(uh,2) \
        +std::pow(th,2))/sh)*Lumi->Lumqqb(x1,x2)/sh;
    d10 += 1./th*2.*log(MUF2/(-th))*qqb0(x1,x2,tm,yh)*Lumi->Lumqqb(x1,x2)/sh;

    // Fix mismatch in plus distribution
    double z2min = x10*(1+lb);
    double m10 = -0.5*a10*std::pow(log(1-z2min),2)-log(1-z2min)*b10;
    m10 *= pre10*(1.-dcut-x20);


    // COMBINED INTEGRALS //
    a1 = log(1-z2)/(1-z2)*(pre1*a1-pre10*a10)*jac;
    b1 = (pre1*b1-pre10*b10)/(1-z2)*jac;
    c1 *= jac*pre1;
    d10 *= pre10*(1-dcut-x20);
    nosingular *= jac*pre1;

    double result = (a1+b1+c1+d10-m10+nosingular)/shad;
    result *= 2.*pt*jacy;

    return result;
}
//==============================================================================================//



//==============================================================================================//
//                                    Non-Delta Terms                                           //
//----------------------------------------------------------------------------------------------//
double HiggsDpT::distrcross(
        double pt,
        double tm,
        double y1,
        double y2,
        double zz1,
        double zz2,
        double zz3
) {
    ////////////////////////////////////////////////////////////////
    // This function computes the remaining terms that are not    //
    // multiplied by either delta(Q^2) (both singular and         //
    // non-Singular)                                              //
    ////////////////////////////////////////////////////////////////

    double yh;
    double jacy = 1;

    if (y1 != y2) {
        jacy *= (y2-y1);
        yh = y1+(y2-y1)*zz3;
    } else {
        yh = y1;
    }

    double x10 = tm/SROOT*exp(yh);
    double x20 = tm/SROOT*exp(-yh);
    double d = x20*std::pow(pt/tm,2)/(1.-x20*(1.-std::pow(pt/tm,2)));
    double z1 = x10+(1.-d-x10)*zz1;
    double jac = (1.-d-x10);
    double la = pt*pt/std::pow(tm,2)*z1/(1.-z1);
    double z2 = x20*(1.+la)+(1.-x20*(1.+la))*zz2;
    jac = jac*(1.-x20*(1.+la));

    double x1 = x10/z1;
    double x2 = x20/z2*(1.+la);
    double sh = std::pow(tm,2)/z1/z2*(1+la);
    double uh = -std::pow(tm,2)/z1*(1-z1)*(1+la);
    double th = q2-std::pow(tm,2)/z2*(1+la);
    double QQ2 = std::pow(tm,2)/z1/z2*(1-z2)*(1-z1)*(1+la);
    double pre1 = std::pow(tm,2)*(1.+la)/std::pow(z1*z2,2);

    // a1:: (log(1-z1)/(1-z1))+ terms
    // b1:: 1/(1-z1))+ terms
    // c1:: regular terms

    // gg-channel
    coeff(pt,uh,th,sh);
    double a1 = (1./(-uh)*pgg(z2)*gg0(x1,x2*z2,tm,yh)+(-z2/uh)*big1) \
                *Lumi->Lumgg(x1,x2)/sh;
    double b1 = (1./uh*pgg(z2)*log(-MUF2*z2/uh)*gg0(x1,x2*z2,tm,yh) \
                +z2/uh*big1*log((QQ2+pt*pt)*z2/(-uh))+z2/uh*big2) \
                *Lumi->Lumgg(x1,x2)/sh;
    double c1 = (1./(-uh)*(-2.*NF*pqg(z2)*log(MUF2/QQ2)+2.*NF*z2*(1.-z2)) \
                *gq0(x1,x2*z2,tm,yh)+0.5*big3)*Lumi->Lumgg(x1,x2)/sh;


    // gq-channel
    a1 += -1./uh*pqq(z2)*gq0(x1,z2*x2,tm,yh)*Lumi->Lumgq(x1,x2)/sh;
    b1 += (1./uh*pqq(z2)*log(-MUF2*z2/uh)*gq0(x1,z2*x2,tm,yh) \
        +z2/uh*8./3*(std::pow(th,2)+std::pow(sh,2))/(-uh)) \
        *Lumi->Lumgq(x1,x2)/sh;
    c1 += (-1./uh*(4./3*(1-z2)*gq0(x1,z2*x2,tm,yh) \
        +(-pgq(z2)*log(MUF2/QQ2)+4./3*z2)*gg0(x1,z2*x2,tm,yh)) \
        +0.5*big5)*Lumi->Lumgq(x1,x2)/sh;

    // qg-channel
    coeff(pt,th,uh,sh);
    a1 += (1./(-uh)*pgg(z2)*qg0(x1,z2*x2,tm,yh)+(-z2/uh)*big4) \
        *Lumi->Lumqg(x1,x2)/sh;
    b1 += (1./uh*pgg(z2)*log(-MUF2*z2/uh)*qg0(x1,z2*x2,tm,yh) \
        +z2/uh*big4*log((QQ2+std::pow(pt,2))*z2/(-uh))) \
        *Lumi->Lumqg(x1,x2)/sh;
    c1 += (1./(-uh)*(-pqg(z2)*log(MUF2/QQ2)+z2*(1-z2)) \
        *qqb0(x1,z2*x2,tm,yh)+0.5*big5)*Lumi->Lumqg(x1,x2)/sh;

    // qqb-channel
    a1 += 1./(-uh)*(pqq(z2)*qqb0(x1,z2*x2,tm,yh)-z2*16./27 \
        *(std::pow(th,2)+std::pow(uh,2)+std::pow(QQ2-th,2) \
        +std::pow(QQ2-uh,2))/sh)*Lumi->Lumqqb(x1,x2)/sh;
    b1 += (1./uh*pqq(z2)*log(-MUF2*z2/uh)*qqb0(x1,z2*x2,tm,yh) \
        -z2/uh*log((QQ2+std::pow(pt,2))*z2/(-uh))*16./27 \
        *(std::pow(th,2)+std::pow(uh,2)+std::pow(QQ2-th,2) \
        +std::pow(QQ2-uh,2))/sh+z2/uh*16./9*beta0*(std::pow(uh,2) \
        +std::pow(th,2))/sh)*Lumi->Lumqqb(x1,x2)/sh;
    c1 += (-1./uh*(4./3*(1-z2)*qqb0(x1,z2*x2,tm,yh)+(-pgq(z2) \
        *log(MUF2/QQ2)+4./3*z2)*qg0(x1,z2*x2,tm,yh))+16./9 \
        *(std::pow(sh-QQ2,2)+std::pow(uh+th-2*QQ2,2))/sh \
        *log(std::pow(pt,2)/(std::pow(pt,2)+QQ2))/std::pow(pt,2)) \
        *Lumi->Lumqqb(x1,x2)/sh;

    // qQ-channel (same flavours)
    c1 += (-1./uh*(-pgq(z2)*log(MUF2/QQ2)+4./3*z2)*qg0(x1,z2*x2,tm,yh) \
        +16./9*(std::pow(sh-QQ2,2)+std::pow(uh+th-2*QQ2,2))/sh \
        *log(std::pow(pt,2)/(std::pow(pt,2)+QQ2))/std::pow(pt,2)) \
        *(Lumi->Lumqq(x1,x2)+Lumi->LumqQb(x1,x2)+Lumi->LumqQ(x1,x2))/sh;

    // Non-Singular terms //
    REG(pt,uh,th,sh);
    double nosingular = 0.5*(REGgg*Lumi->Lumgg(x1,x2) \
                    +REGgq*Lumi->Lumgq(x1,x2) \
                    +REGqg*Lumi->Lumqg(x1,x2) \
                    +REGqq*Lumi->Lumqq(x1,x2) \
                    +REGqqb*Lumi->Lumqqb(x1,x2) \
                    +REGqqpb*(Lumi->LumqQ(x1,x2) \
                    +Lumi->LumqQb(x1,x2)))/sh;

    x1 = x10/z1;
    x2 = x20*(1.+la);
    sh = std::pow(tm,2)/z1*(1.+la);
    uh = -std::pow(tm,2)/z1*(1.-z1)*(1.+la);
    th = q2-std::pow(tm,2)*(1.+la);
    QQ2 = 0.;
    double pre10 = std::pow(tm,2)*(1.+la)/std::pow(z1,2);

    // gg-channel
    coeff(pt,uh,th,sh);
    double a10 = (1./(-uh)*pgg(1.)*gg0(x1,x2,tm,yh)+(-1/uh)*big1) \
                *Lumi->Lumgg(x1,x2)/sh;
    double b10 = (1./uh*pgg(1.)*log(-MUF2/uh)*gg0(x1,x2,tm,yh) \
                +1./uh*big1*log((QQ2+std::pow(pt,2))/(-uh))+1/uh*big2) \
                *Lumi->Lumgg(x1,x2)/sh;
    double d10 = 1./uh*beta0*log(MUF2/(-uh))*gg0(x1,x2,tm,yh) \
                *Lumi->Lumgg(x1,x2)/sh;

    // gq-channel
    a10 += -1./uh*pqq(1.)*gq0(x1,x2,tm,yh)*Lumi->Lumgq(x1,x2)/sh;
    b10 += (1./uh*pqq(1.)*log(-MUF2/uh)*gq0(x1,x2,tm,yh) \
            +1./uh*8./3*(std::pow(th,2)+std::pow(sh,2))/(-uh)) \
            *Lumi->Lumgq(x1,x2)/sh;
    d10 += 1/uh*2*log(MUF2/(-uh))*gq0(x1,x2,tm,yh) \
            *Lumi->Lumgq(x1,x2)/sh;


    // qg-channel
    coeff(pt,th,uh,sh);
    a10 += (1./(-uh)*pgg(1.)*qg0(x1,x2,tm,yh)+(-1./uh)*big4) \
            *Lumi->Lumqg(x1,x2)/sh;
    b10 += (1./uh*pgg(1.)*log(-MUF2/uh)*qg0(x1,x2,tm,yh) \
            +1./uh*big4*log((QQ2+std::pow(pt,2))/(-uh))) \
            *Lumi->Lumqg(x1,x2)/sh;
    d10 += 1./uh*beta0*log(MUF2/(-uh))*qg0(x1,x2,tm,yh)\
            *Lumi->Lumqg(x1,x2)/sh;

    // qqb-channel
    a10 += 1./(-uh)*(pqq(1.)*qqb0(x1,x2,tm,yh)-16./27 \
            *(std::pow(th,2)+std::pow(uh,2)+std::pow(QQ2-th,2) \
            +std::pow(QQ2-uh,2))/sh)*Lumi->Lumqqb(x1,x2)/sh;
    b10 += (1./uh*pqq(1.)*log(-MUF2/uh)*qqb0(x1,x2,tm,yh) \
            -1./uh*log((QQ2+std::pow(pt,2))/(-uh))*16./27 \
            *(std::pow(th,2)+std::pow(uh,2)+std::pow(QQ2-th,2) \
            +std::pow(QQ2-uh,2))/sh+1./uh*16./9*beta0*(std::pow(uh,2) \
            +std::pow(th,2))/sh)*Lumi->Lumqqb(x1,x2)/sh;
    d10 += 1./uh*2*log(MUF2/(-uh))*qqb0(x1,x2,tm,yh) \
            *Lumi->Lumqqb(x1,x2)/sh;

    // Fix mismatch in plus distribution
    double z2min = x20*(1.+la);
    double m10 = -0.5*std::pow(log(1.-z2min),2)*a10-log(1.-z2min)*b10;
    m10 *= pre10*(1.-d-x10);

    // COMBINED INTEGRALS //
    a1 = log(1.-z2)/(1.-z2)*(pre1*a1-pre10*a10)*jac;
    b1 = (pre1*b1-pre10*b10)/(1.-z2)*jac;
    c1 *= pre1*jac;
    d10 = pre10*d10*(1.-d-x10);
    nosingular *= pre1*jac;

    double result = (a1+b1+c1+d10-m10+nosingular)/shad;
    result *= 2.*pt*jacy;

    return result;
}
//**********************************************************************************************//
