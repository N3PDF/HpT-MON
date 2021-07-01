/*
 * =====================================================================================
 *
 *       Filename:  higgsptpartonic.cpp
 *
 *    Description:  This file contains all the terms contribuing to the partonic
 * Higgs fixed order results.
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

#include "./higgsptpartonic.h"

HiggsDpTpartonic::HiggsDpTpartonic(int order, int channel, std::string pdfname,
                                   void *params) {
  PhysParams param = *reinterpret_cast<PhysParams *>(params);

  NF = param.nf;
  MH2 = std::pow(param.mh, 2);
  MUR2 = std::pow(param.mur, 2);
  MUF2 = std::pow(param.muf, 2);
  aass = param.alphas;
  SIGMA0 = param.sigma0;

  ORD = order;
  CHANNEL = channel;

  beta0 = (33 - 2 * NF) / 6.;
}

HiggsDpTpartonic::~HiggsDpTpartonic() {}

//==============================================================================================//
//                          Partonic Functions (Eqs. [2.9]) //
//----------------------------------------------------------------------------------------------//
double HiggsDpTpartonic::gg0(double sh, double th, double uh, double q2) {
  return 3 *
         (std::pow(q2, 4) + std::pow(sh, 4) + std::pow(th, 4) +
          std::pow(uh, 4)) /
         (uh * th * sh);
}

double HiggsDpTpartonic::qg0(double sh, double th, double uh) {
  return 4. / 3. * (std::pow(uh, 2) + std::pow(sh, 2)) / (-th);
}

double HiggsDpTpartonic::gq0(double sh, double th, double uh) {
  return 4. / 3. * (std::pow(th, 2) + std::pow(sh, 2)) / (-uh);
}

double HiggsDpTpartonic::qqb0(double sh, double th, double uh) {
  return 2. * 16. / 9. * (std::pow(th, 2) + std::pow(uh, 2)) / sh;
}

//==============================================================================================//

//==============================================================================================//
//                        Regular Part of Coefficient Functions //
//----------------------------------------------------------------------------------------------//
void HiggsDpTpartonic::REG(double pt, double uh, double th, double sh,
                           double q2) {
  ///////////////////////////////////////////////////////////
  // Below are terms that contribute to the NonSingular    //
  // Real Contribution in Appendix A of G & S. Secifically //
  // these are given by Eqs. (A.1),(A.17),(A.27),(A.35)    //
  ///////////////////////////////////////////////////////////
  double pt2 = std::pow(pt, 2);
  double QQ2 = uh + th + sh - q2;
  double xa = th / (th - QQ2);
  double xb = uh / (uh - QQ2);
  double xQt2 = QQ2 + pt2;

  double E2 = sh * pt2 * QQ2 * (QQ2 - pt2) / std::pow(pt2 + QQ2, 2) *
                  (std::pow(QQ2 - th, 2) + std::pow(QQ2 - uh, 2)) +
              2. * std::pow(sh, 2) * pt2 * QQ2;
  double E4 = 2. * sh * QQ2 * pt2 * (std::pow(sh, 2) + std::pow(QQ2, 2)) /
                  (pt2 + QQ2) * log(pt2 / QQ2) +
              4. * std::pow(sh, 2) * pt2 * QQ2 * log(pt2 / (pt2 + QQ2));

  // eq. (A.4)
  double A0 = (std::pow(th / xa, 4) + std::pow(uh / xb, 4)) * pt2 * QQ2 /
                  std::pow(xQt2, 2) *
                  (5. - 7. * QQ2 / xQt2 + 20. / 3. * std::pow(QQ2 / xQt2, 2)) +
              std::pow(sh, 2) * QQ2 * pt2 * (17. / 3. + 4. * log(pt2 / xQt2));
  double Aep = 4. * std::pow(sh * pt2 * q2 * QQ2, 2) *
               (1. / std::pow(th, 4) + 1. / std::pow(uh, 4));

  // Regular part for specific channel
  REGgg =
      1. / std::pow(sh, 2) / pt2 / QQ2 *
      (9. * (A0 + Aep + A1234(pt, uh, th, sh, q2) + A1234(pt, th, uh, sh, q2) +
             A3412(pt, uh, th, sh, q2) + A3412(pt, th, uh, sh, q2) +
             A1324(pt, uh, th, sh, q2) + A1324(pt, th, uh, sh, q2) +
             A3241(pt, uh, th, sh, q2) + A3241(pt, th, uh, sh, q2)) +
       4. / 3. * NF *
           (B1pmB1pp(pt, uh, th, sh, q2) + B1pmB1pp(pt, th, uh, sh, q2)) +
       3. * NF * (B2pmB2pp(pt, uh, th, sh, q2) + B2pmB2pp(pt, th, uh, sh, q2)));
  REGgq = Rgq(pt, uh, th, sh, q2);
  REGqg = Rgq(pt, th, uh, sh, q2);
  REGqqb = Rqqb(pt, uh, th, sh, q2) + Rqqb(pt, th, uh, sh, q2);
  REGqqpb = 16. / 9. / std::pow(sh, 2) / pt2 / QQ2 * E2;
  REGqq = 16. / 9. / std::pow(sh, 2) / pt2 / QQ2 * (E2 + E4 / 3.);
}

//==============================================================================================//

//==============================================================================================//
//                                      Coefficients //
//----------------------------------------------------------------------------------------------//
void HiggsDpTpartonic::coeff(double pt, double uh, double th, double sh,
                             double q2) {
  double tiny = 1e-4;
  double QQ2 = uh + sh + th - q2;
  double xa = th / (th - QQ2);
  double xb = uh / (uh - QQ2);
  double tmp = tiny * std::pow(pt, 2);

  ////// big1 & big2 & big3 are the terms that contribute to the gg///////
  // Term that gets multiplied by NC^2/2 in Eq. (3,17) (inside the curly
  // brackets)
  big1 = 0.5 * 9. *
         ((std::pow(q2, 4) + std::pow(sh, 4) + std::pow(QQ2, 4) +
           std::pow(uh, 4) + std::pow(th, 4)) +
          xa * xb *
              (std::pow(q2, 4) + std::pow(sh, 4) + std::pow(QQ2, 4) +
               std::pow(uh / xb, 4) + std::pow(th / xa, 4))) /
         sh / uh / th;

  // Term that gets multiplied by Beat0*NC in Eq. (3,17) (inside the curly
  // brackets)
  big2 = beta0 * 3. / 2. *
         (std::pow(q2, 4) + std::pow(sh, 4) +
          xa * xb * (std::pow(uh / xb, 4) + std::pow(th / xa, 4))) /
         sh / uh / th;

  // Terms that get multiplied by NC^2 in the square brackets
  if (QQ2 > tmp) {
    big3 = 9. *
           ((std::pow(q2, 4) + std::pow(sh, 4) + std::pow(QQ2, 4) +
             std::pow(uh / xb, 4) + std::pow(th / xa, 4)) *
                (2 * QQ2 + std::pow(pt, 2)) / std::pow(sh, 2) / QQ2 /
                (QQ2 + std::pow(pt, 2)) +
            2 * std::pow(q2, 2) *
                (std::pow(q2 - th, 4) + std::pow(q2 - uh, 4) + std::pow(uh, 4) +
                 std::pow(th, 4)) /
                sh / uh / th / (q2 - uh) / (q2 - th)) /
           std::pow(pt, 2) * log(std::pow(pt, 2) / (std::pow(pt, 2) + QQ2));
  } else {
    big3 = -9. * ((std::pow(q2, 4) + std::pow(sh, 4) + std::pow(uh, 4) +
                   std::pow(th, 4)) /
                  std::pow(sh, 2) / std::pow(pt, 4));
  }

  ////// big4 & big5 are the terms that contribute to the qg///////
  // Term that gets multiplied by the first NC*CF in Eq. (3.20) (inside square
  // bracket)
  big4 = 4. * ((-std::pow(sh, 3) * th - sh * std::pow(th, 3) +
                std::pow(QQ2, 3) * th + QQ2 * std::pow(th, 3)) /
                   sh / uh / th +
               xa * xb *
                   (-std::pow(sh, 3) * th / xa - sh * std::pow(th / xa, 3) -
                    std::pow(QQ2, 3) * uh / xb - QQ2 * std::pow(uh / xb, 3)) /
                   (sh * uh * th));

  // Term that gets multiplied by the seond NC*CF in Eq. (3.20) (inside square
  // bracket)
  if (QQ2 > tmp) {
    big5 =
        4. *
        ((-std::pow(sh, 3) * th / xa - sh * std::pow(th / xa, 3) -
          std::pow(QQ2, 3) * uh / xb - QQ2 * std::pow(uh / xb, 3)) *
             (2. * QQ2 + std::pow(pt, 2)) / std::pow(sh, 2) /
             (QQ2 + std::pow(pt, 2)) -
         2. * QQ2 * std::pow(q2, 2) * (std::pow(q2 - th, 2) + std::pow(th, 2)) /
             sh / uh / (q2 - uh)) *
        log(std::pow(pt, 2) / (std::pow(pt, 2) + QQ2)) /
        (QQ2 * std::pow(pt, 2));
  } else {
    big5 = -4. * (-std::pow(sh, 3) * th / xa - sh * std::pow(th / xa, 3)) /
           std::pow(sh, 2) / std::pow(pt, 4);
  }
}

//==============================================================================================//
//                                      A Functions //
//----------------------------------------------------------------------------------------------//
double HiggsDpTpartonic::A1234(double pt, double u, double t, double s,
                               double q2) {
  double pt2 = std::pow(pt, 2);
  double q4 = std::pow(q2, 2);
  double QQ2 = u + t + s - q2;
  double x1 = t / (t - QQ2);
  double x2 = u / (u - QQ2);
  double t2 = std::pow(t, 2);
  double t4 = std::pow(t, 4);
  double u2 = std::pow(u, 2);
  double u3 = std::pow(u, 3);
  double u4 = std::pow(u, 4);

  // Eqs. (A.2)-(A.3) G & S
  double A = s + q2 - QQ2;
  double B = sqrt(std::pow(A, 2) - 4. * q2 * s);
  double L1a = log(q2 / s / std::pow(x1, 2));
  double L1b = log(q2 / s / std::pow(x2, 2));
  double L2a = log(q2 * s / std::pow(A - s * x1, 2));
  double L2b = log(q2 * s / std::pow(A - s * x2, 2));
  double L3 = log((A + B) / (A - B));

  double result =
      -0.5 * ((std::pow(s * pt2 / t, 4) + std::pow(q2 * QQ2 / t, 4)) * L1a +
              std::pow(u, 4) * L1b) +
      0.5 * s * q4 * QQ2 * u3 / t / (u - q2) / (t - q2) * (L1a + L1b) +
      0.5 * s * pt2 * QQ2 * u3 / A / (t - q2) * (L2b - L1b) +
      0.5 * s * pt2 * QQ2 * (std::pow(q4, 2) + std::pow(u - q2, 4)) / A / t /
          (u - q2) * (L2a - L1a) +
      s * pt2 * q2 * QQ2 * (u4 / std::pow(B * t, 2) + u2 / std::pow(B, 2)) /
          2. +
      std::pow(s * pt2 * q2 * QQ2, 2) *
          (-6. / std::pow(B, 4) - 4. / t4 + 8. / std::pow(B * t, 2)) +
      L3 * (s * pt2 * u3 * (u + t) / B / t +
            std::pow(s * pt2 * q2 * QQ2, 2) *
                (3. * A / std::pow(B, 5) - 1. / A / std::pow(B, 3)) -
            s * pt2 * q2 * QQ2 *
                ((t2 + t * u + 4 * std::pow(u, 2) - 2. * q2 * QQ2) / B / t +
                 A * (t2 + 3. * t * u + 3. * std::pow(u, 2) - 6. * q2 * QQ2) /
                     2. / std::pow(B, 3) +
                 (t2 + t * u + 7. * u2 - 2. * q2 * QQ2) / (2 * A * B)));

  return result;
}

double HiggsDpTpartonic::A3412(double pt, double u, double t, double s,
                               double q2) {
  double pt2 = std::pow(pt, 2);
  double q4 = std::pow(q2, 2);
  double q6 = std::pow(q2, 3);
  double q8 = std::pow(q2, 4);
  double QQ2 = u + t + s - q2;
  double x1 = t / (t - QQ2);
  double x2 = u / (u - QQ2);
  double xQt2 = QQ2 + pt2;
  double t2 = std::pow(t, 2);
  double t3 = std::pow(t, 3);
  double u2 = std::pow(u, 2);
  double u3 = std::pow(u, 3);
  double u4 = std::pow(u, 4);
  double s2 = std::pow(s, 2);
  double s3 = std::pow(s, 3);
  double s4 = std::pow(s, 4);

  // Eqs. (A.2)-(A.3) G & S
  double A = s + q2 - QQ2;
  double B = sqrt(std::pow(A, 2) - 4. * q2 * s);
  double L1b = log(q2 / s / std::pow(x2, 2));
  double L2a = log(q2 * s / std::pow(A - s * x1, 2));
  double L3 = log((A + B) / (A - B));

  double result =
      s * pt2 * QQ2 * std::pow(A, 3) / 2. / t / (u - q2) * (L2a + L1b) +
      s * pt2 * (u + t) / 16 / u / t / B *
          (std::pow(A, 4) + 6. * std::pow(A, 2) * std::pow(B, 2) +
           std::pow(B, 4)) *
          L3 +
      (-s * pt2 / 2 / u / t *
           (std::pow(s - QQ2, 4) + q8 + 2 * QQ2 * A * std::pow(s - QQ2, 2) -
            2. * QQ2 * q6) -
       std::pow(s * pt2 * q2 * QQ2, 2) / u4 +
       2. * s * pt2 * q2 * QQ2 * (std::pow(A, 2) - s * q2) / u2) *
          L1b +
      s * pt2 / (8 * u * t) *
          (std::pow(QQ2 - u, 3) / (QQ2 - t) * ((QQ2 - t) * s + QQ2 * u) +
           std::pow(QQ2 - t, 3) / (QQ2 - u) * ((QQ2 - u) * s + QQ2 * t)) *
          (4. / 3 + 2 * pt2 / xQt2 + 4. * std::pow(pt2 / xQt2, 2) -
           44. / 3. * std::pow(pt2 / xQt2, 3)) +
      s * pt2 * std::pow(QQ2 - u, 2) / 4. / u / t / (QQ2 - t) *
          (-3. * (t - q2) * ((QQ2 - t) * s + QQ2 * u) -
           QQ2 * (q2 * (t - q2) + QQ2 * (u - q2))) *
          (1. + 2. * pt2 / xQt2 - 6. * std::pow(pt2 / xQt2, 2)) +
      s * pt2 * std::pow(QQ2 - t, 2) / 4. / u / t / (QQ2 - u) *
          (-3. * (u - q2) * ((QQ2 - u) * s + QQ2 * t) +
           3 * QQ2 * (q2 * (u - q2) + QQ2 * (t - q2)) + 4. * u * s2) *
          (1. + 2. * pt2 / xQt2 - 6. * std::pow(pt2 / xQt2, 2)) +
      s * pt2 * (QQ2 - u) / 2. / u / t / (QQ2 - t) *
          (3. * std::pow(t - q2, 2) * ((QQ2 - t) * s + QQ2 * u) +
           3. * (t - q2) * QQ2 * (q2 * (t - q2) + QQ2 * (u - q2)) +
           QQ2 * u * (q2 * (t - QQ2) + QQ2 * (u - q2))) *
          (1. - 2. * pt2 / xQt2) +
      s * pt2 * (QQ2 - t) / 2 / u / t / (QQ2 - u) *
          (3 * std::pow(u - q2, 2) * ((QQ2 - u) * s + QQ2 * t) +
           8. * u * t * s2 + 2. * u * s3 - 2. * QQ2 * u * std::pow(u - QQ2, 2) -
           3. * q2 * QQ2 * std::pow(t - q2, 2) - 3. * QQ2 * (q2 - QQ2) * t2 -
           QQ2 * u * (4. * u * t - u * q2 - QQ2 * t + 2. * t2 - 4. * q4) +
           3. * q2 * std::pow(QQ2, 2) * (t - q2) + q2 * QQ2 * u * (t - QQ2)) *
          (1. - 2. * pt2 / xQt2) -
      4. * std::pow(s * pt2 * q2 * QQ2, 2) / u4 +
      s * pt2 * q2 * QQ2 * (std::pow(B, 2)) / (2. * u2) +
      s2 * pt2 * q4 / 6. * ((s + QQ2) / u / t + QQ2 / u2 + QQ2 / t2) +
      2. * s2 * pt2 * QQ2 * q4 / u2 + s2 * pt2 * q4 / u -
      s2 * pt2 / 12. / u / t *
          (30. * q6 + 54. * (std::pow(QQ2, 2)) * q2 + 8. * std::pow(QQ2, 3)) +
      s * pt2 / 12. / u / t *
          (11. * s4 + 17. * q8 +
           QQ2 * (61. * u2 * t + 17. * u3 + 73. * u * t2 + 29. * t3) +
           q2 * (24. * u2 * t + 6. * u3 + 36. * u * t2 + 18. * t3) +
           std::pow(QQ2, 2) * (-21. * u2 - 33. * t2 - 52. * u * t) +
           q2 * QQ2 * (-73. * u2 - 109. * t2 - 170. * u * t) +
           q4 * (-23. * u2 - 35. * t2 - 52. * u * t) +
           q4 * QQ2 * (134. * t + 110. * u) + 4. * std::pow(QQ2, 4) +
           52. * q2 * std::pow(QQ2, 3) + 20. * q4 * std::pow(QQ2, 2) -
           22. * q6 * QQ2);

  return result;
}

double HiggsDpTpartonic::A1324(double pt, double u, double t, double s,
                               double q2) {
  double pt2 = std::pow(pt, 2);
  double q4 = std::pow(q2, 2);
  double q8 = std::pow(q2, 4);
  double QQ2 = u + t + s - q2;
  double x1 = t / (t - QQ2);
  double x2 = u / (u - QQ2);
  double t2 = std::pow(t, 2);
  double t4 = std::pow(t, 4);
  double u3 = std::pow(u, 3);
  double u4 = std::pow(u, 4);
  double s2 = std::pow(s, 2);

  // Eqs. (A.2)-(A.3) G & S
  double A = s + q2 - QQ2;
  double B = sqrt(std::pow(A, 2) - 4. * q2 * s);
  double L1a = log(q2 / s / std::pow(x1, 2));
  double L1b = log(q2 / s / std::pow(x2, 2));
  double L2a = log(q2 * s / std::pow(A - s * x1, 2));
  double L2b = log(q2 * s / std::pow(A - s * x2, 2));
  double L3 = log((A + B) / (A - B));

  double result =
      -0.5 * ((std::pow(s * pt2 / t, 4) + std::pow(q2 * QQ2 / t, 4)) * L1a +
              u4 * L1b) +
      s2 * pt2 * q4 * QQ2 / t2 * L1a +
      (s * q4 * QQ2 * u3 / t / (u - q2) / (t - q2) + s * pt2 * u3 / t) / 2 *
          (L1a + L1b) +
      s2 * pt2 * (1 - x2) * u3 / A / 2. / (t - q2) * (L2b - L1b) +
      s2 * pt2 * (1 - x1) * (q8 + std::pow(u - q2, 4)) / 2. / A / t / (u - q2) *
          (L2a - L1a) +
      s2 * pt2 * q4 * QQ2 / A / B * L3 +
      s * pt2 * q2 * QQ2 / 2. / t4 *
          (std::pow(s * pt2, 2) - 6. * s * pt2 * q2 * QQ2 +
           q4 * std::pow(QQ2, 2));

  return result;
}

double HiggsDpTpartonic::A3241(double pt, double u, double t, double s,
                               double q2) {
  double pt2 = std::pow(pt, 2);
  double q4 = std::pow(q2, 2);
  double QQ2 = u + t + s - q2;
  double x1 = t / (t - QQ2);
  double xQt2 = QQ2 + pt2;
  double t2 = std::pow(t, 2);
  double t4 = std::pow(t, 4);
  double u2 = std::pow(u, 2);
  double s2 = std::pow(s, 2);

  // Eqs. (A.2)-(A.3) G & S
  double A = s + q2 - QQ2;
  double B = sqrt(std::pow(A, 2) - 4. * q2 * s);
  double L1a = log(q2 / s / std::pow(x1, 2));
  double L2a = log(q2 * s / std::pow(A - s * x1, 2));

  double result =
      s2 * pt2 * std::pow(A, 3) * (1 - x1) / 2. / t / (u - q2) * (L2a - L1a) +
      (-std::pow(s * pt2 * q2 * QQ2, 2) / t4 +
       s * pt2 * q4 * std::pow(QQ2, 2) / u / t -
       s * pt2 * q2 * QQ2 * std::pow(A, 4) / 2. / u / t / (u - q2) / (t - q2) +
       s * pt2 * QQ2 * q2 * (u + t) * (2. * std::pow(A, 2) - s * q2) / u / t2) *
          L1a +
      s2 * pt2 * QQ2 * std::pow(QQ2 - u, 2) / 2. / u / t /
          std::pow(QQ2 - t, 2) * (-u * t - std::pow(QQ2 - t, 2)) *
          (-3. + 10. * QQ2 / xQt2 - 6. * std::pow(QQ2 / xQt2, 2)) +
      s2 * pt2 * QQ2 * (QQ2 - u) / u / t / std::pow(QQ2 - t, 2) *
          (u * t * (QQ2 - u) - std::pow(QQ2 - t, 3) -
           q2 * std::pow(QQ2 - t, 2) - q2 * (QQ2 - t) * (QQ2 - u)) *
          (-1. + 2. * QQ2 / xQt2) +
      s * pt2 * q2 * QQ2 *
          (std::pow(B, 2) / 2. / t2 - 2. * q2 * QQ2 / t2 +
           std::pow(u + t, 2) / 2. / u / t) -
      4. * std::pow(s * pt2 * q2 * QQ2, 2) / t4 +
      s2 * pt2 * QQ2 / 4. / u / t *
          (std::pow(t + u, 2) - (t + u) * (6. * QQ2 + 4. * q2) +
           6. * std::pow(QQ2, 2) + 8. * q2 * QQ2) +
      s2 * pt2 * QQ2 * q4 * std::pow(t + u, 2) / 4. / u2 / t2;

  return result;
}
//==============================================================================================//

//==============================================================================================//
//                                      B Functions //
//----------------------------------------------------------------------------------------------//
double HiggsDpTpartonic::B1pmB1pp(double pt, double u, double t, double s,
                                  double q2) {
  double pt2 = std::pow(pt, 2);
  double QQ2 = u + t + s - q2;
  double x1 = t / (t - QQ2);
  double xQt2 = QQ2 + pt2;
  double t2 = std::pow(t, 2);
  double s2 = std::pow(s, 2);

  double B1pm = std::pow(s2, 2) * pt2 * x1 * std::pow(1 - x1, 3) / t +
                std::pow(s2, 2) * std::pow(pt2, 3) * std::pow(x1, 3) *
                    (1 - x1) / std::pow(t, 3) +
                4 * std::pow(s2, 2) * std::pow(pt2, 2) * std::pow(x1, 2) *
                    std::pow(1 - x1, 2) / t2 -
                s2 * pt2 * QQ2 * (1 + log(pt2 / xQt2));
  double B1pp =
      std::pow(s2, 2) * pt2 * std::pow(x1, 3) * (1 - x1) / t +
      std::pow(s2, 2) * std::pow(pt2, 3) * x1 * std::pow(1 - x1, 3) /
          std::pow(t, 3) +
      4 * std::pow(s2, 2) * std::pow(pt2, 2) * std::pow(x1, 2) *
          std::pow(1 - x1, 2) / t2 -
      s2 * pt2 * QQ2 / std::pow(xQt2, 2) / u / t *
          (std::pow(u * t + pt2 * QQ2, 2) + 2 * s * pt2 * QQ2 * xQt2) +
      s2 * pt2 * QQ2 / u / t * (s2 + std::pow(QQ2, 2)) * log(xQt2 / QQ2);

  double result = B1pm + B1pp;

  return result;
}

double HiggsDpTpartonic::B2pmB2pp(double pt, double u, double t, double s,
                                  double q2) {
  double pt2 = std::pow(pt, 2);
  double q4 = std::pow(q2, 2);
  double q6 = std::pow(q2, 3);
  double QQ2 = u + t + s - q2;
  double x1 = t / (t - QQ2);
  double xQt2 = QQ2 + pt2;
  double t2 = std::pow(t, 2);
  double u2 = std::pow(u, 2);
  double s2 = std::pow(s, 2);
  double s3 = std::pow(s, 3);
  double s4 = std::pow(s, 4);

  double B2pm = 1. / 3. * std::pow(t / x1, 4) * pt2 / xQt2 *
                    (std::pow(pt2, 3) - std::pow(QQ2, 3) - std::pow(xQt2, 3)) /
                    std::pow(xQt2, 3) -
                s2 * pt2 * QQ2 / 3.;

  double B2pp =
      -s2 * pt2 * QQ2 / (2. * u * t) * (s2 + std::pow(QQ2, 2)) *
          log(xQt2 / QQ2) +
      s * pt2 * std::pow(QQ2 - u, 3) / (2. * u * t) / (QQ2 - t) *
          ((QQ2 - t) * s + QQ2 * u) *
          (2. / 3. + QQ2 / xQt2 - 10. / 3. * std::pow(QQ2 / xQt2, 3)) -
      s * pt2 * std::pow(QQ2 - u, 2) / 2. / u / t / std::pow(QQ2 - t, 2) *
          (3. * std::pow(QQ2 - t, 3) * QQ2 +
           (QQ2 - t) * QQ2 * (2. * u * t + q4) +
           std::pow(QQ2 - t, 2) * (s2 + 4. * q2 * QQ2 - u * (QQ2 + q2)) -
           u2 * QQ2 * QQ2 + u * t2 * q2) *
          (1. - 2. * std::pow(QQ2 / xQt2, 2)) +
      s * pt2 * (QQ2 - u) / 2. / u / t / (QQ2 - t) *
          (3. * QQ2 * s * (QQ2 + s) * (QQ2 - t) - t * s3 + q2 * QQ2 * s2 +
           QQ2 * u * std::pow(q2 - QQ2, 2)) *
          (1. - 2. * QQ2 / xQt2) +
      s * pt2 / 12. / u / t *
          (-2. * s4 + 6. * s * q2 * t * (t - q2) + 2. * s * q6 +
           8. * QQ2 * s * std::pow(s - QQ2, 2) - 2. * u * t * s * QQ2 +
           7. * s2 * q2 * QQ2 - 2. * s * std::pow(QQ2, 2) * q2 - q6 * QQ2 +
           3. * q2 * std::pow(QQ2, 3) - 4. * u * t * q2 * QQ2) +
      11. / 6. * s3 * pt2 * std::pow(QQ2, 2) / u / t -
      s2 * pt2 * q4 * QQ2 / 3. / t2;

  double result = B2pm + B2pp;

  return result;
}
//==============================================================================================//

//==============================================================================================//
//                          Regular channel-specific coefficients //
//----------------------------------------------------------------------------------------------//
double HiggsDpTpartonic::Rgq(double pt, double u, double t, double s,
                             double q2) {
  //////////////////////////////////////////////////////////////
  // gq->H+X terms in Eq. (A.17)                              //
  //////////////////////////////////////////////////////////////
  double pt2 = std::pow(pt, 2);
  double q4 = std::pow(q2, 2);
  double q6 = std::pow(q2, 3);
  double QQ2 = u + t + s - q2;
  double x1 = t / (t - QQ2);
  double x2 = u / (u - QQ2);
  double xQt2 = QQ2 + pt2;
  double t2 = std::pow(t, 2);
  double u2 = std::pow(u, 2);
  double u4 = std::pow(u, 4);
  double s2 = std::pow(s, 2);
  double s4 = std::pow(s, 4);

  // Eqs. (A.2)-(A.3) G & S
  double A = s + q2 - QQ2;
  double B = sqrt(std::pow(A, 2) - 4. * q2 * s);
  double L1a = log(q2 / s / std::pow(x1, 2));
  double L1b = log(q2 / s / std::pow(x2, 2));
  double L2a = log(q2 * s / std::pow(A - s * x1, 2));
  double L2b = log(q2 * s / std::pow(A - s * x2, 2));
  double L3 = log((A + B) / (A - B));

  // Eqs. (A.18)-(A.26) G & S
  double C1pm = -2. * s2 * pt2 * QQ2 * log(pt2 / xQt2);
  double C2pm = 0.;
  double C1mp = s2 * pt2 * QQ2 - 1.5 * s2 * pt2 * t2 / u +
                pt2 / 2. / xQt2 *
                    (s * std::pow(QQ2 - t, 3) + QQ2 * std::pow(QQ2 - u, 3)) *
                    (-3. + 10. * QQ2 / xQt2 - 6 * std::pow(QQ2 / xQt2, 2));
  double C2mp = pt2 * QQ2 / std::pow(xQt2, 2) *
                    (s * std::pow(QQ2 - t, 3) + QQ2 * std::pow(QQ2 - u, 3)) *
                    (-2. + 3. * QQ2 / xQt2) +
                2. * s2 * pt2 * QQ2 + 4. * s2 * pt2 * QQ2 * log(pt2 / xQt2);
  double C1pp = -1.5 * s4 * pt2 / u - s * pt2 * QQ2 * A * A / u * L2b +
                s * pt2 / u / t2 * L1a *
                    (std::pow(A - q2, 2) * s * t2 - 2. * QQ2 * q2 * u * t * A -
                     QQ2 * q4 * (QQ2 - t) * u) +
                s * pt2 / u / B * L3 *
                    ((s + QQ2 - q2) * (s * std::pow(A - q2, 2) + QQ2 * A * A) -
                     4. * s * QQ2 * A * (A - q2)) +
                0.5 * s * pt2 *
                    (std::pow(QQ2 - t, 2) * (s / (QQ2 - u) - QQ2 / u) +
                     std::pow(QQ2 - u, 2) * (QQ2 / (QQ2 - t) + s / u)) *
                    (-3. + 10. * QQ2 / xQt2 - 6. * std::pow(QQ2 / xQt2, 2)) +
                s * pt2 * (QQ2 - t) / u / (QQ2 - u) *
                    (2 * s * u * (s + t) -
                     QQ2 * (4. * q2 * QQ2 - QQ2 * t - q2 * u - 2 * u * t)) *
                    (-1 + 2 * QQ2 / xQt2) +
                s * pt2 * (QQ2 - u) / u / (QQ2 - t) *
                    (t * s2 - 2. * u * t * s + 2. * QQ2 * u * (QQ2 - t)) *
                    (-1. + 2. * QQ2 / xQt2) +
                s * pt2 * QQ2 * q2 * (u + t) / t -
                2. * s * pt2 * QQ2 * QQ2 * q4 / t2 +
                s2 * pt2 * QQ2 * q4 / 2. / u2 +
                s * pt2 / 2 / u *
                    (-2. * (QQ2 + q2) * u * s + 2. * q2 * s2 + q4 * s +
                     q2 * QQ2 * (2. * (s - QQ2) + 3. * q2 - u) +
                     5. * QQ2 * s * (s - QQ2));
  double C2pp = 0.5 * s * pt2 * A * A * (1 - x1) * (L2a - L1a) +
                s * pt2 * (q2 - t) * A * A * (1. - x2) / 2. / u * (L1b - L2b) +
                0.5 * s * pt2 / u * (L1b - L1a) * std::pow(s - QQ2, 3) +
                s * pt2 * QQ2 * A * A / (q2 - u) * (L1b - L2a) +
                s * pt2 * QQ2 / u2 * L1b *
                    (2. * u * std::pow(s - QQ2, 2) + 4. * q2 * (q2 - t) * A -
                     2. * q4 * (QQ2 - t) - q6) -
                0.5 * s * pt2 *
                    (std::pow(QQ2 - t, 2) * (s / (QQ2 - u) - QQ2 / u) +
                     std::pow(QQ2 - u, 2) * (QQ2 / (QQ2 - t) + s / u)) *
                    (-3. + 10. * QQ2 / xQt2 - 6. * std::pow(QQ2 / xQt2, 2)) +
                0.5 * s * pt2 * (QQ2 - t) / u / (QQ2 - u) *
                    ((-3. * s * pt2 + u2 - QQ2 * QQ2) * (s + QQ2) - q2 * s * u +
                     2. * QQ2 * QQ2 * (QQ2 - u)) *
                    (-1. + 2. * QQ2 / xQt2) +
                0.5 * s * pt2 * (QQ2 - u) / u / (QQ2 - t) *
                    (3. * s * pt2 * (s + QQ2) - u * QQ2 * (QQ2 - t) +
                     3. * QQ2 * q2 * (QQ2 - u) + s * t * (u + s)) *
                    (-1. + 2. * QQ2 / xQt2) +
                0.5 * s * pt2 * q2 * QQ2 / u2 *
                    (2. * std::pow(s - QQ2, 2) - 2. * q2 * (s - q2) -
                     u * (QQ2 - u) - 4. * q2 * QQ2) +
                s2 * pt2 * (u - t) * (q2 + QQ2) / (2. * u) -
                2. * std::pow(s * pt2 * q2 * QQ2, 2) / u4 * (4. + L1b);
  double C1mm =
      s2 * pt2 * t2 / u * L1a - s * pt2 * QQ2 * std::pow(q2 - t, 2) / u * L2b +
      s * pt2 * q2 * QQ2 / std::pow(B, 2) * (t * (u + t) - 2. * q2 * QQ2) +
      s * pt2 / u / B *
          (t2 * std::pow(B, 2) - q2 * t2 * (u + t) +
           2. * std::pow(QQ2, 2) * q4 + QQ2 * q4 * (3. * t - u) +
           QQ2 * q4 * u / std::pow(B, 2) *
               (-t * (u + t) + 2. * q2 * QQ2 + QQ2 * (t - u))) *
          L3;
  double C2mm =
      s * pt2 * t2 * QQ2 / (2. * u) * (L1a + 3. * L1b) +
      0.5 * s * pt2 * t2 * (1 - x1) * (L2a - L1a) +
      s * pt2 * QQ2 * std::pow(q2 - t, 3) * x2 / (2. * u2) * (L2b - L1b) +
      s * pt2 * t2 * QQ2 / (q2 - u) * (L1b + L2a) +
      s2 * pt2 * t2 / 2. / u * (L1b - L1a) +
      s * pt2 * q2 * QQ2 / u2 * (4. * t * (t - q2) + q4) * L1b -
      2. * std::pow(s * pt2 * q2 * QQ2, 2) / u4 * (L1b + 4.) +
      s * pt2 * t2 * q2 * QQ2 / u2;
  double C2ep = 4. / u4 * std::pow(s * pt2 * q2 * QQ2, 2);

  double result = 1. / s2 / pt2 / QQ2 *
                  (16. / 9. * (C1pm + C1mp + C1pp + C1mm) +
                   4. * (C2pm + C2mp + C2pp + C2mm + C2ep));

  return result;
}

double HiggsDpTpartonic::Rqqb(double pt, double u, double t, double s,
                              double q2) {
  //////////////////////////////////////////////////////////////
  // qqbar->H+X terms in Eq. (A.27)                           //
  //////////////////////////////////////////////////////////////
  double pt2 = std::pow(pt, 2);
  double q4 = std::pow(q2, 2);
  double QQ2 = u + t + s - q2;
  double x1 = t / (t - QQ2);
  double x2 = u / (u - QQ2);
  double xQt2 = QQ2 + pt2;
  double t2 = std::pow(t, 2);
  double s2 = std::pow(s, 2);
  double s3 = std::pow(s, 3);
  double u2 = std::pow(u, 2);

  // Eqs. (A.2)-(A.3) G & S
  double A = s + q2 - QQ2;
  double B = sqrt(std::pow(A, 2) - 4. * q2 * s);
  double L1a = log(q2 / s / std::pow(x1, 2));
  double L1b = log(q2 / s / std::pow(x2, 2));
  double L2a = log(q2 * s / std::pow(A - s * x1, 2));
  double L2b = log(q2 * s / std::pow(A - s * x2, 2));
  double L3 = log((A + B) / (A - B));

  // Eqs. (A.28)-(A.34) G & S
  double D1pm = -s2 * pt2 * QQ2 * (1. + log(pt2 / xQt2)) -
                s3 * std::pow(pt, 2) * x1 * (1. - x1) / t -
                s3 * pt2 * std::pow(1 - x1, 2);
  double D2pm = -s2 * pt2 * QQ2 / 3. -
                s * pt2 * t2 / 6. / std::pow(x1, 2) *
                    (11. - 12. * QQ2 / xQt2 + 3. * std::pow(QQ2 / xQt2, 2)) +
                11. * s * pt2 * t2 / 6.;
  double D1pp = s2 * pt2 * u2 * (1 - x2) * (L2b - L1b) / A +
                s2 * pt2 * std::pow(q2 - u, 2) * (1 - x1) / A * (L2a - L1a) +
                s * pt2 * q2 * QQ2 * (s * pt2 + u * t) / t2 * L1a -
                2. * s2 * pt2 * q4 * QQ2 / A / B * L3 +
                s * pt2 * q2 * QQ2 * (2. * s * pt2 - u * t) / t2;

  double D2pp =
      s * pt2 * u2 * (q2 - t) * (1. - x2) / (2. * A) * (L1b - L2b) +
      s * pt2 * std::pow(q2 - u, 3) * (1 - x1) / (2 * A) * (L1a - L2a) -
      0.5 * s * pt2 * u2 * (L1a + L1b) +
      6. * std::pow(s * pt2 * q2 * QQ2, 2) / std::pow(B, 4) -
      s * pt2 * q2 * QQ2 * u2 / std::pow(B, 2) +
      L3 * (s * pt2 * u2 * (u + t) / B +
            std::pow(s * pt2 * q2 * QQ2, 2) *
                (1. / A / std::pow(B, 3) - 3. * A / std::pow(B, 5)) +
            s * pt2 * q2 * QQ2 *
                ((t - 3. * u) / 2. / B +
                 A * (std::pow(B, 2) + 2 * u2) / 4. / std::pow(B, 3) +
                 (t2 - 6. * u * t + 7. * u2) / 4. / A / B));

  double E1 = 4. / 3. * (2. * s2 * pt2 * QQ2 - s * pt2 * q2 * QQ2);
  double E2 = s * pt2 * QQ2 * (QQ2 - pt2) / std::pow(xQt2, 2) *
                  (std::pow(QQ2 - t, 2) + (std::pow(QQ2 - u, 2))) +
              2. * s2 * pt2 * QQ2;
  double E3 = -2. * s * pt2 * (std::pow(u + t - 2. * QQ2, 2) - 2. * s * pt2) *
                  log(pt2 / xQt2) -
              s * pt2 * QQ2 * (2. * xQt2 + QQ2) / std::pow(xQt2, 2) *
                  (std::pow(QQ2 - t, 2) + (std::pow(QQ2 - u, 2))) -
              6. * s2 * pt2 * QQ2;

  double result =
      (128. / 27. * (D1pm + D1pp) + 32. / 3. * (D2pm + D2pp) +
       0.5 * (16. / 9. * NF * E1 + 16. / 9. * E2 + 16. / 27. * E3)) /
      s2 / pt2 / QQ2;

  return result;
}
//==============================================================================================//

//**********************************************************************************************//
// The below compute the delta and regular terms that enter into the partonic
// cross section     //
//----------------------------------------------------------------------------------------------//

//==============================================================================================//
//                                        Delta Terms //
//----------------------------------------------------------------------------------------------//
double HiggsDpTpartonic::deltapartonic(double pt, double nn, double zz) {
  ///////////////////////////////////////////////////////////////
  // This function computes the terms proportional to delat(Q) //
  // and delta(Q^2) in the SINGULAR part of G &S. These are    //
  // given in Eqs. (3.17), (3.20), (3.24), (3.28).             //
  ///////////////////////////////////////////////////////////////
  double xx = nn;
  double result = 0;
  zz = 1;  // unused integration variable

  // This function calculates the terms of the cross section proportional to
  // delta(Q²)
  double QQ2 = 0;

  double xi = (pt * pt / MH2);
  double tauh = xx * std::pow(sqrt(1 + xi) - sqrt(xi), 2);
  double sh = MH2 / tauh;
  double mt2 = pt * pt + MH2;
  double uh = 0.5 * (QQ2 + MH2 - sh -
                     sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));
  double th = 0.5 * (QQ2 + MH2 - sh +
                     sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));

  // The jacobian is defined  such that the integration variable is Q². This,
  // toghether with setting QQ2=0 is needed to properly handle the delta(Q²)
  // term.
  double jac = 1. / sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2);

  // compute terms proportional to delta(Q^2)
  if (ORD >= 0) {
    switch (CHANNEL) {
      case (0): {
        result += gg0(sh, th, uh, MH2);
        result = 0;
      }  // gg-channel
      break;
      case (1): {
        result += gq0(sh, th, uh);
      }  // gq-channel
      break;
      case (2): {
        result += qg0(sh, th, uh);
      }  // qg-channel
      break;
      case (3): {
        result += qqb0(sh, th, uh);
      }  // qqb-channel
      break;
    }
  }
  if (ORD >= 1) {
    switch (CHANNEL) {
      // We multiply by aass/(2*M_PI) because it is one order higher than the LO
      // part of the xsec that is also calculated in this function.
      case (0):  // gg-channel
      {
        // dominant contribution in small-x are the terms multiplied by uu.
        // specifically the logs (not Li2) dominate uu
        double uu =
            0.5 * std::pow(log(uh / th), 2) + std::pow(M_PI, 2) / 3. -
            log(sh / MH2) * log(-th / MH2) - log(sh / MH2) * log(-uh / MH2) -
            log(-uh / MH2) * log(-th / MH2) + std::pow(log(MH2 / sh), 2) +
            std::pow(log(MH2 / (MH2 - th)), 2) +
            std::pow(log(MH2 / (MH2 - uh)), 2) + 2. * Li2(1. - MH2 / sh) +
            2. * Li2(MH2 / (MH2 - th)) + 2. * Li2(MH2 / (MH2 - uh));
        uu = 0.5 * std::pow(log(uh / th), 2) + std::pow(M_PI, 2) / 3. -
             log(sh / MH2) * log(-th / MH2) - log(sh / MH2) * log(-uh / MH2) -
             log(-uh / MH2) * log(-th / MH2) + std::pow(log(MH2 / sh), 2) +
             std::pow(log(MH2 / (MH2 - th)), 2) +
             std::pow(log(MH2 / (MH2 - uh)), 2) + 2. * Li2(1. - MH2 / sh) +
             2. * Li2(MH2 / (MH2 - th)) + 2. * Li2(MH2 / (MH2 - uh));
        double de = 1.5 * beta0 * (log(-MUR2 / th) + log(-MUR2 / uh)) +
                    67. / 6. - 5. / 9. * NF;
        result += aass / (2. * M_PI) *
                  ((11. + de + 3. * uu) * gg0(sh, th, uh, MH2) +
                   (3. - NF) * (std::pow(MH2, 2) / sh + std::pow(MH2, 2) / th +
                                std::pow(MH2, 2) / uh + MH2));
      } break;
      case (1):  // gq-channel
      {
        double V1 =
            0.5 * (std::pow(log(uh / th), 2) + std::pow(log(-sh / uh), 2) -
                   std::pow(log(-sh / th), 2)) +
            log(sh / MH2) * log(-th / MH2) - log(sh / MH2) * log(-uh / MH2) -
            log(-th / MH2) * log(-uh / MH2) + 2. * Li2(MH2 / (MH2 - uh)) +
            std::pow(log(MH2 / (MH2 - uh)), 2) + std::pow(M_PI, 2);
        double V2 =
            std::pow(log(MH2 / sh), 2) + std::pow(log(MH2 / (MH2 - th)), 2) -
            2. * log(sh / MH2) * log(-th / MH2) + 2. * Li2(1. - MH2 / sh) +
            2. * Li2(MH2 / (MH2 - th)) - 3.5 - 2. * std::pow(M_PI, 2) / 3.;
        double V3 = beta0 * (2 * log(-MUR2 / uh) + log(-MUR2 / th)) + 67. / 3. -
                    10. * NF / 9.;

        result += aass / (2. * M_PI) *
                  ((11. + 3 * V1 + 4. * V2 / 3. + V3) * gq0(sh, th, uh) +
                   20. / 9. *
                       (std::pow(sh, 2) + std::pow(th, 2) + std::pow(uh, 2) -
                        uh * MH2) /
                       (-uh));
      } break;
      case (2):  // qg-channel
      {
        double V1c =
            0.5 * (std::pow(log(th / uh), 2) + std::pow(log(-sh / th), 2) -
                   std::pow(log(-sh / uh), 2)) +
            log(sh / MH2) * log(-uh / MH2) - log(sh / MH2) * log(-th / MH2) -
            log(-uh / MH2) * log(-th / MH2) + 2. * Li2(MH2 / (MH2 - th)) +
            std::pow(log(MH2 / (MH2 - th)), 2) + std::pow(M_PI, 2);
        double V2c =
            std::pow(log(MH2 / sh), 2) + std::pow(log(MH2 / (MH2 - uh)), 2) -
            2 * log(sh / MH2) * log(-uh / MH2) + 2. * Li2(1. - MH2 / sh) +
            2. * Li2(MH2 / (MH2 - uh)) - 3.5 - 2. * std::pow(M_PI, 2) / 3.;
        double V3c = beta0 * (2 * log(-MUR2 / th) + log(-MUR2 / uh)) +
                     67. / 3. - 10. * NF / 9.;
        result += aass / (2. * M_PI) *
                  ((11. + 3. * V1c + 4. * V2c / 3. + V3c) * qg0(sh, th, uh) +
                   20. / 9. *
                       (std::pow(sh, 2) + std::pow(th, 2) + std::pow(uh, 2) -
                        th * MH2) /
                       (-th));
      } break;
      case (3):  // qqbar-channel
      {
        double W1 =
            log(-uh / MH2) * log(-th / MH2) - log(sh / MH2) * log(-uh / MH2) -
            log(sh / MH2) * log(-th / MH2) + 2. * Li2(1. - MH2 / sh) +
            std::pow(log(MH2 / sh), 2) - 0.5 * std::pow(log(uh / th), 2) -
            5. * std::pow(M_PI, 2) / 3.;
        double W2 =
            1.5 * log(std::pow(sh, 2) / th / uh) + std::pow(log(uh / th), 2) -
            2. * log(-uh / MH2) * log(-th / MH2) +
            std::pow(log(MH2 / (MH2 - uh)), 2) +
            std::pow(log(MH2 / (MH2 - th)), 2) + 2. * Li2(MH2 / (MH2 - uh)) +
            2. * Li2(MH2 / (MH2 - th)) - 7. + 2. * std::pow(M_PI, 2);
        double W3 =
            beta0 / 2 *
                (4 * log(MUR2 / sh) + log(-MUR2 / uh) + log(-MUR2 / th)) +
            (67. / 2. - 5. * NF / 3.);
        result += aass / (2. * M_PI) *
                  ((11. + 3. * W1 + 4. * W2 / 3. + W3) * qqb0(sh, th, uh) +
                   160. / 27. *
                       (std::pow(th, 2) + std::pow(uh, 2) + std::pow(sh, 2) -
                        sh * MH2) /
                       sh);
      } break;
    }
  }

  // Jacobian
  result *= jac;

  // 1/sh as in Eq. 2.4 of G&S
  result *= 1. / sh;

  // dsigma/pt² to dsigma/dpt
  result *= 2. * pt;

  // rapidity can be expressed in Q² in two ways
  result *= 2;

  return result;
}
//==============================================================================================//

//==============================================================================================//
//                                    Non-Delta Terms //
//----------------------------------------------------------------------------------------------//

double HiggsDpTpartonic::distrpartonic(double pt, double nn, double zz1,
                                       double zz2) {
  ////////////////////////////////////////////////////////////////
  // This function computes the remaining terms that are not    //
  // multiplied by either delta(Q) or delta(Q^2) (both singular //
  // and non-Singular)                                          //
  ////////////////////////////////////////////////////////////////
  double nonsingular = 0, a1 = 0, b1 = 0, c1 = 0, a10 = 0, b10 = 0, d10 = 0;

  zz2 = 1;          // unused integration variable
  double qq = zz1;  // qq = QQ2/QQ2max is an integration variable used to
                    // integrate out rapidity
  double xx = nn;   // xx = Q²/sh

  double tiny = 1e-12;
  if (xx < tiny || xx > 1. - tiny) {
    return 0.;
  }

  double xi = (pt * pt / MH2);
  double tauh = xx * std::pow(sqrt(1 + xi) - sqrt(xi), 2);
  double sh = MH2 / tauh;
  double mt2 = pt * pt + MH2;
  double QQ2max = (MH2 + sh - 2 * sqrt(sh * mt2));
  double QQ2 = qq * QQ2max;
  // uh and th are defined as per App. B of Ravindran et al. (2002)
  double uh = 0.5 * (QQ2 + MH2 - sh -
                     sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));
  double th = 0.5 * (QQ2 + MH2 - sh +
                     sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));
  double jac1 =
      QQ2max / sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2);  // Jacobian
  double za = -th / (QQ2 - th);

  // `a1factor' and `b1factor' are the functions inside the plus distributions
  // double a1factor = log(1. - za) / (1. - za);
  // double b1factor = 1. / (1. - za);
  // Alternatively one can write `a1factor' and `b1factor' in terms of qq:
  double a1factor =
      (log(qq) / qq + log(QQ2max * za / -th) / qq) / QQ2max * (-th / za);
  double b1factor = 1. / qq / QQ2max * (-th / za);

  coeff(pt, uh, th, sh, MH2);
  REG(pt, uh, th, sh, MH2);

  double shnew = za * sh;
  double thnew = th;
  double uhnew = za * (uh - MH2) + MH2;

  switch (CHANNEL) {
    // a1:: (log(1-za)/(1-za))+ terms
    // b1:: 1/(1-za))+ terms
    // c1:: regular  (non delta or plus distributions) terms
    case (0):  // gg-channel
    {
      a1 += (1. / (-th) * pgg(za) * gg0(shnew, thnew, uhnew, MH2) +
             (-za / th) * big1);
      b1 +=
          (1. / th * pgg(za) * log(-MUF2 * za / th) *
               gg0(shnew, thnew, uhnew, MH2) +
           za / th * big1 * log((QQ2 + pt * pt) * za / (-th)) + za / th * big2);
      c1 +=
          (1 / (-th) *
               (-2 * NF * Pqg(za) * log(MUF2 / QQ2) + 2 * NF * za * (1 - za)) *
               qg0(shnew, thnew, uhnew) +
           0.5 * big3);
      nonsingular += 0.5 * REGgg;
    } break;
    case (1):  // gq-channel
    {
      a1 += (1. / (-th) * pgg(za) * gq0(sh, th, uh) + (-za / th) * big4);
      b1 += (1. / th * pgg(za) * log(-MUF2 * za / th) * gq0(sh, th, uh) +
             za / th * big4 * log((QQ2 + pt * pt) * za / (-th)));
      c1 += (1 / (-th) * (-Pqg(za) * log(MUF2 / QQ2) + za * (1 - za)) *
                 qqb0(sh, th, uh) +
             0.5 * big5);
      nonsingular += 0.5 * REGgq;
    } break;
    case (2):  // qg-channel
    {
      a1 += -1. / th * pqq(za) * qg0(sh, th, uh);
      b1 += (1. / th * pqq(za) * log(-MUF2 * za / th) * qg0(sh, th, uh) +
             za / th * 8. / 3. * (std::pow(uh, 2) + std::pow(sh, 2)) / (-th));
      c1 += (-1 / th *
                 (4. / 3. * (1 - za) * qg0(sh, th, uh) +
                  (-Pgq(za) * log(MUF2 / QQ2) + 4 / 3 * za) *
                      gg0(shnew, th, uhnew, MH2)) +
             0.5 * big5);
      nonsingular += 0.5 * REGqg;
    } break;
    case (3):  // qqb-channel (same flavours)
    {
      a1 += 1. / (-th) *
            (pqq(za) * qqb0(sh, th, uh) -
             za * 16. / 27. *
                 (std::pow(th, 2) + std::pow(uh, 2) + std::pow(QQ2 - th, 2) +
                  std::pow(QQ2 - uh, 2)) /
                 sh);
      b1 += (1. / th * pqq(za) * log(-MUF2 * za / th) * qqb0(sh, th, uh) -
             za / th * log((QQ2 + std::pow(pt, 2)) * za / (-th)) * 16. / 27. *
                 (std::pow(th, 2) + std::pow(uh, 2) + std::pow(QQ2 - th, 2) +
                  std::pow(QQ2 - uh, 2)) /
                 sh +
             za / th * 16. / 9. * beta0 * (std::pow(uh, 2) + std::pow(th, 2)) /
                 sh);
      c1 +=
          (-1 / th *
               (4. / 3 * (1 - za) * qqb0(sh, th, uh) +
                (-Pgq(za) * log(MUF2 / QQ2) + 4. / 3 * za) * gq0(sh, th, uh)) +
           16. / 9 * (std::pow(sh - QQ2, 2) + std::pow(uh + th - 2 * QQ2, 2)) /
               sh * log(std::pow(pt, 2) / (std::pow(pt, 2) + QQ2)) /
               std::pow(pt, 2));
      nonsingular += 0.5 * REGqqb;
    } break;
    case (4):  // qQ, qQb and qq - channel (same flavours)
    {
      c1 += 3. * (-1. / th * (-Pgq(za) * log(MUF2 / QQ2) + 4. / 3 * za) *
                      gq0(sh, th, uh) +
                  16. / 9 *
                      (std::pow(sh - QQ2, 2) + std::pow(uh + th - 2 * QQ2, 2)) /
                      sh * log(std::pow(pt, 2) / (std::pow(pt, 2) + QQ2)) /
                      std::pow(pt, 2));
      nonsingular += 3. * 0.5 * REGqqpb;
    } break;
  }

  // Here we deal with the za=1 part of the plus distribution. za=1 corresponds
  // to setting QQ2=0.
  QQ2 = 0;
  uh = 0.5 *
       (QQ2 + MH2 - sh - sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));
  th = 0.5 *
       (QQ2 + MH2 - sh + sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));
  za = -th / (QQ2 - th);
  double jac10 = QQ2max / sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2);

  double apolefactor =
      (log(qq) / qq + log(QQ2max * za / -th) / qq) / QQ2max * (-th / za);
  double adeltafactor =
      0.5 * std::pow(log(-QQ2max / th), 2) / QQ2max * (-th / za);
  double a10factor = apolefactor - adeltafactor;
  double bdeltafactor = log(-QQ2max / th) / QQ2max * (-th / za);
  double bpolefactor = 1. / qq / QQ2max * (-th / za);
  double b10factor = bpolefactor - bdeltafactor;

  double d10factor = -th / za / QQ2max;

  shnew = za * sh;
  thnew = th;
  uhnew = za * (uh - MH2) + MH2;

  coeff(pt, uh, th, sh, MH2);

  switch (CHANNEL) {
    case (0):  // gg-channel
    {
      a10 += (1. / (-th) * pgg(za) * gg0(shnew, thnew, uhnew, MH2) +
              (-za / th) * big1);
      b10 +=
          (1. / th * pgg(za) * log(-MUF2 * za / th) *
               gg0(shnew, thnew, uhnew, MH2) +
           za / th * big1 * log((QQ2 + pt * pt) * za / (-th)) + za / th * big2);
      d10 += 1. / th * beta0 * log(-MUF2 * za / th) *
             gg0(shnew, thnew, uhnew, MH2);
    } break;
    case (1):  // gq-channel
    {
      a10 += (1. / (-th) * pgg(za) * gq0(sh, th, uh) + (-za / th) * big4);
      b10 += (1. / th * pgg(za) * log(-MUF2 * za / th) * gq0(sh, th, uh) +
              za / th * big4 * log((QQ2 + pt * pt) * za / (-th)));
    } break;
    case (2):  // qg-channel
    {
      a10 += -1. / th * pqq(za) * qg0(sh, th, uh);
      b10 += (1. / th * pqq(za) * log(-MUF2 * za / th) * qg0(sh, th, uh) +
              za / th * 8. / 3. * (std::pow(uh, 2) + std::pow(sh, 2)) / (-th));
    } break;
    case (3):  // qqb-channel
    {
      a10 += 1. / (-th) *
             (pqq(za) * qqb0(sh, th, uh) -
              za * 16. / 27. *
                  (std::pow(th, 2) + std::pow(uh, 2) + std::pow(QQ2 - th, 2) +
                   std::pow(QQ2 - uh, 2)) /
                  sh);
      b10 += (1. / th * pqq(za) * log(-MUF2 * za / th) * qqb0(sh, th, uh) -
              za / th * log((QQ2 + std::pow(pt, 2)) * za / (-th)) * 16. / 27. *
                  (std::pow(th, 2) + std::pow(uh, 2) + std::pow(QQ2 - th, 2) +
                   std::pow(QQ2 - uh, 2)) /
                  sh +
              za / th * 16. / 9. * beta0 * (std::pow(uh, 2) + std::pow(th, 2)) /
                  sh);
    } break;
  }

  double bfinal = b1 * jac1 * b1factor - b10 * jac10 * b10factor;
  double afinal = a1 * jac1 * a1factor - a10 * jac10 * a10factor;
  double cfinal = c1 * jac1;
  double dfinal = d10 * jac10 * d10factor;
  double nonsingularfinal = nonsingular * jac1;

  double result = afinal + bfinal + cfinal + dfinal + nonsingularfinal;

  // dsigma/pt² to dsigma/dpt
  result *= 2. * pt;

  // 1/sh as in Eq. 2.4 of G&S
  result *= 1. / sh;

  // rapidity can be expressed in Q² in two ways
  result *= 2.;

  // result = 0;
  return result;
}

//==============================================================================================//

double HiggsDpTpartonic::distrpartoniccross(double pt, double nn, double zz1,
                                            double zz2) {
  ////////////////////////////////////////////////////////////////
  // This function computes the remaining terms that are not    //
  // multiplied by either delta(Q) or delta(Q^2) (both singular //
  // and non-Singular)                                          //
  ////////////////////////////////////////////////////////////////
  double nonsingular = 0, a1 = 0, b1 = 0, c1 = 0, a10 = 0, b10 = 0, d10 = 0;

  zz2 = 1;          // unused integration variable
  double qq = zz1;  // qq = QQ2/QQ2max is an integration variable used to
                    // integrate out rapidity
  double xx = nn;   // xx = Q²/sh

  double tiny = 1e-12;
  if (xx < tiny || xx > 1. - tiny) {
    return 0.;
  }

  double xi = (pt * pt / MH2);
  double tauh = xx * std::pow(sqrt(1 + xi) - sqrt(xi), 2);
  double sh = MH2 / tauh;
  double mt2 = pt * pt + MH2;
  double QQ2max = (MH2 + sh - 2 * sqrt(sh * mt2));
  double QQ2 = qq * QQ2max;
  // uh and th are defined as per App. B of Ravindran et al. (2002)
  double th = 0.5 * (QQ2 + MH2 - sh -
                     sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));
  double uh = 0.5 * (QQ2 + MH2 - sh +
                     sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));
  double jac1 =
      QQ2max / sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2);  // Jacobian
  double za = -th / (QQ2 - th);

  // `a1factor' and `b1factor' are the functions inside the plus distributions
  // double a1factor = log(1. - za) / (1. - za);
  // double b1factor = 1. / (1. - za);
  // Alternatively one can write `a1factor' and `b1factor' in terms of qq:
  double a1factor =
      (log(qq) / qq + log(QQ2max * za / -th) / qq) / QQ2max * (-th / za);
  double b1factor = 1. / qq / QQ2max * (-th / za);

  coeff(pt, uh, th, sh, MH2);
  REG(pt, uh, th, sh, MH2);

  double shnew = za * sh;
  double thnew = th;
  double uhnew = za * (uh - MH2) + MH2;

  switch (CHANNEL) {
    // a1:: (log(1-za)/(1-za))+ terms
    // b1:: 1/(1-za))+ terms
    // c1:: regular  (non delta or plus distributions) terms
    case (0):  // gg-channel
    {
      a1 += (1. / (-th) * pgg(za) * gg0(shnew, thnew, uhnew, MH2) +
             (-za / th) * big1);
      b1 +=
          (1. / th * pgg(za) * log(-MUF2 * za / th) *
               gg0(shnew, thnew, uhnew, MH2) +
           za / th * big1 * log((QQ2 + pt * pt) * za / (-th)) + za / th * big2);
      c1 +=
          (1 / (-th) *
               (-2 * NF * Pqg(za) * log(MUF2 / QQ2) + 2 * NF * za * (1 - za)) *
               qg0(shnew, thnew, uhnew) +
           0.5 * big3);
      nonsingular += 0.5 * REGgg;
    } break;
    case (1):  // gq-channel
    {
      a1 += (1. / (-th) * pgg(za) * gq0(sh, th, uh) + (-za / th) * big4);
      b1 += (1. / th * pgg(za) * log(-MUF2 * za / th) * gq0(sh, th, uh) +
             za / th * big4 * log((QQ2 + pt * pt) * za / (-th)));
      c1 += (1 / (-th) * (-Pqg(za) * log(MUF2 / QQ2) + za * (1 - za)) *
                 qqb0(sh, th, uh) +
             0.5 * big5);
      nonsingular += 0.5 * REGgq;
    } break;
    case (2):  // qg-channel
    {
      a1 += -1. / th * pqq(za) * qg0(sh, th, uh);
      b1 += (1. / th * pqq(za) * log(-MUF2 * za / th) * qg0(sh, th, uh) +
             za / th * 8. / 3. * (std::pow(uh, 2) + std::pow(sh, 2)) / (-th));
      c1 += (-1 / th *
                 (4. / 3. * (1 - za) * qg0(sh, th, uh) +
                  (-Pgq(za) * log(MUF2 / QQ2) + 4 / 3 * za) *
                      gg0(shnew, th, uhnew, MH2)) +
             0.5 * big5);
      nonsingular += 0.5 * REGqg;
    } break;
    case (3):  // qqb-channel (same flavours)
    {
      a1 += 1. / (-th) *
            (pqq(za) * qqb0(sh, th, uh) -
             za * 16. / 27. *
                 (std::pow(th, 2) + std::pow(uh, 2) + std::pow(QQ2 - th, 2) +
                  std::pow(QQ2 - uh, 2)) /
                 sh);
      b1 += (1. / th * pqq(za) * log(-MUF2 * za / th) * qqb0(sh, th, uh) -
             za / th * log((QQ2 + std::pow(pt, 2)) * za / (-th)) * 16. / 27. *
                 (std::pow(th, 2) + std::pow(uh, 2) + std::pow(QQ2 - th, 2) +
                  std::pow(QQ2 - uh, 2)) /
                 sh +
             za / th * 16. / 9. * beta0 * (std::pow(uh, 2) + std::pow(th, 2)) /
                 sh);
      c1 +=
          (-1 / th *
               (4. / 3 * (1 - za) * qqb0(sh, th, uh) +
                (-Pgq(za) * log(MUF2 / QQ2) + 4. / 3 * za) * gq0(sh, th, uh)) +
           16. / 9 * (std::pow(sh - QQ2, 2) + std::pow(uh + th - 2 * QQ2, 2)) /
               sh * log(std::pow(pt, 2) / (std::pow(pt, 2) + QQ2)) /
               std::pow(pt, 2));
      nonsingular += 0.5 * REGqqb;
    } break;
    case (4):  // qQ, qQb and qq - channel (same flavours)
    {
      c1 += 3. * (-1. / th * (-Pgq(za) * log(MUF2 / QQ2) + 4. / 3 * za) *
                      gq0(sh, th, uh) +
                  16. / 9 *
                      (std::pow(sh - QQ2, 2) + std::pow(uh + th - 2 * QQ2, 2)) /
                      sh * log(std::pow(pt, 2) / (std::pow(pt, 2) + QQ2)) /
                      std::pow(pt, 2));
      nonsingular += 3. * 0.5 * REGqqpb;
    } break;
  }

  // Here we deal with the za=1 part of the plus distribution. za=1 corresponds
  // to setting QQ2=0.
  QQ2 = 0;
  th = 0.5 *
       (QQ2 + MH2 - sh - sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));
  uh = 0.5 *
       (QQ2 + MH2 - sh + sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2));
  za = -th / (QQ2 - th);
  double jac10 = QQ2max / sqrt(std::pow(sh + MH2 - QQ2, 2) - 4. * sh * mt2);

  double apolefactor =
      (log(qq) / qq + log(QQ2max * za / -th) / qq) / QQ2max * (-th / za);
  double adeltafactor =
      0.5 * std::pow(log(-QQ2max / th), 2) / QQ2max * (-th / za);
  double a10factor = apolefactor - adeltafactor;
  double bdeltafactor = log(-QQ2max / th) / QQ2max * (-th / za);
  double bpolefactor = 1. / qq / QQ2max * (-th / za);
  double b10factor = bpolefactor - bdeltafactor;

  double d10factor = -th / za / QQ2max;

  shnew = za * sh;
  thnew = th;
  uhnew = za * (uh - MH2) + MH2;

  coeff(pt, uh, th, sh, MH2);

  switch (CHANNEL) {
    case (0):  // gg-channel
    {
      a10 += (1. / (-th) * pgg(za) * gg0(shnew, thnew, uhnew, MH2) +
              (-za / th) * big1);
      b10 +=
          (1. / th * pgg(za) * log(-MUF2 * za / th) *
               gg0(shnew, thnew, uhnew, MH2) +
           za / th * big1 * log((QQ2 + pt * pt) * za / (-th)) + za / th * big2);
      d10 += 1. / th * beta0 * log(-MUF2 * za / th) *
             gg0(shnew, thnew, uhnew, MH2);
    } break;
    case (1):  // gq-channel
    {
      a10 += (1. / (-th) * pgg(za) * gq0(sh, th, uh) + (-za / th) * big4);
      b10 += (1. / th * pgg(za) * log(-MUF2 * za / th) * gq0(sh, th, uh) +
              za / th * big4 * log((QQ2 + pt * pt) * za / (-th)));
    } break;
    case (2):  // qg-channel
    {
      a10 += -1. / th * pqq(za) * qg0(sh, th, uh);
      b10 += (1. / th * pqq(za) * log(-MUF2 * za / th) * qg0(sh, th, uh) +
              za / th * 8. / 3. * (std::pow(uh, 2) + std::pow(sh, 2)) / (-th));
    } break;
    case (3):  // qqb-channel
    {
      a10 += 1. / (-th) *
             (pqq(za) * qqb0(sh, th, uh) -
              za * 16. / 27. *
                  (std::pow(th, 2) + std::pow(uh, 2) + std::pow(QQ2 - th, 2) +
                   std::pow(QQ2 - uh, 2)) /
                  sh);
      b10 += (1. / th * pqq(za) * log(-MUF2 * za / th) * qqb0(sh, th, uh) -
              za / th * log((QQ2 + std::pow(pt, 2)) * za / (-th)) * 16. / 27. *
                  (std::pow(th, 2) + std::pow(uh, 2) + std::pow(QQ2 - th, 2) +
                   std::pow(QQ2 - uh, 2)) /
                  sh +
              za / th * 16. / 9. * beta0 * (std::pow(uh, 2) + std::pow(th, 2)) /
                  sh);
    } break;
  }

  double bfinal = b1 * jac1 * b1factor - b10 * jac10 * b10factor;
  double afinal = a1 * jac1 * a1factor - a10 * jac10 * a10factor;
  double cfinal = c1 * jac1;
  double dfinal = d10 * jac10 * d10factor;
  double nonsingularfinal = nonsingular * jac1;

  double result = afinal + bfinal + cfinal + dfinal + nonsingularfinal;

  // dsigma/pt² to dsigma/dpt
  result *= 2. * pt;

  // 1/sh as in Eq. 2.4 of G&S
  result *= 1. / sh;

  // rapidity can be expressed in Q² in two ways
  result *= 2.;

  // result = 0;
  return result;
}

//==============================================================================================//
