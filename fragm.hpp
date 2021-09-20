
#pragma once
#ifndef _FRAGM_HPP_
#define _FRAGM_HPP_

#include "photon.hpp"

/**
 *       J0=16==>  PROCESS : G  G   ---> QJ
 *       J0=15==>  PROCESS : G  G   ---> G
 *       J0=14==>  PROCESS : QI G   ---> G
 *       J0=13==>  PROCESS : QI G   ---> QI
 *       J0=12==>  PROCESS : QI QBI ---> G
 *       J0=11==>  PROCESS : QI QBI ---> QI
 *       J0=10==>  PROCESS : QI G   ---> QBI
 *       J0=9 ==>  PROCESS : QI G   ---> QBK
 *       J0=8 ==>  PROCESS : QI G   ---> QK
 *       J0=7 ==>  PROCESS : QI QI  ---> G
 *       J0=6 ==>  PROCESS : QI QI  ---> QI
 *       J0=5 ==>  PROCESS : QI QBI ---> QK
 *       J0=4 ==>  PROCESS : QI QBK ---> G
 *       J0=3 ==>  PROCESS : QI QBK ---> QI
 *       J0=2 ==>  PROCESS : QI QK  ---> G
 *       J0=1 ==>  PROCESS : QI QK  ---> QI
 *
 * */

/**
 * 
 * CONSTANTS
 * 
 * */
constexpr real N = photon_cs_params::nc;
constexpr real N2 = N * N;
constexpr real N3 = N2 * N;
constexpr real N4 = N3 * N;
constexpr real CF = photon_cs_params::cf;
constexpr real AL = 1.0;
constexpr real CQ = 0.0;
constexpr real vC = N2 - 1.0l;
constexpr real v1 = vC * vC / N;
constexpr real v2 = vC / N;
constexpr real v3 = (N4 - 1.0l) / 2.0l / N2;
constexpr real v4 = vC * vC / 2. / N2;
constexpr real GTR = photon_cs_params::nf / 2.0l;
constexpr natural IA1[]{0, 2, 4, 5, 10, 11, 12, 13, 14, 15};
constexpr natural IA2[]{0, 2, 4, 10, 12, 13};
constexpr natural IA3[]{0, 2, 4, 7, 8, 9, 10, 12, 13};
constexpr real CC[]{8.0l * N2, 8.0l * N2, 8.0l * N2, 8.0l * N2, 8.0l * N2, 8.0l * N2, 8.0l * N2, 8.0l * vC *N, 8.0l * vC *N, 8.0l * vC *N, 8.0l * N2, 8.0l * N2, 8.0l * vC *N, 8.0l * vC *N, 8.0l * vC *vC, 8.0l * vC *vC};

void STRU(real const &up1, real const &upb1, real const &do1, real const &dob1, real const &st1, real const &ch1, real const &gl1, 
real const &up2, real const &upb2, real const &do2, real const &dob2, real const &st2, real const &ch2, real const &gl2, 
real const &pff_up, real const &pff_do, real const &pff_st, real const &pff_ch, real const &pff_gl, 
real *GPPV, real *GPPC);

/**
 *  
 * MACROS
 * 
 */

//no X dependece
constexpr real FQQD = 4.0l / 3.0l * (-9.0l / 2.0l + 2.0l * pi2 / 3.0l);

constexpr real FQQW(const real &X);

constexpr real FQQL(const real &X);

//TODO: remove these completely
constexpr real FQGL = 0.0l, FQGD = 0.0l, FQGW = 0.0l, FGGL = 0.0l, FGGD = 0.0l, FGGW = 0.0l, FGQL = 0.0l, FGQD = 0.0l, FGQW = 0.0l;

//Assuming JMAR = 2 for all of the following
constexpr real CGQD = 0.0l;
constexpr real CGQW = 0.0l;
constexpr real CGQL = 0.0l;
constexpr real CQQD = -(9. / 2.0l + pi2 / 3.0l) * CF;
constexpr real CQQW(real const &X);
constexpr real CQQL(real const &X);
constexpr real CQGD = 0.0l;
constexpr real CQGW = 0.0l;
constexpr real CQGL = 0.0l;
constexpr real CGGD = 0.0l;
constexpr real CGGW = 0.0l;
constexpr real CGGL = 0.0l;

constexpr real HQQD(real const &v, natural const &J);

constexpr real HGQD(real const &v, natural const &J);

constexpr real HQGD(real const &v, natural const &J);

constexpr real HGGD(real const &v, natural const &J);

constexpr real HQQW(real const &w, real const &v, natural const &J);

constexpr real HQGW(real const &w, real const &v, natural const &J);

constexpr real HGGW(real const &w, real const &v, natural const &J);

constexpr real HGQW(real const &w, real const &v, natural const &J);

constexpr real HGQL(real const &w, real const &v, natural const &J);

constexpr real HGGL(real const &w, real const &v, natural const &J);

constexpr real HQQL(real const &w, real const &v, natural const &J);

constexpr real HQGL(real const &w, real const &v, natural const &J);

constexpr real HFQQL(real const &w, real const &v);

constexpr real HFQQW(real const &w, real const &v);

constexpr real HFQQD(real const &v);

constexpr real HFQGL(real const &w, real const &v);

constexpr real HFQGW(real const &w, real const &v);

constexpr real HFQGD(real const &v);

constexpr real HFGGL(real const &w, real const &v);

constexpr real HFGGW(real const &w, real const &v);

constexpr real HFGGD(real const &v);

constexpr real HFGQL(real const &w, real const &v);

constexpr real HFGQW(real const &w, real const &v);

constexpr real HFGQD(real const &v);

constexpr real A(real const &S, real const &T, real const &U);

constexpr real B(real const &S, real const &T, real const &U);

constexpr real C(real const &S, real const &T, real const &U);

constexpr real D(real const &S, real const &T, real const &U);

//QJ+QK-->QJ+QK
constexpr real A0(real const &X, real const &S);

//Q+QB --> Q+QB DIFFERENT FLAvORS
constexpr real A2(real const &X, real const &S);

//Q+Q --> Q+Q SAME FLAvOR
constexpr real B0(real const &X, real const &S);

//Q+QB ->Q+QB SAME FLAvOR
constexpr real D0(real const &X, real const &S);

//Q+QB -->G+G    AND ALSO G+G --> Q+QB
constexpr real D1(real const &X, real const &S);

//Q+G -->Q+G
constexpr real E0(real const &X, real const &S);

//G+G -->G+G
constexpr real F2(real const &X, real const &S);


real AvWPL(real const &w, real const &v, real const &S);

real AvDEL(real const &v, real const &S);

real AvLO(real const &w, real const &v, real const &S);

real AvGO(real const &w, real const &v);

real STRUv(real const &w, real const &v, real const &X3, real const &S);

/*
real STRUvC(real const &w, real const &v, real const &X3, real const &S)
{
    real NF = 2.0l * GTR;
    constexpr real DELTA = 0.0l;
    std::cout << "DELTA DELTA DELTA DELTA DELTA DELTA DELTA \n";
    real LTOT = std::log((v * v) * (1.0l - w) * (1.0l - w) * p2 * DELTA * DELTA / q2_fragm / X3 * X3);
    real Z = 1.0l - v + v * w;
    real Y = v * w / (1.0l - v + v * w);
    real T = -S * (1.0l - Y);
    real U = -S * Y;
    real FAC1 = -4.0l * N2 / (1.0l - v) / w / (1.0l - v + v * w);
    real FAC2 = 2.0l * CF * FAC1;
    real FAC3 = 2.0l * CF * FAC2;
    real PQQ = (-LTOT) * CF * (1.0l + Z * Z) / (1.0l - Z) - CF * (1.0l - Z);
    real PQG = (-LTOT) * 1.0l / 2.0l * (Z * Z + (1.0l - Z) * (1.0l - Z)) - Z * (1.0l - Z);
    real PGQ = (-LTOT) * CF * (1.0l + (1.0l - Z) * (1.0l - Z)) / Z - CF * Z;
    real PGG = (-LTOT) * 2.0l * N * ((1.0l - Z) / Z + Z / (1.0l - Z) + Z * (1.0l - Z));

    switch (J0)
    {
    case 0:
        return FAC1 * (CF / N * ((S * S) + (U * U)) / (T * T) * PQQ);
    case 1:
        return FAC1 * (CF / N * ((S * S) + (U * U)) / (T * T) * PGQ + CF / N * ((S * S) + (T * T)) / (U * U) * PGQ);
    case 2:
        return FAC1 * (CF / N * ((S * S) + (U * U)) / (T * T) * PQQ);
    case 3:
        return FAC1 * (CF / N * ((S * S) + (U * U)) / (T * T) * PGQ + CF / N * ((S * S) + (T * T)) / (U * U) * PGQ);
    case 4:
        return FAC1 * (CF / N * ((T * T) + (U * U)) / (S * S) * PQQ + 2.0l * CF * ((T * T) + (U * U)) * (CF / N / T / U - 1.0l / (S * S)) * PQG);
    case 5:
        return FAC1 * (CF / N * ((T * T * T * T) + (U * U * U * U) + (S * S) * ((T * T) + (U * U)) - 2.0l / N * U * T * (S * S)) / (U * U) / (T * T) * PQQ);
    case 6:
        return FAC1 * (CF / N * ((T * T * T * T) + (U * U * U * U) + (S * S) * ((T * T) + (U * U)) - 2.0l / N * U * T * (S * S)) / (U * U) / (T * T) * PGQ);
    case 7:
    case 8:
    case 9:
        return FAC2 * ((1.0l / (U * U) - CF / N / S / T) * ((S * S) + (T * T)) * PQG);
    case 10:
        return FAC1 * (2.0l * CF * ((T * T) + (U * U)) * (CF / N / T / U - 1.0l / (S * S)) * PQG + CF / N * ((S * S * S * S) + (T * T * T * T) + (U * U) * ((S * S) + (T * T)) - 2.0l / N * S * T * (U * U)) / (S * S) / (T * T) * PQQ);
    case 11:
        return FAC1 * (2.0l * CF * ((T * T) + (U * U)) * (CF / N / T / U - 1.0l / (S * S)) * PGG + CF / N * ((S * S * S * S) + (T * T * T * T) + (U * U) * ((S * S) + (T * T)) - 2.0l / N * S * T * (U * U)) / (S * S) / (T * T) * PGQ + CF / N * ((S * S * S * S) + (U * U * U * U) + (T * T) * ((S * S) + (U * U)) - 2.0l / N * S * U * (T * T)) / (S * S) / (U * U) * PGQ + CF / N * ((T * T) + (U * U)) / (S * S) * 2.0l * (NF - 1.0l) * PGQ);
    case 12:
        return FAC2 * ((1.0l / (U * U) - CF / N / S / T) * ((S * S) + (T * T)) * PQG + (1.0l / (T * T) - CF / N / S / U) * ((S * S) + (U * U)) * PQQ);
    case 13:
        return FAC2 * ((1.0l / (U * U) - CF / N / S / T) * ((S * S) + (T * T)) * PGG + (1.0l / (T * T) - CF / N / S / U) * ((S * S) + (U * U)) * PGQ);
    case 14:
        return FAC3 * (N / 2.0l / CF * ((S * S * S * S) + (T * T * T * T) + (U * U * U * U)) * ((S * S) + (T * T) + (U * U)) / (S * S) / (T * T) / (U * U) * PGG + 1.0l / 2.0l / CF * (CF / N / T / U - 1.0l / (S * S)) * ((T * T) + (U * U)) * PGQ * 2.0l * NF);
    case 15:
        return FAC3 * (N / 2.0l / CF * ((S * S * S * S) + (T * T * T * T) + (U * U * U * U)) * ((S * S) + (T * T) + (U * U)) / (S * S) / (T * T) / (U * U) * PQG + 1.0l / 2.0l / CF * (CF / N / T / U - 1.0l / (S * S)) * ((T * T) + (U * U)) * PQQ);
    }
    std::cout << "STRUvC: J0 out of range (" << J0 << "\n";
    return -1.0l;
};
*/

real FDEL1(real const &v, real const &X3);

real FDEL2(real const &v, real const &X3);

real FvWPL1(real const &w, real const &v, real const &X3);

real FvWPL2(real const &w, real const &v, real const &X3);

real FvLO1(real const &w, real const &v, real const &X3);

real FvLO2(real const &w, real const &v, real const &X3);

real FRESC1(real const &w, real const &v, real const &X3);

real FRESC2(real const &w, real const &v, real const &X3);
/*
real FRESCC1(real const &w, real const &v, real const &X3)
{
    real X1 = V * W / v / w / X3;
    real X2 = (1. - V) / (1. - v) / X3;
    real SH = X1 * X2 * P->S;
    real RRESC = STRUvC(w, v, X3, SH);
    return RRESC / SH;
};

real FRESCC2(real const &w, real const &v, real const &X3)
{
    real X1 = V * W / v / w / X3;
    real X2 = (1.0l - V) / (1.0l - v) / X3;
    real vX = 1.0l - v * w;
    real WX = (1.0l - v) / (1.0l - v * w);
    real SH = X1 * X2 * P->S;
    real RRESCC = STRUvC(WX, vX, X3, SH) * v / vX;
    return RRESCC / SH;
};
*/

/**
 * 
 * MODULE STORAGE
 * 
 */
extern natural J0;
extern real V, W, S;
extern real q2_fac, q2_mu, q2_fragm;

#endif