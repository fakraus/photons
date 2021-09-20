#include "photon.hpp"

/**
 * Integrand of direct cross section
 * (formerly 'FUNCHO')
 **/
int direct(
    const int *,
    const real x[],  // x = (w,v,y,p) (y and p optional)
    const int *,
    real f[],        // f = (photon, photon*p)
    void *userdata)
{
    //variables static for better performance
    //no values are actually reused between calls
    static real p, p2, xt, sqy, y_lo, y_up, y, ey, Jy, Jp;
    //real sqy, yddo, yuup;
    static real q2_fac, q2_mu, q2_fragm, alpha_s;
    static real v_min, v_max, Jv, v, v_2, v_3, w_min, Jw, w, w_2, w_3, vw;
    static real V1, vV1, V1_2, V1_3, V1_4, V1_5, V2, V3, V4, V5, W1, X, X_2, Y, Y_2, LV, L1V, LW, L1W, LLVW, L1VVW;
    static real shatd, shatd1, LMS, LMS1, LMU1, LMSS;
    static real sv, sw, svw, sv2, sw2, sv3, sw3, sV1, sV12, sV2, sX, sY, sLV, sL1V, sLW, sL1W, sLLVW, sL1VVW, sshatd, sLMS, sLMSS, sTQG;
    static real UP1, UPB1, DO1, DOB1, ST1, CH1, GL1;
    static real UP2, UPB2, DO2, DOB2, ST2, CH2, GL2;
    static real UP10, UPB10, DO10, DOB10, ST10, CH10, GL10;
    static real SIGQQB, SIGQG, SIGGQ, LAA, LAA2, PQQB0, PQG0, PGQ0, TQQB, TQG;
    static real QQB1, QQB2_W, QQB2_W1, QQB3, GQQB, QG1, QG2_W, QG2_W1, QG3, QG_VW, QG_SWAP, GQ1, GQ2_W, GQ2_W1, GQ3, GQ, GG, QQ, QB, QBS, QS_1_0, QS_0_1, QS_1_1;
    static real FUNCHO1, FUNCHO2;
    static real GQS1, GQS2, GQS3;
    static real GS1, GS2, GS3, GS4, GS5, GS6, GS7, GS8;
    static real PQQB, PQG, PGQ, PGG, PQQ, PQB, PQBS, PQS;
    static real temp1, temp2;

    auto P{reinterpret_cast<const photon_cs_params *const>(userdata)};

    //****************//
    //INTEGRATING P   //
    //****************//
    Jp = P->int_p ? (P->p_max - P->p_min) : 1.0l;
    p = P->int_p ? P->p_min + Jp * x[P->int_y ? 3 : 2] : P->p;
    p2 = p * p;

    //****************//
    //INTEGRATING Y   //
    //****************//

    xt = 2.0l * p / P->sqrtS;
    sqy = std::log((1.0l + std::sqrt(1.0l - xt * xt)) / xt);
    y_lo = std::max(-sqy, P->y_min);
    y_up = std::min(+sqy, P->y_max);

#if _DEBUG_DIRECT_
    std::cout.precision(12);
    std::cout << "\n sqrtS = " << P->sqrtS;
    std::cout << std::scientific << "\np = " << p;
    std::cout << "\n Jp = " << Jp;
    std::cout << "\n xt = " << xt;
    std::cout << "\n sqy = " << sqy;
    std::cout << "\n y_lo = " << y_lo;
    std::cout << "\n y_up = " << y_up;
    std::cout << std::endl;
#endif

    if (y_up < y_lo)
    {
        f[0] = f[1] = 0.0l;
        return 0;
    }
    Jy = P->int_y ? (y_up - y_lo) : 1.0l;
    y = P->int_y ? (y_lo + x[2] * Jy) : P->y;
    ey = std::exp(y);

#if _DEBUG_DIRECT_
    std::cout.precision(12);
    std::cout << std::scientific << "\n y = " << y;
    std::cout << "\n Jy = " << Jy;
    std::cout << "\n ey = " << ey;
    std::cout << std::endl;
#endif

    //******//
    //SCALES//
    // 63   //
    //******//

    q2_fac = P->sf * p2;
    q2_mu = P->sf_mu * p2;
    q2_fragm = P->sf_fragm * p2;

    alpha_s = get_alpha_s(q2_mu); //ALPATOT and ALPHASHO(order = 1 fixed)

#if _DEBUG_DIRECT_
    std::cout << "\n q2_fac = " << q2_fac;
    std::cout << "\n q2_mu = " << q2_mu;
    std::cout << "\n q2_fragm = " << q2_fragm;
    std::cout << "\n alpha_s = " << alpha_s;
    std::cout << std::endl;
#endif
    //skipping isolation stuff(IISOL, EPSILI, EPSILA, DELTAI, DELTA, ...)

    //*************//
    //V-INTEGRATION//
    //     109     //
    //*************//
    v_min = p / P->sqrtS * ey;
    v_max = 1.0l - p / P->sqrtS / ey;
    Jv = v_max - v_min;
    v = v_min + x[1] * Jv;

#if _DEBUG_DIRECT_
    std::cout << "\n v_min = " << v_min;
    std::cout << "\n v_max = " << v_max;
    std::cout << "\n v = " << v;
    std::cout << "\n jV = " << Jv;
    std::cout << std::endl;
#endif

    //*************//
    //W-INTEGRATION//
    //      115    //
    //*************//
    w_min = v_min / v;
    constexpr real w_max = 1.0l;
    Jw = w_max - w_min;
    w = w_min + x[0] * Jw;

#if _DEBUG_DIRECT_
    std::cout << "\n w_min = " << w_min;
    std::cout << "\n w_max = " << w_max;
    std::cout << "\n w = " << w;
    std::cout << "\n Jw = " << Jw;
    std::cout << std::endl;
#endif

    //************************************//
    //STORE SOME VARIABLES FOR CONVENIENCE//
    //  (formerly 'DEFIN')       458      //
    //  need one set with w=1             //
    //************************************//

    v_2 = v * v;
    v_3 = v_2 * v;
    w_2 = w * w;
    w_3 = w_2 * w;
    vw = v * w;

    V1 = 1.0l - v;
    vV1 = v * V1;
    V1_2 = V1 * V1;
    V1_3 = V1_2 * V1;
    V1_4 = V1_3 * V1;
    V1_5 = V1_4 * V1;
    V2 = 2.0l - v;
    V3 = 3.0l - v;
    V4 = 4.0l - v;
    V5 = 5.0l - v;

    W1 = 1.0l - w;

    X = 1.0l - vw;
    X_2 = X * X;
    Y = V1 + vw;
    Y_2 = Y * Y;

    LV = std::log(v);
    L1V = std::log(V1);

    LW = std::log(w);
    L1W = std::log(W1);

    LLVW = std::log(X);
    L1VVW = std::log(Y);

    shatd = p2 / vw / V1;
    shatd1 = p2 / v / V1;
    LMS = std::log(q2_fac / shatd);
    LMS1 = std::log(q2_fac / shatd1);
    LMU1 = std::log(q2_mu / shatd1);
    LMSS = std::log(q2_fragm / shatd);

    //****************************//
    //CALL OF PARTON DISTRIBUTIONS//
    //****************************//

    P->PDF_proj_ptr(v_min / v / w, q2_fac, UP1, UPB1, DO1, DOB1, ST1, CH1, GL1);
    P->PDF_target_ptr((1.0l - v_max) / V1, q2_fac, UP2, UPB2, DO2, DOB2, ST2, CH2, GL2);
    P->PDF_proj_ptr(v_min / v, q2_fac, UP10, UPB10, DO10, DOB10, ST10, CH10, GL10);

    //****************************************************************//
    //LEADING ORDER SUBPROCESS CROSS SECTIONS FOR THE UNPOLARIZED CASE//
    //  (formerly 'SIGLO')          439                               //
    //****************************************************************//

    //skipping IPOL == 1 case
    SIGQQB = 2.0l * P->cf__nc * pi * (v_2 + V1_2);
    SIGQG = 1.0l / P->nc * pi * v * (1.0l + V1_2);
    SIGGQ = 1.0l / P->nc * pi * V1 * (1.0l + v_2);

#if _DEBUG_DIRECT_
    std::cout << "\n SIGQQB = " << SIGQQB;
    std::cout << "\n SIGQG = " << SIGQG;
    std::cout << "\n SIGGQ = " << SIGGQ;
    std::cout << std::endl;
#endif

    //************************************//
    //LOGARITHM TO ACCOUNT FOR W_DOWN != 0//
    //      132                           //
    //************************************//
    LAA = std::log(1.0l - w_min);
    LAA2 = LAA * LAA;

    //************************************//
    //COMBINATIONS OF PARTON DISTRIBUTIONS//
    //      135                           //
    //************************************//

    PQQB0 = 4.0l / 9.0l * (UP10 * UPB2 + UPB10 * UP2) + 1.0l / 9.0l * (DO10 * DOB2 + DOB10 * DO2) + 1.0l / 9.0l * (ST10 * ST2 + ST2 * ST10) + 4.0l / 9.0l * (CH10 * CH2 + CH2 * CH10);

    PQG0 = (4.0l / 9.0l * (UP10 + UPB10) + 1.0l / 9.0l * (DO10 + DOB10) + 1.0l / 9.0l * (ST10 + ST10) + 4.0l / 9.0l * (CH10 + CH10)) * GL2;

    PGQ0 = (4.0l / 9.0l * (UP2 + UPB2) + 1.0l / 9.0l * (DO2 + DOB2) + 1.0l / 9.0l * (ST2 + ST2) + 4.0l / 9.0l * (CH2 + CH2)) * GL10;

#if _DEBUG_DIRECT_
    std::cout << "\n PQQB0 = " << PQQB0;
    std::cout << "\n PQG0 = " << PQG0;
    std::cout << "\n PGQ0 = " << PGQ0;
    std::cout << std::endl;
#endif

    //LEADING ORDER CONTRIB. AND ALL DELTA FUNCTION CONTRIB. OF THE NEXT-TO-LEADING ORDER RESULT
    TQQB = 1.0l - 2.0l * v + 2.0l * v_2;

    //******************************************************//
    //DELTA FUNCTION CONTRIB. TO Q QBAR TO GAMMA GLUON GLUON//
    //  (formerly 'FQQB1(V, W=1))    490                         //
    //******************************************************//
    QQB1 = P->cf__nc * TQQB * (11.0l / 6.0l * P->nc - P->nf / 3.0l) * LMU1 - P->cf__nc * P->nf * TQQB * (5.0l - 3.0l * LV) / 9.0l + P->cf2__nc * (-(7.0l - pi2) * TQQB - 2.0l * TQQB * L1V * LV + L1V * v * (2.0l + v) + LV * V1 * (2.0l + V1) + LV * LV * (3.0l * v_2 + 2.0l * V1) + L1V * L1V * (2.0l * v + 3.0l * V1_2)) + P->cf * ((67.0l - 6.0l * pi2) * TQQB / 18.0l - L1V * vV1 - LV * (11.0l - 16.0l * vV1) / 6.0l - L1V * L1V * (2.0l * v + 3.0l * V1_2) / 2.0l - LV * LV * v * V2 / 2.0l) - (P->cf2__nc * TQQB * LMS1 * (3.0l - 2.0l * L1V + 2.0l * LV));

    //****************************************************//
    //1/(1-W)_PLUS CONTRIB. TO Q QBAR TO GAMMA GLUON GLUON//
    //  (formerly 'FQQB2(V, 1)' and 'FQQB2(V, W))'    522 //
    //****************************************************//
    temp1 = TQQB * (P->cf__nc * P->nf / 3.0l - P->cf * (11.0l - 12.0l * L1V) / 6.0l - 4.0l * P->cf2__nc * (L1V - LV));

    temp2 = -4.0l * P->cf2__nc * (1.0l - 2.0l * vV1); //removed factor LMS=LMS(W)

    //skipping IPOL == 1 case

    QQB2_W = temp1 + temp2 * LMS;
    QQB2_W1 = temp1 + temp2 * LMS1;

    //************************************************************//
    //(LN(1-W)/(1-W))_PLUS CONTRIB. TO Q QBAR TO GAMMA GLUON GLUON//
    //  (formerly 'FQQB3(V, 1)' and 'FQQB3(V, W))'    547         //
    //************************************************************//
    QQB3 = 2.0l * P->cf__nc * (4.0l * P->cf - P->nc) * TQQB;

#if _DEBUG_DIRECT_
    std::cout << "\n QQB1 = " << QQB1;
    std::cout << "\n QQB2_W = " << QQB2_W;
    std::cout << "\n QQB2_W1 = " << QQB2_W1;
    std::cout << "\n QQB3 = " << QQB3;
    std::cout << std::endl;
#endif

    //***************************************************//
    //DELTA FUNCTION CONTRIB. TO Q GLUON TO GAMMA Q GLUON//
    //  (formerly 'FQG1(V, W=1)')        610             //
    //***************************************************//
    TQG = 2.0l - 2.0l * v + v_2;

    //only considering IPOL == 0 case
    QG1 = TQG * (11.0l / 6.0l * P->nc - P->nf / 3.0l) * LMU1 * v / (2.0l * P->nc) - LV * LV * v * (v_2 - 2.0l * V1) / 4.0l + L1V * vV1 / 2.0l - LV * vV1 + pi2 * v * (1.0l + v) * V1 / 4.0l + L1V * L1V * v * (1.0l + v) * V1 / 4.0l - L1V * LV * v * (1.0l + v) * V1 / 2.0l + P->cf__nc * (-7.0l * TQG * v / 4.0l + L1V * v * (1.0l + 2.0l * v) / 2.0l + pi2 * v * (1.0l - 4.0l * v + 5.0l * v_2) / 6.0l - LV * v * (3.0l * v_2 - 2.0l * V1) / 4.0l + LV * LV * v * (3.0l * v_2 + 2.0l * V1) / 2.0l + L1V * L1V * v * (v_2 + V1_2) / 2.0l - L1V * LV * v * (v_2 + V1_2)) + TQG * LMS1 * (-3.0l * P->cf__nc / 4.0l + P->nf / (6.0l * P->nc) - (11.0l - 12.0l * L1V + 12.0l * LV) / 12.0l) * v;

    //*************************************************//
    //1/(1-W)_PLUS CONTRIB. TO Q GLUON TO GAMMA Q GLUON//
    //  (formerly 'FQG2(V,1)' and 'FQG2(V,W)')  657    //
    //*************************************************//

    //only considering IPOL == 0 case
    temp1 = TQG * (-L1V - P->cf__nc * (3.0l - 4.0l * LV) / 4.0l + LV) * v;
    temp2 = -(P->cf + P->nc) * TQG * v / P->nc;

    QG2_W = temp1 + temp2 * LMS;
    QG2_W1 = temp1 + temp2 * LMS1;

#if _DEBUG_DIRECT_
    std::cout << "\n QG1 = " << QG1;
    std::cout << "\n QG2_W = " << QG2_W;
    std::cout << "\n QG2_W1 = " << QG2_W1;
    std::cout << "\n QG3 = " << QG3;
    std::cout << std::endl;
#endif

    //*********************************************************//
    //(LN(1-W)/(1-W))_PLUS CONTRIB. TO Q GLUON TO GAMMA Q GLUON//
    //  (formerly 'FQG3(V,1)' and 'FQG3(V,W)')      684        //
    //*********************************************************//

    //only considering IPOL == 0 case
    QG3 = (P->cf + 2.0l * P->nc) * TQG * v / P->nc;

    //***************************************************//
    //DELTA FUNCTION CONTRIB. TO GLUON Q TO GAMMA Q GLUON//
    //  (formerly FGQ1(V, W=1))     827                  //
    //***************************************************//

    //only considering IPOL == 0 case
    GQ1 = (1.0l + v_2) * (11.0l / 6.0l * P->nc - P->nf / 3.0l) * LMU1 * V1 / (2.0l * P->nc) + (-(L1V * vV1) + LV * vV1 / 2.0l - (L1V * L1V) * (1 - 4.0l * v + v_2) * V1 / 4.0l + (pi2)*v * V2 * V1 / 4.0l - L1V * LV * v * V2 * V1 / 2.0l + (LV * LV) * v * V2 * V1 / 4.0l) + P->cf__nc * (-7.0l * (1.0l + v_2) * V1 / 4.0l + 2.0l * L1V * vV1 + LV * (3.0l - 4.0l * v - 3.0l * v_2) * V1 / 4.0l + (pi2) * (2.0l - 6.0l * v + 5.0l * v_2) * V1 / 6.0l + (L1V * L1V) * V1_3 + (LV * LV) * V1 * (3.0l * v_2 + 2.0l * V1) / 2.0l - L1V * LV * V1 * (v_2 + V1_2)) + (1.0l + v_2) * LMS1 * (-11.0l / 12.0l + P->nf / (6.0l * P->nc) - P->cf__nc * (3.0l - 4.0l * L1V + 4.0l * LV) / 4.0l) * V1;

    //*************************************************//
    //1/(1-W)_PLUS CONTRIB. TO GLUON Q TO GAMMA Q GLUON//
    //  (formerly FGQ2(v,W) and FGW2(V,1))      874    //
    //*************************************************//

    //only considering IPOL == 0 case
    temp1 = (1.0l + v_2) * (-L1V - P->cf__nc * (3.0l - 4.0l * LV) / 4.0l + LV) * V1;
    temp2 = -1.0l * (P->cf + P->nc) * (1.0l + v_2) * V1 / P->nc;

    GQ2_W = temp1 + temp2 * LMS;
    GQ2_W1 = temp1 + temp2 * LMS1;

    //*********************************************************//
    //(LN(1-W)/(1-W))_PLUS CONTRIB. TO GLUON Q TO GAMMA Q GLUON//
    //  (formerly GQ3(V,W) and GQ3(V,1))        901            //
    //*********************************************************//

    //only considering IPOL == 0 case
    GQ3 = (P->cf + 2.0l * P->nc) * (1.0l + v_2) * V1 / P->nc;

#if _DEBUG_DIRECT_
    std::cout << "\n GQ1 = " << GQ1;
    std::cout << "\n GQ2_W = " << GQ2_W;
    std::cout << "\n GQ2_W1 = " << GQ2_W1;
    std::cout << "\n GQ3 = " << GQ3;
    std::cout << std::endl;
#endif

    FUNCHO1 = (SIGQQB + alpha_s * (QQB1 + LAA * QQB2_W1 + 0.5l * LAA2 * QQB3)) * PQQB0 + (SIGQG + alpha_s * (QG1 + LAA * QG2_W1 + 0.5l * LAA2 * QG3)) * PQG0 + (SIGGQ + alpha_s * (GQ1 + LAA * GQ2_W1 + 0.5l * LAA2 * GQ3)) * PGQ0;
    FUNCHO1 *= Jv;

#if _DEBUG_DIRECT_
    std::cout << "\n FUNCHO1 = " << FUNCHO1;
    std::cout << std::endl;
#endif

    //*************************************************//
    //REMAINING CONTRIB. TO Q QBAR TO GAMMA GLUON GLUON//
    //  (formerly 'FQQB(V, W)')         569            //
    //*************************************************//

    //skipping IPOL == 1 case
    GQQB = (P->cf2__nc * LMS * ((1.0l - 2.0l * v) * (1.0l - 2.0l * v) - (1.0l - 2.0l * v) * v / X + vV1 / X_2 + (3.0l * v_2 + V1_2) * w)) * (1.0l - LV / LMS) + P->cf2__nc * L1V * (-2.0l * (2.0l - X) + 2.0l * (1.0l + V1_2) / X) + P->cf * L1V * ((1.0l + v) * V1 - (1.0l + V1_2) / X + vV1 * w) + P->cf * LW * (1.0l - v_2 - (1.0l + V1_2) / X + vV1 * w) + 2.0l * P->cf2__nc * LW * (V1_2 + v_2 * w_2) / X + P->cf * TQQB * LW / W1 + P->cf * L1W * (-1.0l - 2.0l * vV1 + (1.0l + V1_2) / X - (1.0l - 2.0l * v) * vw) + P->cf2__nc * L1W * (1.0l + 8.0l * vV1 - vV1 / X_2 - (3.0l * v + 4.0l * V1_2) / X - (1.0l - 4.0l * v + 8.0l * v_2) * w) - 4.0l * P->cf2__nc * LLVW * v * (V1 - vw) + P->cf * LLVW * v * (V2 - vw) + P->cf__nc * (4.0l * P->cf - P->nc) * TQQB * (L1V - LLVW) / W1 + P->cf * P->nf * v * (V2 - vw) / (3.0l * P->nc) + P->cf2__nc * (-4.0l * vV1 + vV1 / X_2 + vV1 / X + (1.0l + v - 4.0l * v_2) * w) + P->cf * (-(1.0l - 4.0l * v) * v / (2 * X) + vV1 / (2.0l * X * X) - (2.0l - 11.0l * V1_2) / 6.0l + (3.0l - 12.0l * v + 11.0l * v_2) * w / 6.0l);

#if _DEBUG_DIRECT_
    std::cout << "\n GQQB = " << GQQB;
#endif

    //**********************************************//
    //REMAINING CONTRIB. TO Q GLUON TO GAMMA Q GLUON//
    //  (formerly FQG(V,w)      709                 //
    //**********************************************//

    //only considering IPOL == 0 case
    temp1 = P->cf__nc * LV * v * ((3.0l * v_2 - 10.0l * V1) / 2.0l + (3.0l - v) * (3.0l - v) * V1 / (2.0l * Y) + V1_2 / (X_2 * X) - (3.0l - v) * V1_2 / Y_2 + V1_3 / (Y_2 * Y) - V1 * (1.0l + 2.0l * V1) / X_2 + (1.0l + 2.0l * V1) * (1.0l + 2.0l * V1) / (2.0l * X) - (11.0l * v_2 + 2.0l * V1) * w / 2.0l + 4.0l * v_2 * w_2) + LV * v * (-(4.0l - 3.0l * v + v_2) * V1 + (4.0l - v) * V1_2 / Y - V1_3 / (X_2 * X) - V1_3 / Y_2 - V1 * V2 / X + V1_2 * V2 / X_2 + v * (7.0l + v_2) * w / 2.0l - 4.0l * v_2 * w_2 + 2.0l * v_3 * w_3) / V1;

    QG_VW = P->cf__nc * LMS * (-(v * V1_2 / (X_2 * X)) + vV1 * (1.0l + 2.0l * V1) / X_2 - v * (1.0l + 2.0l * V1) * (1.0l + 2.0l * V1) / (2.0l * X) + v * (v_2 + 6.0l * V1) / 2.0l + v * (2.0l + v_2) * w / 2.0l) + LMS * (-v + v * (4.0l - 3.0l * v + v_2) / X + v * V1_2 / (X_2 * X) - vV1 * V2 / X_2 - v_2 * (2.0l * v_2 + 3.0l * V1) * w / V1 + v_3 * (1.0l + v) * w_2 / V1 - (v_3 * v) * w_3 / V1) - LMSS * v * (1.0l + v_2 * W1 * W1) * (P->cf__nc * V1_2 + Y * vw) * (V1_2 + 2.0l * Y * vw) / (2.0l * (Y_2 * Y) * V1) + temp1 - P->cf__nc * L1V * v_3 * (1.0l + w) - L1V * v * ((1.0l + v) * V1_2 - 4.0l * vw * (2.0l * X - vV1 + Y * vw)) / (2.0l * V1) - P->cf__nc * LW * v_3 * (1.0l + 4.0l * w - 2.0l * w_2) + LW * v * (2.0l * TQG / X + v_2 - 2.0l * V1 - vw * (3.0l - 4.0l * v + vw)) / 2.0l + (2.0l * P->cf - P->nc) * LW * v * (3.0l * v_2 + 2.0l * V1) / (2.0l * P->nc * W1) + temp1 * L1W / LV - P->cf__nc * L1W * v_3 * (1.0l - 2.0l * w + 2.0l * w_2) - L1W * v * (TQG / X - (4.0l - 3.0l * v) * vw / 2.0l - v_2 * w_2 / 2.0l) + P->cf__nc * LLVW * v_3 * (3.0l - 2.0l * w) * w + LLVW * v * ((5.0l - 4.0l * v + v_2) * V1 - v * (8.0l - 3.0l * vV1) * w + v_2 * w_2 * (7.0l + v - 4.0l * vw)) / (2.0l * V1) + (L1V - LLVW) / W1 * v * (5.0l - 4.0l * v + v_2) / 2.0l + P->cf__nc * (L1V - LLVW) / W1 * v * (v_2 + V1_2) + P->cf__nc * L1VVW * v * (-4.0l * V1 + (3.0l - v) * (3.0l - v) * V1 / Y - 2.0l * (3.0l - v) * V1_2 / Y_2 + 2.0l * V1_3 / (Y_2 * Y) + v * (2.0l * X + v) * w) + L1VVW * v * (-11.0l / 2.0l + v + 2.0l * (4.0l - v) * V1 / Y - 2.0l * V1_2 / Y_2 + 3.0l * vw / 2.0l + v_2 * w_2 / 2.0l) - L1VVW / W1 * (1.0l - 2.0l * v) * v / 2.0l - P->cf__nc * L1VVW / W1 * v * (1.0l + v_2) + vw * (-TQG * V1 / 2.0l - 2.0l * v * V1_2 / Y_2 - (9.0l - 4.0l * v) * v * V1_2 / (2.0l * X_2) + 3.0l * v * V1_3 / (X_2 * X) + vV1 * (1.0l + 3.0l * V1) / (2.0l * X) + vV1 * V2 / Y + Y * v_2 * w) / V1 + P->cf__nc * (-v * (22.0l - 35.0l * v + 11.0l * v_2) / (2.0l * X) + (23.0l - 16.0l * v) * vV1 / (2.0l * X_2) - v * (13.0l - 3.0l * v - 3.0l * v_2) * V1 / (2.0l * Y) - 4.0l * v * V1_2 / (X_2 * X) + (17.0l - 5.0l * v) * v * V1_2 / (2.0l * Y_2) - 3.0l * v * V1_3 / (Y_2 * Y) + 3.0l * v * (v_2 + 10.0l * V1) / 4.0l + v * (4.0l - 7.0l * v_2) * w / 4.0l);

#if _DEBUG_DIRECT_
    std::cout << "\n QG_VW = " << QG_VW;
#endif

    //**********************************************//
    //REMAINING CONTRIB. TO Q GLUON TO GAMMA Q GLUON//
    //  (formerly FQG(1-V*W,(1-V)/(1-V*W)))  in 823 //
    //**********************************************//

    sv = X;
    sw = V1 / X;

    svw = sv * sw;

    sv2 = sv * sv;
    sw2 = sw * sw;
    sv3 = sv2 * sv;
    sw3 = sw2 * sw;

    sV1 = 1.0l - sv;
    sV12 = sV1 * sV1;
    sV2 = 2.0l - sv;

    sX = 1.0l - svw;
    sY = sV1 + svw;

    sLV = std::log(sv);
    sL1V = std::log(sV1);

    sLW = std::log(sw);
    sL1W = std::log(1.0l - sw);
    //L1W1 = -inf

    sLLVW = std::log(sX);
    sL1VVW = std::log(sY);

    sshatd = p2 / svw / sV1;
    sLMS = std::log(q2_fac / sshatd);
    sLMSS = std::log(q2_fragm / sshatd);

    sTQG = 2.0l - 2.0l * sv + sv2;

    //only considering IPOL == 0 case
    temp1 = P->cf__nc * sLV * sv * ((3.0l * sv2 - 10.0l * sV1) / 2.0l + (3.0l - sv) * (3.0l - sv) * sV1 / (2.0l * sY) + sV12 / (sX * sX * sX) - (3.0l - sv) * sV12 / (sY * sY) + (sV12 * sV1) / (sY * sY * sY) - sV1 * (1.0l + 2.0l * sV1) / (sX * sX) + (1.0l + 2.0l * sV1) * (1.0l + 2.0l * sV1) / (2.0l * sX) - (11.0l * sv2 + 2.0l * sV1) * sw / 2.0l + 4.0l * sv2 * sw2) + sLV * sv * (-(4.0l - 3.0l * sv + sv2) * sV1 + (4.0l - sv) * sV12 / sY - (sV12 * sV1) / (sX * sX * sX) - (sV12 * sV1) / (sY * sY) - sV1 * sV2 / sX + sV12 * sV2 / (sX * sX) + sv * (7.0l + sv2) * sw / 2.0l - 4.0l * sv2 * sw2 + 2.0l * sv3 * sw3) / sV1;

    QG_SWAP = P->cf__nc * sLMS * (-(sv * sV12 / (sX * sX * sX)) + sv * sV1 * (1.0l + 2.0l * sV1) / (sX * sX) - sv * (1.0l + 2.0l * sV1) * (1.0l + 2.0l * sV1) / (2.0l * sX) + sv * (sv2 + 6.0l * sV1) / 2.0l + sv * (2.0l + sv2) * sw / 2.0l) + sLMS * (-sv + sv * (4.0l - 3.0l * sv + sv2) / sX + sv * sV12 / (sX * sX * sX) - sv * sV1 * sV2 / (sX * sX) - sv2 * (2.0l * sv2 + 3.0l * sV1) * sw / sV1 + sv3 * (1.0l + sv) * sw2 / sV1 - sv * sv3 * sw3 / sV1) - sLMSS * sv * (1.0l + sv2 * (1.0l - sw) * (1.0l - sw)) * (P->cf__nc * sV12 + sY * svw) * (sV12 + 2.0l * sY * svw) / (2.0l * (sY * sY * sY) * sV1) + temp1 - P->cf__nc * sL1V * sv3 * (1.0l + sw) - sL1V * sv * ((1.0l + sv) * sV12 - 4.0l * svw * (2.0l * sX - sv * sV1 + sY * svw)) / (2.0l * sV1) - P->cf__nc * sLW * sv3 * (1.0l + 4.0l * sw - 2.0l * sw2) + sLW * sv * (2.0l * sTQG / sX + sv2 - 2.0l * sV1 - svw * (3.0l - 4.0l * sv + svw)) / 2.0l + (2.0l * P->cf - P->nc) * sLW * sv * (3.0l * sv2 + 2.0l * sV1) / (2.0l * P->nc * (1.0l - sw)) + temp1 * sL1W / sLV - P->cf__nc * sL1W * sv3 * (1.0l - 2.0l * sw + 2.0l * sw2) - sL1W * sv * (sTQG / sX - (4.0l - 3.0l * sv) * svw / 2.0l - sv2 * sw2 / 2.0l) + P->cf__nc * sLLVW * sv3 * (3.0l - 2.0l * sw) * sw + sLLVW * sv * ((5.0l - 4.0l * sv + sv2) * sV1 - sv * (8.0l - 3.0l * sv * sV1) * sw + sv2 * sw2 * (7.0l + sv - 4.0l * svw)) / (2.0l * sV1) + (sL1V - sLLVW) / (1.0l - sw) * sv * (5.0l - 4.0l * sv + sv2) / 2.0l + P->cf__nc * (sL1V - sLLVW) / (1.0l - sw) * sv * (sv2 + sV12) + P->cf__nc * sL1VVW * sv * (-4.0l * sV1 + (3.0l - sv) * (3.0l - sv) * sV1 / sY - 2.0l * (3.0l - sv) * sV12 / (sY * sY) + 2.0l * (sV12 * sV1) / (sY * sY * sY) + sv * (2.0l * sX + sv) * sw) + sL1VVW * sv * (-11.0l / 2.0l + sv + 2.0l * (4.0l - sv) * sV1 / sY - 2.0l * sV12 / (sY * sY) + 3.0l * svw / 2.0l + sv2 * sw2 / 2.0l) - sL1VVW / (1.0l - sw) * (1.0l - 2.0l * sv) * sv / 2.0l - P->cf__nc * sL1VVW / (1.0l - sw) * sv * (1.0l + sv2) + svw * (-sTQG * sV1 / 2.0l - 2.0l * sv * sV12 / (sY * sY) - (9.0l - 4.0l * sv) * sv * sV12 / (2.0l * (sX * sX)) + 3.0l * sv * (sV12 * sV1) / (sX * sX * sX) + sv * sV1 * (1.0l + 3.0l * sV1) / (2.0l * sX) + sv * sV1 * sV2 / sY + sY * sv2 * sw) / sV1 + P->cf__nc * (-sv * (22.0l - 35.0l * sv + 11.0l * sv2) / (2.0l * sX) + (23.0l - 16.0l * sv) * sv * sV1 / (2.0l * (sX * sX)) - sv * (13.0l - 3.0l * sv - 3.0l * sv2) * sV1 / (2.0l * sY) - 4.0l * sv * sV12 / (sX * sX * sX) + (17.0l - 5.0l * sv) * sv * sV12 / (2.0l * (sY * sY)) - 3.0l * sv * (sV12 * sV1) / (sY * sY * sY) + 3.0l * sv * (sv2 + 10.0l * sV1) / 4.0l + sv * (4.0l - 7.0l * sv2) * sw / 4.0l);

#if _DEBUG_DIRECT_
    std::cout << "\n QG_SWAP = " << QG_SWAP;
#endif

    //**********************************************//
    //REMAINING CONTRIB. TO GLUON Q TO GAMMA Q GLUON//
    //  (formerly FGQ(v,W))         926             //
    //**********************************************//

    //only considering IPOL == 0 case
    GQ = ((1.0l + v_2) * (L1V - LLVW) / W1 * V1 - (1.0l + v_2) * LW * V1 / W1 + (-3.0l * P->cf + 4.0l * P->cf * L1W + 8.0l * P->nc * L1W - 4.0l * P->nc * LLVW + 4.0l * P->cf * LV + 4.0l * P->nc * LV - 4.0l * P->nc * LW) * v * (V1 + v_2 - vw + v_2 * w + v_2 * w_2) / (4.0l * P->nc)) + (-(P->cf + P->nc) * LMS * v * (V1 + v_2 - vw + v_2 * w + v_2 * w_2) / P->nc) + v / (1.0l - v * w) * QG_SWAP;

#if _DEBUG_DIRECT_
    std::cout << "\n GQ = " << GQ;
#endif

    //***************************************//
    //CONTRIB. TO GLUON GLUON TO GAMMA Q QBAR//
    //  (formerly GG(V,W))      963          //
    //***************************************//

    //only considering IPOL == 0 case

    GG = 1.0l / P->nc * LMS * v * (-(7.0l - 8.0l * v + 3.0l * v_2) / 2.0l - 2.0l * V1_2 / X_2 + 2.0l * V1 * V2 / X + (2.0l - 3.0l * v + 2.0l * v_2) * w - (3.0l * v_2 + 4.0l * V1) * w_2 / 2.0l) - LMSS * v * (X_2 + (V2 * V2) - 4.0l * V2 * V1 / Y + 4.0l * V1_2 / Y_2) * (1.0l / P->nc - V1 * vw / (P->cf * Y_2)) / 2.0l + 1.0l / P->cf * LV * v_2 * V1 * (3.0l - 4.0l * V1_2 / (Y_2 * Y_2) + 4.0l * V1 / (Y * V2) - 2.0l * (1.0l + V1_2) / (X * V2) + 4.0l * V1 * V2 / (Y_2 * Y) - 2.0l * (V2 * V2) / Y_2) * w / 2.0l + 2.0l * 1.0l / P->nc * LV * v * (V1_2 / X_2 + V1_2 / Y_2 - 2.0l * V1_2 / (X * V2) - 2.0l * V1_2 / (Y * V2) + (v_2 + V1) * (W1 + w_2)) + 1.0l / P->cf * L1V * v_2 * V1 * w * (1.0l - 2.0l * w + 2.0l * (1.0l + V1_2) * w / (X * Y)) / 2.0l + 1.0l / P->nc * L1V * v * ((4.0l - 3.0l * v + v_3) / v - 2.0l * V1 * (1.0l + V1_2) / (X * Y * v) - 2.0l * v * W1 * w) + 1.0l / P->cf * LW * vV1 * (1.0l - 3.0l * V1 + 2.0l * V1 * (1.0l + V1_2) / (X * Y)) * w / 2.0l + 1.0l / P->nc * LW * v * (-2.0l * V1 * W1 + 2.0l * V1 * (1.0l + V1_2) * W1 / (X * Y) + v_2 * w_2) + 1.0l / P->cf * L1W * vV1 * (1.0l - 2.0l * v * V1_2 / (Y_2 * Y_2) - (1.0l + V1_2) / (X * V2) + 2.0l * vV1 * V2 / (Y_2 * Y) - v * (V2 * V2) / Y_2 - V1 * (2.0l * V1 - v * V2) / (Y * V2)) * w + 1.0l / P->nc * L1W * v * (2.0l * V1_2 / X_2 + 2.0l * V1_2 / Y_2 + 3.0l * (1.0l + V1_2) - 2.0l * (2.0l - v_2) * V1 / (Y * v * V2) + 2.0l * V1_2 * (2.0l * V1 - v * V2) / (X * v * V2) - 4.0l * (v_2 + V1) * w + (3.0l * v_2 + 2.0l * V1) * w_2) - 1.0l / P->nc * LLVW * v * (v_2 + V1_2 + vw * (2.0l * V1 + vw)) + 1.0l / P->cf * L1VVW * v_2 * V1 * (-3.0l - 4.0l * V1_2 / (Y_2 * Y_2) + 2.0l * V2 / Y + 4.0l * V1 * V2 / (Y_2 * Y) - 2.0l * (V2 * V2) / Y_2) * w - 4.0l * 1.0l / P->nc * L1VVW * v_3 * V1 * W1 * w / Y_2 + 1.0l / P->cf * vV1 * w * (-2.0l * v / Y + 2.0l * (7.0l - 6.0l * v) * v / Y_2 - 4.0l * vV1 / X_2 + 12.0l * v * V1_2 / (Y_2 * Y_2) + V2 + 2.0l * v * V2 / X - 12.0l * vV1 * V2 / (Y_2 * Y) - 2.0l * V2 * w) / 2.0l + 1.0l / P->nc * v * (-(3.0l - v) * V1 - 4.0l * V1_2 / X_2 + 4.0l * V1 * V2 / X + 2.0l * (2.0l - 3.0l * v + 2.0l * v_2) * w - (V2 * V2) * w_2) / 2.0l;

#if _DEBUG_DIRECT_
    std::cout << "\nGG = " << GG;
#endif

    //****************************//
    //CONTRIB. TO Q Q TO GAMMA Q Q//
    //  (formerly FQQ(V,W)) 1058  //
    //****************************//

    //only considering IPOL == 0 case
    temp1 = P->cf / P->nc2 * LV * v * (-(1.0l + v) * (1.0l + v) + 4.0l * V1 / (Y * V2) - 2.0l * (1.0l + V1_2) / (X * V2) + 2.0l * v * (1.0l + v) * w - v_2 * w_2) + P->cf__nc * LV * ((43.0l - 33.0l * v + 12.0l * v_2 - 4.0l * v_3) * V1 + V1_2 / X_2 - 8.0l * (5.0l - 4.0l * v + v_2) * V1_2 / Y + 4.0l * V1_4 / Y_2 - 2.0l * V1 * V2 / X - 2.0l * V1 * (1.0l + 11.0l * v + v_2 * (1.0l + 2.0l * V1)) * w + (22.0l * v_2 - 18.0l * v_3 + 4.0l * (v_3 * v) + V1_2) * w_2 - 4.0l * (2.0l - Y) * v_3 * w_3) / (2.0l * V1 * w);

    QQ = P->cf__nc * LMS * (-(V1_2 / X_2) - V1 * (5.0l - 3.0l * v + 2.0l * v_2 * V1) + 2.0l * V1 * V2 / X - 2.0l * V1 * (v_3 - V1_2) * w - (4.0l * v_2 - 2.0l * v_3 * V1 + V1_2) * w_2 + 2.0l * (2.0l - Y) * v_3 * w_3) / (2.0l * V1 * w) + LMSS * (1.0l + v_2 * W1 * W1) * (P->cf / P->nc2 * Y_2 * vV1 * w - P->cf__nc * (V1_4 + v * V1_3 * w + v_2 * V1_2 * w_2 + Y * v_3 * w_3)) / (Y_2 * V1 * w) + temp1 + P->cf / P->nc2 * L1V * (V1 - 4.0l * V1 / (Y * V2) + v * (1.0l + V1_2) / (X * V2) - v_3 * W1 * W1 - vw) + P->cf__nc * L1V * v * (2.0l * (1.0l - 2.0l / Y) * V1_2 + v * (3.0l + v_2) * w - 2.0l * (2.0l - Y) * v_2 * w_2) / V1 + P->cf__nc * LW * (2.0l * (3.0l + v_2) * V1 - 4.0l * V1_2 / Y - 2.0l * v * (1.0l + v) * V2 * w + (3.0l - v) * v_2 * w_2) / w + P->cf / P->nc2 * LW * (-(1.0l + v + v_3) + 4.0l * V1 / (Y * V2) - v * (1.0l + V1_2) / (X * V2) + v * (1.0l + 2.0l * v_2) * w - v_3 * w_2) + temp1 * L1W / LV - P->cf__nc * L1W * (X - v) * v * (X + v) / Y - P->cf / P->nc2 * L1W * (1.0l + v_2 * W1 * W1) * (V1_2 - v * (1.0l + v_2 * W1) * w) / (X * Y) + P->cf__nc * LLVW * v * ((1.0l + v) * V1 - v * (4.0l - vV1) * w + 2.0l * (2.0l - Y) * v_2 * w_2) / V1 + 4.0l * P->cf__nc * L1VVW * v_3 * V1 * W1 * w / Y_2 - P->cf / P->nc2 * Y_2 * v + P->cf__nc * (-(27.0l - 30.0l * v + 2.0l * v_2 + 2.0l * v_3) * V1 - V1_2 / X_2 + 8.0l * (5.0l - v) * V1_3 / Y - 16.0l * V1_4 / Y_2 + V1 * (1.0l + 5.0l * V1) / X + V1 * (1.0l + 2.0l * (3.0l - v) * vV1) * w - V1 * (1.0l - 4.0l * v * (1.0l + v) * V1) * w_2 + 2.0l * Y * v_3 * w_3) / (2.0l * V1 * w);

#if _DEBUG_DIRECT_
    std::cout << "\nQQ = " << QQ;
#endif

    //**********************************//
    //CONTRIB. TO Q QBAR TO GAMMA Q QBAR//
    //  (formerly FQB(V,W))     1153    //
    //**********************************//

    //only considering IPOL == 0 case
    temp1 = P->cf / P->nc2 * LV * v * (-2.0l - 2.0l * (1.0l + v - v_2) * V1 / Y - V1_2 + 6.0l * V1_2 / Y_2 + v * (5.0l + v) * w - 3.0l * v_2 * w_2) + P->cf__nc * LV * ((11.0l - 53.0l * v + 40.0l * v_2 - 12.0l * v_3) * V1 + V1_2 / X_2 - 4.0l * (10.0l - 27.0l * v + 16.0l * v_2 - 3.0l * v_3) * V1_2 / Y + 4.0l * V3 * V5 * V1_4 / Y_2 - 8.0l * V4 * V1_5 / (Y_2 * Y) + 8.0l * (V1_5 * V1) / (Y_2 * Y_2) - 2.0l * V1 * V2 / X - 2.0l * (1.0l - 5.0l * v - 7.0l * v_2) * V1 * w + (2.0l * v_3 + 6.0l * (v_3 * v) + V1_2) * w_2 - 4.0l * (2.0l - Y) * v_3 * w_3) / (2.0l * V1 * w);

    QB = P->cf__nc * LMS * (-V1_2 / X_2 - V1 * (5.0l - 3.0l * v + 2.0l * v_2 * V1) + 2.0l * V1 * V2 / X - 2.0l * V1 * (v_3 - V1_2) * w - (4.0l * v_2 - 2.0l * v_3 * V1 + V1_2) * w_2 + 2.0l * (2.0l - Y) * v_3 * w_3) / (2.0l * V1 * w) + LMSS * (1.0l + v_2 * W1 * W1) * (-P->cf / P->nc2 * Y_2 * vV1 * w * (V1_2 - vV1 * w + v_2 * w_2) - P->cf__nc * ((V1_5 * V1) + 3.0l * v * V1_5 * w + 5.0l * v_2 * V1_4 * w_2 + 4.0l * v_3 * V1_3 * w_3 + 5.0l * (v_3 * v) * V1_2 * w * w_3 + 3.0l * v_3 * v * vV1 * w_2 * w_3 + v_3 * v_3 * (w_3 * w_3))) / ((Y_2 * Y_2) * V1 * w) + temp1 + P->cf / P->nc2 * L1V * v * (4.0l - 3.0l * v - 2.0l * V1 * (1.0l + v_2 + 2.0l * V1) / Y - v * (X + 2.0l * V1) * w) + P->cf__nc * L1V * v * (-2.0l * V1_2 + 4.0l * V1_2 / Y + v * (3.0l + v_2) * w - 2.0l * (2.0l - Y) * v_2 * w_2) / V1 + P->cf / P->nc2 * LW * v * (-(5.0l - 4.0l * v + v_2) + 2.0l * V1 * (1.0l + v_2 + 2.0l * V1) / Y + 2.0l * (1.0l + X) * vw) + P->cf__nc * LW * (4.0l * V1_2 / Y - 2.0l * (1.0l + v) * V1_2 + 2.0l * v * (v_2 + V2) * w - v_2 * (1.0l + v) * w_2) / w - P->cf / P->nc2 * LW * (1.0l + v_2) / W1 + temp1 * L1W / LV + L1W * (P->cf__nc * (2.0l - Y) * v * (Y - 2.0l * vw) / Y - P->cf / P->nc2 * v_2 * (Y - 2.0l * vw) * (2.0l * v - w - 3.0l * vw + 2.0l * v * w_2) / Y) + P->cf / P->nc2 * LLVW * v * (V2 - (3.0l - Y) * vw) - P->cf__nc * LLVW * v * ((1.0l + v) * V1 + 2.0l * X * vw + v_2 * w * (1.0l + X - Y + 2.0l * v * w_2)) / V1 - P->cf / P->nc2 * (L1V - LLVW) / W1 * (1.0l + V1_2) + 2.0l * P->cf__nc * L1VVW * v * (-2.0l * (4.0l - 3.0l * v + v_2) - 2.0l * V3 * V3 * V1_2 / Y_2 + 4.0l * V3 * V1_3 / (Y_2 * Y) - 4.0l * V1_4 / (Y_2 * Y_2) + 4.0l * V1 * (V2 * V2) / Y + V3 * vw) + P->cf / P->nc2 * L1VVW * v * (15.0l - 13.0l * v + 3.0l * v_2 + 12.0l * V1_2 / Y_2 - 12.0l * V1 * V2 / Y - 7.0l * vw + 3.0l * v_2 * w_2) + P->cf / P->nc2 * L1VVW / W1 * (3.0l - 2.0l * vV1) + P->cf / P->nc2 * (-3.0l - (5.0l - 2.0l * v) * v / X + vV1 / X_2 + 4.0l * V3 * vV1 / Y + 2.0l * v * V1_2 - 8.0l * v * V1_2 / Y_2 + (v_2 + 2.0l * v_3 + V1_2) * w + 2.0l * v_3 * w_2) / 2.0l + P->cf__nc * (-(31.0l - 54.0l * v + 22.0l * v_2 + 2.0l * v_3) * V1 - V1_2 / X_2 + 4.0l * (27.0l - 19.0l * v) * V1_3 / Y - 4.0l * (37.0l - 18.0l * v) * V1_4 / Y_2 + 24.0l * V4 * V1_5 / (Y_2 * Y) - 24.0l * (V1_5 * V1) / (Y_2 * Y_2) + V1 * V2 / X + V1 * (1.0l + 2.0l * V3 * vV1) * w - V1 * (1.0l - 4.0l * v_2 * V1) * w_2 + 2.0l * Y * v_3 * w_3) / (2.0l * V1 * w);

#if _DEBUG_DIRECT_
    std::cout << "\n QB = " << QB;
#endif

    //********************************************//
    //CONTRIB. TO Q QBAR TO GAMMA QPRIME QBARPRIME//
    //  (formerly FQBS(v,W))        1281          //
    //********************************************//

    //only considering IPOL == 0 case
    QBS = (P->cf__nc * v_2 * V1 * (-1.0l - 4.0l * V1_2 / (Y_2 * Y_2) + 2.0l * V2 / Y + 4.0l * V1 * V2 / (Y_2 * Y) - 2.0l * (V2 * V2) / Y_2) * w) * (LMSS - LV - L1W - 2.0l * L1VVW) + 2.0l * P->cf__nc * v_3 * V1 * W1 * w * (-1.0l + 6.0l * vV1 * w / Y_2) / Y_2;

#if _DEBUG_DIRECT_
    std::cout << "\n QBS = " << QBS;
    std::cout << std::endl;
#endif

    //**************************************//
    //CONTRIB. TO Q QPRIME TO GAMMA Q QPRIME//
    // (formerly FQS(V,W,CQ=0/1,CQS=1/0)) 1309//
    //**************************************//

    //only considering IPOL == 0 case
    temp1 = 2.0l * P->cf__nc * LV * v * (-(5.0l + v) + (7.0l - 4.0l * v + v_2) / Y + 3.0l * vw) + P->cf__nc * LV * v_2 * w * (-2.0l * V1 / X - 2.0l * (3.0l - v) * V1 / Y + V1_2 / X_2 + 2.0l * V1_2 / Y_2 + 2.0l * (4.0l - vV1) - 4.0l * (2.0l - Y) * vw) / (2.0l * V1) + P->cf__nc * LV * V1 * (2.0l * (2.0l * v_2 + V2) + 2.0l * v * V1_2 * W1 / Y_2 - (1.0l + 2.0l * v_2) * (2.0l - w) * w) / (2.0l * w);

    QS_1_1 = P->cf__nc * LMS * v_2 * w * (-(3.0l + v_2) / 2.0l + V1 / X - V1_2 / (2.0l * X_2) + (2.0l - Y) * vw) / V1 - P->cf__nc * LMS * (1.0l + v_2) * V1 * (2.0l - 2.0l * w + w_2) / (2.0l * w) + P->cf__nc * LMSS * (1.0l + v_2 * W1 * W1) * (-(v_2 * w_2 * (V1_2 + 2.0l * Y * vw)) - V1_2 * (2.0l * V1_2 + vw * (2.0l * Y - vw))) / (2.0l * Y_2 * V1 * w) + temp1 + 2.0l * P->cf__nc * (1.0l - 2.0l / Y) * L1V * vV1 + P->cf__nc * L1V * v_2 * w * (3.0l + v_2 - 2.0l * (2.0l - Y) * vw) / V1 - 2.0l * P->cf__nc * (2.0l - Y) * LW * v_2 * w / Y + P->cf__nc * LW * V1 * (2.0l + v_2 * (2.0l - 2.0l * w + w_2)) / w + temp1 * L1W / LV - P->cf__nc * (2.0l - Y) * L1W * v * (V1 - vw) / Y + P->cf__nc * (2.0l - Y) * LLVW * v - P->cf__nc * LLVW * v_2 * w * (3.0l + v_2 - 2.0l * (2.0l - Y) * vw) / V1 + 2.0l * P->cf__nc * (2.0l - Y) * L1VVW * v + 2.0l * P->cf__nc * L1VVW * vV1 * (-v / Y - V1 / Y_2) + 2.0l * P->cf__nc * L1VVW * v_2 * (1.0l - V3 / Y + V1 / Y_2) * w + P->cf__nc * (2.0l - Y) * v * (1.0l / X - 2.0l * v / Y_2) * V1 * w + P->cf__nc * v_2 * w * ((1.0l - 2.0l * v) * V1 - V1_2 / X_2 - 4.0l * V1_2 / Y_2 + V1 * V2 / X + 2.0l * V1 * V2 / Y + 2.0l * Y * vw) / (2.0l * V1) + P->cf__nc * V1 * (2.0l * V4 * V1 / Y - 4.0l * V1_2 / Y_2 - 2.0l * (1.0l + vV1) + (1.0l + 3.0l * v - 2.0l * v_2) * w - (1.0l - 2.0l * v) * (1.0l + v) * w_2) / (2.0l * w);

#if _DEBUG_DIRECT_
    std::cout << "\n QS_1_1 = " << QS_1_1;
#endif

    //only considering IPOL == 0 case

    temp1 = P->cf__nc * LV * V1 * (2.0l * (2.0l * v_2 + V2) + 2.0l * v * V1_2 * W1 / Y_2 - (1.0l + 2.0l * v_2) * (2.0l - w) * w) / (2.0l * w);

    QS_0_1 = -P->cf__nc * LMS * (1.0l + v_2) * V1 * (2.0l - 2.0l * w + w_2) / (2.0l * w) + P->cf__nc * LMSS * (1.0l + v_2 * W1 * W1) * (-V1_2 * (2.0l * V1_2 + vw * (2.0l * Y - vw))) / (2.0l * Y_2 * V1 * w) + temp1 + P->cf__nc * LW * V1 * (2.0l + v_2 * (2.0l - 2.0l * w + w_2)) / w + temp1 * L1W / LV + 2.0l * P->cf__nc * L1VVW * vV1 * (-v / Y - V1 / Y_2) + P->cf__nc * V1 * (2.0l * V4 * V1 / Y - 4.0l * V1_2 / Y_2 - 2.0l * (1.0l + vV1) + (1.0l + 3.0l * v - 2.0l * v_2) * w - (1.0l - 2.0l * v) * (1.0l + v) * w_2) / (2.0l * w);

#if _DEBUG_DIRECT_
    std::cout << "\n QS_0_1 = " << QS_0_1;
#endif

    //only considering IPOL == 0 case
    temp1 = P->cf__nc * LV * v_2 * w * (-2.0l * V1 / X - 2.0l * (3.0l - v) * V1 / Y + V1_2 / X_2 + 2.0l * V1_2 / Y_2 + 2.0l * (4.0l - vV1) - 4.0l * (2.0l - Y) * vw) / (2.0l * V1);

    QS_1_0 = P->cf__nc * LMS * v_2 * w * (-(3.0l + v_2) / 2.0l + V1 / X - V1_2 / (2.0l * X_2) + (2.0l - Y) * vw) / V1 + P->cf__nc * LMSS * (1.0l + v_2 * W1 * W1) * (-(v_2 * w_2 * (V1_2 + 2.0l * Y * vw))) / (2.0l * Y_2 * V1 * w) + temp1 + P->cf__nc * L1V * v_2 * w * (3.0l + v_2 - 2.0l * (2.0l - Y) * vw) / V1 + temp1 * L1W / LV - P->cf__nc * LLVW * v_2 * w * (3.0l + v_2 - 2.0l * (2.0l - Y) * vw) / V1 + 2.0l * P->cf__nc * L1VVW * v_2 * (1.0l - V3 / Y + V1 / Y_2) * w + P->cf__nc * v_2 * w * ((1.0l - 2.0l * v) * V1 - V1_2 / X_2 - 4.0l * V1_2 / Y_2 + V1 * V2 / X + 2.0l * V1 * V2 / Y + 2.0l * Y * vw) / (2.0l * V1);

#if _DEBUG_DIRECT_
    std::cout << "\n QS_1_0 = " << QS_1_0;
#endif

    GQS1 = QS_1_0;
    GQS2 = QS_0_1;
    GQS3 = QS_1_1 - GQS1 - GQS2;

    GS1 = 4.0l / 9.0l * GQS1 - 2.0l / 9.0l * GQS3 + 1.0l / 9.0l * GQS2;
    GS2 = 4.0l / 9.0l * GQS1 + 2.0l / 9.0l * GQS3 + 1.0l / 9.0l * GQS2;
    GS3 = 1.0l / 9.0l * GQS1 - 2.0l / 9.0l * GQS3 + 4.0l / 9.0l * GQS2;
    GS4 = 1.0l / 9.0l * GQS1 + 2.0l / 9.0l * GQS3 + 4.0l / 9.0l * GQS2;
    GS5 = 1.0l / 9.0l * GQS1 + 1.0l / 9.0l * GQS3 + 1.0l / 9.0l * GQS2;
    GS6 = 1.0l / 9.0l * GQS1 - 1.0l / 9.0l * GQS3 + 1.0l / 9.0l * GQS2;
    GS7 = 4.0l / 9.0l * GQS1 + 4.0l / 9.0l * GQS3 + 4.0l / 9.0l * GQS2;
    GS8 = 4.0l / 9.0l * GQS1 - 4.0l / 9.0l * GQS3 + 4.0l / 9.0l * GQS2;

#if _DEBUG_DIRECT_
    std::cout << "\n GS1 = " << GS1;
    std::cout << "\n GS2 = " << GS2;
    std::cout << "\n GS3 = " << GS3;
    std::cout << "\n GS4 = " << GS4;
    std::cout << "\n GS5 = " << GS5;
    std::cout << "\n GS6 = " << GS6;
    std::cout << "\n GS7 = " << GS7;
    std::cout << "\n GS8 = " << GS8;
    std::cout << std::endl;
#endif
    //************************************//
    //COMBINATIONS OF PARTON DISTRIBUTIONS//
    //************************************//

    //Q QBAR TO GAMMA GLUON GLUON
    PQQB = 4.0l / 9.0l * (UP1 * UPB2 + UPB1 * UP2) + 1.0l / 9.0l * (DO1 * DOB2 + DOB1 * DO2) + 1.0l / 9.0l * (ST1 * ST2 + ST2 * ST1) + 4.0l / 9.0l * (CH1 * CH2 + CH2 * CH1);
    //Q GLUON TO GAMMA Q GLUON
    PQG = (4.0l / 9.0l * (UP1 + UPB1) + 1.0l / 9.0l * (DO1 + DOB1) + 1.0l / 9.0l * (ST1 + ST1) + 4.0l / 9.0l * (CH1 + CH1)) * GL2;
    //GLUON Q TO GAMMA Q GLUON
    PGQ = (4.0l / 9.0l * (UP2 + UPB2) + 1.0l / 9.0l * (DO2 + DOB2) + 1.0l / 9.0l * (ST2 + ST2) + 4.0l / 9.0l * (CH2 + CH2)) * GL1;
    //GLUON GLUON TO GAMMA Q QBAR (INCLUDING CHARM IN THE FINAL STATE)
    PGG = GL1 * GL2 * (4.0l / 9.0l + 1.0l / 9.0l + 1.0l / 9.0l + 4.0l / 9.0l);
    //C+1.0l/9.0l) //?? looks lost
    //Q Q TO GAMMA Q Q
    PQQ = 4.0l / 9.0l * (UP1 * UP2 + UPB1 * UPB2) + 1.0l / 9.0l * (DO1 * DO2 + DOB1 * DOB2) + 1.0l / 9.0l * (ST1 * ST2 + ST1 * ST2) + 4.0l / 9.0l * (CH1 * CH2 + CH2 * CH1);
    //Q QBAR TO GAMMA Q QBAR
    PQB = 4.0l / 9.0l * (UP1 * UPB2 + UPB1 * UP2) + 1.0l / 9.0l * (DO1 * DOB2 + DOB1 * DO2) + 1.0l / 9.0l * (ST1 * ST2 + ST2 * ST1) + 4.0l / 9.0l * (CH1 * CH2 + CH2 * CH1);
    //Q QBAR TO GAMMA QPRIME QBARPRIME (INCL. CHARM IN THE FINAL STATE)
    PQBS = (UP1 * UPB2 + UPB1 * UP2) * (1.0l / 9.0l + 1.0l / 9.0l + 4.0l / 9.0l) + (DO1 * DOB2 + DOB1 * DO2) * (1.0l / 9.0l + 4.0l / 9.0l + 4.0l / 9.0l) + (ST1 * ST2 + ST2 * ST1) * (4.0l / 9.0l + 1.0l / 9.0l + 4.0l / 9.0l) + (CH1 * CH2 + CH2 * CH1) * (4.0l / 9.0l + 1.0l / 9.0l + 1.0l / 9.0l);
    //PQBS = (UP1*UPB2+UPB1*UP2)*(1.0l/9.0l+1.0l/9.0l+4.0l/9.0l+1.0l/9.0l)+            (DO1*DOB2+DOB1*DO2)*(1.0l/9.0l+4.0l/9.0l+4.0l/9.0l+1.0l/9.0l)+            (ST1*ST2+ST2*ST1)*(4.0l/9.0l+1.0l/9.0l+4.0l/9.0l+1.0l/9.0l)+            (CH1*CH2+CH2*CH1)*(4.0l/9.0l+1.0l/9.0l+1.0l/9.0l+1.0l/9.0l);
    //Q Q(BAR)PRIME TO GAMMA Q Q(BAR)PRIME
    PQS = (UP1 * DO2 + UPB1 * DOB2) * GS1 + (UPB1 * DO2 + UP1 * DOB2) * GS2 + (DO1 * UP2 + DOB1 * UPB2) * GS3 + (DOB1 * UP2 + DO1 * UPB2) * GS4 + (UP1 * ST2 + UPB1 * ST2 + CH1 * DO2 + CH1 * DOB2 + 2.0l * CH1 * ST2) * (GS1 + GS2) + (ST1 * UP2 + ST1 * UPB2 + DO1 * CH2 + DOB1 * CH2 + 2.0l * ST1 * CH2) * (GS3 + GS4) + (ST1 * DO2 + ST1 * DOB2 + DO1 * ST2 + DOB1 * ST2) * (GS5 + GS6) + (UP1 * CH2 + UPB1 * CH2 + CH1 * UP2 + CH1 * UPB2) * (GS7 + GS8);

#if _DEBUG_DIRECT_
    std::cout << "\n PQQB = " << PQQB;
    std::cout << "\n PQG = " << PQG;
    std::cout << "\n PGQ = " << PGQ;
    std::cout << "\n PGG = " << PGG;
    std::cout << "\n PQQ = " << PQQ;
    std::cout << "\n PQB = " << PQB;
    std::cout << "\n PQBS = " << PQBS;
    std::cout << "\n PQS = " << PQS;
    std::cout << std::endl;
#endif
    //*******************************//
    //INTEGRAND FOR W-DEPENDENT TERMS//
    //*******************************//
    FUNCHO2 = (QQB2_W * PQQB - QQB2_W1 * PQQB0) / W1 + QQB3 * (PQQB - PQQB0) * L1W / W1 + GQQB * PQQB + (QG2_W * PQG - QG2_W1 * PQG0) / W1 + (QG3 * PQG - QG3 * PQG0) * L1W / W1 + QG_VW * PQG + (GQ2_W * PGQ - GQ2_W1 * PGQ0) / W1 + GQ3 * (PGQ - PGQ0) * L1W / W1 + GQ * PGQ + GG * PGG + QQ * PQQ + QB * PQB + QBS * PQBS + PQS;
    FUNCHO2 *= alpha_s * Jv * Jw;

#if _DEBUG_DIRECT_
    std::cout << "\nFUNCHO2 = " << FUNCHO2;
    std::cout << std::endl;
#endif
    //******************************//
    // - - - - FINAL RESULT - - - - //
    //******************************//

    //factors 2*hbarc2 * alpha_em pulled out of integral
    //this gives d2sigma_dydp in f[0] & the first moment of it in f[1]!
    f[1] = alpha_s * (FUNCHO1 + FUNCHO2) * Jp * Jy / p2;
    if (isnanl(f[1]) or isinfl(f[1]))
        f[1] = 0.0l;
    f[0] = f[1] / p;

#if _DEBUG_DIRECT_
    std::cout << "\nFUNCHO = " << 2.0l * f[0] * hbarc2 * alpha_em;
    std::cout << std::endl;
#endif

    return 0;
}
