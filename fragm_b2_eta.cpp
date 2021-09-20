#include "photon.hpp"
#include "fragm.hpp"

/**
  * Integrand of eta decay background B2 cross sections
  * */
int fragm_b2_eta(
    const int *,
    const real x[], // x = (w,v,z,x,y,p) (y and p optional) (p->sqrt(p2 + m_t2) for mesons)
    const int *,
    real f[],       // f = (eta, eta*p)
    void *userdata)
{
    auto P{reinterpret_cast<const photon_cs_params *const>(userdata)};

    static real p, Jp, p2, xt;
    static real sqy, y_lo, y_up, y, Jy, ey;
    static real alpha_s;
    static real z, z_min, Jz;
    static real v_min, v_max, Jv, v;
    static real w_min, Jw, w;
    static real x_min, x_max, Jx, xx, Egamma, Eeta, beta;
    static real v_2, v_3, v_4, _vm1, _vm1_2, _vm1_3, _vm1_4, _1mv, _1mv_2, _1mv_3, _1mv_4, _1mw;
    static real X1, X2, BX1, BX2, BXJAC, XJAC;
    static real DPU, DPD, DPS, DPC, DPG;
    static real UP1, UPB1, DO1, DOB1, ST1, CH1, GL1;
    static real UP2, UPB2, DO2, DOB2, ST2, CH2, GL2;
    static real SHD;
    static real BORN[16]{};
    static real DMU[16]{};
    static real help1, help2, help3, help4;
    static real F01[16]{}, F02[16]{}, FHODEL1[16]{}, FHODEL2[16]{}, FHOREST1[16]{}, FHOREST2[16]{}, FDELMU1[16], FDELMU2[16], GRRT[16], GRRC[16], GPPT[16], GPPC[16];
    static real GHD, GHE;
    static unsigned short n;

    //*************//
    //INTEGRATING P//
    //*************//
    Jp = P->int_p ? (P->p_max - P->p_min) : 1.0l;
    p = P->int_p ? P->p_min + Jp * x[P->int_y ? 5 : 4] : P->p;
    p2 = p * p;

    //*************//
    //INTEGRATING Y//
    //*************//

    xt = 2.0l * p / P->sqrtS;
    sqy = std::log((1.0l+std::sqrt(1.0l-xt*xt))/xt);
    y_lo = std::max(-sqy, P->y_min);
    y_up = std::min(+sqy, P->y_max);

#if _DEBUG_FRAGM_
    std::cout.precision(12);
    std::cout << std::scientific << "\n SQS = " << P->sqrtS;
    std::cout << "\n p = " << p;
    std::cout << "\n Jp = " << Jp;
    std::cout << "\n xt = " << xt;
    std::cout << "\n sqy = " << sqy;
    std::cout << "\n y_lo = " << y_lo;
    std::cout << "\n y_up = " << y_up;
    std::cout << std::endl;
#endif

    if(y_up < y_lo)
    {
        f[0] = f[1] = 0.0l;
        return 0;
    }

    Jy = P->int_y ? (P->y_max - P->y_min) : 1.0l;
    y = P->int_y ? (P->y_min + x[4] * Jy) : P->y;
    ey = std::exp(y);

#if _DEBUG_FRAGM_
    std::cout << std::scientific << "\n y = " << y;
    std::cout << "\n Jy = " << Jy;
    std::cout << "\n ey = " << ey;
    std::cout << std::endl;
#endif

    //*************//
    //INTEGRATING X//
    //*************//

    Egamma = p * cosh(y);
    x_min = std::max(1.0l/(1.0l + eps/Egamma), 
                        1.0l/(1.0l + std::pow(mass_eta*std::cosh(y)/Egamma, 2)/rho2(y, P->energy)));
    x_max = 4.0l/(4.0l + std::pow(mass_eta/Egamma, 2));
    xx = x_min + x[3]*(x_max - x_min);
    Jx = std::max(0.0l, x_max - x_min);
    p /= xx;
    Eeta = p * cosh(y);
    beta = std::sqrt(1.0l - std::pow(mass_eta/Eeta, 2));

    //******//
    //SCALES//
    //******//
    q2_fac = P->sf * p2;
    q2_mu = P->sf_mu * p2;
    q2_fragm = P->sf_fragm * p2;
    alpha_s = get_alpha_s(q2_mu);

    V = 1.0l - p / P->sqrtS / ey;
    W = p2 / P->S / V / (1.0l - V);
    S = P->S;

#if _DEBUG_FRAGM_
    std::cout << "\n q2_fac = " << q2_fac;
    std::cout << "\n q2_mu = " << q2_mu;
    std::cout << "\n q2_fragm = " << q2_fragm;
    std::cout << "\n alpha_s = " << alpha_s;
    std::cout << "\n V = " << V;
    std::cout << "\n W = " << W;
    std::cout << "\n S = " << S;
    std::cout << std::endl;
#endif
    //skipping isolation stuff(IISOL, EPSILI, EPSILA, DELTAI, DELTA, ...)


    //**************************************//
    //Z - INTEGRATION FOR NON-ISOLATED CASE//
    //**************************************//
    z_min = 1.0l - V + V * W;
    constexpr real z_max = 1.0l;
    Jz = z_max - z_min;
    z = z_min + Jz * x[2];

#if _DEBUG_FRAGM_
    std::cout << "\n z_min = " << z_min;
    std::cout << "\n z_max = " << z_max;
    std::cout << "\n z = " << z;
    std::cout << "\n Jz = " << Jz;
    std::cout << std::endl;
#endif

    //*************//
    //v-INTEGRATION//
    //*************//
    v_min = V * W / z;
    v_max = 1.0l - (1.0l - V) / z;
    Jv = v_max - v_min;
    v = v_min + Jv * x[1];

#if _DEBUG_FRAGM_
    std::cout << "\n v_min = " << v_min;
    std::cout << "\n v_max = " << v_max;
    std::cout << "\n v = " << v;
    std::cout << "\n Jv = " << Jv;
    std::cout << std::endl;
#endif

    //*************//
    //w-INTEGRATION//
    //*************//
    w_min = V * W / z / v;
    constexpr real w_max = 1.0l;
    Jw = w_max - w_min;
    w = w_min + Jw * x[0];

#if _DEBUG_FRAGM_
    std::cout << "\n w_min = " << w_min;
    std::cout << "\n w_max = " << w_max;
    std::cout << "\n w = " << w;
    std::cout << "\n Jw = " << Jw;
    std::cout << std::endl;
#endif
    //**********//
    //SHORTHANDS//
    //**********//
    v_2 = v * v;
    v_3 = v_2 * v;
    v_4 = v_3 * v;
    _vm1 = v - 1.0l;
    _vm1_2 = _vm1 * _vm1;
    _vm1_3 = _vm1_2 * _vm1;
    _vm1_4 = _vm1_3 * _vm1;
    _1mv = 1.0l - v;
    _1mv_2 = _1mv * _1mv;
    _1mv_3 = _1mv_2 * _1mv;
    _1mv_4 = _1mv_3 * _1mv;
    _1mw = 1.0l - w;

    X1 = V * W / v / w / z;
    X2 = (1.0l - V) / _1mv / z;
    BX1 = V * W / v / z;
    BX2 = X2;
    BXJAC = Jz * Jv / BX1 / BX2 / (z * z);
    XJAC = Jz * Jv / X1 / X2 / (z * z);

#if _DEBUG_FRAGM_
    std::cout << "\n X1 = " << X1;
    std::cout << "\n X2 = " << X2;
    std::cout << "\n BX1 = " << BX1;
    std::cout << "\n BX2 = " << BX2;
    std::cout << "\n BXJAC = " << BXJAC;
    std::cout << "\n XJAC = " << XJAC;
    std::cout << std::endl;
#endif

    //*****************************************************//
    //CALL OF FRAGMENTATION FUNCTIONS - RETURN D, NOT z*D !//
    //*****************************************************//

    P->etaFF_ptr(z, q2_fragm, DPU, DPD, DPS, DPC, DPG);


#if _DEBUG_FRAGM_
    /*owens data*/
    std::cout << "\n\tPHOTONS:";
    std::cout << "\n DPU = " << DPU[0];
    std::cout << "\n DPD = " << DPD[0];
    std::cout << "\n DPS = " << DPS[0];
    std::cout << "\n DPC = " << DPC[0];
    std::cout << "\n DPG = " << DPG[0];
    std::cout << "\n";
    std::cout << "\n\tetaS:";
    std::cout << "\n DPU = " << DPU[1];
    std::cout << "\n DPD = " << DPD[1];
    std::cout << "\n DPS = " << DPS[1];
    std::cout << "\n DPC = " << DPC[1];
    std::cout << "\n DPG = " << DPG[1];
    std::cout << "\n";
#endif

    P->PDF_proj_ptr(BX1, q2_fac, UP1, UPB1, DO1, DOB1, ST1, CH1, GL1);

    P->PDF_target_ptr(BX2, q2_fac, UP2, UPB2, DO2, DOB2, ST2, CH2, GL2);

#if _DEBUG_FRAGM_
    std::cout << "\n UP1 = " << UP1;
    std::cout << "\n UPB1 = " << UPB1;
    std::cout << "\n DO1 = " << DO1;
    std::cout << "\n DOB1 = " << DOB1;
    std::cout << "\n ST1 = " << ST1;
    std::cout << "\n CH1 = " << CH1;
    std::cout << "\n GL1 = " << GL1;
    std::cout << "\n";
#endif
#if _DEBUG_FRAGM_
    std::cout << "\n UP2 = " << UP2;
    std::cout << "\n UPB2 = " << UPB2;
    std::cout << "\n DO2 = " << DO2;
    std::cout << "\n DOB2 = " << DOB2;
    std::cout << "\n ST2 = " << ST2;
    std::cout << "\n CH2 = " << CH2;
    std::cout << "\n GL2 = " << GL2;
    std::cout << "\n";
#endif

    STRU(UP1, UPB1, DO1, DOB1, ST1, CH1, GL1,
        UP2, UPB2, DO2, DOB2, ST2, CH2, GL2,
        DPU, DPD, DPS, DPC, DPG,
        GRRT, GRRC);

    SHD = BX1 * BX2 * P->S;

#if _DEBUG_FRAGM_
    std::cout << "\n\tPHOTONS:";
    std::cout << "\n GRRT:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << GRRT[0][i];

    std::cout << "\n\n GRRC:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << GRRC[0][i];

    std::cout << "\n";
    std::cout << "\n\tetaS:";

    std::cout << "\n GRRT:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << GRRT[1][i];

    std::cout << "\n\n GRRC:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << GRRC[1][i];

    std::cout << "\n";
    std::cout << "\n SHD = " << SHD;
#endif

    //*************************************************************//
    //SQUARE OF THE MATRIX ELEMENTS (INTEGRATED ON THE PHASE SPACE)//
    //  (formerly 'FBOR(v, SHD, F01/2)')                           //
    //*************************************************************//
    
    help1 = pi * P->cf / SHD / P->nc;
    F01[0] = help1 * (v_2 + 1.0l) / v / _1mv_3;
    F01[2] = F01[0];
    F01[4] = help1 * (2.0l * v_2 - 2.0l * v + 1.0l) / v / _1mv;
    F01[5] = 2.0l * help1 * (v - 2.0l + (4 + 1.0l / P->nc) / v - (3.0l + 1.0l / P->nc) / v_2 + 1.0l / v_3) / _1mv_3;
    F01[10] = 2.0l * help1 / P->nc * (P->nc * v_4 - (3. * P->nc + 1.) * v_3 + (4. * P->nc + 1.) * v_2 - 2.0l * P->nc * v + P->nc) / _1mv_3 / v;
    F01[11] = help1 / P->nc * (2. * v_2 - 2. * v + 1.) * (2. * N2 * v_2 - 2. * N2 * v + N2 - 1.) / v_2 / _1mv_2;
    F01[12] = pi / SHD / v / _1mv / (2. * N2) * (v_2 + 1.) * ((N2 - 1.) * v_2 + 2. * v + (N2 - 1.)) / v / _1mv_2;
    F01[13] = pi / (2. * N2) / SHD / v / _1mv * (v_2 - 2. * v + 2.) * ((N2 - 1.) * v_2 - 2. * N2 * v + 2. * N2) / v_2 / _1mv;
    F01[14] = pi * (4. * N2) / vC / SHD / v / _1mv * (3. - v * _1mv + v / _1mv_2 + _1mv / v_2);
    F01[15] = pi / (2. * P->nc) / vC / SHD / v / _1mv * (v_2 + _1mv_2) * (2. * N2 * (v_2 - v) + N2 - 1.) / v / _1mv;

    //*************************************************************//

    //FBOR(1.0l - _1mv, SHD, F02);
    F02[0] = pi * P->cf / P->nc / SHD / _1mv * (_1mv_2 + 1.) / v_3;
    F02[2] = F02[0];
    F02[4] = F01[4];
    F02[5] = F01[5];
    F02[10] = pi * 2. * P->cf / N2 / SHD / _1mv / v * (P->nc * _1mv_4 - (3. * P->nc + 1.) * _1mv_3 + (4. * P->nc + 1.) * _1mv_2 - 2. * P->nc * _1mv + P->nc) / (v_2);
    F02[11] = F01[11];
    F02[12] = F01[13];
    F02[13] = F01[12];
    F02[14] = F01[14];
    F02[15] = F01[15];

    //*************************************************************//

#if _DEBUG_FRAGM_
    std::cout << "\nF01:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl << i << " -> " << F01[i];

    std::cout << "\n\nF02:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl << i << " -> " << F02[i];
#endif

    //*********************************//
    //formerly AvDELMU(v, SHD, FDELMU1)//
    //*********************************//

    //FMU(v, SHD, FDELMU1)//
    FDELMU1[0] = (12 * (v_2 + 1) * (4 * GTR * v4 + 4 * GTR * v3 - 11 * v2 - 11 * v1) / ((_vm1 * _vm1) * v) / 18.0l / (8.0l * N2)) / _1mv / SHD;
    FDELMU1[2] = FDELMU1[0];
    FDELMU1[4] = (12 * (2 * v_2 - 2 * v + 1) * (4 * GTR * v4 + 4 * GTR * v3 - 11 * v2 - 11 * v1) / v / 18.0l / (8.0l * N2)) / _1mv / SHD;
    FDELMU1[5] = (12 * (4 * GTR * (v_4 - 2 * v_3 + 4 * v_2 - 3 * v + 1) * v4 - 11 * _vm1 * v * v4 + 4 * GTR * (v_4 - 2 * v_3 + 4 * v_2 - 3 * v + 1) * v3 - 11 * _vm1 * v * v3 - 11 * (v_4 - 2 * v_3 + 4 * v_2 - 3 * v + 1) * v2 + 4 * GTR * _vm1 * v * v2 - 11 * (v_4 - 2 * v_3 + 4 * v_2 - 3 * v + 1) * v1) / ((_vm1 * _vm1) * v_3) / 9.0l / (8.0l * N2)) / _1mv / SHD;
    FDELMU1[10] = (12 * (4 * GTR * (v_4 - 3 * v_3 + 4 * v_2 - 2 * v + 1) * v4 + 11 * _vm1 * v_2 * v4 + 4 * GTR * (v_4 - 3 * v_3 + 4 * v_2 - 2 * v + 1) * v3 + 11 * _vm1 * v_2 * v3 - 11 * (v_4 - 3 * v_3 + 4 * v_2 - 2 * v + 1) * v2 - 4 * GTR * _vm1 * v_2 * v2 - 11 * (v_4 - 3 * v_3 + 4 * v_2 - 2 * v + 1) * v1) / ((_vm1 * _vm1) * v) / 9.0l / (8.0l * N2)) / _1mv / SHD;
    FDELMU1[11] = ((132 * N4 * _vm1_4 * ((2 * v_2 - 2 * v + 1) * (2 * v_2 - 2 * v + 1)) * vC - 48 * GTR * N3 * _vm1_4 * ((2 * v_2 - 2 * v + 1) * (2 * v_2 - 2 * v + 1)) * vC - 132 * N2 * _vm1_4 * (2 * v_2 - 2 * v + 1) * vC + 48 * GTR * P->nc * _vm1_4 * (2 * v_2 - 2 * v + 1) * vC) / (N2 * (_vm1 * _vm1 * _vm1 * _vm1 * _vm1) * v_2) / 18.0l / (8.0l * N2)) / _1mv / SHD;
    FDELMU1[12] = ((-132 * N4 * _vm1_4 * ((v_2 + 1) * (v_2 + 1)) * vC + 48 * GTR * N3 * _vm1_4 * ((v_2 + 1) * (v_2 + 1)) * vC + 132 * N2 * (_vm1 * _vm1 * _vm1 * _vm1 * _vm1 * _vm1) * (v_2 + 1) * vC - 48 * GTR * P->nc * (_vm1 * _vm1 * _vm1 * _vm1 * _vm1 * _vm1) * (v_2 + 1) * vC) / (N2 * (_vm1 * _vm1 * _vm1 * _vm1 * _vm1 * _vm1) * v_2) / 18.0l / (8.0l * P->nc * vC)) / _1mv / SHD;
    FDELMU1[13] = ((44 * N4 * _vm1_4 * ((v_2 - 2 * v + 2) * (v_2 - 2 * v + 2)) * vC - 16 * GTR * N3 * _vm1_4 * ((v_2 - 2 * v + 2) * (v_2 - 2 * v + 2)) * vC - 44 * N2 * _vm1_4 * v_2 * (v_2 - 2 * v + 2) * vC + 16 * GTR * P->nc * _vm1_4 * v_2 * (v_2 - 2 * v + 2) * vC) / (N2 * (_vm1 * _vm1 * _vm1 * _vm1 * _vm1) * v_3) / 6.0l / (8.0l * P->nc * vC)) / _1mv / SHD;
    FDELMU1[14] = ((192 * GTR * N2 * _vm1_4 * ((v_2 - v + 1) * (v_2 - v + 1) * (v_2 - v + 1)) * vC - 528 * N3 * _vm1_4 * ((v_2 - v + 1) * (v_2 - v + 1) * (v_2 - v + 1)) * vC) / ((_vm1 * _vm1 * _vm1 * _vm1 * _vm1 * _vm1) * v_3) / 9.0l / (8.0l * (vC * vC))) / _1mv / SHD;
    FDELMU1[15] = ((44 * N4 * _vm1_4 * ((2 * v_2 - 2 * v + 1) * (2 * v_2 - 2 * v + 1)) * vC - 16 * GTR * N3 * _vm1_4 * ((2 * v_2 - 2 * v + 1) * (2 * v_2 - 2 * v + 1)) * vC - 44 * N2 * _vm1_4 * (2 * v_2 - 2 * v + 1) * vC + 16 * GTR * P->nc * _vm1_4 * (2 * v_2 - 2 * v + 1) * vC) / (N2 * (_vm1 * _vm1 * _vm1 * _vm1 * _vm1) * v_2) / 6.0l / (8.0l * (vC * vC))) / _1mv / SHD;

    //*********************************//

    //FMU(1.0l - v, SHD, FDELMU2)
    FDELMU2[0] = (12.0l * (_1mv_2 + 1.0l) * (4.0l * GTR * v4 + 4.0l * GTR * v3 - 11.0l * v2 - 11.0l * v1) / (v_2 * _1mv) / 18.0l / (8.0l * N2)) / v / SHD;
    FDELMU2[2] = FDELMU2[0];
    FDELMU2[4] = FDELMU1[4];
    FDELMU2[5] = FDELMU1[5];
    FDELMU2[10] = ((12 * (4 * GTR * (_1mv_4 - 3 * _1mv_3 + 4 * _1mv_2 - 2 * _1mv + 1) * v4 + 11 * (-v) * _1mv_2 * v4 + 4 * GTR * (_1mv_4 - 3 * _1mv_3 + 4 * _1mv_2 - 2 * _1mv + 1) * v3 + 11 * (-v) * _1mv_2 * v3 - 11 * (_1mv_4 - 3 * _1mv_3 + 4 * _1mv_2 - 2 * _1mv + 1) * v2 - 4 * GTR * (-v) * _1mv_2 * v2 - 11 * (_1mv_4 - 3 * _1mv_3 + 4 * _1mv_2 - 2 * _1mv + 1) * v1) / (v_2 * _1mv) / 9.0 / (8.0l * N2)) / v / SHD);
    FDELMU2[11] = FDELMU1[11];
    FDELMU2[12] = FDELMU1[13];
    FDELMU2[13] = FDELMU1[12];
    FDELMU2[14] = FDELMU1[14];
    FDELMU2[15] = FDELMU1[15];

#if _DEBUG_FRAGM_
    std::cout << "\n\nFDELMU1:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << FDELMU1[i];

    std::cout << "\n\nFDELMU2:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << FDELMU2[i];
#endif
    BORN[0] = (F01[0] * GRRT[0] + F02[0] * GRRC[0]) * BXJAC,
    //BORN[1] = 0.0l;
    BORN[2] = (F01[2] * GRRT[2] + F02[2] * GRRC[2]) * BXJAC;
    // BORN[3] = 0.0l;
    BORN[4] = (F01[4] * GRRT[4] + F02[4] * GRRC[4]) * BXJAC;
    BORN[5] = (F01[5] * GRRT[5] + F02[5] * GRRC[5]) * BXJAC;
    // BORN[6] = 0.0l;
    // BORN[7] = 0.0l;
    // BORN[8] = 0.0l;
    // BORN[9] = 0.0l;
    BORN[10] = (F01[10] * GRRT[10] + F02[10] * GRRC[10]) * BXJAC;
    BORN[11] = (F01[11] * GRRT[11] + F02[11] * GRRC[11]) * BXJAC;
    BORN[12] = (F01[12] * GRRT[12] + F02[12] * GRRC[12]) * BXJAC;
    BORN[13] = (F01[13] * GRRT[13] + F02[13] * GRRC[13]) * BXJAC;
    BORN[14] = (F01[14] * GRRT[14] + F02[14] * GRRC[14]) * BXJAC;
    BORN[15] = (F01[15] * GRRT[15] + F02[15] * GRRC[15]) * BXJAC;

    help1 = std::log(SHD / q2_mu);

    DMU[0] = (FDELMU1[0] * GRRT[0] + FDELMU2[0] * GRRC[0]) * help1 * BXJAC;
    // DMU[1] = 0.0l;
    DMU[2] = (FDELMU1[2] * GRRT[2] + FDELMU2[2] * GRRC[2]) * help1 * BXJAC;
    // DMU[3] = 0.0l;
    DMU[4] = (FDELMU1[4] * GRRT[4] + FDELMU2[4] * GRRC[4]) * help1 * BXJAC;
    DMU[5] = (FDELMU1[5] * GRRT[5] + FDELMU2[5] * GRRC[5]) * help1 * BXJAC;
    // DMU[6] = 0.0l;
    // DMU[7] = 0.0l;
    // DMU[8] = 0.0l;
    // DMU[9] = 0.0l;
    DMU[10] = (FDELMU1[10] * GRRT[10] + FDELMU2[10] * GRRC[10]) * help1 * BXJAC;
    DMU[11] = (FDELMU1[11] * GRRT[11] + FDELMU2[11] * GRRC[11]) * help1 * BXJAC;
    DMU[12] = (FDELMU1[12] * GRRT[12] + FDELMU2[12] * GRRC[12]) * help1 * BXJAC;
    DMU[13] = (FDELMU1[13] * GRRT[13] + FDELMU2[13] * GRRC[13]) * help1 * BXJAC;
    DMU[14] = (FDELMU1[14] * GRRT[14] + FDELMU2[14] * GRRC[14]) * help1 * BXJAC;
    DMU[15] = (FDELMU1[15] * GRRT[15] + FDELMU2[15] * GRRC[15]) * help1 * BXJAC;
        

#if _DEBUG_FRAGM_
    std::cout << "\nBORN:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << BORN[0][i];

    std::cout << "\n\nDMU:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << DMU[0][i];
#endif

    P->PDF_proj_ptr(X1, q2_fac, UP1, UPB1, DO1, DOB1, ST1, CH1, GL1);

#if _DEBUG_FRAGM_
    std::cout << "\n UP1 = " << UP1;
    std::cout << "\n UPB1 = " << UPB1;
    std::cout << "\n DO1 = " << DO1;
    std::cout << "\n DOB1 = " << DOB1;
    std::cout << "\n ST1 = " << ST1;
    std::cout << "\n CH1 = " << CH1;
    std::cout << "\n GL1 = " << GL1;
    std::cout << "\n";
#endif

    STRU(UP1, UPB1, DO1, DOB1, ST1, CH1, GL1,
        UP2, UPB2, DO2, DOB2, ST2, CH2, GL2,
        DPU, DPD, DPS, DPC, DPG,
        GPPT, GPPC);

#if _DEBUG_FRAGM_
    std::cout << "\n\tPHOTONS:";
    std::cout << "\nGPPT:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << GPPT[0][i];

    std::cout << "\n\nGPPC:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << GPPC[0][i];

    std::cout << "\n\tetaS:";
    std::cout << "\nGPPT:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << GPPT[1][i];

    std::cout << "\n\nGPPC:\n";
    for (auto i{0u}; i < 16; ++i)
        std::cout << std::endl
                  << i << " -> " << GPPC[1][i];
#endif

    help1 = std::log(1. - w_min);
    help4 = std::log(_1mw);

    for (n = 0; n < 10u; ++n)
    {
        J0 = IA1[n];
        help2 = FvWPL1(1.0l, v, z);
        help3 = FvLO1(1.0l, v, z);
        FHODEL1[J0] = (FDEL1(v, z) + help1 * (help2 + help3 * help1  / 2.0l)) / _1mv - Jw / _1mv * (help2 + help3 * help4) / _1mw;
        FHODEL1[J0] = FHODEL1[J0] / CC[J0] * BXJAC;
    }

#if _DEBUG_FRAGM_
    for (J0 = 0u; J0 < 16u; ++J0)
        std::cout << "\nFHODEL1[" << J0 << "] = " << FHODEL1[J0];

#endif

    for (J0 = 0u; J0 < 16u; ++J0)
        FHODEL2[J0] = FHODEL1[J0];

    for (n = 0; n < 6u; ++n)
    {
        J0 = IA2[n];
        help2 = FvWPL2(1.0l, v, z);
        help3 = FvLO2(1.0l, v, z);
        FHODEL2[J0] = (FDEL2(v, z) + (help2  + help3 * help1  / 2.0) * help1) / _1mv - Jw / _1mv * (help2 + help3 * help4) / _1mw;
        FHODEL2[J0] = FHODEL2[J0] / CC[J0] * BXJAC;
    }

#if _DEBUG_FRAGM_
    for (J0 = 0u; J0 < 16u; ++J0)
        std::cout << "\nFHODEL2[" << J0 << "] = " << FHODEL2[J0];
#endif

    for (J0 = 0u; J0 < 16u; ++J0)
    {
        FHOREST1[J0] = (w_max - w_min) / _1mv * (FvWPL1(w, v, z) / w / _1mw + FvLO1(w, v, z) / w * help4 / _1mw + FRESC1(w, v, z) / w);
        FHOREST1[J0] = FHOREST1[J0] / CC[J0] * XJAC;
    }

#if _DEBUG_FRAGM_
    for (J0 = 0u; J0 < 16u; ++J0)
        std::cout << "\nFHOREST1[" << J0 << "] = " << FHOREST1[J0];
#endif
    for (J0 = 0u; J0 < 16u; ++J0)
    {
        FHOREST2[J0] = FHOREST1[J0];
    }

    for (n=0; n < 9u; ++n)
    {
        J0 = IA3[n];
        FHOREST2[J0] = Jw / _1mv * (FvWPL2(w, v, z) / w / _1mw + FvLO2(w, v, z) / w * help4 / _1mw + FRESC2(w, v, z) / w);

        FHOREST2[J0] = FHOREST2[J0] / CC[J0] * XJAC;
    }

    GHD = GHE = 0.0l;
    for (J0 = 0u; J0 < 16u; ++J0)
    {
#if _DEBUG_FRAGM_
        std::cout << "\nFHOREST2[" << J0 << "] = " << FHOREST2[J0];
#endif
        GHD += (FHODEL1[J0] * GRRT[J0] + FHODEL2[J0] * GRRC[J0] + DMU[J0]) * alpha_s + BORN[J0];
        GHE += (FHOREST1[J0] * GPPT[J0] + FHOREST2[J0] * GPPC[J0]) * alpha_s;
    }

#if _DEBUG_FRAGM_
    std::cout << "\n\tPHOTONS:";
    std::cout << "\nGHD = " << GHD[0];
    std::cout << "\nGHE = " << GHE[0];

    std::cout << "\n\tetaS:";
    std::cout << "\nGHD = " << GHD[1];
    std::cout << "\nGHE = " << GHE[1];
#endif




    //*******************************//
    // - - - - FINAL RESULTS - - - - //
    //*******************************//

    //Factor 2*hbarc2/S pulled out of integral
    //this gives d2sigma_dydp for etas in f[0]

    f[0] = (GHD + GHE) * alpha_s * alpha_s * Jp * Jy * p;
    
#if _DEBUG_FRAGM_
    std::cout << "\n\n FRAGM PHOTONS = " << f[0] * hbarc2 * 2.0l / P->S;´
    std::cout << "\n FRAGM etaS = " << f[2] * hbarc2 * 2.0l / P->S;
    std::cout << "\n FRAGM ETAS = " << f[4] * hbarc2 * 2.0l / P->S;
#endif

#if not ALTMODE

    //**********************************//
    // - - - - MESON DECAY (B2) - - - - //
    //**********************************//    

    //factor 2*b pulled out of integral
    f[0] *= Jx/(xx*beta);

    //IF BOTH B1 and B2 contribute, dont count twice!

#endif //ALTMODE

    if(isnanl(f[0]) or isinfl(f[0])) f[0] = 0.0l;
    f[1] = f[0] * p;

    return 0;
}
