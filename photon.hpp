// Fabian Kraus, 2021
// Proton+Proton->Photon NLO cross section
// Direct photons, fragmentation photons + pion & eta decay background for E706



/*  compile with:

c++ photon.cpp direct.cpp fragm.cpp fragm_b2_pi.cpp fragm_b2_eta.cpp -o photoncs -lLHAPDF -lcubal `root-config --cflags --glibs` -O3

or (experimental)

g++-10 photon.cpp direct.cpp fragm.cpp fragm_b2_pi.cpp fragm_b2_eta.cpp -o photoncs -lLHAPDF -lcubal -lm `root-config --cflags --glibs` -O3 -std=c++17 -Wall -ffast-math

*/



/*  example calls (Fermilab E706)

./photoncs -v -y=0 -E=800 -N=100000 -p=8 -cs=E
./photoncs -v -E=800 -N=100000 -p_bins={3.5,4} -y_bins={-1,0.5} -avg_y

./photoncs -v -y_bins={-0.75,0.75} -avg_y -E=530 -N=200000 -p_bins={3.5,4,4.5,5,5.5,6,7,8,10,12}

./photoncs -v -y_bins={-1.00,0.50} -avg_y -E=800 -N=200000 -p_bins={3.5,4,4.5,5,5.5,6,7,8,10,12}

./photoncs -v -y_bins={-0.75,0.75} -avg_y -E=530 -N=200000 -p={2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14}
./photoncs -v -y_bins={-1.00,0.50} -avg_y -E=800 -N=200000 -p={2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14}

./photoncs -v -y={-1.,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.} -E=530 -N=200000 -p={4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14} -cs=E
./photoncs -v -y={-1.,-.9,-.8,-.7,-.6,-.5,-.4,-.3,-.2,-.1,0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.} -E=800 -N=200000 -p={4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14} -cs=E
./photoncs -v -y={-1.,-.5,0.0,.5,1.} -E=800 -N=200000 -p={4,6,8,10,12,14} -cs=E

./photoncs -v -pBe -y_bins={-0.75,0.75} -E=530 -N=200000 -p_bins={3.5,3.75,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.5,7,7.5,8,9.0,10,12} 
./photoncs -v -pBe -y_bins={-1.00,0.50} -E=800 -N=200000 -p_bins={3.5,3.75,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.5,7,7.5,8,9.0,10,12} 

*/



#pragma once
#ifndef _PHOTON_HPP_
#define _PHOTON_HPP_


/**
 * 
 * DEPENDENCIES
 * 
 * */

//for PDFs
//https://lhapdf.hepforge.org/downloads/, 6.3.0
#include <LHAPDF/LHAPDF.h>

//for numerical integration (VEGAS)
//http://www.feynarts.de/cuba/, 4.2.1
// #include <cuba.h>
//http://www.feynarts.de/cuba/, 4.2.1 with long double (--with-real=10)
#include <cubal.h>

//for interpolation algorithm
//https://root.cern/, 6.24/02
#include <TGraph2D.h>



/**
 *
 * COMPILATION FLAGS
 *
 * */

#define ALTMODE false   //for using different main fcts
#define _DEBUG_DIRECT_ false
#define _DEBUG_FRAGM_  false



/**
 *
 * INCLUDES
 *
 * */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
//using namespace std;

/**
 *
 * TYPE ALIASES
 *
 * */

typedef unsigned int natural;
typedef int integer;
typedef long double real; //remember to change cuba<->cubal!

/**
 *
 * CONSTANTS & MACROS
 *
 * */

// math
constexpr real pi   {3.14159265358979323846264338l};
constexpr real pi2  {9.86960440108935861883449099l};

// physics
constexpr real alpha_em {0.007297352569l};
constexpr real hbarc2   {3.8937937217185937212565161306E8l}; //(hbar*c)^2/(GeV^2*pb) (picobarn)

// masses in GeV/c^2
constexpr real mass_pion    {0.13497685l};
constexpr real mass_eta     {0.547862l};
constexpr real mass_proton  {0.9382720813l};


// E706 detector resolutions & measurements
constexpr real delta_r  {0.55l};             //radial resolution in cm
constexpr real d        {900.0l};            //distance target to detector in cm
constexpr real delta_phi_inner{pi / 192.0l}; //azimuthal resolution for r < r_phi
constexpr real delta_phi_outer{pi / 384.0l}; //azimuthal resolution for r > r_phi
constexpr real r_phi    {40.0l};             //in cm
constexpr real eps      {0.1l};              //energy detection threshold in GeV

// rapidity of CMS in lab frame for beam momentum p
constexpr real y_cms_lab(const real &p)
{
    const real E{std::sqrt(p * p + mass_proton * mass_proton)};
    return -0.25l * std::log((E + p) / (E - p));
}

// CMS rapidity from lab radius r
constexpr real y_cms_from_r(const real &r, const real p)
{
    return std::log((std::sqrt(d*d+r*r) + d)/r) + y_cms_lab(p);
}

// lab radius r from CMS rapidity
constexpr real r_from_y_cms(const real &y, const real p)
{
    const real yprime{y-y_cms_lab(p)};
    return 2.0l*d*std::exp(yprime)/(std::exp(2.0l*yprime)-1.0l);
}

// CMS rapidity that corresponds to the difference between inner and outer phi-board for beam momentum p
constexpr real y_phi(real const &p)
{
    return y_cms_from_r(r_phi, p);
}

// azimuthal resolution at given rapidity
constexpr real delta_phi(const real &y, real const &p)
{
    return ((y < y_phi(p)) ? delta_phi_outer : delta_phi_inner);
}

// rapidity resolution at given rapidity
constexpr real delta_y(const real &y, real const &p)
{
    const real r{r_from_y_cms(y,p)};
    return d*delta_r/(r*sqrt(d*d+r*r));
}

// total resolution
constexpr real rho2(const real &y, real const &p)
{
    const real dy{delta_y(y,p)};
    const real dphi{delta_phi(y,p)};
    return dy*dy + dphi*dphi;
}

/**
 *
 * FUNCTION DECLARATIONS
 *
 * */

//integrands
int direct(const int *ndim, const real x[], const int *ncomp, real f[], void *userdata);
int fragm(const int *ndim, const real x[], const int *ncomp, real f[], void *userdata);
int fragm_b2_pi(const int *ndim, const real x[], const int *ncomp, real f[], void *userdata);
int fragm_b2_eta(const int *ndim, const real x[], const int *ncomp, real f[], void *userdata);

//use LHAPDF to look up PDFs, returns x*f(x)!
void PDF_p(real const &x, real const &Q2, real &UP, real &UPB, real &DO, real &DOB, real &ST, real &CH, real &GL);
void PDF_Be(real const &x, real const &Q2, real &UP, real &UPB, real &DO, real &DOB, real &ST, real &CH, real &GL);

//uses LHAPDF to look up alpha_strong at given scale Q^2
real get_alpha_s(real const &Q2);

//look up selected fragmentation functions, returns D(x) NOT x*D(x)!
void gammaFF_owens(real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL);

void gammaFF_RGV(real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL);
void initgammaFF_RGV();
void freegammaFF_RGV();

void pionFF_DSS14(real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL);
void initpionFF_DSS14();
void freepionFF_DSS14();

void etaFF_Aidala11(real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL);
void initetaFF_Aidala11();
void freeetaFF_Aidala11();

/**
 *
 * SETTINGS
 *
 * */

//determines units, format of cross sections calculated
enum class e_cs_form
{
    //d3s_dy_dp2,
    E_d3s_dp3,
    d2s_dydp,
    //p3_d2s_dy_dp
};

//determines the target, i.e. what to do with the PDFs
enum class e_target
{
    HYDROGEN,
    BERYLLIUM
};

//determines the fragmentation functions to use
enum class e_gammaFF
{
    RGV,
    owens,
};

enum class e_pionFF
{
    DSS14,
};

enum class e_etaFF
{
    Aidala11,
};

//collection of all parameters with some default values
struct photon_cs_params
{
    //Main functionality
    void calc();

    /**
     * Main parameters
     * Most can be changed through the command line
     * No recompile necessary!
     * */
    real energy{800.0l};                           //beam energy in GeV
    std::string pdfname{"NNPDF31_nlo_as_0118"}; //name of parton distribution to use
    e_target target{e_target::HYDROGEN};
    e_gammaFF gammaFF{e_gammaFF::RGV};
    e_pionFF  pionFF{e_pionFF::DSS14};
    e_etaFF   etaFF{e_etaFF::Aidala11};
    e_cs_form cs_form{e_cs_form::E_d3s_dp3}; //what form the resulting diff. cross sections should be in
    
    //transverse momentum of the photon
    real p{};               //photon transverese momentum, not used if in bin mode
    real p_min{};           //lower momentum bin end
    real p_max{};           //upper momentum bin end
    bool int_p{false};      //integrate over momentum bin
    bool avg_p{false};      //average over momentum bin by adding factor 1/(delta p_t), not used if int_p false

    //rapidity of the photon (->angle)
    real y{};               //photon (pseudo-) rapidity, not used if in bin mode
    real y_min{};           //lower rapidity bin end
    real y_max{};           //upper rapidity bin end
    bool int_y{false};      //integrate over rapidity bin
    bool avg_y{false};      //average over rapiity range by adding factor 1/(delta y), not used if int_y false

    //scales
    real sf{1.0l};
    real sf_mu{1.0l};
    real sf_fragm{1.0l};

    //Numerical integration settings
    real eps_rel{1.e-6l};
    real eps_abs{0.0l}; //unused

    natural neval_min{100000u};
    natural neval_max{200000u};

    natural N_start{1000u}; //how many function calls in the first iteration
    natural N_incr{1000u};  //how much to increase function calls per iteration
    natural N_batch{1000u}; //how many function calls in a batch

    bool random_seed{true}; //will always use seed 0 if turned off
    natural seed{0};

    //results storage
    real out_cs_direct{0.0};     //direct cross section contribution
    real out_cs_direct_err{0.0}; //direct cross section contribution error
    real out_cs_direct_xi2{0.0}; //direct cross section contribution error probability

    real out_cs_fragm{0.0};      //photon fragmentation cross section contribution
    real out_cs_fragm_err{0.0};
    real out_cs_fragm_xi2{0.0};

    real out_cs_pion_B1{0.0};         //pion decay background B1 (both photons hit the same cell)
    real out_cs_pion_B1_err{0.0};
    real out_cs_pion_B1_xi2{0.0};

    real out_cs_pion_B2{0.0};         //pion decay background B2 (one photon is below energy detection threshold)
    real out_cs_pion_B2_err{0.0};
    real out_cs_pion_B2_xi2{0.0};

    real out_cs_eta_B1{0.0};         //eta decay background B1 (both photons hit the same cell)
    real out_cs_eta_B1_err{0.0};
    real out_cs_eta_B1_xi2{0.0};

    real out_cs_eta_B2{0.0};         //eta decay background B2 (one photon is below energy detection threshold)
    real out_cs_eta_B2_err{0.0};
    real out_cs_eta_B2_xi2{0.0};

    real out_cs{0.0};            //signal - background
    real out_cs_err{0.0};

    integer out_evaluations_dir{0};
    integer out_evaluations_fragm{0};
    integer out_evaluations_b2_pi{0};
    integer out_evaluations_b2_eta{0};

    //etc
    constexpr static real nc{3.0l};//number of QCD colors
    constexpr static real nf{5.0l};//number of active flavors
    constexpr static real cf{4.0l / 3.0l};
    constexpr static real nc2{nc * nc}, cf__nc{cf / nc}, cf2__nc{cf * cf / nc};

    //Set by Init(), better not touch them
    real S;   //set from 'energy'
    real sqrtS; //set from 'energy'

    //called by "calc" before integration, sets parameters that dont depend on p_t, y
    void Init()
    {
        S = 2.0l * mass_proton * energy;
        sqrtS = std::sqrt(S);
    }

    //Prints settings to stream
    void PrintDetails(std::ostream &out);

    //Internal
    //pointer to PDF function based on target class
    void (*PDF_target_ptr)(real const &x, real const &Q2, real &UP, real &UPB, real &DO, real &DOB, real &ST, real &CH, real &GL){PDF_p};
    //pointer to PDF function based on projectile class
    void (*PDF_proj_ptr)(real const &x, real const &Q2, real &UP, real &UPB, real &DO, real &DOB, real &ST, real &CH, real &GL){PDF_p};
    //pointer to FF functions
    void (*gammaFF_ptr)(real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL){gammaFF_RGV};
    void (*pionFF_ptr) (real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL){pionFF_DSS14};
    void (*etaFF_ptr)  (real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL){etaFF_Aidala11};
};

struct photon_cs_params;
std::ostream &operator<<(std::ostream &out, const photon_cs_params &P);

#endif // _PHOTON_HPP_
