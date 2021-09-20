#include "photon.hpp"

//pointer to loaded parton distributions & some additional data storage
LHAPDF::PDF *PDFptr{nullptr};

//pointers to ROOT graphs used for interpolation
TGraph2D *gr_up{nullptr};
TGraph2D *gr_down{nullptr};
TGraph2D *gr_strange{nullptr};
TGraph2D *gr_charm{nullptr};
TGraph2D *gr_gluon{nullptr};

TGraph2D *XUTOTF{nullptr};
TGraph2D *XDTOTF{nullptr};
TGraph2D *XSTOTF{nullptr};
TGraph2D *XCTOTF{nullptr};
//TGraph2D *XBTOTF{nullptr};
TGraph2D *XGF{nullptr};
TGraph2D *XUVALF{nullptr};
TGraph2D *XDVALF{nullptr};
TGraph2D *XSVALF{nullptr};

TGraph2D *eXUTOTF{nullptr};
TGraph2D *eXDTOTF{nullptr};
TGraph2D *eXSTOTF{nullptr};
TGraph2D *eXCTOTF{nullptr};
TGraph2D *eXBTOTF{nullptr};
TGraph2D *eXGF{nullptr};

#if not ALTMODE
int main(
    int argc,
    char **argv)
{
    using namespace std;
    using namespace chrono;

    photon_cs_params P;

    vector<real> p_vec;
    vector<real> y_vec;
    int cores {4}; 

#if not _DEBUG_DIRECT_ and not _DEBUG_FRAGM_
    bool verbose{false};
    //parse command line arguments into 'P'
    for (auto i = 0; i < argc; ++i)
    {
        if (string("-E=").compare(0, 3, argv[i], 0, 3) == 0)
        {
            P.energy = stold(argv[i] + 3);
            if (not(P.energy > 0))
            {
                std::cout << "-E: Invalid energy\n";
                return -1;
            }
        }

        else if (string("-PDF=").compare(0, 5, argv[i], 0, 5) == 0)
        {
            if (strlen(argv[i] + 5) > 100)
            {
                std::cout << "-PDF: String too long\n";
                return -1;
            }
            P.pdfname = argv[i] + 5;
            //LHAPDF will stop the execution if this PDF doesnt exist
        }

        else if (string("-pp").compare(0, 3, argv[i], 0, 3) == 0)
        {
            P.target = e_target::HYDROGEN;
            P.PDF_target_ptr = PDF_p;
        }

        else if (string("-pBe").compare(0, 4, argv[i], 0, 4) == 0)
        {
            P.target = e_target::BERYLLIUM;
            P.PDF_target_ptr = PDF_Be;
        }

        //use Reya-Glueck-Vogt photon fragmentation functions
        else if (string("-RGV").compare(0, 4, argv[i], 0, 4) == 0)
        {
            P.gammaFF = e_gammaFF::RGV;
            P.gammaFF_ptr = gammaFF_RGV;
        }
        //use owens photon fragmentation functions
        else if (string("-owens").compare(0, 6, argv[i], 0, 6) == 0)
        {
            P.gammaFF = e_gammaFF::owens;
            P.gammaFF_ptr = gammaFF_owens;
        }
        //TODO: enable if multiple piFFs used, same for etaFFs
        // //use DSS14 pion fragmentation functions
        // else if (string("-DSS14").compare(0, 6, argv[i], 0, 6) == 0)
        // {
        //     P.pionFF = e_pionFF::DSS14;
        //     P.pionFF_ptr = pionFF_DSS14;
        // }

        //output cross section format
        //might require additional integration
        else if (string("-cs=").compare(0, 4, argv[i], 0, 4) == 0)
        {
            if (string("E_d3s_dp3").compare(0, 9, argv[i], 4, 4 + 9) == 0 or
                string("E").compare(0, 1, argv[i], 4, 4 + 1) == 0)
            {
                P.cs_form = e_cs_form::E_d3s_dp3;
            }
            else if (string("d2s_dydp").compare(0, 8, argv[i], 4, 4 + 8) == 0 or
                     string("d2s").compare(0, 3, argv[i], 4, 4 + 3) == 0)
            {
                P.cs_form = e_cs_form::d2s_dydp;
            }
            else
            {
                cout << "\t-cs: Unrecognized cross section format."
                        "Try E_d3s_dp3 ('E') or d2s_dydp ('d2s')"
                     << endl;
                return -1;
            }
        }
        //p_t at points
        else if (string("-p=").compare(0, 3, argv[i], 0, 3) == 0)
        {
            p_vec.push_back(stold(argv[i] + 3));
            if (not(p_vec.back() > 0))
            {
                cout << "\t-pt: Invalid transverse momentum " << p_vec.back() << "\n";
                return -1;
            }
            P.int_p = false;
        }
        //p_t integrated over bins
        else if (string("-p_bins=").compare(0, 8, argv[i], 0, 8) == 0)
        {
            p_vec.push_back(stold(argv[i] + 8));
            if (not(p_vec.back() > 0))
            {
                cout << "\t-ptbins: Invalid transverse momentum " << p_vec.back() << "\n";
                return -1;
            }
            P.int_p = true;
        }
        //divide result by pt bin width
        else if (string("-avg_p").compare(0, 6, argv[i], 0, 6) == 0)
        {
            P.avg_p = true;
        }

        //y at points
        else if (string("-y=").compare(0, 3, argv[i], 0, 3) == 0)
        {
            y_vec.push_back(stold(argv[i] + 3));
            P.int_y = false;
        }
        //y integrated over bins
        else if (string("-y_bins=").compare(0, 8, argv[i], 0, 8) == 0)
        {
            y_vec.push_back(stold(argv[i] + 8));
            P.int_y = true;
        }
        //divide result by y bin width
        else if (string("-avg_y").compare(0, 6, argv[i], 0, 6) == 0)
        {
            P.avg_y = true;
        }

        //scales
        else if (string("-scale=").compare(0, 7, argv[i], 0, 7) == 0)
        {
            P.sf = P.sf_fragm = P.sf_mu = stold(argv[i] + 7);
            if (not(P.sf > 0))
            {
                cout << "\t-scale: Invalid scale factor" << endl;
                return -1;
            }
        }

        /**
         * 
         * INTEGRATOR OPTIONS
         * 
         * */

        else if (string("-eps=").compare(0, 5, argv[i], 0, 5) == 0)
        {
            P.eps_rel = stold(argv[i] + 5);
            if (not(P.eps_rel < 1.0l and P.eps_rel > 0.0l))
            {
                cout << "\t-eps: Choose a number between 0 and 1, " << endl;
                return -1;
            }
        }

        else if (string("-N=").compare(0, 3, argv[i], 0, 3) == 0)
        {
            P.neval_max = stoi(argv[i] + 3);
            if (not(P.neval_max >= 1000))
            {
                cout << "\t-N: Invalid integer" << endl;
                return -1;
            }
        }

        else if (string("-Nmin=").compare(0, 6, argv[i], 0, 6) == 0)
        {
            P.neval_min = stoi(argv[i] + 6);
            if (not(P.neval_min >= 1000))
            {
                cout << "\t-Nmin: Invalid integer" << endl;
                return -1;
            }
        }

        else if (string("-noran").compare(0, 6, argv[i], 0, 6) == 0)
        {
            P.random_seed = false;
        }

        /**
         * 
         * OPTIONS NOT IN P
         * 
         * */
        else if (string("-cores=").compare(0, 7, argv[i], 0, 7) == 0)
        {
            cores = stoi(argv[i] + 7);
            if (not(cores > 0))
            {
                cout << "\t-cores: Choose a number greater than 0, " << endl;
                return -1;
            }
        }
        else if (string("-v").compare(0, 2, argv[i], 0, 2) == 0)
        {
            verbose = true;
        }
        else if (i != 0)
        {
            cout << "Unrecognized Command '" << argv[i] << "', aborting\n";
            return -1;
        }
    } //end argv loop

    /**
     * 
     * SANITY CHECKS
     * 
     * */

    if (p_vec.size() == 0 or
        (p_vec.size() == 1 and P.int_p))
    {
        cout << "\tUse '-pt={x,y,z,...}' or '-p_bins={x,y,z,...}' to select photon transverse momentum points or ranges\n";
        return -1;
    }
    if (y_vec.size() == 0 or
        (y_vec.size() == 1 and P.int_y))
    {
        cout << "\tUse '-y={x,y,z,...}' or '-y_bins={x,y,z,...}' to select photon rapidity points or ranges\n";
        return -1;
    }
    if (P.neval_max < P.neval_min)
    {
        cout << "\t-N needs to be larger than -Nmin!\n";
        return -1;
    }

#else //set pt, y vals for debugging
    p_vec.resize(2);
    p_vec[0] = 13.5l;
    p_vec[1] = 14.0l;
    P.int_p = true;
    y_vec.resize(2);
    y_vec[0] = -0.l;
    y_vec[1] = +1.l;
    P.int_y = true;
    P.energy = 800.0l;
    P.PDF_target_ptr = PDF_Be;
#endif

    /**
     * 
     * PDFs
     * 
     * */

    //if (!verbose)
    LHAPDF::setVerbosity(0); //LHAPDF causes some weird output with CUBA multithreading
    PDFptr = LHAPDF::mkPDF(P.pdfname, 0);
    PDFptr->setForcePositive(1); //will make PDF positive semi definite
    /**
     * 
     * integrator setup
     * 
     * */

    int p = 10000;
    cubacores(&cores, &p);

    /**
     * 
     * FFs
     * 
     * */

    //TODO: replace with classes -> ctors/dtors
    if (P.gammaFF == e_gammaFF::RGV)
    {
        initgammaFF_RGV();
    }
    if (P.pionFF == e_pionFF::DSS14)
    {
        initpionFF_DSS14();
    }
    if (P.etaFF == e_etaFF::Aidala11)
    {
        initetaFF_Aidala11();
    }

#if !_DEBUG_DIRECT_ and !_DEBUG_FRAGM_
    /**
     * 
     * OUTPUT FILE
     * 
     * */

    time_t tt;
    time(&tt);
    auto ti{localtime(&tt)};
    if (P.random_seed)
        P.seed = tt;

    string filename{"logs/photoncs_" + to_string(tt) + ".log"};
    ofstream out(filename);

    out << "# PHOTON CROSS SECTION\n# "
        << asctime(ti);

    P.PrintDetails(out);
#endif
    /**
     * 
     * MAIN LOOP
     * 
     * */

    const auto p_it_max{p_vec.size() - (P.int_p ? 1ul : 0ul)};
    for (auto p_it{0ul}; p_it < p_it_max; ++p_it)
    {
        if (P.int_p) //bin mode
        {
            P.p_min = p_vec.at(p_it);
            P.p_max = p_vec.at(p_it + 1);
        }
        else //single point mode
        {
            P.p = p_vec.at(p_it);
        }

        const auto y_it_max{y_vec.size() - (P.int_y ? 1ul : 0ul)};
        for (auto y_it{0ul}; y_it < y_it_max; ++y_it)
        {
            if (P.int_y) //bin mode
            {
                P.y_min = y_vec.at(y_it);
                P.y_max = y_vec.at(y_it + 1);
            }
            else //single point mode
            {
                P.y = y_vec.at(y_it);
            }

            /**
             * 
             * CALCULATION CALL
             * 
             * */
            
            auto start{high_resolution_clock::now()};

            P.calc();

            auto stop{high_resolution_clock::now()};
            auto duration{duration_cast<milliseconds>(stop - start)};

#if not _DEBUG_DIRECT_ and not _DEBUG_FRAGM_
            //save results to file
            out << P << endl;
            // out.flush();

            //also write results to caller if -v option used
            if (verbose)
            {
                //print kinematic range into console
                cout.precision(2);
                cout << fixed << "E = " << P.energy << ", ";

                if (P.int_p)
                    cout << "p = [" << P.p_min << ", " << P.p_max << "], ";
                else
                    cout << "p = " << P.p << ", ";

                if (P.int_y)
                    cout << "y = [" << P.y_min << ", " << P.y_max << "] => ";
                else
                    cout << "y = " << P.y << " => ";
                cout.precision(4);
                cout << scientific << "cs = "
                     << P.out_cs << " +- " << P.out_cs_err
                     << "\n\txi2 direct = " << P.out_cs_direct_xi2
                     << ", xi2 fragm = " << P.out_cs_fragm_xi2
                     << ", xi2 B1 = " << P.out_cs_pion_B1_xi2
                     << ", xi2 B2 = " << P.out_cs_pion_B2_xi2
                     << "\n\tT_cpu = " << duration.count() << "ms for "
                     << P.out_evaluations_dir << " + " << P.out_evaluations_fragm << " + " << P.out_evaluations_b2_pi << " + " << P.out_evaluations_b2_eta
                     << " integrand evaluations\n\n";
            }
#endif
        } //end y loop

    } //end p loop

    delete PDFptr;

    //TODO: replace with a class
    if (P.gammaFF == e_gammaFF::RGV)
    {
        freegammaFF_RGV();
    }
    if (P.pionFF == e_pionFF::DSS14)
    {
        freepionFF_DSS14();
    }
    if (P.etaFF == e_etaFF::Aidala11)
    {
        freeetaFF_Aidala11();
    }

    return 0;
}

void photon_cs_params::calc()
{
    /**
     * 
     * DEBUG SECTION
     * TODO: remove
     * */

    Init();
    real xt = 2.0l * p / sqrtS;
    real sqy = std::log((1.0l+std::sqrt(1.0l-xt*xt))/xt);
    if (not int_y and std::abs(y) > sqy) 
    {
        std::cout << "\ny = "  << y << " is outside the allowed kinematic region (given by |y| < " << sqy << ")\n";
    }

#if _DEBUG_FRAGM_
    /**
     * 
     * DEBUG SECTION
     * 
     * */

    //just call the function at some random point to see what it does
    real x[]{0.5l, 0.5l, 0.5l, 0.5l, 0.5l};
    real test[4]{};
    fragm(nullptr, x, nullptr, test, this);
    std::cout << std::endl << "PHOTONS  ->  "
              << test[0] * 2 * hbarc2/S;
    std::cout << std::endl << "PIONS (B1)  ->  "
              << test[2] * 2 * hbarc2/S  << std::endl;
    return;
#endif

#if _DEBUG_DIRECT_
    //just call the function at some random point to see what it does
    real x[]{0.9991l, 0.5l, 0.5l, 0.5};
    real test{-1.23456789};
    direct(nullptr, x, nullptr, &test, this);
    std::cout << std::endl
              << test*2.0l * hbarc2 * alpha_em << std::endl;
    return;
#endif

    /**
     * 
     * INTEGRATION:
     * DIRECT PART
     * ONLY EXISTS FOR PHOTONS
     * */
    {
        const int dim{4 - (int_p ? 0 : 1) - (int_y ? 0 : 1)};
        int fail{};

        //always calculate the integral and its first moment for convenience
        real I[2];
        real Err[2];
        real Xi2[2];

        Vegas(dim, 2, direct, this, 1,
              eps_rel, eps_abs, 0, seed,
              neval_min, neval_max, N_start, N_incr, N_batch,
              0, nullptr, nullptr,
              &out_evaluations_dir, &fail,
              I, Err, Xi2);

        if (fail == 0)
        {
            std::cout << "\n ! FAILURE IN DIRECT INTEGRATION ! \n";
            // still use output
            // out_cs_direct = 0.0;
            // out_cs_direct_err = 0.0;
        }
        //else
        std::cout << "direct:  I[0]=" << I[0] << " +- " << Err[0] << "\n";
        std::cout << "direct:  I[1]=" << I[1] << " +- " << Err[1] << "\n";
        if (cs_form == e_cs_form::E_d3s_dp3)
        {
            if (!int_p)
            {
                out_cs_direct = I[0] / (2.0l * pi * p);
                out_cs_direct_err = Err[0] / (2.0l * pi * p);
            }
            else
            {
                if (std::abs(I[1]) > 0.0l)
                {
                    out_cs_direct = I[0] * I[0] / I[1] / (2 * pi);
                    out_cs_direct_err = (2.0l * Err[0] * I[0] * I[0] + Err[1]) * out_cs_direct;
                }
                else out_cs_direct = out_cs_direct_err = 0.0l;
            }
        }
        else if (cs_form == e_cs_form::d2s_dydp)
        {
            //dont need first momement here
            out_cs_direct = I[0];
            out_cs_direct_err = Err[0];
        }

        //constant factors pulled out of the integrand
        out_cs_direct *= 2.0l * hbarc2 * alpha_em;
        out_cs_direct_err *= 2.0l * hbarc2 * alpha_em;
        //TODO: is this correct?
        out_cs_direct_xi2 = std::sqrt(Xi2[0] * Xi2[0] + Xi2[1] * Xi2[1]);
    }

    /**
     * 
     * INTEGRATION:
     * FRAGMENTATION PART
     * CALCULATES PHOTONS, PIONS, ETAS SIMULTANEOUSLY FOR B1
     * 
     * */
    {
        int dim{5 - (int_p ? 0 : 1) - (int_y ? 0 : 1)};
        int fail{};
        real I[6];
        real Err[6];
        real Xi2[6];

        Vegas(dim, 6, fragm, this, 1,
              eps_rel, eps_abs, 0, seed,
              neval_min, neval_max, N_start, N_incr, N_batch,
              0, nullptr, nullptr,
              &out_evaluations_fragm, &fail,
              I, Err, Xi2);

        if (fail == 0)
        {
            std::cout << "\n ! FAILURE IN FRAGMENTATION INTEGRATION ! \n";
            // still use output
            // out_cs_fragm = 0.0;
            // out_cs_fragm_err = 0.0;
            // out_cs_pion = 0.0;
            // out_cs_pion_err = 0.0;
        }
        //else
        std::cout << "fragm:   I[0]=" << I[0] << " +- " << Err[0] << "\n";
        std::cout << "fragm:   I[1]=" << I[1] << " +- " << Err[1] << "\n";
        std::cout << "pion B1: I[2]=" << I[2] << " +- " << Err[2] << "\n";
        std::cout << "pion B1: I[3]=" << I[3] << " +- " << Err[3] << "\n";
        std::cout << "eta  B1: I[4]=" << I[4] << " +- " << Err[4] << "\n";
        std::cout << "eta  B1: I[5]=" << I[5] << " +- " << Err[5] << "\n";

        if (cs_form == e_cs_form::E_d3s_dp3)
        {
            if (!int_p)
            {
                out_cs_fragm = I[0] / (2.0l * pi * p);
                out_cs_fragm_err = Err[0] / (2.0l * pi * p);
                out_cs_pion_B1 = I[2] / (2.0l * pi * p);
                out_cs_pion_B1_err = Err[2] / (2.0l * pi * p);
                out_cs_eta_B1 = I[4] / (2.0l * pi * p);
                out_cs_eta_B1_err = Err[4] / (2.0l * pi * p);
            }
            else 
            {
                if (std::abs(I[1]) > 0.0l)
                {
                    out_cs_fragm = I[0] * I[0] / I[1] / (2 * pi);
                    out_cs_fragm_err = (2.0l * Err[0] * I[0] * I[0] + Err[1]) * out_cs_fragm;
                }
                else out_cs_fragm = out_cs_fragm_err = 0.0l;

                if (std::abs(I[3]) > 0.0l)
                {
                    out_cs_pion_B1 = I[2] * I[2] / I[3] / (2 * pi);
                    out_cs_pion_B1_err = (2.0l * Err[2] * I[2] * I[2] + Err[3]) * out_cs_pion_B1;
                }
                else out_cs_pion_B1 = out_cs_pion_B1_err = 0.0l;

                if (std::abs(I[5]) > 0.0l)
                {
                    out_cs_eta_B1 = I[4] * I[4] / I[5] / (2 * pi);
                    out_cs_eta_B1_err = (2.0l * Err[4] * I[4] * I[4] + Err[5]) * out_cs_eta_B1;
                }
                else out_cs_eta_B1 = out_cs_eta_B1_err = 0.0l;
            }
        }
        else if (cs_form == e_cs_form::d2s_dydp)
        {
            //dont need first momements here
            out_cs_fragm = I[0];
            out_cs_fragm_err = Err[0];
            out_cs_pion_B1 = I[2];
            out_cs_pion_B1_err = Err[2];
            out_cs_eta_B1 = I[4];
            out_cs_eta_B1_err = Err[4];
        }

        //constant factors pulled out of the integrand
        out_cs_fragm *= 2.0l * hbarc2 / S;
        out_cs_fragm_err *= 2.0l * hbarc2 / S;

        //last factor is the pi->2 gamma branching ratio
        out_cs_pion_B1 *= 2.0l * hbarc2 / S * 0.98823l;
        out_cs_pion_B1_err *= 2.0l * hbarc2 / S * 0.98823l;

        //last factor is the eta->2 gamma branching ratio
        out_cs_eta_B1 *= 2.0l * hbarc2 / S * 0.3938l;
        out_cs_eta_B1_err *= 2.0l * hbarc2 / S * 0.3938l;

        //TODO: is this correct? probably not
        out_cs_fragm_xi2 = std::sqrt(Xi2[0] * Xi2[0] + Xi2[1] * Xi2[1]);
        out_cs_pion_B1_xi2 = std::sqrt(Xi2[2] * Xi2[2] + Xi2[3] * Xi2[3]);
        out_cs_eta_B1_xi2 = std::sqrt(Xi2[4] * Xi2[4] + Xi2[5] * Xi2[5]);
    }
    
    /**
     * 
     * INTEGRATION:
     * PION and ETA DECAY BACKGROUND B2 -> need one additional integration!
     * 
     * */
    {
        int dim{6 - (int_p ? 0 : 1) - (int_y ? 0 : 1)};
        int fail{};
        real I[2];
        real Err[2];
        real Xi2[2];

        Vegas(dim, 2, fragm_b2_pi, this, 1,
              eps_rel, eps_abs, 0, seed,
              neval_min, neval_max, N_start, N_incr, N_batch,
              0, nullptr, nullptr,
              &out_evaluations_b2_pi, &fail,
              I, Err, Xi2);

        if (fail == 0)
        {
            std::cout << "\n ! FAILURE IN PION B2 INTEGRATION ! \n";
            // still use output
            // out_cs_fragm = 0.0;
            // out_cs_fragm_err = 0.0;
            // out_cs_pion = 0.0;
            // out_cs_pion_err = 0.0;
        }
        //else
        std::cout << "pion B2: I[0]=" << I[0] << " +- " << Err[0] << "\n";
        std::cout << "pion B2: I[1]=" << I[1] << " +- " << Err[1] << "\n";

        if (cs_form == e_cs_form::E_d3s_dp3)
        {
            if (!int_p)
            {
                out_cs_pion_B2 = I[0] / (2.0l * pi * p);
                out_cs_pion_B2_err = Err[0] / (2.0l * pi * p);
            }
            else 
            {
                if (std::abs(I[1]) > 0.0l)
                {
                    out_cs_pion_B2 = I[0] * I[0] / I[1] / (2 * pi);
                    out_cs_pion_B2_err = (2.0l * Err[0] * I[0] * I[0] + Err[1]) * out_cs_pion_B2;
                }
                else out_cs_pion_B2 = out_cs_pion_B2_err = 0.0l;
            }
        }
        else if (cs_form == e_cs_form::d2s_dydp)
        {
            //dont need first momements here
            out_cs_pion_B2 = I[0];
            out_cs_pion_B2_err = Err[0];
        }

        out_cs_pion_B2 *= 2.0l * hbarc2 / S * 0.98823l * 2.0l;
        out_cs_pion_B2_err *= 2.0l * hbarc2 / S * 0.98823l * 2.0l;

        //TODO: is this correct? probably not
        out_cs_pion_B2_xi2 = std::sqrt(Xi2[0] * Xi2[0] + Xi2[1] * Xi2[1]);
    }
    {
        int dim{6 - (int_p ? 0 : 1) - (int_y ? 0 : 1)};
        int fail{};
        real I[2];
        real Err[2];
        real Xi2[2];

        Vegas(dim, 2, fragm_b2_eta, this, 1,
              eps_rel, eps_abs, 0, seed,
              neval_min, neval_max, N_start, N_incr, N_batch,
              0, nullptr, nullptr,
              &out_evaluations_b2_eta, &fail,
              I, Err, Xi2);

        if (fail == 0)
        {
            std::cout << "\n ! FAILURE IN ETA B2 INTEGRATION ! \n";
            // still use output
            // out_cs_fragm = 0.0;
            // out_cs_fragm_err = 0.0;
            // out_cs_pion = 0.0;
            // out_cs_pion_err = 0.0;
        }
        //else
        std::cout << "eta  B2: I[0]=" << I[0] << " +- " << Err[0] << "\n";
        std::cout << "eta  B2: I[1]=" << I[1] << " +- " << Err[1] << "\n";

        if (cs_form == e_cs_form::E_d3s_dp3)
        {
            if (!int_p)
            {
                out_cs_eta_B2 = I[0] / (2.0l * pi * p);
                out_cs_eta_B2_err = Err[0] / (2.0l * pi * p);
            }
            else 
            {
                if (std::abs(I[1]) > 0.0l)
                {
                    out_cs_eta_B2 = I[0] * I[0] / I[1] / (2 * pi);
                    out_cs_eta_B2_err = (2.0l * Err[0] * I[0] * I[0] + Err[1]) * out_cs_eta_B2;
                }
                else out_cs_eta_B2 = out_cs_eta_B2_err = 0.0l;
            }
        }
        else if (cs_form == e_cs_form::d2s_dydp)
        {
            //dont need first momements here
            out_cs_eta_B2 = I[0];
            out_cs_eta_B2_err = Err[0];
        }

        out_cs_eta_B2 *= 2.0l * hbarc2 / S * 0.3938l * 2.0l;
        out_cs_eta_B2_err *= 2.0l * hbarc2 / S * 0.3938l * 2.0l;

        //TODO: is this correct? probably not
        out_cs_eta_B2_xi2 = std::sqrt(Xi2[0] * Xi2[0] + Xi2[1] * Xi2[1]);
    }

    /**
     * 
     * RESULTS PROCESSING
     * 
     * */
    //return average values if so specified in the settings
    if (int_y and avg_y)
    {
        out_cs_direct /= (y_max - y_min);
        out_cs_direct_err /= (y_max - y_min);
        out_cs_fragm /= (y_max - y_min);
        out_cs_fragm_err /= (y_max - y_min);
        out_cs_pion_B1 /= (y_max - y_min);
        out_cs_pion_B1_err /= (y_max - y_min);
        out_cs_pion_B2 /= (y_max - y_min);
        out_cs_pion_B2_err /= (y_max - y_min);
        out_cs_eta_B1 /= (y_max - y_min);
        out_cs_eta_B1_err /= (y_max - y_min);
        out_cs_eta_B2 /= (y_max - y_min);
        out_cs_eta_B2_err /= (y_max - y_min);
    }

    if (int_p and avg_p)
    {
        out_cs_direct /= (p_max - p_min);
        out_cs_direct_err /= (p_max - p_min);
        out_cs_fragm /= (p_max - p_min);
        out_cs_fragm_err /= (p_max - p_min);
        out_cs_pion_B1 /= (p_max - p_min);
        out_cs_pion_B1_err /= (p_max - p_min);
        out_cs_pion_B2 /= (p_max - p_min);
        out_cs_pion_B2_err /= (p_max - p_min);
        out_cs_eta_B1 /= (p_max - p_min);
        out_cs_eta_B1_err /= (p_max - p_min);
        out_cs_eta_B2 /= (p_max - p_min);
        out_cs_eta_B2_err /= (p_max - p_min);
    }

    out_cs = out_cs_direct + out_cs_fragm + out_cs_pion_B1 + out_cs_pion_B2 + out_cs_eta_B1 + out_cs_eta_B2;
    out_cs_err = std::sqrt(  std::pow(out_cs_direct_err, 2) 
                           + std::pow(out_cs_fragm_err, 2)
                           + std::pow(out_cs_pion_B1_err, 2)
                           + std::pow(out_cs_pion_B2_err, 2)
                           + std::pow(out_cs_eta_B1_err, 2)
                           + std::pow(out_cs_eta_B2_err, 2)
                           );
}
#endif //ALTMODE
/**
 * 
 * I/O
 * 
 * */
std::ostream &operator<<(std::ostream &out, const photon_cs_params &P)
{
    out.precision(6);

    if (P.int_p)
        out << std::scientific << P.p_min << " \t" << P.p_max << " \t";
    else
        out << std::scientific << P.p << " \t";

    if (P.int_y)
        out << P.y_min << " \t" << P.y_max << " \t";
    else
        out << P.y << " \t";

    out << P.out_cs_direct << " \t"
        << P.out_cs_direct_err << " \t\t"
        << P.out_cs_fragm << " \t"
        << P.out_cs_fragm_err << " \t\t"
        << P.out_cs_pion_B1 << " \t"
        << P.out_cs_pion_B1_err << " \t\t"
        << P.out_cs_pion_B2 << " \t"
        << P.out_cs_pion_B2_err << " \t\t"
        << P.out_cs_eta_B1 << " \t"
        << P.out_cs_eta_B1_err << " \t\t"
        << P.out_cs_eta_B2 << " \t"
        << P.out_cs_eta_B2_err << " \t\t"
        << P.out_cs << " \t"
        << P.out_cs_err;

    return out;
}

void photon_cs_params::PrintDetails(std::ostream &out)
{
    out << "\n#Settings used:" << std::scientific;
    out << "\nEnergy[GeV] = " << energy;
    out << "\nPDFSet = " << pdfname;
    out << "\nTarget = " << (target == e_target::HYDROGEN ? "H2" : "Be");
    out << "\ngammaFF = ";
    switch (gammaFF)
    {
    case e_gammaFF::owens:
        out << "Owens";
        break;
    case e_gammaFF::RGV:
        out << "RGV";
        break;
    }
    out << "\npionFF = ";
    switch (pionFF)
    {
    case e_pionFF::DSS14:
        out << "DSS14";
        break;
    }
    out << "\netaFF = ";
    switch (etaFF)
    {
    case e_etaFF::Aidala11:
        out << "Aidala11";
        break;
    }
    out << "\nCS format/units = ";
    switch (cs_form)
    {
    case e_cs_form::E_d3s_dp3:
        out << "E: E*d^3(sigma)/d^3P [pb*c^3/GeV^2]";
        break;
    case e_cs_form::d2s_dydp:
        out << "d2S: d^2(sigma)/dydp_t [pb/(GeV/c)]";
        break;
    }
    out << "\nScale = " << sf;

    out << "\n\nIntegrating momentum = \t" << int_p;
    out << "\nAveraging momentum = \t" << avg_p;

    out << "\n\nIntegrating rapidity = \t" << int_y;
    out << "\nAveraging rapidity = \t" << avg_y;

    out << "\n\nIntegrator = \tCUBA Vegas";
    out << "\nEPSREL = \t\t" << eps_rel;
    out << "\nEPSABS = \t\t" << eps_abs;
    out << "\nMINEVAL = \t\t" << neval_min;
    out << "\nMAXNEVAL = \t\t" << neval_max;
    out << "\nNSTART = \t\t" << N_start;
    out << "\nNINCREASE = \t" << N_incr;
    out << "\nNBATCH = \t\t" << N_batch;
    out << "\nSEED = \t\t\t" << seed;

    //print legend
    out << "\n\n# ";
    if (int_p)
        out << "p_min \t\tp_max   \t\t";
    else
        out << "p     \t\t";

    if (int_y)
        out << "y_min \t\t\ty_max   \t\t";
    else
        out << "y     \t\t\t";

    out << "direct    \t\tdirect err    \t\t";
    out << "fragm    \t\tfragm err    \t\t";
    out << "pion B1 \t\tpion B1 err \t\t";
    out << "pion B2 \t\tpion B2 err \t\t";
    out << "eta B1  \t\teta B1 err  \t\t";
    out << "eta B2  \t\teta B2 err  \t\t";
    out << "total \t\t\ttotal err\n";
}

/**
 * 
 * Looks up PDFs at given x, Q2
 * IORD = 1, ITAR = 2, IIP = 0
 * 
 * */
void PDF_Be(
    real const &x,
    real const &Q2,
    real &UP,
    real &UPB,
    real &DO,
    real &DOB,
    real &ST,
    real &CH,
    real &GL)
{
    static std::map<int, double> xfx_store;
    if (x <= 0.0l or x >= 1.0l)
    {
        GL = DO = DOB = UP = UPB = ST = CH = 0.0l;
        return;
    }
    PDFptr->xfxQ2(x, Q2, xfx_store);

    GL = xfx_store[21]; //fixed, gluons at #21 instead of #0

    DO = xfx_store[1];
    DOB = xfx_store[-1];
    UP = xfx_store[2];
    UPB = xfx_store[-2];

    ST = xfx_store[3];
    //STB = xfx_store[-3];      //considered equal to ST
    CH = xfx_store[4];
    //CHB = xfx_store[-4];      //considered equal to CH

    //BO = xfx_store[5];        //ignored completely
    //BOB = xfx_store[-5];      //ignored completely
    //TO = xfx_store[6];        //ignored completely
    //TOB = xfx_store[-6];      //ignored completely

    auto UDN{(UP * 4.0l + DO * 5.0l) / 9.0l};
    auto DDN{(UP * 5.0l + DO * 4.0l) / 9.0l};
    auto UBN{(UPB * 4.0l + DOB * 5.0l) / 9.0l};
    auto DBN{(UPB * 5.0l + DOB * 4.0l) / 9.0l};
    UP = UDN;
    DO = DDN;
    UPB = UBN;
    DOB = DBN;
}

/**
 * 
 * Looks up PDFs at given x, Q2
 * IORD = 1, ITAR = 0, IIP = 0
 * 
 * */
void PDF_p(
    real const &x,
    real const &Q2,
    real &UP,
    real &UPB,
    real &DO,
    real &DOB,
    real &ST,
    real &CH,
    real &GL)
{
    static std::map<int, double> xfx_store;
    if (x <= 0.0l or x >= 1.0l)
    {
        GL = DO = DOB = UP = UPB = ST = CH = 0.0l;
        return;
    }
    PDFptr->xfxQ2(x, Q2, xfx_store);

    GL = xfx_store[21]; //fixed, gluons at #21 instead of #0

    DO = xfx_store[1];
    DOB = xfx_store[-1];
    UP = xfx_store[2];
    UPB = xfx_store[-2];

    ST = xfx_store[3];
    //STB = xfx_store[-3];      //considered equal to ST
    CH = xfx_store[4];
    //CHB = xfx_store[-4];      //considered equal to CH

    //BO = xfx_store[5];        //ignored completely
    //BOB = xfx_store[-5];      //ignored completely
    //TO = xfx_store[6];        //ignored completely
    //TOB = xfx_store[-6];      //ignored completely
}

real get_alpha_s(real const &Q2)
{
    return PDFptr->alphasQ2(Q2);
}

void gammaFF_owens(real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL)
{
    auto F = alpha_em / 2.0l / pi * std::log(Q2 / (0.2l * 0.2l)) / z;
    auto help1 = (2.21l - 1.28l * z + 1.29l * z * z) * std::pow(z, 0.049l) / (1.0l - 1.63l * std::log(1.0l - z));
    help1 *= F / 9.0l;
    auto help2 = 0.002l * (1.0l - z) * (1.0l - z) * std::pow(z, -1.54l);
    help2 *= F;
    UP = CH = 4.0l * help1 + help2;
    DO = ST = help1 + help2;
    GL = 0.0243l * F * std::pow(1.0l - z, 1.03l) * std::pow(z, -0.97l);
}

void initgammaFF_RGV()
{
    gr_up = new TGraph2D();
    gr_down = new TGraph2D();
    gr_strange = new TGraph2D();
    gr_charm = new TGraph2D();
    gr_gluon = new TGraph2D();
    constexpr auto Nx{36u};
    constexpr auto NQ2{11};
    constexpr double x[Nx + 2]{0.0, 1.E-3, 1.5E-3,
                               2.E-3, 3.E-3, 5.E-3, 7.E-3, 0.01, 0.015, 0.02, 0.03,
                               0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                               0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
                               0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 0.995, 0.999, 1.0};
    constexpr double Q2[NQ2 + 1]{0.0, 4.001, 8., 15., 30., 50., 100., 300., 1.E3, 3.E3, 1.E4, 1.E5};
    std::ifstream f("fraghoxx.dat");

    auto N{0u};
    real up, down, strange, charm, gluon;
    for (auto iQ2{0u}; iQ2 < NQ2 + 1; ++iQ2)
    {
        for (auto ix{0u}; ix < Nx + 2; ++ix)
        {
            if ((ix == 0u) or (ix == Nx + 1) or (iQ2 == 0u))
                continue;
            f >> up >> down >> strange >> charm >> gluon;
            // up /= (std::pow(1.0l-x[ix], 5.0l)*std::pow(x[ix], 0.3l));
            // down /= (std::pow(1.0l-x[ix], 5.0l)*std::pow(x[ix], 0.3l));
            // strange /= (std::pow(1.0l-x[ix], 5.0l)*std::pow(x[ix], 0.3l));
            // charm /= (std::pow(1.0l-x[ix], 5.0l)*std::pow(x[ix], 0.3l));
            // gluon /= (std::pow(1.0l-x[ix], 5.0l)*std::pow(x[ix], 0.3l));
            if (iQ2 == 1u)
            {
                gr_up->SetPoint(N, x[ix], Q2[0], up);
                gr_down->SetPoint(N, x[ix], Q2[0], down);
                gr_strange->SetPoint(N, x[ix], Q2[0], strange);
                gr_charm->SetPoint(N, x[ix], Q2[0], charm);
                gr_gluon->SetPoint(N, x[ix], Q2[0], gluon);
                N++;
            }

            if (ix == 1u)
            {
                gr_up->SetPoint(N, x[0], Q2[iQ2], 0);
                gr_down->SetPoint(N, x[0], Q2[iQ2], 0);
                gr_strange->SetPoint(N, x[0], Q2[iQ2], 0);
                gr_charm->SetPoint(N, x[0], Q2[iQ2], 0);
                gr_gluon->SetPoint(N, x[0], Q2[iQ2], 0);
                N++;
            }

            gr_up->SetPoint(N, x[ix], Q2[iQ2], up);
            gr_down->SetPoint(N, x[ix], Q2[iQ2], down);
            gr_strange->SetPoint(N, x[ix], Q2[iQ2], strange);
            gr_charm->SetPoint(N, x[ix], Q2[iQ2], charm);
            gr_gluon->SetPoint(N, x[ix], Q2[iQ2], gluon);
            N++;

            if (ix == Nx)
            {
                gr_up->SetPoint(N, x[Nx + 1], Q2[iQ2], up);
                gr_down->SetPoint(N, x[Nx + 1], Q2[iQ2], down);
                gr_strange->SetPoint(N, x[Nx + 1], Q2[iQ2], strange);
                gr_charm->SetPoint(N, x[Nx + 1], Q2[iQ2], charm);
                gr_gluon->SetPoint(N, x[Nx + 1], Q2[iQ2], gluon);
                N++;
            }
        }
    }
}

void freegammaFF_RGV()
{
    delete gr_up;
    delete gr_down;
    delete gr_strange;
    delete gr_charm;
    delete gr_gluon;
}

void gammaFF_RGV(real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL)
{
    if (z <= 0.0l or z >= 1.0l)
    {
        //std::cout << "WARNING: gammaFF_RGV: z out of range (" << z << ")\n";
        UP = DO = ST = CH = GL = 0.0l;
        return;
    }
    if (Q2 <= 0.0l or Q2 >= 1.E5l)
    {
        //std::cout << "WARNING: gammaFF_RGV: Q2 out of range (" << Q2 << ")\n";
        UP = DO = ST = CH = GL = 0.0l;
        return;
    }
    const real fac = alpha_em / z;
    UP = fac * gr_up->Interpolate(z, Q2);      //*(std::pow(1.0l-z, 5.0l)*std::pow(z, 0.3l));
    DO = fac * gr_down->Interpolate(z, Q2);    //*(std::pow(1.0l-z, 5.0l)*std::pow(z, 0.3l));
    ST = fac * gr_strange->Interpolate(z, Q2); //*(std::pow(1.0l-z, 5.0l)*std::pow(z, 0.3l));
    CH = fac * gr_charm->Interpolate(z, Q2);   //*(std::pow(1.0l-z, 5.0l)*std::pow(z, 0.3l));
    GL = fac * gr_gluon->Interpolate(z, Q2);   //*(std::pow(1.0l-z, 5.0l)*std::pow(z, 0.3l));
}

void initpionFF_DSS14()
{
    XUTOTF = new TGraph2D();
    XDTOTF = new TGraph2D();
    XSTOTF = new TGraph2D();
    XCTOTF = new TGraph2D();
    //XBTOTF = new TGraph2D();
    XGF = new TGraph2D();
    XUVALF = new TGraph2D();
    XDVALF = new TGraph2D();
    XSVALF = new TGraph2D();

    constexpr auto Nx{47u};
    constexpr auto NQ2{24};
    constexpr double x[Nx]{
        0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5,
        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
        0.925, 0.95, 0.975, 1.0};

    constexpr double Q2[NQ2]{
        1.0, 1.25, 1.5, 2.5,
        4.0, 6.4, 1.0e1, 1.5e1, 2.5e1, 4.0e1, 6.4e1,
        1.0e2, 1.8e2, 3.2e2, 5.8e2, 1.0e3, 1.8e3,
        3.2e3, 5.8e3, 1.0e4, 1.8e4, 3.2e4, 5.8e4, 1.0e5};

    double Parton[9][Nx][NQ2];

    std::ifstream f("PI_DSS14.GRID");

    for (auto ix{0u}; ix < Nx - 1; ++ix)
    {
        for (auto iQ2{0u}; iQ2 < NQ2; ++iQ2)
        {
            //if(not (ix == 0u) and not (ix == Nx+1) and not (iQ2 == 0u))
            f >> Parton[0][ix][iQ2] >> Parton[1][ix][iQ2] >> Parton[2][ix][iQ2] >> Parton[3][ix][iQ2] >> Parton[4][ix][iQ2] >> Parton[5][ix][iQ2] >> Parton[6][ix][iQ2] >> Parton[7][ix][iQ2] >> Parton[8][ix][iQ2];
        }
    }

    int N[9]{};
    for (auto iQ2{0u}; iQ2 < NQ2; ++iQ2)
    {
        XUTOTF->SetPoint(N[0]++, -1e12, std::log(Q2[iQ2]), 0.0);
        XDTOTF->SetPoint(N[1]++, -1e12, std::log(Q2[iQ2]), 0.0);
        XSTOTF->SetPoint(N[2]++, -1e12, std::log(Q2[iQ2]), 0.0);
        XCTOTF->SetPoint(N[3]++, -1e12, std::log(Q2[iQ2]), 0.0);
        //XBTOTF->SetPoint(N[4]++, -1e5, std::log(Q2[iQ2]), 0.0);
        XGF->SetPoint(N[5]++, -1e12, std::log(Q2[iQ2]), 0.0);
        XUVALF->SetPoint(N[6]++, -1e12, std::log(Q2[iQ2]), 0.0);
        XDVALF->SetPoint(N[7]++, -1e12, std::log(Q2[iQ2]), 0.0);
        XSVALF->SetPoint(N[8]++, -1e12, std::log(Q2[iQ2]), 0.0);
        for (auto ix{0u}; ix < Nx - 1; ++ix)
        {
            XUTOTF->SetPoint(N[0]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[0][ix][iQ2] / ((1.0l - x[ix]) * (1.0l - x[ix]) * std::sqrt(x[ix])));
            XDTOTF->SetPoint(N[1]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[1][ix][iQ2] / ((1.0l - x[ix]) * (1.0l - x[ix]) * std::sqrt(x[ix])));
            XSTOTF->SetPoint(N[2]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[2][ix][iQ2] / ((1.0l - x[ix]) * (1.0l - x[ix]) * std::sqrt(x[ix])));
            XCTOTF->SetPoint(N[3]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[3][ix][iQ2] / (std::pow(1.0l - x[ix], 5.l) * std::pow(x[ix], 0.3l)));
            //XBTOTF->SetPoint(N[4]++, std::log(x[ix]), std::log(Q2[iQ2]),
            //Parton[4][ix][iQ2]/(std::pow(1.0-x[ix], 5)*std::pow(x[ix], 0.3l)));
            XGF->SetPoint(N[5]++, std::log(x[ix]), std::log(Q2[iQ2]),
                          Parton[5][ix][iQ2] / ((1.0l - x[ix]) * (1.0l - x[ix]) * std::pow(x[ix], 0.3l)));
            //std::cout << x[ix] << ", " << Q2[iQ2] << " = " << Parton[5][ix][iQ2]/((1.0l-x[ix])*(1.0l-x[ix])*std::pow(x[ix], 0.3l)) << std::endl;
            XUVALF->SetPoint(N[6]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[6][ix][iQ2] / ((1.0l - x[ix]) * (1.0l - x[ix]) * std::sqrt(x[ix])));
            XDVALF->SetPoint(N[7]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[7][ix][iQ2] / ((1.0l - x[ix]) * (1.0l - x[ix]) * std::sqrt(x[ix])));
            XSVALF->SetPoint(N[8]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[8][ix][iQ2] / ((1.0l - x[ix]) * (1.0l - x[ix]) * std::sqrt(x[ix])));
        }
        XUTOTF->SetPoint(N[0]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        XDTOTF->SetPoint(N[1]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        XSTOTF->SetPoint(N[2]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        XCTOTF->SetPoint(N[3]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        //XBTOTF->SetPoint(N[4]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        XGF->SetPoint(N[5]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        XUVALF->SetPoint(N[6]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        XDVALF->SetPoint(N[7]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        XSVALF->SetPoint(N[8]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
    }
}

void freepionFF_DSS14()
{
    delete XUTOTF;
    delete XDTOTF;
    delete XSTOTF;
    delete XCTOTF;
    //delete XBTOTF;
    delete XGF;
    delete XUVALF;
    delete XDVALF;
    delete XSVALF;
}

void pionFF_DSS14(real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL)
{
    if (z <= 0.0l or z >= 1.0l)
    {
        //std::cout << "WARNING: pionFF_DSS14: z out of range (" << z << ")\n";
        UP = DO = ST = CH = GL = 0.0l;
        return;
    }
    if (Q2 <= 0.0l or Q2 >= 1.E5l)
    {
        //std::cout << "WARNING: pionFF_DSS14: Q2 out of range (" << Q2 << ")\n";
        UP = DO = ST = CH = GL = 0.0l;
        return;
    }
    const auto lz{std::log(z)};
    const auto lQ2{std::log(Q2)};
    const auto z12{(1.0l - z) * (1.0l - z)};
    const auto z13{std::pow(z, 0.3l)};
    const auto sqrtz{std::sqrt(z)};
    const auto fac{(z12 * sqrtz) / z};

    real UTOT = XUTOTF->Interpolate(lz, lQ2) * fac;
    real DTOT = XDTOTF->Interpolate(lz, lQ2) * fac;
    real STOT = XSTOTF->Interpolate(lz, lQ2) * fac;
    real CTOT = XCTOTF->Interpolate(lz, lQ2) * (z12 * z12 * (1.0l - z) * z13) / z;
    //real BTOT = XBTOTF->Interpolate(lz, lQ2)*(std::pow(1.0-z, 5)*z13));
    GL = std::max(XGF->Interpolate(lz, lQ2) * (z12 * z13) / z, 0.0l);
    real UVAL = XUVALF->Interpolate(lz, lQ2) * fac;
    real DVAL = XDVALF->Interpolate(lz, lQ2) * fac;
    real SVAL = XSVALF->Interpolate(lz, lQ2) * fac;

    real Up = (UTOT + UVAL) * 0.5l;
    real UBp = (UTOT - UVAL) * 0.5l;
    real Dp = (DTOT + DVAL) * 0.5l;
    real DBp = (DTOT - DVAL) * 0.5l;
    real Sp = (STOT + SVAL) * 0.5l;
    real SBp = (STOT - SVAL) * 0.5l;
    real Cp = CTOT * 0.5l;
    //real Bp  =  BTOT*0.5l;

    UP = std::max((UBp + Up) * 0.5l, 0.0l);
    // UB =  U
    DO = std::max((DBp + Dp) * 0.5l, 0.0l);
    // DB =  D
    ST = std::max((SBp + Sp) * 0.5l, 0.0l);
    // SB =  S
    CH = std::max(Cp, 0.0l);
    // B  =  Bp
}

void initetaFF_Aidala11()
{
    eXUTOTF = new TGraph2D();
    eXDTOTF = new TGraph2D();
    eXSTOTF = new TGraph2D();
    eXCTOTF = new TGraph2D();
    //eXBTOTF = new TGraph2D();
    eXGF = new TGraph2D();

    constexpr auto Nx{47u};
    constexpr auto NQ2{24};
    constexpr double x[Nx]{
            0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
             0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
             0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5, 
             0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
             0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 
             0.925, 0.95, 0.975, 1.0
             };

    constexpr double Q2[NQ2]{
                1., 1.25, 1.5, 2.5, 
                4.0, 6.4, 1.0e1, 1.5e1, 2.5e1, 4.0e1, 6.4e1,
                1.0e2, 1.8e2, 3.2e2, 5.8e2, 1.0e3, 1.8e3,
                3.2e3, 5.8e3, 1.0e4, 1.8e4, 3.2e4, 5.8e4, 1.0e5
                };

    double Parton[9][Nx][NQ2];

    std::ifstream f("ETANLO.GRID");

    for (auto ix{0u}; ix < Nx - 1; ++ix)
    {
        for (auto iQ2{0u}; iQ2 < NQ2; ++iQ2)
        {
            //if(not (ix == 0u) and not (ix == Nx+1) and not (iQ2 == 0u))
            f >> Parton[0][ix][iQ2] >> Parton[1][ix][iQ2] >> Parton[2][ix][iQ2] >> Parton[3][ix][iQ2] >> Parton[4][ix][iQ2] >> Parton[5][ix][iQ2] >> Parton[6][ix][iQ2] >> Parton[7][ix][iQ2] >> Parton[8][ix][iQ2];
        }
    }

    int N[9]{};
    for (auto iQ2{0u}; iQ2 < NQ2; ++iQ2)
    {
        eXUTOTF->SetPoint(N[0]++, -1e12, std::log(Q2[iQ2]), 0.0);
        eXDTOTF->SetPoint(N[1]++, -1e12, std::log(Q2[iQ2]), 0.0);
        eXSTOTF->SetPoint(N[2]++, -1e12, std::log(Q2[iQ2]), 0.0);
        eXCTOTF->SetPoint(N[3]++, -1e12, std::log(Q2[iQ2]), 0.0);
        //eXBTOTF->SetPoint(N[4]++, -1e5, std::log(Q2[iQ2]), 0.0);
        eXGF->SetPoint(N[5]++,    -1e12, std::log(Q2[iQ2]), 0.0);
        for (auto ix{0u}; ix < Nx - 1; ++ix)
        {
            eXUTOTF->SetPoint(N[0]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[0][ix][iQ2] / (std::pow(1.0l - x[ix], 2) * std::pow(x[ix], 0.5l)));
            eXDTOTF->SetPoint(N[1]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[1][ix][iQ2] / (std::pow(1.0l - x[ix], 2) * std::pow(x[ix], 0.5l)));
            eXSTOTF->SetPoint(N[2]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[2][ix][iQ2] / (std::pow(1.0l - x[ix], 2) * std::pow(x[ix], 0.5l)));
            eXCTOTF->SetPoint(N[3]++, std::log(x[ix]), std::log(Q2[iQ2]),
                             Parton[3][ix][iQ2] / (std::pow(1.0l - x[ix], 7) * std::pow(x[ix], 0.3l)));
            //eXBTOTF->SetPoint(N[4]++, std::log(x[ix]), std::log(Q2[iQ2]),
            //                 Parton[4][ix][iQ2] / (std::pow(1.0l - x[ix], 7) * std::pow(x[ix], 0.3l)));
            eXGF->SetPoint(N[5]++, std::log(x[ix]), std::log(Q2[iQ2]),
                          Parton[5][ix][iQ2] / (std::pow(1.0l - x[ix], 5) * std::pow(x[ix], 0.3l)));
        }
        eXUTOTF->SetPoint(N[0]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        eXDTOTF->SetPoint(N[1]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        eXSTOTF->SetPoint(N[2]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        eXCTOTF->SetPoint(N[3]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        //eXBTOTF->SetPoint(N[4]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
        eXGF->SetPoint(N[5]++, std::log(1.0), std::log(Q2[iQ2]), 0.0);
    }
}

void freeetaFF_Aidala11()
{
    delete eXUTOTF;
    delete eXDTOTF;
    delete eXSTOTF;
    delete eXCTOTF;
    //delete eXBTOTF;
    delete eXGF;
}

void etaFF_Aidala11(real const &z, real const &Q2, real &UP, real &DO, real &ST, real &CH, real &GL)
{
    if (z <= 0.0l or z >= 1.0l)
    {
        //std::cout << "WARNING: etaFF: z out of range (" << z << ")\n";
        UP = DO = ST = CH = GL = 0.0l;
        return;
    }
    if (Q2 <= 0.0l or Q2 >= 1.E5l)
    {
        //std::cout << "WARNING: etaFF: Q2 out of range (" << Q2 << ")\n";
        UP = DO = ST = CH = GL = 0.0l;
        return;
    }
    const auto lz{std::log(z)};
    const auto lQ2{std::log(Q2)};
    const auto z12{(1.0l - z) * (1.0l - z)};
    const auto z17{std::pow(1.0l - z, 7)};
    const auto z15{std::pow(1.0l - z, 5)};
    const auto z03{std::pow(z, 0.3l)};
    const auto sqrtz{std::sqrt(z)};

    UP = std::max(0.5l * eXUTOTF->Interpolate(lz, lQ2) * z12 * sqrtz / z, 0.0l);
    DO = std::max(0.5l * eXDTOTF->Interpolate(lz, lQ2) * z12 * sqrtz / z, 0.0l);
    ST = std::max(0.5l * eXSTOTF->Interpolate(lz, lQ2) * z12 * sqrtz / z, 0.0l);
    CH = std::max(0.5l * eXCTOTF->Interpolate(lz, lQ2) * z17 * z03 / z, 0.0l);
    //real BTOT = 0.5l * eXBTOTF->Interpolate(lz, lQ2)*(std::pow(1.0-z, 5)*z13));
    GL = std::max(eXGF->Interpolate(lz, lQ2) * z15 * z03 / z, 0.0l);
}
