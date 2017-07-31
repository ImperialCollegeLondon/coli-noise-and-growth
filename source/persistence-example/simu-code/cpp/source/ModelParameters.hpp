/*
__ Francois Bertaux, Imperial College
__ f.bertaux@imperial.ac.uk
__ December 2015
*/


#ifndef MODEL_PARAMETERS
#define MODEL_PARAMETERS

#include "libs/nr3.h"

struct ModelParameters
{
    Int mf_NumSpecies;
    Int mf_NumReacs;
    VecDoub mf_ReacRates;

    ModelParameters();

    // for stochastic gene expression
    inline void set_km(Doub km) { mf_ReacRates[0] = km; }
    inline void set_rm(Doub rm) { mf_ReacRates[1] = rm; }
    inline void set_kp(Doub kp) { mf_ReacRates[2] = kp; }
    inline void set_rp(Doub rp) { mf_ReacRates[3] = rp; }

    // for birth and division size (LNM)
    Doub lnm_a;
    Doub lnm_b;
    Doub lnm_sigma_1;
    Doub lnm_sigma_2;

    // for growth and toxin inhibition
    Doub regulation_conc;
    Doub regulation_hill;
    Doub mu_max;
    Doub alpha_max;
    Doub alpha_CV;

    // for dependent GE
    Doub km_0;
    Doub km_per_size;
    Doub km_per_alpha;
    Doub kp_0;
    Doub kp_per_size;

};

#endif
