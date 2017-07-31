#ifndef GEMODELFUNCTIONS
#define GEMODELFUNCTIONS

// hard-coded fits to klumpp 'data point' for gene expression params of constitutive / P-like proteins
// results are given in some 'relative' scale except for rm (but we will all use them as relative, with 2 dblgs/hr as ref)
Doub klumpp_gd (Doub mu_per_hr) { return 0.0006 + 1.1846 * exp ( 0.4403 * mu_per_hr ) ; }
Doub klumpp_kmg (Doub mu_per_hr) { return  max ( 0. , -1.7678 + 3.378 * exp (1.5568 * mu_per_hr) / ( 1 + exp (1.5568 * mu_per_hr) ) ) ; }
Doub klumpp_km (Doub mu_per_hr) { return klumpp_gd (mu_per_hr) * klumpp_kmg (mu_per_hr) ; }
Doub klumpp_vol (Doub mu_per_hr) { return max ( -5.2691 + 5.0621 * exp ( 0.2186 * mu_per_hr ) , 0.4) ; }
Doub klumpp_rm (Doub mu_per_hr) { return 0.5591 - 0.0443 * mu_per_hr ; }
Doub klumpp_kp (Doub mu_per_hr) { return 1.0144 - 0.0406 * mu_per_hr ; }

#include "CellState.hpp"
#include "ModelParameters.hpp"

void update_gene_expression_rates ( CellState* cell , ModelParameters* modelParameters , string GE_model )
{
    if ( GE_model == "constant_rates" ) { return ; }

    else if ( GE_model == "dependent_rates" )
    {
        Doub km = modelParameters->km_0 + modelParameters->km_per_alpha * cell->alpha + modelParameters->km_per_size * cell->L ;
        modelParameters->set_km( max(km,0.) );

        Doub kp = modelParameters->kp_0 + modelParameters->kp_per_size * cell->L;
        modelParameters->set_kp( max (kp,0.) );
    }

    else if ( GE_model == "klumpp_constitutive" )
    {
        cout << "error: GE_model unsupported now (" << GE_model << ")" << endl ; exit(1) ;
    }

    else { cout << "error: unknown GE_model (" << GE_model << ")" << endl ; exit(2) ; }
}


#endif // GEMODELFUNCTIONS
