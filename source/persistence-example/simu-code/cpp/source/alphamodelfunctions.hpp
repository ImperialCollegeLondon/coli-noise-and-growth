#ifndef ALPHAMODELFUNCTIONS
#define ALPHAMODELFUNCTIONS

Doub give_repression_factor ( Doub regulation_conc , Doub regulation_hill , Doub conc )
{
    return 1. / ( 1 + pow ( conc / regulation_conc , regulation_hill ) ) ;
}

#include "CellState.hpp"
#include "StochSimulator.hpp"

void update_alpha ( CellState* cell , ModelParameters* modelParameters , string alpha_model , StochSimulator* stochSimulator )
{
    if ( alpha_model == "instantaneous" )
    {
        cell->alpha = modelParameters->alpha_max * give_repression_factor ( modelParameters->regulation_conc , modelParameters->regulation_hill , cell->get_prot_Level() / cell->L ) ;
    }
    else if ( alpha_model == "cell_cycle_average" )
    {
        if ( cell->L == cell->Lb ) // update only if just born
        {
            Doub duration_cell_cycle = log(cell->Ld/cell->previous_Lb) / cell->alpha ;
            Doub temp_average_toxin_conc = cell->integrated_toxin_conc / duration_cell_cycle ;
            cell->alpha = modelParameters->alpha_max * give_repression_factor ( modelParameters->regulation_conc , modelParameters->regulation_hill , temp_average_toxin_conc ) ;
            cell->integrated_toxin_conc = 0. ;
        }
    }
    else if ( alpha_model == "norm_distrib" )
    {
        if ( cell->L == cell->Lb ) // update only if just born
        {
            Doub alpha = 0 ; while ( alpha <= 0 ) alpha = stochSimulator->ran.norm ( modelParameters->alpha_max , modelParameters->alpha_max * modelParameters->alpha_CV ) ;
            cell->alpha = alpha;
        }
    }
    else
    {
        cout << "error: unknown alpha-model (" << alpha_model << ")." << endl ;
        exit (10) ;
    }
}

#endif // ALPHAMODELFUNCTIONS

