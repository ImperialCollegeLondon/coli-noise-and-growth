#ifndef SIZEMODELFUNCTIONS
#define SIZEMODELFUNCTIONS

#include "CellState.hpp"
#include "ModelParameters.hpp"
#include "StochSimulator.hpp"

void choose_division_size ( CellState* cell , ModelParameters* modelParameters , StochSimulator* stochSimulator )
{
    cell->previous_Ld = cell->Ld ;
   cell->Ld = cell->Lb ;
    while ( cell->Ld <= cell->Lb )
    {
        cell->Ld = modelParameters->lnm_a * cell->Lb + modelParameters->lnm_b + stochSimulator->ran.norm ( 0 , modelParameters->lnm_sigma_1 ) ;
    }
}

void choose_birth_size ( CellState* cell , ModelParameters* modelParameters , StochSimulator* stochSimulator )
{
    cell->previous_Lb = cell->Lb ;
    cell->Lb = 0 ; while (cell->Lb <= 0 || cell->Lb >= cell->Ld) cell->Lb = cell->Ld * stochSimulator->ran.norm ( 0.5 , modelParameters->lnm_sigma_2 ) ;
    cell->L = cell->Lb ;
}


#endif // SIZEMODELFUNCTIONS

