/*
__ Francois Bertaux, Imperial College
__ f.bertaux@imperial.ac.uk
__ December 2015
*/


#include "CellState.hpp"

CellState::CellState ( ModelParameters* modelParameters ) : mf_ModelParameters (modelParameters)
{
	mf_SpecieCounts = VecDoub ( mf_ModelParameters->mf_NumSpecies , 0. ) ;
    Lb = mf_ModelParameters->lnm_b / (2. - mf_ModelParameters->lnm_a) ;
    L = Lb ;
    alpha = mf_ModelParameters->alpha_max ;
    integrated_toxin_conc = 0. ;
}
