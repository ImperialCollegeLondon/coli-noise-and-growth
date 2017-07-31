/*
__ Francois Bertaux, Imperial College
__ f.bertaux@imperial.ac.uk
__ December 2015
*/


#include "ModelParameters.hpp"

ModelParameters::ModelParameters ()
{
    // init parameter values
	mf_NumSpecies = 8 ;
	mf_NumReacs = 10 ;

	mf_ReacRates = VecDoub ( mf_NumReacs , 0. ) ;

	mf_ReacRates[0] = 0.0 ; //transcription_unrepressed
	mf_ReacRates[1] = 0.0 ; //mRNA_degradation
	mf_ReacRates[2] = 0.0 ; //translation
	mf_ReacRates[3] = 0.0 ; //transcription_unrepressed
	mf_ReacRates[4] = 0.0 ; //mRNA_degradation
	mf_ReacRates[5] = 0.0 ; //translation
	mf_ReacRates[6] = 0.0 ; //repression
	mf_ReacRates[7] = 0.0 ; //derepression
	mf_ReacRates[8] = 0.0 ; //repression
	mf_ReacRates[9] = 0.0 ; //derepression
}
