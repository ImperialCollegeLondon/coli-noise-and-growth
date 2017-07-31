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
	Int mf_NumSpecies ;
	Int mf_NumReacs ;
	VecDoub mf_ReacRates ;

	ModelParameters () ;
};

#endif
