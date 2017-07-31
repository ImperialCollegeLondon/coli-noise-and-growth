/*
__ Francois Bertaux, Imperial College
__ f.bertaux@imperial.ac.uk
__ December 2015
*/

#ifndef CELL_STATE
#define CELL_STATE

#include "libs/nr3.h"

#include "ModelParameters.hpp"

struct CellState
{
	ModelParameters* mf_ModelParameters ;
	VecDoub mf_SpecieCounts ;

	CellState ( ModelParameters* modelParameters ) ;

    // fields
    Doub alpha ;
    Doub L ;
    Doub Lb ;
    Doub Ld ;
    Doub previous_Lb ;
    Doub previous_Ld ;

    // for alpha model cell-cycle average
    Doub integrated_toxin_conc ;

    // compute next div time
    inline Doub get_time_to_div () { return log (Ld/L) / alpha ; }

	// methods for name access to species
	inline Doub get_mRNA_Level () { return mf_SpecieCounts[0] ; }
	inline Doub get_prot_Level () { return mf_SpecieCounts[1] ; }
    inline void set_mRNA_Level ( Doub value ) { mf_SpecieCounts[0] = value ; }
    inline void set_prot_Level ( Doub value ) { mf_SpecieCounts[1] = value ; }
} ;


#endif
