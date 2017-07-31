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

	// methods for name access to species
	inline Doub get_A_gene_Level () { return mf_SpecieCounts[0] ; }
	inline Doub get_A_mRNA_Level () { return mf_SpecieCounts[1] ; }
	inline Doub get_A_prot_Level () { return mf_SpecieCounts[2] ; }
	inline Doub get_B_gene_Level () { return mf_SpecieCounts[3] ; }
	inline Doub get_B_mRNA_Level () { return mf_SpecieCounts[4] ; }
	inline Doub get_B_prot_Level () { return mf_SpecieCounts[5] ; }
	inline Doub get_A_gene_B_bound_Level () { return mf_SpecieCounts[6] ; }
	inline Doub get_B_gene_A_bound_Level () { return mf_SpecieCounts[7] ; }
};


#endif