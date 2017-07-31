/*
__ Francois Bertaux, Imperial College
__ f.bertaux@imperial.ac.uk
__ December 2015
*/

#ifndef STOCH_SIMULATOR
#define STOCH_SIMULATOR


#include "CellState.hpp"

#include "libs/ran.h"
#include "libs/sparse.h"

struct StochSimulator
{

        // main fields and methods
        ModelParameters* mf_ModelParameters ;
        StochSimulator ( ModelParameters* modelParameters , Int seed = 1 ) ;
        void simulate ( CellState* state , Doub duration ) ;
        void prepareForSteps ( CellState* cellState ) ;
        Doub doStep ( CellState* cellState , Doub targetTime ) ;

        // other fields and methods
        Ran ran; // random generator
        Int mm, nn;
        VecDoub a;  // reaction rates
        MatDoub instate, outstate; // reactants matrix, change state matrix
        NRvector<NRsparseCol> outchg, depend; // change state sparse matrix, reaction dependancy sparse matrix
        VecInt pr;  // priority list for reaction
        Doub t; // time
        Doub asum; // sum of all reactions rates
        typedef Doub(StochSimulator::*rateptr)(VecDoub &s);
        rateptr *dispatch;
        void describereactions();
        ~StochSimulator() {delete [] dispatch;}


	// rate functions
	Doub rate0(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[0] * s[0] ; } // transcription_unrepressed
	Doub rate1(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[1] * s[1] ; } // mRNA_degradation
	Doub rate2(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[2] * s[1] ; } // translation
	Doub rate3(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[3] * s[3] ; } // transcription_unrepressed
	Doub rate4(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[4] * s[4] ; } // mRNA_degradation
	Doub rate5(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[5] * s[4] ; } // translation
	Doub rate6(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[6] * s[0] * s[5] ; } // repression
	Doub rate7(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[7] * s[6] ; } // derepression
	Doub rate8(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[8] * s[3] * s[2] ; } // repression
	Doub rate9(VecDoub &s) {return mf_ModelParameters->mf_ReacRates[9] * s[7] ; } // derepression
};


#endif
