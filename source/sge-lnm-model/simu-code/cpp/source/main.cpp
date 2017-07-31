
#include "StochSimulator.hpp"
#include "libs/common.h"
#include "simparams.h"


int main ( int argc, char *argv[] )
{

    // parameters
    SimParams sim_params ( argc , argv ) ;

    Doub km = sim_params.get_double_param ( "km" ) ;
    Doub rm = sim_params.get_double_param ( "rm" ) ;
    Doub kp = sim_params.get_double_param ( "kp" ) ;

    Doub lnm_a = sim_params.get_double_param ( "lnm_a" ) ;
    Doub lnm_b = sim_params.get_double_param ( "lnm_b" ) ;
    Doub lnm_sigma_1 = sim_params.get_double_param ( "lnm_sigma_1" ) ;
    Doub lnm_sigma_2 = sim_params.get_double_param ( "lnm_sigma_2" ) ;
    Doub elong_rate_CV = sim_params.get_double_param ( "elong_rate_CV" ) ;

    Int num_generations = sim_params.get_int_param ( "num_generations" ) ;
    Int num_lineages = sim_params.get_int_param ( "num_lineages" ) ;

    Doub doublings_per_time_unit = sim_params.get_double_param ( "mu" ) ;
    Int num_time_points_per_cc = sim_params.get_int_param ( "num_timepoints_per_cc" ) ;

    string out_folder = sim_params.get_string_param ( "out_folder" ) ;
    sim_params.display_params () ;
    sim_params.write_params () ;

    // init param and sim structs
    ModelParameters* modelParameters = new ModelParameters () ;
    StochSimulator* stochSimulator = new StochSimulator ( modelParameters ) ;

    // set model parameters
    modelParameters->mf_ReacRates[0] = km ;
    modelParameters->mf_ReacRates[1] = rm ;
    modelParameters->mf_ReacRates[2] = kp ;
    modelParameters->mf_ReacRates[3] = 0. ; // no degradation (dilution only)

    // Matrice storing data
    NRvector<MatDoub> mRNA_data(num_time_points_per_cc) ;
    for (int tp=0 ; tp<num_time_points_per_cc ; tp++) mRNA_data[tp] = MatDoub ( num_lineages , num_generations ) ;
    NRvector<MatDoub> prot_data(num_time_points_per_cc) ;
    for (int tp=0 ; tp<num_time_points_per_cc ; tp++) prot_data[tp] = MatDoub ( num_lineages , num_generations ) ;
    NRvector<MatDoub> size_data(num_time_points_per_cc) ;
    for (int tp=0 ; tp<num_time_points_per_cc ; tp++) size_data[tp] = MatDoub ( num_lineages , num_generations ) ;


    // iterate on cells
//    cout << "starting lineages loop" << endl ;
    for (Int c=0;c<num_lineages;c++)
    {
        // construction of a cell
        CellState *cell = new CellState (modelParameters) ;
        // first birth size is average one at steady-state
        Doub Lb = lnm_b / ( 2. - lnm_a ) ;
//        cout << "\taverage birth size = " << Lb << endl;

        // iterate on the number of generation asked
        for ( Int gen=0 ; gen<num_generations ; gen++ )
        {
            // choose division size
            Doub Ld = Lb ; while ( Ld <= Lb ) Ld = lnm_a * Lb + lnm_b + stochSimulator->ran.norm ( 0 , lnm_sigma_1 ) ;
            size_data[0][c][gen] = Lb ;

            // choose elongation rate ( and thus division time with size increment )
            Doub alpha = 0 ; while ( alpha <= 0 ) alpha = stochSimulator->ran.norm ( log(2.) * doublings_per_time_unit , elong_rate_CV * log(2.) * doublings_per_time_unit ) ;
            Doub cc_duration = log ( Ld / Lb ) / alpha ;
//            cout << "\t\talpha = " << alpha << " // size growth ratio = " << Ld / Lb << " // " << "cycle time = " << cc_duration << endl;

            // fluctuations during cell cycle
            mRNA_data[0][c][gen] = cell->get_mRNA_Level () ;
            prot_data[0][c][gen] = cell->get_prot_Level () ;
            Doub delta_t = cc_duration / (Doub) ( num_time_points_per_cc - 1.0 ) ;
            for (Int tp = 1 ; tp<num_time_points_per_cc ; tp++ )
            {
                stochSimulator->simulate ( cell , delta_t ) ;
                mRNA_data[tp][c][gen] = cell->get_mRNA_Level () ;
                prot_data[tp][c][gen] = cell->get_prot_Level () ;
                size_data[tp][c][gen] = size_data[tp-1][c][gen] * exp ( alpha * delta_t ) ;
            }

            // do imperfect cell splitting
            Lb = 0 ; while (Lb <= 0 || Lb >= Ld) Lb = Ld * stochSimulator->ran.norm ( 0.5 , lnm_sigma_2 ) ;

            // division: binomial partitioning of mRNA and prot
            Int new_mrna = 0 ;
            for (Int m=0;m<cell->get_mRNA_Level ();m++)
            {
                if ( stochSimulator->ran.doub () < Lb/Ld ) new_mrna ++ ;
            }
            cell->set_mRNA_Level (new_mrna) ;
            Int new_prot = 0 ;
            for (Int p=0;p<cell->get_prot_Level ();p++)
            {
                if ( stochSimulator->ran.doub () < Lb/Ld ) new_prot ++ ;
            }
            cell->set_prot_Level (new_prot) ;
        }

        // free memory used by the cell if not needed
        delete cell ;
    }

    // write data on file
    for (Int tp=0 ; tp<num_time_points_per_cc ; tp++)
    {
        ostringstream filename_mRNA ;
        filename_mRNA << "mRNA_data_" << tp << ".dat" ;
        writeMatrixInTextFile ( out_folder , filename_mRNA.str() , mRNA_data[tp] ) ;
        ostringstream filename_prot ;
        filename_prot << "prot_data_" << tp << ".dat" ;
        writeMatrixInTextFile ( out_folder , filename_prot.str() , prot_data[tp] ) ;
        ostringstream filename_size ;
        filename_size << "size_data_" << tp << ".dat" ;
        writeMatrixInTextFile ( out_folder , filename_size.str() , size_data[tp] ) ;
    }
}
