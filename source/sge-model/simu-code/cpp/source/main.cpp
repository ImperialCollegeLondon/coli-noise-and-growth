
#include "StochSimulator.hpp"
#include "libs/common.h"
#include "simparams.h"

Doub do_normal_partitioning_single ( Doub level , Doub proba , StochSimulator* stoch_simulator )
{
    Doub new_level = 0. ;
    for (Int m=0;m<level;m++)
    {
        if ( stoch_simulator->ran.doub () < proba ) new_level ++ ;
    }
    return new_level ;
}

void do_normal_partitioning ( Doub proba , CellState* cell , StochSimulator* stoch_simulator )
{
    cell->set_mRNA_Level ( do_normal_partitioning_single ( cell->get_mRNA_Level() , proba , stoch_simulator ) ) ;
    cell->set_prot_Level ( do_normal_partitioning_single ( cell->get_prot_Level() , proba , stoch_simulator ) ) ;
}

Doub do_equal_partitioning_single ( Doub before_level , Doub proba , StochSimulator* stoch_simulator )
{
    Doub after_level = (Doub) ( (Int) ( before_level * proba ) ) ;
    if ( after_level - before_level * proba != 0. )
    {
        if ( stoch_simulator->ran.doub() < proba ) after_level += 1.0 ;
    }
    return after_level ;
}

void do_equal_partitioning ( Doub proba , CellState* cell , StochSimulator* stoch_simulator )
{
    cell->set_mRNA_Level ( do_equal_partitioning_single ( cell->get_mRNA_Level() , proba , stoch_simulator ) ) ;
    cell->set_prot_Level ( do_equal_partitioning_single ( cell->get_prot_Level() , proba , stoch_simulator ) ) ;
}

void do_partitioning ( string partitioning_type , Doub proba , CellState* cell , StochSimulator* stoch_simulator )
{
    if ( partitioning_type == "normal" )
    {
        do_normal_partitioning (proba,cell,stoch_simulator) ;
    }
    else if ( partitioning_type == "equal" )
    {
        do_equal_partitioning (proba,cell,stoch_simulator) ;
    }
}

int main ( int argc, char *argv[] )
{

    // parameters
    SimParams sim_params ( argc , argv ) ;

    Doub km = sim_params.get_double_param ( "km" ) ;
    Doub rm = sim_params.get_double_param ( "rm" ) ;
    Doub kp = sim_params.get_double_param ( "kp" ) ;

    Int num_generations = sim_params.get_int_param ( "num_generations" ) ;
    Int num_lineages = sim_params.get_int_param ( "num_lineages" ) ;

    Doub doublings_per_time_unit = sim_params.get_double_param ( "mu" ) ;
    Doub cc_duration = 1.0 / doublings_per_time_unit ;

    string partitioning_type = sim_params.get_string_param ( "partitioning_type" ) ;

    Int num_time_points_per_cc = sim_params.get_int_param ( "num_timepoints_per_cc" ) ;
    Doub delta_t = cc_duration / (Doub) ( num_time_points_per_cc - 1.0 ) ;

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

    // iterate on cells
//    cout << "starting lineages loop" << endl ;
    for (Int c=0;c<num_lineages;c++)
    {
        // construction of a cell
        CellState *cell = new CellState (modelParameters) ;

        // iterate on the number of generation asked
//        cout << "starting generations loop" << endl ;
        for ( Int gen=0 ; gen<num_generations ; gen++ )
        {
            // fluctuations during cell cycle
            mRNA_data[0][c][gen] = cell->get_mRNA_Level () ;
            prot_data[0][c][gen] = cell->get_prot_Level () ;
//            cout << "starting cell sim" << endl ;
            for (Int tp = 1 ; tp<num_time_points_per_cc ; tp++ )
            {
                stochSimulator->simulate ( cell , delta_t ) ;
                mRNA_data[tp][c][gen] = cell->get_mRNA_Level () ;
                prot_data[tp][c][gen] = cell->get_prot_Level () ;
            }
//            cout << "finished cell sim" << endl ;

            // division: binomial partitioning of mRNA and prot
            do_partitioning ( partitioning_type , 0.5 , cell , stochSimulator ) ;
        }

        // free memory used by the cell if not needed
        delete cell ;
    }

    // write data on file
//    cout << "about to write..." << endl ;
    for (Int tp=0 ; tp<num_time_points_per_cc ; tp++)
    {
        ostringstream filename_mRNA ;
        filename_mRNA << "mRNA_data_" << tp << ".dat" ;
        writeMatrixInTextFile ( out_folder , filename_mRNA.str() , mRNA_data[tp] ) ;
        ostringstream filename_prot ;
        filename_prot << "prot_data_" << tp << ".dat" ;
        writeMatrixInTextFile ( out_folder , filename_prot.str() , prot_data[tp] ) ;
    }
}
