
#include "StochSimulator.hpp"
#include "common.h"
#include "simparams.h"

void doBinomialPartioning (CellState* cell,StochSimulator* stochSimulator, Doub q)
{
    // mRNA A
    Int count = 0 ;
    for (Int m=0;m<cell->get_A_mRNA_Level();m++)
    {
        if ( stochSimulator->ran.doub () < q ) count ++ ;
    }
    cell->mf_SpecieCounts[1] = count ;
    // prot A
    count = 0 ;
    for (Int p=0;p<cell->get_A_prot_Level();p++)
    {
        if ( stochSimulator->ran.doub () < q ) count ++ ;
    }
    cell->mf_SpecieCounts[2] = count ;
    // mRNA B
    count = 0 ;
    for (Int p=0;p<cell->get_B_mRNA_Level();p++)
    {
        if ( stochSimulator->ran.doub () < q ) count ++ ;
    }
    cell->mf_SpecieCounts[4] = count ;
    // prot B
    count = 0 ;
    for (Int p=0;p<cell->get_B_prot_Level();p++)
    {
        if ( stochSimulator->ran.doub () < q ) count ++ ;
    }
    cell->mf_SpecieCounts[5] = count ;
}

int main ( int argc, char *argv[] )
{

    // parameters
    SimParams sim_params ( argc , argv ) ;

    Doub km = sim_params.get_double_param ( "km" ) ;
    Doub kp = sim_params.get_double_param ( "kp" ) ;
    Doub rm = sim_params.get_double_param ( "rm" ) ;
    Doub koff = sim_params.get_double_param ( "koff" ) ;
    Doub kon = sim_params.get_double_param ( "kon" ) ;

    Doub lnm_a = sim_params.get_double_param ( "lnm_a" ) ;
    Doub lnm_b = sim_params.get_double_param ( "lnm_b" ) ;
    Doub lnm_sigma_1 = sim_params.get_double_param ( "lnm_sigma_1" ) ;
    Doub lnm_sigma_2 = sim_params.get_double_param ( "lnm_sigma_2" ) ;
    Doub elong_rate_CV = sim_params.get_double_param ( "elong_rate_CV" ) ;

    Doub doublings_per_time_unit = sim_params.get_double_param ( "mu" ) ;

    Int num_lineages = sim_params.get_int_param ( "num_lineages" ) ;
    Doub sim_duration = sim_params.get_double_param ( "sim_duration" ) ;
    Doub update_period = sim_params.get_double_param ("update_period" ) ;
    Int num_updates_per_output = sim_params.get_double_param ("num_updates_per_output" ) ;
    Int num_timepoints = (Int) ( sim_duration / update_period / (Doub) num_updates_per_output ) ;

    string out_folder = sim_params.get_string_param ( "out_folder" ) ;
    sim_params.display_params () ;
    sim_params.write_params () ;

    // init param and sim structs
    ModelParameters* modelParameters = new ModelParameters () ;
    StochSimulator* stochSimulator = new StochSimulator ( modelParameters ) ;

    // set model parameters
    modelParameters->mf_ReacRates[0] = km ; //transcription_unrepressed
    modelParameters->mf_ReacRates[1] = rm ; //mRNA_degradation
    modelParameters->mf_ReacRates[2] = kp ; //translation
    modelParameters->mf_ReacRates[3] = km ; //transcription_unrepressed
    modelParameters->mf_ReacRates[4] = rm ; //mRNA_degradation
    modelParameters->mf_ReacRates[5] = kp ; //translation
    modelParameters->mf_ReacRates[6] = kon ; //repression
    modelParameters->mf_ReacRates[7] = koff ; //derepression
    modelParameters->mf_ReacRates[8] = kon ; //repression
    modelParameters->mf_ReacRates[9] = koff ; //derepression

    // Matrice storing data
    MatDoub lineage_traj_mRNA_data ( num_timepoints , num_lineages , 0. ) ;
    MatDoub lineage_traj_prot_data ( num_timepoints , num_lineages , 0. ) ;
    MatDoub lineage_traj_size_data ( num_timepoints , num_lineages , 0. ) ;
    Int num_timepoints_done = 0 ;

    // iterate on cells
    for (Int c=0;c<num_lineages;c++)
    {
        // construction of a cell
        CellState *cell = new CellState (modelParameters) ;
        // gene A and B copy number
        cell->mf_SpecieCounts[0] = 1 ;
        cell->mf_SpecieCounts[3] = 1 ;

        // first birth size is average one at steady-state
        Doub Lb = lnm_b / ( 2. - lnm_a ) ;
        Doub L = Lb ;

        // choose division size
        Doub Ld = Lb ; while ( Ld <= Lb ) Ld = lnm_a * Lb + lnm_b + stochSimulator->ran.norm ( 0 , lnm_sigma_1 ) ;

        // choose elongation rate ( and thus division time with size increment )
        Doub alpha = 0 ; while ( alpha <= 0 ) alpha = stochSimulator->ran.norm ( log(2.) * doublings_per_time_unit , elong_rate_CV * log(2.) * doublings_per_time_unit ) ;

        // set the rates of bi-molecular reactions, that depends on cell size !
        modelParameters->mf_ReacRates[6] = kon / L ;
        modelParameters->mf_ReacRates[8] = kon / L ;

        // timers (all in absolute time)
        Doub t_birth = 0 ;
        Doub t = 0 ;
        Doub next_update_t = 0 ;
        Doub next_div_t = log ( Ld / Lb ) / alpha ;
        num_timepoints_done = 0 ;
        Int num_updates_since_last_output = 0 ;

        // simulate the whole lineage
        while ( t < sim_duration )
        {
            // if next div before next update
            if ( next_div_t < next_update_t )
            {
                // simulate until div
                stochSimulator->simulate ( cell , next_div_t - t ) ; // remaining time before div
                // update time
                t = next_div_t ;
                // do imperfect cell splitting
                Lb = 0 ; while (Lb <= 0 || Lb >= Ld) Lb = Ld * stochSimulator->ran.norm ( 0.5 , lnm_sigma_2 ) ;
                L = Lb ;
                // division: binomial partitioning of molecular species
                doBinomialPartioning (cell,stochSimulator,Lb/Ld) ;
                // next division size and elongation rate
                Doub Ld = Lb ; while ( Ld <= Lb ) Ld = lnm_a * Lb + lnm_b + stochSimulator->ran.norm ( 0 , lnm_sigma_1 ) ;
                Doub alpha = 0 ; while ( alpha <= 0 ) alpha = stochSimulator->ran.norm ( log(2.) * doublings_per_time_unit , elong_rate_CV * log(2.) * doublings_per_time_unit ) ;                // update timers
                t_birth = t ;
                next_div_t = t + log ( Ld / Lb ) / alpha ;
            }
            else
            {
                if ( num_updates_since_last_output == num_updates_per_output - 1 )
                {

                    lineage_traj_mRNA_data[num_timepoints_done][c] = cell->get_A_mRNA_Level() ;
                    lineage_traj_prot_data[num_timepoints_done][c] = cell->get_A_prot_Level() ;
                    lineage_traj_size_data[num_timepoints_done][c] = L ;
                    num_timepoints_done++ ;
                    num_updates_since_last_output = 0 ;
                }
                else
                {
                    num_updates_since_last_output++ ;
                }
                // simulate
                stochSimulator->simulate ( cell , next_update_t - t ) ; // remaining time before update
                // update cell variables
                L *= exp ( alpha * (next_update_t-t) ) ;
                // update the rates of bi-molecular reactions !
                modelParameters->mf_ReacRates[6] = kon / L ;
                modelParameters->mf_ReacRates[8] = kon / L ;
                // update time
                t = next_update_t ;
                // update timers
                next_update_t = t + update_period ;
                next_div_t = t + log ( Ld / L ) / alpha ;
            }
        }

        // free memory used by the cell if not needed
        delete cell ;
    }

    // write data on file
    ostringstream suffix ; suffix << ".dat" ;
    writeMatrixInTextFile ( out_folder , "lineage_traj_mRNA_data"+suffix.str ()  , lineage_traj_mRNA_data ) ;
    writeMatrixInTextFile ( out_folder , "lineage_traj_prot_data"+suffix.str ()  , lineage_traj_prot_data ) ;
    writeMatrixInTextFile ( out_folder , "lineage_traj_size_data"+suffix.str ()  , lineage_traj_size_data ) ;
    return 0 ;

}

