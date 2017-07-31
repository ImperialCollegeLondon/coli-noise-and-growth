
#include "StochSimulator.hpp"
#include "libs/common.h"
#include "simparams.h"

#include "alphamodelfunctions.hpp"
#include "gemodelfunctions.hpp"
#include "partitioningfunctions.hpp"
#include "sizemodelfunctions.hpp"

int main( int argc, char *argv[] )
{

    // struct for parsing parameters
    SimParams sim_params( argc , argv );

    // init model param and sim structs
    ModelParameters* modelParameters = new ModelParameters();
    Int seed = sim_params.get_int_param("random_seed");
    StochSimulator* stochSimulator = new StochSimulator( modelParameters , seed );

    // parsing model variants
    string partitioning_type = sim_params.get_string_param( "partitioning_type" );
    string GE_model = sim_params.get_string_param( "GE_model" );
    string alpha_model = sim_params.get_string_param( "alpha_model" );

    // parsing alpha model parameters
    modelParameters->mu_max = sim_params.get_double_param( "mu_max" ); // given in per hour
    modelParameters->alpha_max = modelParameters->mu_max * log(2.) / 60.; // we use minutes as the time-scale
    if(alpha_model == "norm_distrib" ) {
        modelParameters->alpha_CV = sim_params.get_double_param( "alpha_CV" );
    }
    else {
        modelParameters->regulation_hill = sim_params.get_double_param( "regulation_hill" );
        modelParameters->regulation_conc = sim_params.get_double_param( "regulation_conc" );
    }

    //  parsing GE model parameters
    if( GE_model == "constant_rates" )
    {
        modelParameters->set_km( sim_params.get_double_param( "km" ) );
        modelParameters->set_rm( sim_params.get_double_param( "rm" ) );
        modelParameters->set_kp( sim_params.get_double_param( "kp" ) );
        modelParameters->set_rp( sim_params.get_double_param( "rp" ) );
    }
    if( GE_model == "dependent_rates" )
    {
        modelParameters->km_0 = sim_params.get_double_param(("km_0") );
        modelParameters->km_per_size = sim_params.get_double_param(("km_per_size") );
        modelParameters->km_per_alpha = sim_params.get_double_param(("km_per_alpha") );

        modelParameters->kp_0 = sim_params.get_double_param(("kp_0") );
        modelParameters->kp_per_size = sim_params.get_double_param(("kp_per_size") );

        modelParameters->set_rm( sim_params.get_double_param( "rm" ) );
        modelParameters->set_rp( sim_params.get_double_param( "rp" ) );
    }

    // parsing size model parameters
    modelParameters->lnm_a = sim_params.get_double_param("lnm_a");
    modelParameters->lnm_b = sim_params.get_double_param("lnm_b");
    modelParameters->lnm_sigma_1 = sim_params.get_double_param("lnm_sigma_1");
    modelParameters->lnm_sigma_2 = sim_params.get_double_param("lnm_sigma_2");

    // simulation parameters(always the same, do not depend on model variant)
    Int num_lineages = sim_params.get_int_param( "num_lineages" );
    Doub sim_duration = sim_params.get_double_param( "sim_duration" );
    Doub update_period = sim_params.get_double_param("update_period" );
    Int num_updates_per_output = sim_params.get_double_param("num_updates_per_output" );
    Int num_timepoints =(Int)( sim_duration / update_period /(Doub) num_updates_per_output ) + 1;
    string out_folder = sim_params.get_string_param( "out_folder" );
    sim_params.display_params();
    sim_params.write_params();

    // for debugging
    Bool verbose = false;

    // Matrice storing data
    MatDoub lineage_traj_mRNA_data( num_timepoints , num_lineages , 0. );
    MatDoub lineage_traj_prot_data( num_timepoints , num_lineages , 0. );
    MatDoub lineage_traj_size_data( num_timepoints , num_lineages , 0. );
    MatDoub lineage_traj_alpha_data( num_timepoints , num_lineages , 0. );
    MatDoub lineage_traj_time_data( num_timepoints , num_lineages , 0. );
    Int num_timepoints_done = 0;

    // iterate on cells
    for(Int c=0;c<num_lineages;c++)
    {
        // construction of a cell(look constructor for init state)
        CellState *cell = new CellState(modelParameters);

        // update gene expression rates
        update_gene_expression_rates( cell , modelParameters , GE_model );

        // choose division size
        choose_division_size(cell,modelParameters,stochSimulator);

        // verbose
        if(verbose)
        {
            cout << endl << "Lineage start." << endl;
            cout << "L_birth = " << cell->Lb << " , L_div = " << cell->Ld << endl << endl;
        }

        // timers(all in absolute time)
        Doub t_birth = 0;
        Doub t = 0;
        Doub next_update_t = update_period;
        Doub next_div_t = t + cell->get_time_to_div();
        num_timepoints_done = 0;
        Int num_updates_since_last_output = 0;

        while( t < sim_duration )
        {
            if(verbose) cout << endl << "__ t = " << t << ", next_div_t = " << next_div_t << " __" << endl;

            // if next div before next update
            if( next_div_t < next_update_t )
            {
                // simulate until div
                stochSimulator->simulate( cell , next_div_t - t ); // remaining time before div

                // update time
                t = next_div_t;

                //verbose
                if(verbose) cout << "Division reached." << endl;
                //                if(verbose) cout << "Time average toxin conc = " << time_average_toxin_conc_previous_cycle << endl;
                if(verbose) cout << "Toxin concentration just before div = " << cell->get_prot_Level() / cell->L << endl;

                // choose birth size(splitting)
                if(verbose) cout << "about to choose birth size from div size = " << cell->Ld << endl;
                choose_birth_size(cell,modelParameters,stochSimulator);
                if(verbose) cout << "L_birth chosen = " << cell->Lb << endl;

                // binomial partitioning of mRNA and prot
                do_partitioning( partitioning_type , cell->Lb/cell->Ld , cell , stochSimulator );

                if(verbose) cout << "Toxin concentration just after div = " << cell->get_prot_Level() / cell->Lb << endl;

                // update alpha(how its done depends on the alpha_model)
                update_alpha( cell , modelParameters , alpha_model , stochSimulator );

                if(verbose) cout << "Next alpha = " << cell->alpha * 60. << endl;

                // update gene expression rates
                update_gene_expression_rates( cell , modelParameters , GE_model );

                // choose division size
                choose_division_size( cell , modelParameters , stochSimulator );

                // update timers
                t_birth = t;
                next_div_t = t + cell->get_time_to_div();
            }

            else // just simulate until next update
            {
                // if output needed
                if((num_updates_since_last_output == num_updates_per_output) ||(t == 0)  )
                {
                    if(num_timepoints_done == num_timepoints) { cout << "error output table" << endl; exit(5); }
                    lineage_traj_mRNA_data[num_timepoints_done][c] = cell->get_mRNA_Level();
                    lineage_traj_prot_data[num_timepoints_done][c] = cell->get_prot_Level();
                    lineage_traj_size_data[num_timepoints_done][c] = cell->L;
                    lineage_traj_alpha_data[num_timepoints_done][c] = cell->alpha;
                    lineage_traj_time_data[num_timepoints_done][c] = t;
                    num_timepoints_done++;
                    num_updates_since_last_output = 0;
                    if(verbose) cout << "outputs done = " << num_timepoints_done << endl;
                }

                // simulate
                stochSimulator->simulate( cell , next_update_t - t ); // remaining time before update

                // update alpha(will be done only if required by alpha_model)
                update_alpha( cell , modelParameters , alpha_model , stochSimulator );

                // update gene expression rates(if needed by GE_model)
                update_gene_expression_rates(cell , modelParameters , GE_model );

                // update time
                t = next_update_t;

                // update timers
                next_update_t = t + update_period;
                next_div_t = t + cell->get_time_to_div();
                num_updates_since_last_output++;
            }
        }

        // free memory used by the cell if not needed
        delete cell;
    }

    // verify all points in the tables have been filled, if needed, remove last point
    if(num_timepoints_done < num_timepoints)
    {
        //        cout << "timepoints not filled" << endl;
        if( num_timepoints_done < num_timepoints - 1 )
        {
            cout << "error: more than one timepoint missing ?!" << endl; exit(4);
        }
        lineage_traj_time_data[num_timepoints-1][0] = -1; // flag that will be detected by matlab
    }

    // write data on file
    ostringstream suffix; suffix << ".dat";
    if(verbose) cout << endl << endl << "writing tables" << endl;
    writeMatrixInTextFile( out_folder , "lineage_traj_mRNA_data"+suffix.str()  , lineage_traj_mRNA_data );
    writeMatrixInTextFile( out_folder , "lineage_traj_prot_data"+suffix.str()  , lineage_traj_prot_data );
    writeMatrixInTextFile( out_folder , "lineage_traj_size_data"+suffix.str()  , lineage_traj_size_data );
    writeMatrixInTextFile( out_folder , "lineage_traj_alpha_data"+suffix.str()  , lineage_traj_alpha_data );
    writeMatrixInTextFile( out_folder , "lineage_traj_time_data"+suffix.str() , lineage_traj_time_data );
    return 0;
}
