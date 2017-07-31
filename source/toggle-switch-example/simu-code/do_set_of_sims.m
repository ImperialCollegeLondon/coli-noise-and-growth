function do_set_of_sims ( par_file )
%DO_SET_OF_SIMS Simulate all param set in a csv file
%   par_file in input-data

% read the parameter data
par_data = readtable ( par_file ) ;
lnm_trends = readtable ( '../input-data/linear_trends_lnm_params.csv' , 'ReadRowNames' , true ) ;

% sim params
sim_pars.sim_duration = 30000000 ;
sim_pars.update_period = 0.1 ;
sim_pars.num_updates_per_output = 150 ;
sim_pars.num_lineages = 1 ;

% create the output folder
folder_name_splits = strsplit(par_file,'/') ;
folder_name = folder_name_splits{end} ;
folder_name = strsplit(folder_name,'.csv') ;
folder_name = folder_name{1} ;

% iterate on parameter sets
for i_sim = 1:size(par_data,1)

    
    % build the struct to pass to the simulation function
    par_struct = table2struct ( par_data(i_sim,:) ) ;
    names = [ fieldnames(par_struct) ; fieldnames(sim_pars) ] ;
    pars = cell2struct ( [struct2cell(par_struct); struct2cell(sim_pars)] , names , 1 ) ;
    
    % compute the lnm parameters
    pars.lnm_a = lnm_trends{'a_LNM','intercept'} + pars.mu * 60 * lnm_trends{'a_LNM','slope'} ;
    pars.lnm_b = pars.V_birth * ( 2 - pars.lnm_a ) ; % b such that average size at birth respected
    lnm_b_orig_unit = lnm_trends{'b_LNM','intercept'} + pars.mu * 60 * lnm_trends{'b_LNM','slope'} ;
    lnm_sigma_1_orig_unit = lnm_trends{'sig1_LNM','intercept'} + pars.mu * 60 * lnm_trends{'sig1_LNM','slope'} ;
    pars.lnm_sigma_1 = lnm_sigma_1_orig_unit * pars.lnm_b / lnm_b_orig_unit ; 
    pars.lnm_sigma_2 = lnm_trends{'sig2_LNM','intercept'} + pars.mu * 60 * lnm_trends{'sig2_LNM','slope'} ;
    pars.elong_rate_CV = lnm_trends{'alpha_CV','intercept'} + pars.mu * 60 * lnm_trends{'alpha_CV','slope'} ;    
        
    % do the simulation
    sim_data{i_sim} = do_single_sim ( pars ) ;
    
    % save pars 
    pars_data{i_sim} = pars ;
    
end
% create and write the table
save ( [ '../output-data/sim-data_' folder_name '.mat' ] , 'sim_data' , 'pars_data' ) ;

end
