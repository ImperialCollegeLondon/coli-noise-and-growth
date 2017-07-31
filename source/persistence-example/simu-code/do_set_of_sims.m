function do_set_of_sims ( par_file )
%DO_SET_OF_SIMS Simulate all param set in a csv file
%   par_file in input-data

% read the parameter data
par_data = readtable ( par_file ) ;

% create the output folder
folder_name_splits = strsplit(par_file,'/') ;
folder_name = folder_name_splits{end} ;
folder_name = strsplit(folder_name,'.csv') ;
folder_name = folder_name{1} ;
output_parent_folder = folder_name_splits{end-1} ;
mkdir ( [ '../output-data/' output_parent_folder '/' folder_name ] ) ;

% iterate on parameter sets
for i_sim = 1:size(par_data,1)

    % build the struct to pass to the simulation function
    par_struct = table2struct ( par_data(i_sim,:) ) ;
    
    % do the simulation
    sim_data = do_single_sim ( par_struct ) ;

    % construct the table lineage trajs
    size_fast = sim_data.traj_size(:,1) ;
    toxin_number_fast = sim_data.traj_prot(:,1) ;
    toxin_conc_fast = toxin_number_fast ./ size_fast ;
    alpha_fast = sim_data.traj_alpha(:,1) ;
    size_slow = sim_data.traj_size(:,2) ;
    toxin_number_slow = sim_data.traj_prot(:,2) ;
    toxin_conc_slow = toxin_number_slow ./ size_slow ;
    alpha_slow = sim_data.traj_alpha(:,2) ;
    
    % create and write the table
    result = table ( size_fast , toxin_number_fast , toxin_conc_fast , alpha_fast , size_slow , toxin_number_slow , toxin_conc_slow , alpha_slow ) ;
    writetable( result , [ '../output-data/' output_parent_folder '/' folder_name '/parset-' num2str(i_sim) '_lineage-trajs.csv' ] ) ;
    
end

end
