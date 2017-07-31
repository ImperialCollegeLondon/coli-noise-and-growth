function do_all_sets_of_sims ( file_pattern )

% find all par files
par_files = dir ( [ '../input-data/' file_pattern '.csv' ] ) ;

% simulate all the corresponding sets
for i_set = 1:length(par_files)
    par_data = [ '../input-data/' par_files(i_set).name ] ;
    result_table = do_set_of_sims ( par_data ) ;
    writetable ( result_table , [ '../output-data/sim-result_' par_files(i_set).name ] ) ;
end

end