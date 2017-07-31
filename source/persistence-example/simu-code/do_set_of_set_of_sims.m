function do_set_of_set_of_sims ( input_folder )

par_files = dir ( [ '../input-data/' input_folder '/*.csv' ] ) ;

% simulate all the corresponding sets
for i_set = 1:length(par_files)
    do_set_of_sims ( [ '../input-data/' input_folder '/' par_files(i_set).name ] ) ;
end

end

