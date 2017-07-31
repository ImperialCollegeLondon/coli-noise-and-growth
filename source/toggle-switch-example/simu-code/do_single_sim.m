function sim_data = do_single_sim ( pars )

% cpp whould be called
path2cpp = 'cpp/build/simulator.app/Contents/MacOS/simulator' ;
pars.out_folder = 'cpp/sim-data' ;

% call the cpp code
call_cpp_program ( path2cpp , pars ) ;

% read results
sim_data = read_sim_data ( pars.out_folder ) ;

end

function sim_data = read_sim_data (folder)
% load the data
sim_data.params = parse_metadata ( [ folder '/params.txt' ] ) ;
free_A_traj = load ( [ folder '/lineage_traj_prot_data.dat' ] ) ;
% mRNA_traj =  load ( [ folder '/lineage_traj_mRNA_data.dat' ] ) ;
size_traj =  load ( [ folder '/lineage_traj_size_data.dat' ] ) ;
sim_data.free_A_conc_traj = free_A_traj ./ size_traj ;
sim_data.output_period = str2double ( sim_data.params.update_period ) * str2double ( sim_data.params.num_updates_per_output ) ;
sim_data.time = linspace ( 0 , sim_data.output_period * (length(sim_data.free_A_conc_traj)-1) , length(sim_data.free_A_conc_traj) ) ;
% delete it to avoid problems with next calls
system ( ['rm ' folder '/*' ] ) ;
end

