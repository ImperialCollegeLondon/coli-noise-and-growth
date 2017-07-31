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
traj_mRNA_data = load ( [folder '/lineage_traj_mRNA_data.dat'] ) ;
traj_prot_data = load ( [folder '/lineage_traj_prot_data.dat'] ) ;
traj_size_data = load ( [folder '/lineage_traj_size_data.dat'] ) ;
traj_alpha_data = load ( [folder '/lineage_traj_alpha_data.dat'] ) ;
traj_time_data = load ( [folder '/lineage_traj_time_data.dat'] ) ;
if ( traj_time_data(end,1) == -1 ) % remove empty timepoint if needed
    traj_mRNA_data = traj_mRNA_data(1:end-1,:) ;
    traj_prot_data = traj_prot_data(1:end-1,:) ;
    traj_size_data = traj_size_data(1:end-1,:) ;
    traj_alpha_data = traj_alpha_data(1:end-1,:) ;
    traj_time_data = traj_time_data(1:end-1,:) ;
end
sim_data.traj_mRNA = traj_mRNA_data ;
sim_data.traj_prot = traj_prot_data ;
sim_data.traj_size = traj_size_data ;
sim_data.traj_alpha = traj_alpha_data ;
sim_data.traj_time = traj_time_data ;
% delete it to avoid problems with next calls
system ( ['rm ' folder '/*' ] ) ;
end

