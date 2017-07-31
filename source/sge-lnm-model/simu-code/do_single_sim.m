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
num_tps = str2double ( sim_data.params.num_timepoints_per_cc ) ;
sim_data.prot_avg = zeros ( num_tps , 1 ) ;
sim_data.prot_CV = zeros ( num_tps , 1 ) ;
sim_data.mRNA_avg = zeros ( num_tps , 1 ) ;
sim_data.mRNA_CV = zeros ( num_tps , 1 ) ;
for tp=1:num_tps
    prot = load ( [ folder '/prot_data_' num2str(tp-1) '.dat' ] ) ;
    mRNA = load ( [ folder '/mRNA_data_' num2str(tp-1) '.dat' ] ) ;
    size = load ( [ folder '/size_data_' num2str(tp-1) '.dat' ] ) ;
    sim_data.prot_avg(tp) = mean ( prot(:,end) ) ;
    sim_data.prot_CV(tp) = compute_CV ( prot(:,end) ) ;
    sim_data.mRNA_avg(tp) = mean ( mRNA(:,end) ) ;
    sim_data.mRNA_CV(tp) = compute_CV ( mRNA(:,end) ) ;
    sim_data.size_avg(tp) = mean ( size(:,end) ) ;
    sim_data.size_CV(tp) = compute_CV ( size(:,end) ) ;
    sim_data.prot_conc_avg(tp) = mean ( prot(:,end) ./ size(:,end) ) ;
    sim_data.prot_conc_CV(tp) = compute_CV ( prot(:,end) ./ size(:,end) ) ;
end

% delete it to avoid problems with next calls
system ( ['rm ' folder '/*' ] ) ;
end
