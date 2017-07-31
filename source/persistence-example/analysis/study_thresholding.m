
%% load reference trace
sig = 0.325 ; rp = 1.3 ; i_sim = 5 ;
trace = readtable ( [ '../output-data/prot_mRNA_noise_exploration/sig-' num2str(sig) '_rp-' num2str(rp) '/parset-' num2str(i_sim) '_lineage-trajs.csv' ] ) ;
trace = trace.alpha_fast ;

%% study the threshold dependence of waiting times and number of changes
num_stays_min = 1 ;
thresh_vec = 0:0.0001:1 ;
avg_waiting_vec = zeros ( size(thresh_vec) ) ;
num_waitings_vec = zeros ( size(thresh_vec) ) ;
for i_t = 1:length(thresh_vec)
    waitings = compute_exit_times ( trace < thresh_vec(i_t) ) ;
    avg_waiting_vec(i_t) = mean ( waitings ) ;
    num_waitings_vec(i_t) = length ( waitings ) ;
end
I_num_wait_OK = find ( num_waitings_vec >= num_stays_min ) ;

%% plot stuff
subplot ( 2 , 2 , [1 2] ) ;
plot ( trace , 'b' ) ;
subplot ( 2 , 2 , 3 ) ;
semilogx ( thresh_vec(I_num_wait_OK) , num_waitings_vec(I_num_wait_OK) , 'b' ) ;
subplot ( 2 , 2 , 4 ) ;
semilogx ( thresh_vec(I_num_wait_OK) , avg_waiting_vec(I_num_wait_OK) , '-r' ) ;