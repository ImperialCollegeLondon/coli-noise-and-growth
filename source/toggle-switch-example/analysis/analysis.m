
%% define data sets to analyse
data_sets = { 'Q-response-km' , 'Q-response-kp' , 'klumpp-response' } ;
meta_data = cell ( size(data_sets) ) ;


%% do the analysis on all data sets
num_bins = 50 ;
for i_d = 1:length(data_sets)

    % load data
    clear sim_data pars_data ;
    load ( [ '/Users/fbertaux/BoxSync/stuff-too-heavy-for-dropbox/coli-paper/toggle-switch/sim-data_' ...
        data_sets{i_d} '.mat' ] ) ;

    % analyse modality for all sims
    modality_res = cell ( size(sim_data) ) ;
    for i_s = 1:length(sim_data)
        modality_res{i_s} = analyze_modality ( sim_data{i_s}.free_A_conc_traj , ...
            num_bins , sim_data{i_s}.time(2) - sim_data{i_s}.time(1) ) ;
    end

    % construct vectors of mus, P_ON and rates ..
    mu_vec = cellfun ( @(x)(x.mu * 60) , pars_data ) ;
    P_on_vec = cellfun ( @(x)(x.P_on_from_rates) , modality_res ) ;
    rate_on_vec = cellfun ( @(x)(x.rate_OFF_ON * 60) , modality_res ) ;
    rate_off_vec = cellfun ( @(x)(x.rate_ON_OFF * 60) , modality_res ) ;

    % save data
    meta_data{i_d}.mu_vec = mu_vec ;
    meta_data{i_d}.P_on_vec = P_on_vec ;
    meta_data{i_d}.rate_on_vec = rate_on_vec ;
    meta_data{i_d}.rate_off_vec = rate_off_vec ;

end


%% make meta plot
marker_size = 6 ; axes_size = 20 ; axes_lw = 1.5 ;
% P_on
subplot ( 1 , 3 , 1 ) ;
plot ( meta_data{1}.mu_vec , meta_data{1}.P_on_vec , '-ko' , 'MarkerFaceColor' , 'k' , 'MarkerSize' , marker_size ) ; hold on ;
plot ( meta_data{2}.mu_vec , meta_data{2}.P_on_vec , '-bo' , 'MarkerFaceColor' , 'b' , 'MarkerSize' , marker_size ) ; hold on ;
plot ( meta_data{3}.mu_vec , meta_data{3}.P_on_vec , '-ro' , 'MarkerFaceColor' , 'r' , 'MarkerSize' , marker_size ) ; hold on ;
xlabel ( 'Division rate (dblgs/hr)' ) ; ylabel ( 'ON state occupancy (%)' ) ; set ( gca , 'FontSize' , axes_size , 'LineWidth' , axes_lw ) ;
xlim ( [0 4] ) ;  ylim ( [0 0.7] ) ;
legend ( { 'P expression' , 'Q expression - k_m' , 'Q expression - k_p' } , 'FontSize' , 12 ) ;

% ON rate
subplot ( 1 , 3 , 2 ) ;
plot ( meta_data{1}.mu_vec , meta_data{1}.rate_on_vec , '-ko' , 'MarkerFaceColor' , 'k' , 'MarkerSize' , marker_size ) ; hold on ;
plot ( meta_data{2}.mu_vec , meta_data{2}.rate_on_vec , '-bo' , 'MarkerFaceColor' , 'b' , 'MarkerSize' , marker_size ) ; hold on ;
plot ( meta_data{3}.mu_vec , meta_data{3}.rate_on_vec , '-ro' , 'MarkerFaceColor' , 'r' , 'MarkerSize' , marker_size ) ; hold on ;
xlabel ( 'Division rate (dblgs/hr)' ) ; ylabel ( 'ON switch rate (hr^{-1})' ) ;  set ( gca , 'FontSize' , axes_size , 'LineWidth' , axes_lw ) ;
xlim ( [0 4] ) ;  ylim ( [0. 0.15] ) ;

% OFF rate
subplot ( 1 , 3 , 3 ) ;
plot ( meta_data{1}.mu_vec , meta_data{1}.rate_off_vec , '-ko' , 'MarkerFaceColor' , 'k' , 'MarkerSize' , marker_size ) ; hold on ;
plot ( meta_data{2}.mu_vec , meta_data{2}.rate_off_vec , '-bo' , 'MarkerFaceColor' , 'b' , 'MarkerSize' , marker_size ) ; hold on ;
plot ( meta_data{3}.mu_vec , meta_data{3}.rate_off_vec , '-ro' , 'MarkerFaceColor' , 'r' , 'MarkerSize' , marker_size ) ; hold on ;
xlim ( [0 4] ) ;  ylim ( [0. 0.5] ) ;
xlabel ( 'Division rate (dblgs/hr)' ) ; ylabel ( 'OFF switch rate (hr^{-1})' ) ; set ( gca , 'FontSize' , axes_size , 'LineWidth' , axes_lw ) ;


% write fig
set ( gcf , 'Color' , 'None' ) ; set ( gcf , 'Position' , [200         585        1582         415] ) ;
export_fig ( '../plots/meta_analysis.pdf' ) ; close ;


%% just do two traces of the ref plot
clear sim_data pars_data ;
load ( [ '/Users/fbertaux/BoxSync/stuff-too-heavy-for-dropbox/coli-paper/toggle-switch/sim-data_' ...
    data_sets{1} '.mat' ] ) ;
i_sim_ref = 6 ; num_tps = 800 ; burning_time = 100000 ; % in mins
i_burn = find ( sim_data{i_sim_ref}.time > burning_time , 1 , 'first' ) ;
plot ( (sim_data{i_sim_ref}.time(i_burn:i_burn+num_tps) - sim_data{i_sim_ref}.time(i_burn))./60 , sim_data{i_sim_ref}.free_A_conc_traj(i_burn:i_burn+num_tps) , 'b' , ...
    'LineWidth' , 3.0 ) ; hold on ;

% plot threshold
modality_res = analyze_modality ( sim_data{i_sim_ref}.free_A_conc_traj , ...
            num_bins , sim_data{i_sim_ref}.time(2) - sim_data{i_sim_ref}.time(1) ) ;
plot ( [0 200] , [1 1] .* modality_res.threshold , '--k' ) ;

xlabel ( 'Time (hours)' ) ; ylabel ( 'Free A concentration (#/\mum^3)' ) ;
set ( gca , 'FontSize' , 45 ) ;
set ( gca , 'LineWidth' , 3 ) ;
set ( gcf , 'Position' , [0 0 2000 1600] ) ;
set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/reference_traj.pdf' ) ; close ;
