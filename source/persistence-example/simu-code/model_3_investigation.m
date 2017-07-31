
% vary km slope and look effect on SS distrib

% model type definition
pars.size_model = 'constant_lnm' ;
pars.GE_model = 'dependent_rates' ;
pars.alpha_model = 'instantaneous' ;
pars.partitioning_type = 'normal' ;

% Dependent GE parameters
pars.km_per_size = 0 ;
pars.rm = 0.139 ;
pars.kp = 0.939 ;
pars.rp = 0.001 ;

% lnm parameters
pars.lnm_a = 1.0 ;
pars.lnm_b = 1.0 ;
pars.lnm_sigma_1 = 0.2 ;
pars.lnm_sigma_2 = 0.05 ;

% growth media and toxin parameters
pars.mu_max = 2.0 ;
pars.regulation_hill = 2.0 ;
pars.regulation_conc = 140 ;

% simulation parameters
pars.num_lineages = 1 ;
pars.sim_duration = 50000. * 60 ;
pars.update_period = 0.1 ;
pars.num_updates_per_output = 100 ;

% vector to vary
km_mu_middle = 0.279 ;
km_0_vec = [0.9 1 1.1] .* km_mu_middle ;
km_per_alpha_vec = ( km_mu_middle - km_0_vec ) ./ ( pars.mu_max * log(2) / 60 / 2 ) ;

% create the structure getting the alpha distrib
edges = linspace ( 0 , 1. * pars.mu_max , 10 ) ;
centers = ( edges(1:end-1) + edges(2:end) ) ./ 2 ;
distribs = zeros ( length(edges)-1 , length(km_per_alpha_vec) ) ;
alpha_trajs = cell ( size(km_per_alpha_vec) ) ;


% vary the par and do sim each time
for i=1:length(km_per_alpha_vec)
    pars.km_0 = km_0_vec(i) ;
    pars.km_per_alpha = km_per_alpha_vec(i) ;
    sim_data = do_single_sim ( pars ) ;
    distribs(:,i) = histcounts ( sim_data.traj_alpha .* 60 ./ log(2) , edges , 'Normalization' , 'probability' ) ;
    alpha_trajs{i} = sim_data.traj_alpha(1:end) .* 60 ./ log(2) ;
    if i==1
        time_traj = sim_data.traj_time(1:end) ;
    end
end

% plot
for i=1:length(km_per_alpha_vec)
    
    subplot ( 3 , length(km_per_alpha_vec) , i ) ;
    pars.km_per_alpha = km_per_alpha_vec(i) ;
    pars.km_0 = km_0_vec(i) ;
    plot ( [0 , pars.mu_max] , [pars.km_0 , pars.km_0 + pars.mu_max * log(2) / 60 * pars.km_per_alpha] ) ; hold on ;
    plot ( [0 , pars.mu_max] , [1 1] .* km_mu_middle , '--k' ) ;
    ylim ( [0 1.2 * max(km_0_vec)] ) ;
    xlabel ( '\mu_{cell} (dblg/hr)' ) ; ylabel ( 'km_{cell}' ) ;
    
    subplot ( 3 , length(km_per_alpha_vec) , i + length(km_per_alpha_vec) ) ;
    bar ( centers , distribs(:,i) ) ; ylim ( [0 , 1] ) ;
    xlabel ( '\mu_{cell} (dblg/hr)' ) ; ylabel ( 'Frequency' ) ;
    
    subplot ( 3 , length(km_per_alpha_vec) , i + 2*length(km_per_alpha_vec) ) ;
    plot ( time_traj(1:30:end) ./ 60 , alpha_trajs{i}(1:30:end) ) ;
    ylabel ( '\mu_{cell} (dblg/hr)' ) ; xlabel ( 'Time (hrs)' ) ;

end
set ( gcf , 'Color' , 'w' ) ;




