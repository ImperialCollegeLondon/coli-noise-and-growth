% model type definition
pars.size_model = 'constant_lnm' ;
pars.GE_model = 'constant_rates' ;
pars.alpha_model = 'instantaneous' ;
pars.partitioning_type = 'normal' ;

% GE rates (per min)
pars.km = 0.279 ;
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
pars.regulation_conc = 140 ; % # per size

% simulation parameters
pars.num_lineages = 1 ;
pars.sim_duration = 10000. * 60 ;
pars.update_period = 0.1 ;
pars.num_updates_per_output = 100 ;

% varying param
rp_vec = [0.0006 0.0008 0.001] ;

% create the structure getting the alpha distrib
edges = linspace ( 0 , 1. * pars.mu_max , 10 ) ;
centers = ( edges(1:end-1) + edges(2:end) ) ./ 2 ;
distribs = zeros ( length(edges)-1 , length(rp_vec) ) ;
alpha_trajs = cell ( size(rp_vec) ) ;


% vary the par and do sim each time
for i=1:length(rp_vec)
    pars.rp = rp_vec(i) ;
    sim_data = do_single_sim ( pars ) ;
    distribs(:,i) = histcounts ( sim_data.traj_alpha .* 60 ./ log(2) , edges , 'Normalization' , 'probability' ) ;
    alpha_trajs{i} = sim_data.traj_alpha(1:end) .* 60 ./ log(2) ;
    if i==1
        time_traj = sim_data.traj_time(1:end) ;
    end
end

% plot
for i=1:length(rp_vec)
    
    subplot ( 2 , length(rp_vec) , i ) ;
    bar ( centers , distribs(:,i) ) ; ylim ( [0 , 1] ) ;
    text ( 0.8 , 0.8 , ['rp = ' num2str(rp_vec(i)) ' min^{-1}' ] , 'FontSize' , 20 ) ;
    xlabel ( '\mu_{cell} (dblg/hr)' ) ; ylabel ( 'Frequency' ) ;
    
    subplot ( 2 , length(rp_vec) , i + length(rp_vec) ) ;
    plot ( time_traj(1:30:end) ./ 60 , alpha_trajs{i}(1:30:end) ) ;
    ylabel ( '\mu_{cell} (dblg/hr)' ) ; xlabel ( 'Time (hrs)' ) ;

end
set ( gcf , 'Color' , 'w' ) ;


