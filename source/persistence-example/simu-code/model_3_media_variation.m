
% model type definition
pars.size_model = 'constant_lnm' ;
pars.GE_model = 'dependent_rates' ;
pars.alpha_model = 'instantaneous' ;
pars.partitioning_type = 'normal' ;

% parameters that varies: mu_max
mu_max_vec = [ 1.75 2.0 2.25 ] ;

km_0_vec = [1 1 1] .* 0.3 ;
km_per_alpha_vec = [1 1 1] .* (- 5) ;

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

pars.regulation_hill = 2.0 ;
pars.regulation_conc = 140 ;

% simulation parameters
pars.num_lineages = 1 ;
pars.sim_duration = 20000. * 60 ;
pars.update_period = 0.1 ;
pars.num_updates_per_output = 100 ;

% create the structure getting the alpha distrib
edges = linspace ( 0 , 1 , 10 ) ;
centers = ( edges(1:end-1) + edges(2:end) ) ./ 2 ;
distribs = zeros ( length(edges)-1 , length(mu_max_vec) ) ;
alpha_trajs = cell ( size(mu_max_vec) ) ;

% vary the par and do sim each time
for i=1:length(mu_max_vec)
    pars.mu_max = mu_max_vec(i) ;
    pars.km_0 = km_0_vec(i) ;
    pars.km_per_alpha = km_per_alpha_vec(i) ;
    sim_data = do_single_sim ( pars ) ;
    distribs(:,i) = histcounts ( sim_data.traj_alpha .* 60 ./ log(2) ./ pars.mu_max , edges , 'Normalization' , 'probability' ) ;
    alpha_trajs{i} = sim_data.traj_alpha(1:end) .* 60 ./ log(2) ;
    if i==1
        time_traj = sim_data.traj_time(1:end) ;
    end
end

%% plot
for i=1:length(mu_max_vec)
    
    subplot ( 3 , length(mu_max_vec) , i ) ;
    plot ( [0 , mu_max_vec(end)] , [km_0_vec(i) , km_0_vec(i) + mu_max_vec(end) * log(2) / 60 * km_per_alpha_vec(i)] ) ; hold on ;
    h = plot ( mu_max_vec(i) , km_0_vec(i) + km_per_alpha_vec(i) * mu_max_vec(i) * log(2) / 60 , 'ro' , 'MarkerFaceColor' , 'r' , 'MarkerSize' , 15 ) ;
%     ylim ( [0 1.2 * ( max(km_0_vec) + max(km_per_alpha_vec .* mu_max_vec) * log(2) / 60 )] ) ; xlim ( [0 max(mu_max_vec)] ) ;
    xlabel ( '\mu_{cell} (dblg/hr)' ) ; 
    if (i==1) ; ylabel ( 'km_{cell}' ) ; end
    if (i==1) ; legend ( h , {'\mu_{max} (media)'} , 'Location' , 'SouthEast' , 'FontSize' , 15 ) ; end
    
    subplot ( 3 , length(mu_max_vec) , i + length(mu_max_vec) ) ;
    bar ( centers , distribs(:,i) ) ; ylim ( [0 , 1] ) ;
    xlabel ( '\mu_{cell} / \mu_{max}' ) ; 
    if (i==1) ; ylabel ( 'Frequency' ) ; end
    
    subplot ( 3 , length(mu_max_vec) , i + 2*length(mu_max_vec) ) ;
    plot ( time_traj(1:30:end) ./ 60 , alpha_trajs{i}(1:30:end) ./ mu_max_vec(i) ) ;
    xlabel ( 'Time (hrs)' ) ; ylim ([0 1]) ;
    if (i==1) ; ylabel ( '\mu_{cell} / \mu_{max}' ) ; end

end
subplot ( 3 , length(mu_max_vec) , 1 ) ;
set ( gcf , 'Color' , 'w' ) ;

