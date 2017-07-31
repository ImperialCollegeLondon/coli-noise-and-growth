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
pars.sim_duration = 50000. * 60 ;
pars.update_period = 0.1 ;
pars.num_updates_per_output = 100 ;

% do the simulation
sim_data_1 = do_single_sim ( pars ) ;

% change the alpha model
pars.alpha_model = 'cell_cycle_average' ;
sim_data_2 = do_single_sim ( pars ) ;

subplot ( 2 , 2 , 1 ) ;
plot ( sim_data_1.traj_time(:,1) ./ 60 , sim_data_1.traj_alpha(:,1) .* 60 / log(2) , 'Color' , [1 1 1] .* 0.5 ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( '\mu_{cell} (dblgs/hr)' ) ;
subplot ( 2 , 2 , 2 ) ;
plot ( sim_data_2.traj_time(:,1) ./ 60 , sim_data_2.traj_alpha(:,1) .* 60 / log(2) , 'Color' , [1 1 1] .* 0.5 ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( '\mu_{cell} (dblgs/hr)' ) ;
subplot ( 2 , 2 , 3 ) ;
h1 = histogram ( sim_data_1.traj_alpha .* 60 / log(2) , 'Normalization' , 'probability' ) ;
xlabel ( '\mu_{cell} (dblgs/hr)' ) ; ylabel ( 'Frequency' ) ; xlim ( [0 2]) ;
subplot ( 2 , 2 , 4 ) ;
h2 = histogram ( sim_data_2.traj_alpha .* 60 / log(2) , 'Normalization' , 'probability' ) ;
h2.BinWidth = h1.BinWidth ;
xlabel ( '\mu_{cell} (dblgs/hr)' ) ; ylabel ( 'Frequency' ) ; xlim ( [0 2]) ;
set ( gcf , 'Color' , 'w' ) ;


