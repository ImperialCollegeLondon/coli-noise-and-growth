
% model 1 (rocco like) <-> dilution only (constant GE rates, constant lnm parameters, instantaneous effect)

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
pars.sim_duration = 20000. * 60 ;
pars.update_period = 0.1 ;
pars.num_updates_per_output = 100 ;

% do the simulation
sim_data = do_single_sim ( pars ) ;

% plot result
subplot ( 4 , 1 , 1 ) ;
plot ( sim_data.traj_time(:,1) ./ 60 , sim_data.traj_alpha(:,1) .* 60 / log(2) , 'b' ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( '\mu_{cell} (dblgs/hr)' ) ; xlim ( [0 8000] ) ;
ylim ( [0 3.0] ) ;
subplot ( 4 , 1 , 2 ) ;
plot ( sim_data.traj_time(:,1) ./ 60 , sim_data.traj_prot(:,1) ./ sim_data.traj_size(:,1) , 'r' ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( '[Toxin] (#/\mum^3)' ) ; xlim ( [0 8000] ) ; 
subplot ( 4 , 1 , 3 ) ;
plot ( sim_data.traj_time(:,1) ./ 60 , sim_data.traj_mRNA(:,1) , 'g' ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( 'Toxin mRNA' ) ; xlim ( [0 8000] ) ;
subplot ( 4 , 1 , 4 ) ;
plot ( sim_data.traj_time(:,1) ./ 60 , sim_data.traj_size(:,1) , 'k' ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( 'size' ) ; xlim ( [0 8000] ) ;

% save plot
mkdir ( '../output-data/model-1/' ) ;
set ( gcf , 'Color' , 'None' ) ;
export_fig ( gcf , '../output-data/model-1/lineage-traj' ) ; close ;

% plot and save hist of alpha values
histogram ( sim_data.traj_alpha .* 60 / log(2) , 'Normalization' , 'probability' ) ;
xlim ( [-0.1 2.1] ) ;
xlabel ( '\mu_{cell} (dblgs/hr)' ) ; ylabel ( 'Frequency' ) ;
set ( gcf , 'Color' , 'None' ) ;
export_fig ( gcf , '../output-data/model-1/distrib-mu-cell' ) ; close ;


