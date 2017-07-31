% model 2 <-> dilution only, Matthew like average of [toxin] (constant GE rates, constant lnm parameters, cell-cycle average effect)

% model type definition
pars.size_model = 'constant_lnm' ;
pars.GE_model = 'constant_rates' ;
pars.alpha_model = 'cell_cycle_average' ;
pars.partitioning_type = 'normal' ;

% GE rates
pars.km = 0.3 ;
pars.rm = 0.139 ;
pars.kp = 1.0 ;
pars.rp = 0.002 ;

% lnm parameters
pars.lnm_a = 1.0 ;
pars.lnm_b = 1.0 ;
pars.lnm_sigma_1 = 0.2 ;
pars.lnm_sigma_2 = 0.05 ;

% growth media and toxin parameters
pars.mu_max = 2.0 ;
pars.regulation_hill = 2.0 ;
pars.regulation_conc = 115 ;

% simulation parameters
pars.num_lineages = 1 ;
pars.sim_duration = 20000. * 60 ;
pars.update_period = 0.1 ;
pars.num_updates_per_output = 10 ;

% do the simulation
sim_data = do_single_sim ( pars ) ;

% plot result
subplot ( 4 , 1 , 1 ) ;
plot ( sim_data.traj_time(:,1) ./ 60 , sim_data.traj_alpha(:,1) .* 60 / log(2) , 'b' ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( '\mu_{cell} (dblgs/hr)' ) ;
ylim ( [0 3.0] ) ;
subplot ( 4 , 1 , 2 ) ;
plot ( sim_data.traj_time(:,1) ./ 60 , sim_data.traj_prot(:,1) ./ sim_data.traj_size(:,1) , 'r' ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( '[Toxin] (#/\mum^3)' ) ;
subplot ( 4 , 1 , 3 ) ;
plot ( sim_data.traj_time(:,1) ./ 60 , sim_data.traj_mRNA(:,1) , 'g' ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( 'Toxin mRNA' ) ;
subplot ( 4 , 1 , 4 ) ;
plot ( sim_data.traj_time(:,1) ./ 60 , sim_data.traj_size(:,1) , 'k' ) ;
xlabel ( 'Time (hours)' ) ; ylabel ( 'size' ) ;

% save plot
mkdir ( '../output-data/model-2/' ) ;
set ( gcf , 'Color' , 'None' ) ;
export_fig ( gcf , '../output-data/model-2/lineage-traj' ) ; close ;

% plot and save hist of alpha values
histogram ( sim_data.traj_alpha .* 60 / log(2) , 'Normalization' , 'probability' ) ;
xlim ( [-0.1 2.1] ) ;
xlabel ( '\mu_{cell} (dblgs/hr)' ) ; ylabel ( 'Frequency' ) ;
set ( gcf , 'Color' , 'None' ) ;
export_fig ( gcf , '../output-data/model-2/distrib-mu-cell' ) ; close ;
