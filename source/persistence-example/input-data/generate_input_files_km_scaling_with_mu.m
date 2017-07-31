% sigma and rp
det.sigma = 10^(-0.25) ;
det.rp = 10^(-1) ;

% vector of division rates
mu_vec_params = readtable ( '../../deterministic-model/input-data/mu_vector_params.csv' ) ;
mu = linspace(mu_vec_params{1,'mu_min'},mu_vec_params{1,'mu_max'},mu_vec_params{1,'num_points'})' ; % in doublings per hr
num_sets = length(mu) ;
i_ref = round ( num_sets/2 ) ;

% 'ref' Vb and EM
V_birth = 50 .* ones ( num_sets , 1 ) ;
EM_birth = 1 .* ones ( num_sets , 1 ) ;

% alpha_max, mRNA, prot degradation rate
alpha_max = mu ./ mu(i_ref) ;
rm = 10 .* ones ( num_sets , 1 ) ;
rp = det.rp .* ones ( num_sets , 1 ) ;

% regulation params
regulation_hill = 2 .* ones ( num_sets , 1 ) ;
regulation_conc = ones ( num_sets , 1 ) ; % by definition of V_birth unit

% other model parameters that vary
km_intercept = zeros ( num_sets , 1 ) ;
km_slope_alpha = zeros ( num_sets , 1 ) ;

km_slope_size = ones ( num_sets , 1) .* ( EM_birth(i_ref) * rm(i_ref) / V_birth(i_ref) ) .* alpha_max ;
kp = ones ( num_sets , 1) .* ( rm(i_ref) * det.sigma / km_slope_size(i_ref) ) ;

% partitioning parameter
partitioning_type = cell ( num_sets , 1 ) ;
partitioning_type(:) = {'normal'} ;

% simulation parameters
num_lineages = 2 .* ones ( num_sets , 1 ) ;
sim_duration = 3000 .* ones ( num_sets , 1 ) ;
update_period = 0.05 .* ones ( num_sets , 1 ) ;
num_updates_per_output = 10 .* ones ( num_sets , 1 ) ;

% create and write the table with params
param_table = table ( alpha_max , km_intercept , km_slope_size , km_slope_alpha , rm , kp , rp , V_birth , partitioning_type , regulation_hill , regulation_conc , num_lineages , sim_duration , update_period , num_updates_per_output , EM_birth ) ;
mkdir ( 'km_scaling_with_mu' ) ;
writetable ( param_table , [ 'km_scaling_with_mu/set_1.csv' ] ) ;
