function demands = compute_demands ()
%COMPUTE_DEMANDS Compute all the demands (i.e. P_conc_birth vs mu)

% vector of division rates
mu_vec_params = readtable ( '../input-data/mu_vector_params.csv' ) ;
mu_vec = linspace(mu_vec_params{1,'mu_min'},mu_vec_params{1,'mu_max'},mu_vec_params{1,'num_points'}) ; % in doublings per hr

% ref parameters
def_pars = give_ref_pars () ;

% P,R params
P_R_params = readtable ( '../input-data/P_R_demand_params.csv' ) ;

% Q demand
demands.Q = give_P_conc_birth(def_pars) .* ones(size(mu_vec)) ;

% P and R demand
for i_mu=1:length(mu_vec)
    mu = mu_vec(i_mu) / 60 ;
    demands.P(i_mu) = give_P_conc_birth(def_pars) * ( 1 - P_R_params{1,'slope'} * (mu-def_pars.mu) ) ;
    demands.R(i_mu) = give_P_conc_birth(def_pars) * ( 1 + P_R_params{1,'slope'} * (mu-def_pars.mu) ) ;
end

end

