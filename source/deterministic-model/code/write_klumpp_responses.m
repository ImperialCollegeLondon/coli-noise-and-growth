
%% vector of division rates
mu_vec_params = readtable ( '../input-data/mu_vector_params.csv' ) ;
mu_vec = linspace(mu_vec_params{1,'mu_min'},mu_vec_params{1,'mu_max'},mu_vec_params{1,'num_points'}) ; % in doublings per hr

%% ref parameters
def_pars = give_ref_pars () ;

%% compute responses
for i_mu = 1:length(mu_vec)
    responses(i_mu) = compute_klumpp_response ( def_pars , mu_vec(i_mu) / 60 ) ;
end

%% write the responses
mu = [responses.mu]' ;
km = [responses.km]' ;
rm = [responses.rm]' ;
kp = [responses.kp]' ;
V_birth = [responses.V_birth]' ;
result_table = table ( mu , km , rm , kp , V_birth ) ;
writetable ( result_table , '../output-data/klumpp-response.csv' ) ;