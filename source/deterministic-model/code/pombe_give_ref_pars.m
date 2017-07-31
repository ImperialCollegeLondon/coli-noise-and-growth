function ref_pars = pombe_give_ref_pars ()
%GIVE_DEFAULT_PARS Gives the refault set of pars
%   Reads from input data
%   Time unit = minutes

% div rate (middle of the mu vector defined in input data)
mu_vec_params = readtable ( '../input-data/pombe_mu_vector_params.csv' ) ;
mu_vec = linspace(mu_vec_params{1,'mu_min'},mu_vec_params{1,'mu_max'},mu_vec_params{1,'num_points'}) ; % in doublings per hr
ref_pars.mu = mu_vec(round(length(mu_vec)/2)) / 60 ;

% read the parameters from the input data
ref_params_data = readtable ( '../input-data/pombe_reference_params.csv') ;

% mRNA params
ref_pars.rm = log(2) / ref_params_data{1,'mRNA_HL_mins'} ;
ref_pars.km = ref_params_data{1,'avg_mRNA_number_at_birth'} * ref_pars.rm * (2-exp(-ref_pars.rm/ref_pars.mu)) / (1-exp(-ref_pars.rm/ref_pars.mu)) ;

% protein synthesis
ref_pars.V_birth = 1 ; % just for having conc = copy number when computing response
ref_pars = pombe_compute_response ( ref_params_data{1,'avg_prot_number_at_birth'} , ref_pars , 'kp' , 'no' ) ;

% volume at birth
ref_pars.V_birth = 1 ;

end

