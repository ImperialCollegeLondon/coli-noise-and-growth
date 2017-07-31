function pars = give_ref_pars ()

% ref params gene expression
klumpp = readtable ( '../input-data/klumpp-response.csv' ) ;
pars.km = klumpp.km(6) ;
pars.kp = klumpp.kp(6) ;
pars.rm = klumpp.rm(6) ;
pars.V_birth = klumpp.V_birth(6) ;

% ref params for growth and division
pars.mu = 2 / 60 ;
lnm_trends = readtable ( '../input-data/linear_trends_lnm_params.csv' , 'ReadRowNames' , true ) ;
pars.lnm_a = lnm_trends{'a_LNM','intercept'} + pars.mu * 60 * lnm_trends{'a_LNM','slope'} ;
pars.lnm_b = pars.V_birth * ( 2 - pars.lnm_a ) ; % b such that average size at birth respected
lnm_b_orig_unit = lnm_trends{'b_LNM','intercept'} + pars.mu * 60 * lnm_trends{'b_LNM','slope'} ;
lnm_sigma_1_orig_unit = lnm_trends{'sig1_LNM','intercept'} + pars.mu * 60 * lnm_trends{'sig1_LNM','slope'} ;
pars.lnm_sigma_1 = lnm_sigma_1_orig_unit * pars.lnm_b / lnm_b_orig_unit ;
pars.lnm_sigma_2 = lnm_trends{'sig2_LNM','intercept'} + pars.mu * 60 * lnm_trends{'sig2_LNM','slope'} ;
pars.elong_rate_CV = lnm_trends{'alpha_CV','intercept'} + pars.mu * 60 * lnm_trends{'alpha_CV','slope'} ;
 
% params repression
pars.kon = 1.0 ;
pars.koff = 0.25 ;

% simulation params
pars.sim_duration = 50000 ;
pars.update_period = 0.1 ;
pars.num_updates_per_output = 100 ;
pars.num_lineages = 1 ;

end

