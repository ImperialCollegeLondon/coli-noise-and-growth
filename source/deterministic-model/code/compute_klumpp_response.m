function pars = compute_klumpp_response ( ref_pars , mu )

pars = ref_pars ;
pars.mu = mu ;

pars.V_birth = give_V_birth_from_SJ (pars.mu*60) ;

addpath ( '../input-data' ) ;

pars.rm =  ref_pars.rm * give_klumpp_fitted_pars ( pars.mu*60 , 'rm' ) / give_klumpp_fitted_pars ( ref_pars.mu*60 , 'rm' ) ;

pars.km = ref_pars.km * give_klumpp_fitted_pars ( pars.mu*60 , 'km' ) / give_klumpp_fitted_pars ( ref_pars.mu*60 , 'km' ) ;

pars.kp = ref_pars.kp * give_klumpp_fitted_pars ( pars.mu*60 , 'kp' ) / give_klumpp_fitted_pars ( ref_pars.mu*60 , 'kp' ) ;

end

