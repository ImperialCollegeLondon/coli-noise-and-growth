function P_conc_birth = give_P_conc_birth ( pars )
%GIVE_P_CONC_BIRTH gives concentration of P at birth for given pars
%   For the deterministic model

% compute P_birth
P_birth = pars.km * pars.kp / pars.mu / pars.rm * ( 1 - pars.mu/pars.rm*(1-exp(-pars.rm/pars.mu))/(2-exp(-pars.rm/pars.mu)) ) ;

% compute P_conc_birth
P_conc_birth = P_birth / pars.V_birth ;

end

