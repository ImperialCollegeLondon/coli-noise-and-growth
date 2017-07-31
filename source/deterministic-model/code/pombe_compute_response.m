function pars = pombe_compute_response ( P_conc_birth_demand , ref_pars , par_to_change , size_change  )
%COMPUTE_RESPONSE Compute parameter change to match a new demand

pars = ref_pars ;

% if size change with div rate
if ( strcmp ( size_change , 'yes' ) )
    error ( 'not yet implemented for pombe' ) ;
end

% km response
if ( strcmp ( par_to_change , 'km' ) )
    pars.km = P_conc_birth_demand * pars.V_birth / pars.kp * pars.mu * pars.rm ... 
        / ( 1 - pars.mu/pars.rm*(1-exp(-pars.rm/pars.mu))/(2-exp(-pars.rm/pars.mu)) ) ;
% kp resonse
elseif ( strcmp ( par_to_change , 'kp' ) )
    pars.kp = P_conc_birth_demand * pars.V_birth / pars.km * pars.mu * pars.rm ...
        / ( 1 - pars.mu/pars.rm*(1-exp(-pars.rm/pars.mu))/(2-exp(-pars.rm/pars.mu)) ) ;
end

end

