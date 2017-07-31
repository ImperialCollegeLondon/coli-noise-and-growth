function pars = compute_response ( P_conc_birth_demand , ref_pars , par_to_change , size_change  )
%COMPUTE_RESPONSE Compute parameter change to match a new demand

pars = ref_pars ;

% if size change with div rate
if ( strcmp ( size_change , 'yes' ) )
    pars.V_birth = give_V_birth_from_SJ (pars.mu*60) ;
end

% km response
if ( strcmp ( par_to_change , 'km' ) )
    pars.km = P_conc_birth_demand * pars.V_birth / pars.kp * pars.mu * pars.rm ... 
        / ( 1 - pars.mu/pars.rm*(1-exp(-pars.rm/pars.mu))/(2-exp(-pars.rm/pars.mu)) ) ;
end

% kp response
if ( strcmp ( par_to_change , 'kp' ) )
    pars.kp = P_conc_birth_demand * pars.V_birth / pars.km * pars.mu * pars.rm ...
        / ( 1 - pars.mu/pars.rm*(1-exp(-pars.rm/pars.mu))/(2-exp(-pars.rm/pars.mu)) ) ;
end

end

