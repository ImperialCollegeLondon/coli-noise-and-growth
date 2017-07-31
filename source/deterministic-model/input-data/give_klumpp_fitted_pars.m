function par_val = give_klumpp_fitted_pars ( doublings_per_hr , par_name )

% rm is return in per min !
% km and kp are relative (so compare to a reference)

if ( strcmp ( par_name , 'rm' ) )
    par_val = 0.5591 - 0.0443 * doublings_per_hr ;
    return ;
end

if ( strcmp ( par_name , 'kp' ) )
    par_val = 1.0144 - 0.0406 * doublings_per_hr ;
    return ;
end

if ( strcmp ( par_name , 'km' ) )
    g = 0.0006 + 1.1846 * exp ( 0.4403 * doublings_per_hr ) ;
    ksm_no_g = -1.7678 + 3.378 * exp ( 1.5568 * doublings_per_hr ) / ( 1 + exp ( 1.5568 * doublings_per_hr ) ) ;
    par_val = g * ksm_no_g ;
    return ;
end



end

