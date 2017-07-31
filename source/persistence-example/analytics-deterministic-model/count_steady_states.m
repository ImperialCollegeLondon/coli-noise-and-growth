function roots = count_steady_states ( pars )

RE = 1e-6 ;
ABS = 1e-10 ;

P_conc_min = pars.sigma / ( pars.rp + pars.alpha_max ) ;
P_conc_max = pars.sigma / pars.rp ;

root_1 = fzero ( @(P_conc)diff_nullclines(P_conc,pars) , P_conc_min ) ;
root_2 = fzero ( @(P_conc)diff_nullclines(P_conc,pars) , P_conc_max ) ;

if ( compare2roots(root_1,root_2,RE,ABS) )
    roots.num_roots = 1 ;
    roots.values = root_1 ;
    return ;
end

if ( isnan(root_1) || isnan(root_2) )
    roots.num_roots = 1 ;
    roots.values = [root_1 root_2] ;
    roots.values(isnan(roots.values)) = [] ;
    return ;
end

if ( diff_nullclines(root_1*(1+RE),pars) * diff_nullclines(root_2*(1-RE),pars) > 0 )
    roots.num_roots = 2 ;
    roots.values = [root_1 , root_2] ;
    return ;
end

root_3 = fzero ( @(P_conc)diff_nullclines(P_conc,pars) , [root_1*(1+RE) , root_2*(1-RE)] ) ;
if ( compare2roots(root_1,root_3,RE,ABS) )
    roots.num_roots = 2 ;
    roots.values = [root_1 , root_2] ;
    return ;
end

if ( compare2roots(root_2,root_3,RE,ABS) )
    roots.num_roots = 2 ;
    roots.values = [root_1 , root_2] ;
    return ;
end

roots.num_roots = 3 ;
roots.values = sort ( [root_1 , root_2 , root_3] ) ;

end

function comparison = compare2roots ( r1 , r2 , RE , ABS )
    comparison = abs(r1-r2) < max ( RE * max(r1,r2) , ABS ) ;
end

function y = diff_nullclines(P_conc,pars)
    alpha_toxin = pars.alpha_max / ( 1 + P_conc ^ pars.hill ./ pars.P_toxin ^ pars.hill ) ;
    alpha_expression = ( pars.sigma - pars.rp * P_conc ) ./ P_conc ;
    y = alpha_toxin - alpha_expression ;
end

