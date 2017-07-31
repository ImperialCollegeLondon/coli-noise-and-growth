function draw_nullclines ( pars , number_sets , current_set )

P_max = pars.sigma / pars.rp ;
P_vec = logspace ( -2 , log10(P_max) , 1000 ) ;
alpha_toxin = pars.alpha_max ./ ( 1 + P_vec .^ pars.hill ./ pars.P_toxin ^ pars.hill ) ;
alpha_expression = ( pars.sigma - pars.rp .* P_vec ) ./ P_vec ;


subplot ( number_sets , 2 , 2 * (current_set - 1) + 1 ) ;
semilogx ( P_vec , alpha_toxin , 'r' ) ; hold on ;
semilogx ( P_vec , alpha_expression , 'g' ) ; hold on ;
xlabel ( 'P_{toxin}' ) ; ylabel ( '\alpha' ) ;
ylim ( [0 pars.alpha_max] ) ;
subplot ( number_sets , 2 , 2 * current_set ) ;
semilogx ( P_vec , alpha_expression - alpha_toxin , 'b' ) ; hold on ;
plot ( P_vec , zeros ( size(P_vec) ) , 'k' ) ;
xlabel ( 'P_{toxin}' ) ; ylabel ( '\alpha_{expression}-\alpha_{toxin}' ) ;
ylim ( [ -0.05 0.05 ] ) ;

end

