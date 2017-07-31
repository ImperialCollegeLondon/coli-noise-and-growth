


pars(1).alpha_max = 1 ;
pars(1).hill = 2 ;
pars(1).P_toxin = 1 ;

i_pars = 1 ;
pars(i_pars).sigma = 10^(-0.35) ;
pars(i_pars).rp = 10^(-2.1) ;
i_pars = i_pars + 1 ;

pars(i_pars) = pars(i_pars-1) ;
pars(i_pars).sigma = 10^(-0.25) ;
pars(i_pars).rp = 10^(-1.1) ;
i_pars = i_pars + 1 ;

pars(i_pars) = pars(i_pars-1) ;
pars(i_pars).sigma = 10^(-0.25) ;
pars(i_pars).rp = 10^(-1.1) ;
i_pars = i_pars + 1 ;

pars(i_pars) = pars(i_pars-1) ;
pars(i_pars).sigma = 10^(-0.25) ;
pars(i_pars).rp = 10^(-1.1) ;
i_pars = i_pars + 1 ;

pars(i_pars) = pars(i_pars-1) ;
pars(i_pars).sigma = 10^(-0.25) ;
pars(i_pars).rp = 10^(-1.1) ;
i_pars = i_pars + 1 ;

% sims done
pars(3) = pars(1) ;
pars(3).sigma = 10^(-0.4) ;
pars(3).rp = 10^(-1.8) ;

% sims done
pars(3) = pars(1) ;
pars(3).sigma = 10^(-0.4) ;
pars(3).rp = 10^(-1.8) ;
