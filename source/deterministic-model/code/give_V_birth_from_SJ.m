function V_birth = give_V_birth_from_SJ ( mu )
%GIVE_ Simply parse the fitted trend parameters and apply them
%   mu in doublings per hour
%   V_birth in um^3

trend_pars = readtable ( '../input-data/SJ-average-data.txt' , 'ReadRowNames' , true ) ;
V_birth = trend_pars{'newborn_volume','factor'} * 2 ^ ( trend_pars{'newborn_volume','exponent'} * mu ) ;

end

