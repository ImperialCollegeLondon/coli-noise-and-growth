
% load orig klumpp response
orig_klumpp = readtable ( 'klumpp-response.csv' ) ;

% compute
addpath ( '../../deterministic-model/code' ) ;
P_conc_birth = zeros( size(orig_klumpp,1) , 1 ) ;
for i=1:size(orig_klumpp,1)
    pars.mu = orig_klumpp.mu(i) ;
    pars.km = orig_klumpp.km(i) ;
    pars.rm = orig_klumpp.rm(i) ;
    pars.kp = orig_klumpp.kp(i) ;
    pars.V_birth = orig_klumpp.V_birth(i) ;
    % get the ref params for the middle mu value
    if (i==round(size(orig_klumpp,1)/2))
        ref_pars = pars ;
    end
    P_conc_birth(i) = give_P_conc_birth (pars) ;
end
clear pars orig_pars;

% create new response
klumpp_no_size_change_km_affected = orig_klumpp ;
klumpp_no_size_change_kp_affected = orig_klumpp ;

% compute it
for i=1:size(orig_klumpp,1)
    
    % % % km
    
    % equalize to ref_pars except mu
    klumpp_no_size_change_km_affected.V_birth(i) = ref_pars.V_birth ;
    klumpp_no_size_change_km_affected.km(i) = ref_pars.km ;
    klumpp_no_size_change_km_affected.kp(i) = ref_pars.kp ;
    klumpp_no_size_change_km_affected.rm(i) = ref_pars.rm ;
    pars = klumpp_no_size_change_km_affected(i,:);
    new_pars = compute_response (P_conc_birth(i),pars,'km','no') ;    
    klumpp_no_size_change_km_affected.km(i) = new_pars.km ;
    clear pars ;
    
    % % % kp
    
    % equalize to ref_pars except mu
    klumpp_no_size_change_kp_affected.V_birth(i) = ref_pars.V_birth ;
    klumpp_no_size_change_kp_affected.km(i) = ref_pars.km ;
    klumpp_no_size_change_kp_affected.kp(i) = ref_pars.kp ;
    klumpp_no_size_change_kp_affected.rm(i) = ref_pars.rm ;
    pars = klumpp_no_size_change_kp_affected(i,:) ;
    new_pars = compute_response (P_conc_birth(i),pars,'kp','no') ;
    klumpp_no_size_change_kp_affected.kp(i) = new_pars.kp ;
    clear pars ;
    
end

% verify P_conc birt for km
P_conc_birth_no_size_change_km_affected = zeros( size(orig_klumpp,1) , 1 ) ;
for i=1:size(orig_klumpp,1)
    pars = klumpp_no_size_change_km_affected(i,:) ;
    P_conc_birth_no_size_change_km_affected(i) = give_P_conc_birth (pars) ;
    clear pars ;
end

% verify P_conc birt for kp
P_conc_birth_no_size_change_kp_affected = zeros( size(orig_klumpp,1) , 1 ) ;
for i=1:size(orig_klumpp,1)
    pars = klumpp_no_size_change_kp_affected(i,:) ;
    P_conc_birth_no_size_change_kp_affected(i) = give_P_conc_birth (pars) ;
    clear pars ;
end

% write the table
writetable(klumpp_no_size_change_km_affected,'klumpp_no_size_change_km_affected.csv') ;
writetable(klumpp_no_size_change_kp_affected,'klumpp_no_size_change_kp_affected.csv') ;