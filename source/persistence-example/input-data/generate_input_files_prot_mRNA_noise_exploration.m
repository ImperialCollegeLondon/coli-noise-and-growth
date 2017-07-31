
% vectors of sigma and rp
sigmas = [ 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 ] ;
rps = [ 1 1.2 1.4 1.6 1.8 2.0 2.2 ] ;

% vectors of V_birth and EM_birth
V_birth_vec =  [10 50 100] ;
EM_birth_vec = [0.2 1 5] ;
num_sets = length ( V_birth_vec ) * length ( EM_birth_vec ) ;

for i_sig = 1:length(sigmas)
    for i_rp = 1:length(rps)
        
        % sigma and rp
        det.sigma = 10^(-sigmas(i_sig)) ;
        det.rp = 10^(-rps(i_rp)) ;
                
        % alpha_max, mRNA, prot degradation rate
        alpha_max = ones ( num_sets , 1 ) ;
        rm = 10 .* ones ( num_sets , 1 ) ;
        rp = det.rp .* ones ( num_sets , 1 ) ;
        
        % regulation params
        regulation_hill = 2 .* ones ( num_sets , 1 ) ;
        regulation_conc = ones ( num_sets , 1 ) ; % by definition of V_birth unit
        
        % other model parameters that vary
        kp = zeros ( num_sets , 1 ) ;
        km_intercept = zeros ( num_sets , 1 ) ;
        km_slope_alpha = zeros ( num_sets , 1 ) ;
        km_slope_size = zeros ( num_sets , 1 ) ;
        V_birth = ones ( num_sets , 1 ) ;
        EM_birth = ones ( num_sets , 1 ) ;
        for i_V = 1:length(V_birth_vec)
            for i_EM = 1:length(EM_birth_vec)
                I = sub2ind ( [length(V_birth_vec) length(EM_birth_vec)] , i_V , i_EM ) ;
                V_birth(I) = V_birth_vec (i_V) ;
                EM_birth(I) = EM_birth_vec(i_EM) ;
                km_slope_size(I) = EM_birth(I) * rm(I) / V_birth(I) ;
                kp(I) = rm(I) * det.sigma / km_slope_size(I) ;
            end
        end
        
        % partitioning parameter
        partitioning_type = cell ( num_sets , 1 ) ;
        partitioning_type(:) = {'normal'} ;
        
        % simulation parameters
        num_lineages = 2 .* ones ( num_sets , 1 ) ;
        sim_duration = 500 .* ones ( num_sets , 1 ) ;
        update_period = 0.05 .* ones ( num_sets , 1 ) ;
        num_updates_per_output = 10 .* ones ( num_sets , 1 ) ;
        
        % create and write the table with params
        param_table = table ( alpha_max , km_intercept , km_slope_size , km_slope_alpha , rm , kp , rp , V_birth , partitioning_type , regulation_hill , regulation_conc , num_lineages , sim_duration , update_period , num_updates_per_output , EM_birth ) ;
        mkdir ( 'prot_mRNA_noise_exploration' ) ;
        writetable ( param_table , [ 'prot_mRNA_noise_exploration/sig-' num2str(-log10(det.sigma)) '_rp-' num2str(-log10(det.rp)) '.csv' ] ) ;
        
        
    end
end





