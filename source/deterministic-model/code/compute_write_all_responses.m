
%% vector of division rates
mu_vec_params = readtable ( '../input-data/mu_vector_params.csv' ) ;
mu_vec = linspace(mu_vec_params{1,'mu_min'},mu_vec_params{1,'mu_max'},mu_vec_params{1,'num_points'}) ; % in doublings per hr

%% ref parameters
def_pars = give_ref_pars () ;

%% compute demands
demands = compute_demands () ;

%% compute responses
response_types = { 'km' , 'kp' } ;
demand_types = fieldnames(demands) ;
% iterate on demand types
for i_demand = 1:length(demand_types)
    % iterate on response types
    for i_response = 1:length(response_types)
        demand_vector = demands.(demand_types{i_demand}) ;
        % iterate on division rate values
        for i_mu = 1:length(mu_vec)
            pars = def_pars ;
            pars.mu = mu_vec(i_mu) / 60 ;
            % compute and store the response
            responses.(demand_types{i_demand}).([response_types{i_response} '_no_size_change'])(i_mu) = ...
                compute_response ( demand_vector(i_mu) , pars , response_types{i_response} , 'no' ) ;
            responses.(demand_types{i_demand}).(response_types{i_response})(i_mu) = ...
                compute_response ( demand_vector(i_mu) , pars , response_types{i_response} , 'yes' ) ;
        end
    end
end

%% write the responses
% iterate on demand types
for i_demand = 1:length(demand_types)
    % iterate on response types
    response_types = fieldnames ( responses.(demand_types{i_demand}) ) ;
    for i_response = 1:length(response_types)
        mu = [responses.(demand_types{i_demand}).(response_types{i_response}).mu]' ;
        km = [responses.(demand_types{i_demand}).(response_types{i_response}).km]' ;
        rm = [responses.(demand_types{i_demand}).(response_types{i_response}).rm]' ;
        kp = [responses.(demand_types{i_demand}).(response_types{i_response}).kp]' ;
        V_birth = [responses.(demand_types{i_demand}).(response_types{i_response}).V_birth]' ;
        prot_conc_demand = demands.(demand_types{i_demand})' ;
        result_table = table ( mu , km , rm , kp , V_birth , prot_conc_demand ) ;
        writetable ( result_table , ['../output-data/demand-' demand_types{i_demand} ...
            '_response-' response_types{i_response} '.csv' ] ) ;
    end
end

%% make the plots

% demand functions
for i=1:3
    subplot(2,3,i) ; 
    plot ( mu_vec , demands.(demand_types{i}) , '-ko' ) ;
    ylim ( [ 0 120 ] ) ;
    set ( gca , 'FontSize' , 20 ) ;
    title ( [ demand_types{i} ' demand' ] , 'FontSize' , 20 ) ;
    xlabel ( 'Doublings per hour' ) ; 
    ylabel ( 'Prot. conc. at birth (#/\mum^3)' ) ;
end

% response: no size change
for i=1:3
    subplot(2,3,3+i) ;
    semilogy ( mu_vec , [responses.(demand_types{i}).km_no_size_change.km] ./ def_pars.km , '-bo' ) ; hold on ;
    semilogy ( mu_vec , [responses.(demand_types{i}).km.km] ./ def_pars.km , '-ro' ) ; hold on ;
    ax = axis () ; semilogy ( ax(1:2) , [1 1] , 'k' ) ;
    ylim ( [ 0.01 100 ] ) ;
    set ( gca , 'FontSize' , 20 ) ;
    title ( [ 'km or kp response to ' demand_types{i} ' demand' ] , 'FontSize' , 20 ) ;
    xlabel ( 'Doublings per hour' ) ; 
    ylabel ( 'Relative change of k_m (or k_p)' ) ;
end
subplot ( 2 , 3 , 4 ) ; 
legend ( { 'no size change' , 'empirical size change' } , 'FontSize' , 15 , 'Location' , 'NorthWest' ) ;
set ( gcf , 'Position' , [0 0 1600 1000] ) ; set ( gcf , 'Color' , 'w' ) ;
export_fig ( gcf , '../plots/Q_P_R_demands_km_kp_responses.pdf' ) ; close ;


