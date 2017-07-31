
%% path (for analyze_modality function)
addpath ( '../analysis/') ;

%% model type definition
pars.size_model = 'constant_lnm' ;
pars.GE_model = 'dependent_rates' ;
pars.alpha_model = 'instantaneous' ;
pars.partitioning_type = 'normal' ;

%% model parameters that don't change
% Dependent GE parameters
pars.km_per_size = 0 ;
pars.rm = 0.139 ;
pars.kp = 0.939 ;
pars.rp = 0.001 ;
% lnm parameters
pars.lnm_a = 1.0 ;
pars.lnm_b = 1.0 ;
pars.lnm_sigma_1 = 0.2 ;
pars.lnm_sigma_2 = 0.05 ;
% growth media and toxin parameters
pars.regulation_hill = 2.0 ;
pars.regulation_conc = 140 ;
% simulation parameters
pars.num_lineages = 3 ;
sim_duration_ref = 50000. * 60 ;
pars.sim_duration = sim_duration_ref ;
pars.update_period = 0.1 ;
pars.num_updates_per_output = 100 ;

%% parameters that varies: mu_max and km_slope
mu_max_vec = linspace ( 1.5 , 2.5 , 25 ) ;
mu_max_ref = mu_max_vec ( round(length(mu_max_vec)/2) ) ;
km_ref = 0.279 ;
km_slope_vec = linspace ( -4 / log(2) , 4 / log(2) , 25 ) ;
km_0_vec = km_ref - km_slope_vec .* mu2alpha(mu_max_ref) ;

%% parameters of the exploration
n_bins = 20 ;

%% to store during explo
mod_res = cell ( length(mu_max_vec) , length(km_slope_vec) , pars.num_lineages ) ;

%% do all sims
% num_done = 0 ;
% for i=1:length(mu_max_vec)
%     for j=1:length(km_slope_vec)
%         % params
%         pars.mu_max = mu_max_vec(i) ;
%         pars.km_0 = km_0_vec(j) ;
%         pars.km_per_alpha = km_slope_vec(j) ;
%         % do sim, get alpha traj
%         sim_data = do_single_sim ( pars ) ;
%         % iterate on lineages
%         for k=1:pars.num_lineages
%             mu_traj = alpha2mu ( sim_data.traj_alpha(1:end,k) ) ;
%             if (i==1 && j==1 && k==1)
%                 time_traj = sim_data.traj_time(1:end) ./ 60 ;
%                 dt_hrs = time_traj(2) - time_traj(1) ;
%                 clear time_traj ;
%             end
%             % analyse modality
%             mod_res{i,j,k} = analyze_modality ( mu_traj , n_bins , dt_hrs ) ;
%         end
%         % for display progress
%         num_done = num_done + 1 ;
%         disp ( ['sims done = ' num2str(num_done) ' / ' num2str( length(km_slope_vec) * length(mu_max_vec) )] ) ;
%     end
% end

%% save the data
% clear sim_data mu_traj ;
load ( '../output-data/m3_kms_media_first_stage' ) ;

%% next stages: simulate longer when needed
should_continue = true ; round = 1 ;
while ( should_continue )
    pars.sim_duration = 2 * pars.sim_duration ;
    round = round + 1 ;
    num_redone = 0 ; num_done = 0 ;
    for i=1:length(mu_max_vec)
        for j=1:length(km_slope_vec)
            % do something only if at least one is bistable
            is_bimodals = cellfun ( @(x)(x.is_bimodal) , mod_res(i,j,:) ) ;
            [mb,I] = max ( is_bimodals ) ;
            if mb == 1 && mod_res{i,j,I}.is_bistable
                % compute CV is occupancy estimate
                rate_ONs = cellfun ( @(x)(x.rate_OFF_ON) , mod_res(i,j,:) ) ;
                rate_OFFs = cellfun ( @(x)(x.rate_ON_OFF) , mod_res(i,j,:) ) ;
                on_occupancies = rate_ONs ./ ( rate_OFFs + rate_ONs ) ;
                CV_on_occupancy = compute_CV ( on_occupancies ) ;
                % if too high, simulate new traces for longer
                if CV_on_occupancy > 0.3
                    CV_on_occupancy
                    pars.mu_max = mu_max_vec(i) ;
                    pars.km_0 = km_0_vec(j) ;
                    pars.km_per_alpha = km_slope_vec(j) ;
                    % do sim, get alpha traj
                    sim_data = do_single_sim ( pars ) ;
                    % iterate on lineages
                    for k=1:pars.num_lineages
                        mu_traj = alpha2mu ( sim_data.traj_alpha(1:end,k) ) ;
                        if (i==1 && j==1 && k==1)
                            time_traj = sim_data.traj_time(1:end) ./ 60 ;
                            dt_hrs = time_traj(2) - time_traj(1) ;
                            clear time_traj ;
                        end
                        % analyse modality
                        mod_res{i,j,k} = analyze_modality ( mu_traj , n_bins , dt_hrs ) ;
                    end
                    num_redone = num_redone + 1 ;
                end
            end
            % for display progress
            num_done = num_done + 1 ;
            disp ( [ 'Round ' num2str(round) ': redone/looked/tot = ' num2str(num_redone) ' / ' num2str(num_done) ' / ' num2str( length(km_slope_vec) * length(mu_max_vec) )] ) ;
        end
    end
    % decide if other round
    if ( num_redone == 0 || round > 4 )
        should_continue = false ;
    end
    % save after this round
    clear sim_data mu_traj ;
    save ( '../output-data/m3_kms_media' ) ;
end


%% plotting: ON occupancy heatmap
color_map = jet ( 100 ) ; marker_size = 25 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        
        % check bimodality
        is_bimodals = cellfun ( @(x)(x.is_bimodal) , mod_res(i,j,:) ) ;
        % if one replicate not detected as bimodal, color as white
        if min ( is_bimodals ) == 0
            color = [1 1 1] ;
        else
            % check rates
            rate_ONs = cellfun ( @(x)(x.rate_OFF_ON) , mod_res(i,j,:) ) ;
            rate_OFFs = cellfun ( @(x)(x.rate_ON_OFF) , mod_res(i,j,:) ) ;
            on_occupancies = rate_ONs ./ ( rate_OFFs + rate_ONs ) ;
            on_occupancy = mean ( on_occupancies ) ;
            CV_on_occupancy = compute_CV ( on_occupancies ) ;
            % if not a number (could not be detected) or CV > 30% -> grey
            if ( isnan(on_occupancy) || CV_on_occupancy > 0.3 )
                on_occupancies
                CV_on_occupancy
                color = 0.5 .* [1 1 1] ;
            else
                % color using colormap
                I_color = ceil ( on_occupancy * 100 ) ;
                color = color_map ( I_color , : ) ;
            end
        end
        
        % plot
        plot ( mu_max_vec(i) , km_slope_vec(j) * log(2) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
% xlim ( [1.4 2.6] ) ; ylim ( [-2.2 2.2] ) ;
set ( gcf , 'Position' , [1077         206         573         555] ) ;
% xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ; set ( gca , 'FontSize' , 15 ) ;  set ( gcf , 'Color' , 'None' ) ;
% export_fig ( '../plots/m3_kms_media_4_ON_occupancy' ) ; close ;

