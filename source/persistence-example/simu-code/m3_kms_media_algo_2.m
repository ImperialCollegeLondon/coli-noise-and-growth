
%% path (for analyze_modality function)
addpath ( '../analysis/') ;

%% model type definition
pars.size_model = 'constant_lnm' ;
pars.GE_model = 'dependent_rates' ;
pars.alpha_model = 'instantaneous' ;
pars.partitioning_type = 'normal' ;

%% exploration parameters
T_sim_hrs = 60000 ;
CV_limit = 0.1 ;
min_switches_first_sim = 10 ;
n_bins = 20 ;

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
pars.num_lineages = 1 ;
sim_duration_ref = T_sim_hrs * 60 ;
pars.sim_duration = sim_duration_ref ;
pars.update_period = 0.1 ;
pars.num_updates_per_output = 100 ;
dt_hrs = pars.update_period * pars.num_updates_per_output / 60. ;

%% parameters that varies: mu_max and km_slope
mu_max_vec = linspace ( 1.5 , 2.5 , 30 ) ;
mu_max_ref = mu_max_vec ( round(length(mu_max_vec)/2) ) ;
km_ref = 0.279 ;
km_slope_vec = linspace ( -3 / log(2) , 3 / log(2) , 30 ) ;
km_0_vec = km_ref - km_slope_vec .* mu2alpha(mu_max_ref) ;

%% do all sims
tic ;
if ( ~exist('rates_ON','var') )
    rates_ON =  zeros ( length(mu_max_vec) , length(km_slope_vec) ) ;
    rates_OFF =  zeros ( length(mu_max_vec) , length(km_slope_vec) ) ;
end
num_done = I_stop ;
done_since_last_save = 0 ;
for I = I_stop:length(mu_max_vec)*length(km_slope_vec)
    % transform in an index
    [i,j] = ind2sub(size(rates_ON),I)
    % do only if not done already
    if ( rates_ON(i,j) == 0 )
        % params
        pars.mu_max = mu_max_vec(i) ;
        pars.km_0 = km_0_vec(j) ;
        pars.km_per_alpha = km_slope_vec(j) ;
        pars.random_seed = 0 ;
        % do sim, get alpha traj
        sim_data = do_single_sim ( pars ) ;
        mu_traj = alpha2mu ( sim_data.traj_alpha(1:end) ) ;
        % analyse modality
        mod_res = analyze_modality ( mu_traj , n_bins , dt_hrs ) ;
        % continue only if bimodal and at least n_min switches
        if mod_res.n_switches > min_switches_first_sim
            waitings_ON = mod_res.waitings_ON ;
            waitings_OFF = mod_res.waitings_OFF ;
            threshold = mod_res.threshold ;
            % do more lineages until the CV of both waiting times is low
            should_continue = true ;
            while ( should_continue )
                % do additional lineage
                pars.random_seed = pars.random_seed + 1 ;
                sim_data = do_single_sim ( pars ) ;
                mu_traj = alpha2mu ( sim_data.traj_alpha(1:end) ) ;
                % compute waiting times (use same threshold as first sim !)
                [new_waitings_ON,new_waitings_OFF] = compute_waiting_times ( mu_traj,threshold ) ;
                waitings_ON = [ waitings_ON ; new_waitings_ON .* dt_hrs ] ;
                waitings_OFF = [ waitings_OFF ; new_waitings_OFF .* dt_hrs ] ;
                % compute CVs, decide to continue or not, save result
                CV_ON = compute_CV ( 1 ./ [mean(waitings_ON(1:3:end)) mean(waitings_ON(2:3:end)) mean(waitings_ON(3:3:end))] ) ;
                CV_OFF = compute_CV ( 1 ./ [mean(waitings_OFF(1:3:end)) mean(waitings_OFF(2:3:end)) mean(waitings_OFF(3:3:end))] ) ;
                disp ( [ 'CV_ON = ' num2str(CV_ON) ' , CV_OFF = ' num2str(CV_OFF) , ' threshold = ' num2str(threshold) ] ) ;
                if ( CV_ON > 0 && CV_OFF > 0 && CV_ON < CV_limit && CV_OFF < CV_limit )
                    should_continue = false ;
                    rates_ON(i,j) = 1 / mean(waitings_OFF) ;
                    rates_OFF(i,j) = 1 / mean(waitings_ON) ;
                end
            end
        else
            rates_ON(i,j) = -1 ;
            rates_OFF(i,j) = -1 ;
        end
    end
    % for display progress
    num_done = num_done + 1 ; done_since_last_save = done_since_last_save + 1 ;
    disp ( ['sims done = ' num2str(num_done) ' / ' num2str( length(km_slope_vec) * length(mu_max_vec) )] ) ;
    if ( done_since_last_save == 10 )
        I_stop = I ;
        save ( '../output-data/m3_kms_media_algo_2_ongoing' ) ;
        done_since_last_save = 0 ;
    end
end
toc ;

%% save the data
save ( '../output-data/m3_kms_media_algo_2' ) ;

%% plotting: ON rate heatmap
color_map = jet ( 100 ) ; marker_size = 20 ; scale =  0.15 ; 
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        if ( rates_ON(i,j) > 0 )
            I_color = min ( 100 , ceil ( rates_ON(i,j) * 100 / scale ) ) ;
            color = color_map ( I_color , : ) ;
        else
            color = 0.5 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) * log(2) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
set ( gcf , 'Position' , [2237 459 497 466] ) ;
xlabel ( '\mu_{max} (dblgs/hr)') ; ylabel ( 'km_{slope}' ) ;
ylim ( [ -3 3] ) ; xlim ( [1.5 2.5] ) ; 
set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/m3_kms_algo2_ON_rate.pdf' ) ; close ;
% make colorbar
scale_max = scale ; scale_min = 0 ;
for i=1:length(color_map)
    plot ( scale_min + i / 100 * (scale_max-scale_min)  , 0.5 , 's' , 'Color' , 'None' , 'MarkerSize' , 40 , 'MarkerFaceColor' , color_map(i,:) ) ; hold on ;
end
set ( gca , 'YTickLabel' , [] ) ; set ( gcf , 'Position' , [157 403 853 53] , 'Color' , 'None' ) ;
export_fig ( '../plots/colormap_ON_rate' ) ; close ;

%% plotting: OFF rate heatmap
color_map = jet ( 100 ) ; marker_size = 20 ; scale = 0.1 ; 
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        if ( rates_OFF(i,j) > 0 )
            I_color = min ( 100 , ceil ( rates_OFF(i,j) * 100 / scale ) ) ;
            color = color_map ( I_color , : ) ;
        else
            color = 0.5 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) * log(2) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
set ( gcf , 'Position' , [2237 459 497 466] ) ;
xlabel ( '\mu_{max} (dblgs/hr)') ; ylabel ( 'km_{slope}' ) ;
ylim ( [ -3 3] ) ; xlim ( [1.5 2.5] ) ; 
set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/m3_kms_algo2_OFF_rate.pdf' ) ; close ;
% make colorbar
scale_max = scale ; scale_min = 0 ;
for i=1:length(color_map)
    plot ( scale_min + i / 100 * (scale_max-scale_min)  , 0.5 , 's' , 'Color' , 'None' , 'MarkerSize' , 40 , 'MarkerFaceColor' , color_map(i,:) ) ; hold on ;
end
set ( gca , 'YTickLabel' , [] ) ; set ( gcf , 'Position' , [157 403 853 53] , 'Color' , 'None' ) ;
export_fig ( '../plots/colormap_OFF_rate' ) ; close ;


%% plotting: ON occupancy heatmap
color_map = jet ( 100 ) ; marker_size = 20 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        if ( rates_OFF(i,j) > 0 && rates_ON(i,j) > 0 )
            I_color = ceil ( rates_ON(i,j) * 100 / (rates_ON(i,j)+rates_OFF(i,j)) ) ;
            color = color_map ( I_color , : ) ;
        else
            color = 0.5 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) * log(2) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
set ( gcf , 'Position' , [2237 459 497 466] ) ;
xlabel ( '\mu_{max} (dblgs/hr)') ; ylabel ( 'km_{slope}' ) ;
ylim ( [ -3 3] ) ; xlim ( [1.5 2.5] ) ; 
set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/m3_kms_algo2_on_occupancy.pdf' ) ; close ;

