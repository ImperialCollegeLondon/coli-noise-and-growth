

%% load data
clear ;
load ( '/Users/fbertaux/BoxSync/stuff-too-heavy-for-dropbox/coli-paper/persistence/model_3_km_slope_variation_data_large_scale_zoom' ) ;
clear i j ;

%% parameters for analysis
num_bins = 20 ;

%% analyse modality and rates
mod_res = cell ( size(mu_trajs) ) ;
num_done = 0 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        mod_res{i,j} = analyze_modality ( mu_trajs{i,j} , num_bins , time_traj(2) - time_traj(1) ) ;
        num_done = num_done + 1 ;
        disp ( ['sims done = ' num2str(num_done) ' / ' num2str( length(km_slope_vec) * length(mu_max_vec) )] ) ;
    end
end

%% plot ON occupancy heatmap
% heatmap
color_map = jet ( 100 ) ; marker_size = 20 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        if mod_res{i,j}.is_bimodal && mod_res{i,j}.is_bistable
            I_color = ceil ( mod_res{i,j}.P_on_from_rates * 100 ) ;
            color = color_map ( I_color , : ) ;
        else
            color = 0.4 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
set ( gcf , 'Position' , [ 273   143   672   559 ] ) ;
xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ;  set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/model_3_km_slope_variation_data_large_scale_zoom_ON_occupancy' ) ; close ;

%% plot ON rate heatmap
% heatmap
scale_max_on_rate = 0.1 ;
color_map = jet ( 100 ) ; marker_size = 20 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        if mod_res{i,j}.is_bimodal && mod_res{i,j}.is_bistable
            I_color = ceil ( mod_res{i,j}.rate_OFF_ON * 100 / scale_max_on_rate ) ;
            if I_color <= 100
            color = color_map ( I_color , : ) ;
            else
                color = [1 1 1] ;
            end
        else
            color = 0.4 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
set ( gcf , 'Position' , [ 273   143   672   559 ] ) ;
xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ; set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/model_3_km_slope_variation_data_large_scale_zoom_ON_rate' ) ; close ;


%% plot OFF rate heatmap
% heatmap
scale_max_off_rate = 0.25 ;
color_map = jet ( 100 ) ; marker_size = 20 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        if mod_res{i,j}.is_bimodal && mod_res{i,j}.is_bistable
            I_color = ceil ( mod_res{i,j}.rate_ON_OFF * 100 / scale_max_off_rate ) ;
            if I_color <= 100
            color = color_map ( I_color , : ) ;
            else
                color = [1 1 1] ;
            end
        else
            color = 0.4 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
set ( gcf , 'Position' , [ 273   143   672   559 ] ) ;
xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ;  set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/model_3_km_slope_variation_data_large_scale_zoom_OFF_rate' ) ; close ;

%% make a plot illustration what is km_slope
addpath ( '../simu-code' ) ;
num_slopes = 5 ;
color_map = jet ( num_slopes ) ;
small_km_slope_vec = linspace ( min(km_slope_vec) , max(km_slope_vec) , num_slopes ) ;
small_km_0_vec = km_ref - small_km_slope_vec .* mu2alpha(mu_max_ref) ;
for i=1:num_slopes
    plot ( [0 4] , small_km_0_vec(i) + small_km_slope_vec(i) .* mu2alpha ([0 4]) , 'Color' , color_map(i,:) ) ;
    hold on ;
end
xlim ( [0 4] ) ; xlabel ( '\mu_{cell} (dblgs/hr)' ) ; ylabel ( 'km_{cell} (min^{-1})' ) ;
set ( gca , 'FontSize' , 30 , 'LineWidth' , 2 ) ; set ( gcf , 'Position' , [ 273   143   672   559 ] , 'Color' , 'None' ) ;
export_fig ( '../plots/illustr_km_slope' ) ; close ;
