

%% load data
clear ;
load ( '../output-data/model_3_km_slope_variation_data_large_scale' ) ;
clear i j ;

%% parameters for analysis
mu_edges = linspace ( 0 , mu_max_vec(end) , 20 ) ;
mu_centers = ( mu_edges(1:end-1)+mu_edges(2:end) ) ./ 2 ;

%% analyse modality and rates
mod_res = cell ( size(mu_trajs) ) ;
num_done = 0 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        mod_res{i,j} = analyze_modality ( mu_trajs{i,j} , mu_edges , time_traj(2) - time_traj(1) ) ;
        num_done = num_done + 1 ;
        disp ( ['sims done = ' num2str(num_done) ' / ' num2str( length(km_slope_vec) * length(mu_max_vec) )] ) ;
    end
end

%% plot ON occupancy heatmap
% heatmap
color_map = jet ( 100 ) ; marker_size = 40 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        if mod_res{i,j}.is_bistable
            I_color = ceil ( mod_res{i,j}.P_on_from_rates * 100 ) ;
            color = color_map ( I_color , : ) ;
        else
            color = 0.4 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
xlim ( [-0.2 0.2] + minmax(mu_max_vec) ) ; ylim ( [-0.3 0.3] + minmax(km_slope_vec) ) ;
set ( gcf , 'Position' , [ 1 1 900 710 ] ) ;
xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ;
export_fig ( '../plots/model_3_km_slope_variation_data_large_scale_1_ON_occupancy' ) ; close ;
%% make colorbar
color_map = jet ( 100 ) ;
scale_max = 3 ; scale_min = -3 ;
for i=1:length(color_map)
    plot ( scale_min + i / 100 * (scale_max-scale_min)  , 0.5 , 's' , 'Color' , 'None' , 'MarkerSize' , 40 , 'MarkerFaceColor' , color_map(i,:) ) ; hold on ;
end
set ( gca , 'YTickLabel' , [] ) ;
set ( gcf , 'Position' , [157 403 853 53] , 'Color' , 'None' ) ;
export_fig ( '../plots/m3_3_colormap' ) ; close ;


