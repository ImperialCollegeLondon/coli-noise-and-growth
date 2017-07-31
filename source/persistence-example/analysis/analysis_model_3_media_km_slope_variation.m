

%% load data
clear ;
load ( '../output-data/model_3_km_slope_variation_data' ) ;
clear i j ;

%% parameters for analysis
num_bins = 20 ;
% just for plotting
mu_edges = linspace ( 0 , mu_max_vec(end) , 20 ) ;
mu_centers = ( mu_edges(1:end-1)+mu_edges(2:end) ) ./ 2 ;


%% show subplots of histograms
I = 1 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        subplot ( length(mu_max_vec) , length(km_slope_vec) , I ) ; I = I + 1 ;
        counts = histcounts ( mu_trajs{i,j} , mu_edges, 'Normalization' , 'probability' ) ;
                
        modality_result = analyze_modality ( mu_trajs{i,j} , num_bins , time_traj(2) - time_traj(1) ) ;
        if modality_result.is_bimodal
            bar ( mu_centers , counts , 'r' ) ; hold on ;
            plot ( [1 1] .* modality_result.threshold , [0 1] , 'k' ) ;
            text ( 0.2 , 0.7 , [ 'P_{ON} = ' num2str(modality_result.P_on_from_distrib) ' (' num2str(modality_result.P_on_from_rates) ') from rates' ] ) ;
            text ( 0.2 , 0.5 , [ '<T_{OFF}> = ' num2str(modality_result.mean_OFF_time) ' hrs' ] ) ;
            text ( 0.2 , 0.3 , [ '<T_{ON}> = ' num2str(modality_result.mean_ON_time) ' hrs' ] ) ;
        else
            bar ( mu_centers , counts ) ;
        end
        xlim ( [0 1.1 .* mu_edges(end)] ) ; ylim ( [0 1.0] ) ;
        
        text ( 0.2 , 0.9 , [ '\mu_{max} = ' num2str(mu_max_vec(i)) ' , km_{\alpha} = ' num2str(km_slope_vec(j)) ] ) ;
        set ( gca , 'FontSize' , 10 ) ;
    end
end


%% write plot
set ( gcf , 'Color' , 'None' ) ;
export_fig ( gcf , '../plots/model_3_media_km_slope_variation_low_scale.png' ) ;
close ;
