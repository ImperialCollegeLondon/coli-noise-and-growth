%% load
load ( '../output-data/m3_kms_media_4.mat' ) ;
clear i j pars ;


%% plot ON occupancy heatmap
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
            % if not a number (could not be detected) or CV > 20% -> grey
            if ( isnan(on_occupancy) || CV_on_occupancy > 0.2 )
                color = 0.5 .* [1 1 1] ;
            else
                % color using colormap
                I_color = ceil ( on_occupancy * 100 ) ;
                color = color_map ( I_color , : ) ;
            end
        end
        
        % plot
        plot ( mu_max_vec(i) , km_slope_vec(j) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
xlim ( [1.4 2.6] ) ; ylim ( [-2.2 2.2] ) ;
set ( gcf , 'Position' , [ 363   371   352   326 ] ) ;
xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ; set ( gca , 'FontSize' , 15 ) ;  set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/m3_kms_media_4_ON_occupancy' ) ; close ;


%% plot OFF rate heatmap
scale_max_off_rate = 0.07 ;
color_map = jet ( 100 ) ; marker_size = 25 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        
        % check bimodality
        is_bimodals = cellfun ( @(x)(x.is_bimodal) , mod_res(i,j,:) ) ;
        % if one replicate not detected as bimodal, color as white
        if min ( is_bimodals ) == 0
            color = [1 1 1] ;
        else
            % check rate OFF
            rate_OFFs = cellfun ( @(x)(x.rate_ON_OFF) , mod_res(i,j,:) ) ;
            rate_OFF = mean ( rate_OFFs ) ;
            CV_rate_OFF = compute_CV ( rate_OFFs ) ;
            % if not a number (could not be detected) or CV > 20% -> grey
            if ( isnan(rate_OFF) || CV_rate_OFF > 0.2 )
                color = 0.5 .* [1 1 1] ;
            else
                % color using colormap
                I_color = ceil ( rate_OFF * 100 / scale_max_off_rate ) ;
                % if we go over the limit, black
                if I_color > 100
                    color = [0 0 0] ;
                else
                    color = color_map ( I_color , : ) ;
                end
            end
        end
        
        % plot
        plot ( mu_max_vec(i) , km_slope_vec(j) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
xlim ( [1.4 2.6] ) ; ylim ( [-2.2 2.2] ) ;
set ( gcf , 'Position' , [ 363   371   352   326 ] ) ;
xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ; set ( gca , 'FontSize' , 15 ) ;  set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/m3_kms_media_4_OFF_rate' ) ; close ;