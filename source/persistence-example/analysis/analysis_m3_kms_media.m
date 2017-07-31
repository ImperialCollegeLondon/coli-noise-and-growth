
%% load
load ( '../output-data/m3_kms_media.mat' ) ;
clear i j pars ;


%% plot ON occupancy heatmap
% heatmap
color_map = jet ( 100 ) ; marker_size = 20 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        on_occupancy = rate_ON(i,j)/ ( rate_OFF(i,j) + rate_ON(i,j) ) ;
        if ~ isnan ( on_occupancy )
            I_color = ceil ( on_occupancy * 100 ) ;
            color = color_map ( I_color , : ) ;
        else
            color = 0.4 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
set ( gcf , 'Position' , [ 391   227   528   497 ] ) ;
xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ;  set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/m3_kms_media_ON_occupancy' ) ; close ;


%% plot ON rate heatmap
% heatmap
scale_max_on_rate = 0.1 ;
color_map = jet ( 100 ) ; marker_size = 20 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        I_color = ceil ( rate_ON(i,j) * 100 / scale_max_on_rate ) ;
        if I_color > 0 && I_color <= 100
            color = color_map ( I_color , : ) ;
        else
            color = [1 1 1] ;
        end
        if rate_ON(i,j) == 0
            color = 0.4 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
set ( gcf , 'Position' , [ 391   227   528   497 ] ) ;
xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ; set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/m3_kms_media_ON_rate' ) ; close ;


%% plot OFF rate heatmap
% heatmap
scale_max_off_rate = 0.1 ;
color_map = jet ( 100 ) ; marker_size = 20 ;
for i=1:length(mu_max_vec)
    for j=1:length(km_slope_vec)
        I_color = ceil ( rate_OFF(i,j) * 100 / scale_max_off_rate ) ;
        if I_color > 0 && I_color <= 100
            color = color_map ( I_color , : ) ;
        else
            color = [1 1 1] ;
        end
        if rate_OFF(i,j) == 0
            color = 0.4 .* [1 1 1] ;
        end
        plot ( mu_max_vec(i) , km_slope_vec(j) , 's' , 'Color' , 'None' , ...
            'MarkerFaceColor' , color , 'MarkerSize' , marker_size ) ; hold on ;
    end
end
set ( gcf , 'Position' , [ 391   227   528   497 ] ) ;
xlabel ( '\mu_{max} (dblgs/hr)' ) ; ylabel ( 'k_m slope' ) ;  set ( gcf , 'Color' , 'None' ) ;
export_fig ( '../plots/m3_kms_media_OFF_rate' ) ; close ;
