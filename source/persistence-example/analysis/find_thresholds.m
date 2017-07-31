function [thresholds,low_avg_T,high_avg_T] = find_thresholds ( trace , max_fraction_out_states )

low_ts = linspace ( 0 , 1 , 50 ) ;
high_ts = low_ts ;
low_avg_T_matrix = zeros ( length(low_ts) , length(high_ts) ) ;
high_avg_T_matrix = zeros ( length(low_ts) , length(high_ts) ) ;
cost_matrix = zeros ( length(low_ts) , length(high_ts) ) ;


for i_low = 1:length(low_ts)
    for i_high = 1:length(high_ts)
        if ( low_ts(i_low) < high_ts(i_high) )
            [~,low_avg_T,high_avg_T,fraction_in_states] = compute_switching(trace,low_ts(i_low),high_ts(i_high)) ;
            if ( fraction_in_states > 1 - max_fraction_out_states )
                cost_matrix(i_low,i_high) = min ( [low_avg_T,high_avg_T] ) ;
            end
            low_avg_T_matrix(i_low,i_high) = low_avg_T ;
            high_avg_T_matrix(i_low,i_high) = high_avg_T ;
        end
    end
end

[~,I_max] = max ( cost_matrix(:) ) ;
[i_low,i_high] = ind2sub(size(cost_matrix),I_max) ;
thresholds = [low_ts(i_low) high_ts(i_high)] ;
low_avg_T = low_avg_T_matrix(i_low,i_high) ;
high_avg_T = high_avg_T_matrix(i_low,i_high) ;

plot ( trace , 'b' ) ; hold on ;
plot ( [1 length(trace)] , [thresholds(1) thresholds(1)] , 'r' ) ; hold on ;
plot ( [1 length(trace)] , [thresholds(2) thresholds(2)] , 'g' ) ; hold on ;

end






% function [low_avg_T,high_avg_T] = compute_avg_residency_times ( trace , thresholds )
%
%     above = (trace > thresholds(2)) ;
%     spanLocs = logical (above) ;   %identify contiguous ones
%     spanLength = regionprops (spanLocs, 'area') ;  %length of each span
%     spanLength = [ spanLength.Area ] ;
%     if (isempty(spanLength))
%         high_avg_T = 0 ;
%     else
%         high_avg_T = mean ( spanLength ) ;
%     end
%     below = (trace < thresholds(1)) ;
%     spanLocs = logical (below) ;   %identify contiguous ones
%     spanLength = regionprops (spanLocs, 'area') ;  %length of each span
%     spanLength = [ spanLength.Area ] ;
%     if (isempty(spanLength))
%         low_avg_T = 0 ;
%     else
%         low_avg_T = mean ( spanLength ) ;
%     end
%
% end


