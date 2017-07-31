function modality_result = analyze_modality ( traj , num_bins , dt  )

% bin data
edges = linspace ( 0 , 1.01 * max(traj) , num_bins+1 ) ;
disc_data = discretize ( traj , edges ) ;
disc_data_full = disc_data ;
data_size = length(disc_data) ;

% for return the edges vector
modality_result.edges = edges ;

% by default, rates are zeros
modality_result.rate_ON_OFF = 0 ;
modality_result.rate_OFF_ON = 0 ;

% find main mode
first_mode = mode ( disc_data ) ;

% delete it in the data to find other modes
disc_data(disc_data==first_mode)=[] ;
second_mode = mode ( disc_data ) ;
modes_ordered = minmax ( [first_mode second_mode] ) ;
mode_low = modes_ordered(1) ; mode_high = modes_ordered(2) ;

% compute how much points are in-between
in_between_modes = sum ( disc_data>mode_low & disc_data<mode_high ) ;

% iterate this until non-zero or no more data
while ( ~isempty(disc_data) && in_between_modes == 0 && ~isnan(second_mode) )

    disc_data(disc_data==second_mode)=[] ;
    second_mode = mode (disc_data) ;
    modes_ordered = minmax ( [first_mode second_mode] ) ;
    mode_low = modes_ordered(1) ; mode_high = modes_ordered(2) ;
    in_between_modes = sum ( disc_data>mode_low & disc_data<mode_high ) ;

end

% bimodality if two modes separated by some points...
modality_result.is_bimodal = false ;
if in_between_modes <= 0
    return ;
end

% and  if more that 5% of expected points for uniform
if ( min ( [sum(disc_data_full==first_mode) sum(disc_data_full==second_mode)] ) > 0.05 * data_size / num_bins )
    modality_result.is_bimodal = true ;
else
    return ;
end

% compute the threshold between the two modes (i.e. the bin in-between
% with the lowest count)
counts = zeros(mode_high-mode_low+1,1) ;
for i=1:length(counts)
    counts(i) = sum(disc_data_full==mode_low+i-1) ;
end
[~,modality_result.threshold_bin] = min (counts) ;
modality_result.threshold_bin = modality_result.threshold_bin + mode_low - 1 ;
centers = ( edges(1:end-1) + edges(2:end) ) ./ 2 ;
modality_result.threshold = centers(modality_result.threshold_bin) ;

% compute P_on (from distrib)
modality_result.P_on_from_distrib = sum ( traj > modality_result.threshold ) / data_size ;

% compute the rates OFF -> ON and ON -> OFF
modality_result.is_bistable = false ;
binary_traj = traj > modality_result.threshold ;
conn_comp = bwconncomp ( binary_traj ) ;
num_OFF_ON = conn_comp.NumObjects - 1 ;
num_ON_OFF = num_OFF_ON ; % equal by construction of conn_comp (draw it to be convinced)
waiting_times_OFF_ON = zeros ( num_OFF_ON , 1 ) ;
waiting_times_ON_OFF = zeros ( num_ON_OFF , 1 ) ;
for i=1:num_OFF_ON
    waiting_times_OFF_ON(i) = min(conn_comp.PixelIdxList{i+1}) - max(conn_comp.PixelIdxList{i}) - 1 ;
    waiting_times_ON_OFF(i) = max(conn_comp.PixelIdxList{i}) - min(conn_comp.PixelIdxList{i}) + 1 ;
end
modality_result.mean_OFF_time = mean ( waiting_times_OFF_ON ) * dt ;
modality_result.mean_ON_time = mean ( waiting_times_ON_OFF ) * dt ;
modality_result.rate_OFF_ON = 1 / modality_result.mean_OFF_time ;
modality_result.rate_ON_OFF = 1 / modality_result.mean_ON_time ;
modality_result.n_switches = num_OFF_ON ;
if num_OFF_ON < 2
    modality_result.is_bistable = false ;
else
    modality_result.is_bistable = true ;
    if num_OFF_ON < 5
        modality_result.low_switches_number = true ;
    else
        modality_result.low_switches_number = false ;
    end
end
modality_result.P_on_from_rates = modality_result.rate_OFF_ON / (modality_result.rate_OFF_ON + modality_result.rate_ON_OFF) ;


end
