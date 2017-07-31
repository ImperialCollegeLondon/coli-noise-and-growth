function [ waitings_ON , waitings_OFF ] = compute_waiting_times ( traj , threshold )

binary_traj = traj > threshold ;
conn_comp = bwconncomp ( binary_traj ) ;
num_switches = conn_comp.NumObjects - 1 ;
waitings_OFF = zeros ( num_switches , 1 ) ;
waitings_ON = zeros ( num_switches , 1 ) ;
for i=1:num_switches
    waitings_OFF(i) = min(conn_comp.PixelIdxList{i+1}) - max(conn_comp.PixelIdxList{i}) - 1 ;
    waitings_ON(i) = max(conn_comp.PixelIdxList{i}) - min(conn_comp.PixelIdxList{i}) + 1 ;
end


end

