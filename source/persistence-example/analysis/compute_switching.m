function [num_switches,avg_residency_time_low,avg_residency_time_high,time_fraction_in_states] = compute_switching ( trace , low_thresh , high_thresh )

if ( low_thresh >= high_thresh )
    error ( 'wrong thresholds : lt >= ht ?' ) ;
end

% find first state
i_t = 1 ;
is_above = trace(i_t) > high_thresh ;
is_below = trace(i_t) < low_thresh ;
while ( (~is_above) && (~is_below) && (i_t<length(trace)) )
    i_t = i_t + 1 ;
    is_above = trace(i_t) > high_thresh ;
    is_below = trace(i_t) < low_thresh ;
end
was_high_state = is_above ;
was_low_state = is_below ;
if ( was_high_state + was_low_state == 0 )
    num_switches = 0 ;
    avg_residency_time_low = 0 ;
    avg_residency_time_high = 0 ;
    time_fraction_in_states = 0 ;
    return ;
end

residency_time = 1 ;
num_switches = 0 ;
num_high_state = 0 ;
num_low_state = 0 ;
sum_residency_time_high = 0 ;
sum_residency_time_low = 0 ;

% go over the trace
while ( i_t < length(trace) )
    
    % compute new state
    i_t = i_t + 1 ;
    is_above = trace(i_t) > high_thresh ;
    is_below = trace(i_t) < low_thresh ;
    
    % switch from high to low ?
    if ( was_high_state )
        if (is_below)
            was_high_state = false ;
            was_low_state = true ;
            num_switches = num_switches + 1 ;
            num_high_state = num_high_state + 1 ;
            last_residency_time_high = residency_time ;
            sum_residency_time_high = sum_residency_time_high + last_residency_time_high ;
            residency_time = 1 ;
        else
            residency_time = residency_time + 1 ;
        end
    end
    
    % switch from low to high ?
    if ( was_low_state )
        if (is_above)
            was_low_state = false ;
            was_high_state = true ;
            num_switches = num_switches + 1 ;
            num_low_state = num_low_state + 1 ;
            last_residency_time_low = residency_time ;
            sum_residency_time_low = sum_residency_time_low + last_residency_time_low ;
            residency_time = 1 ;
        else
            residency_time = residency_time + 1 ;
        end
    end
    
end

% compute avg res time
avg_residency_time_high = sum_residency_time_high / num_high_state ;
avg_residency_time_low = sum_residency_time_low / num_low_state ;
time_fraction_in_states = (sum_residency_time_high + sum_residency_time_low ) / length(trace) ;

end

