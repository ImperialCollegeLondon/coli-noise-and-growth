function waiting_times = compute_exit_times( bin_trace )

% if always true, no waiting times
if ( all(bin_trace) )
    waiting_times = [] ;
    return ;
end

% remove the first 1111's because we don't know the start
if ( bin_trace(1) )
    bin_trace(1:find(~bin_trace,1,'first')-1) = [] ;
end

% remove the end 1111's because we don't know the end
if ( bin_trace(end) )
    bin_trace(find(~bin_trace,1,'last')+1:end) = [] ;
end

d = diff([0 ; bin_trace ; 0]);
waiting_times = find(d<0)-find(d>0) ;

end

