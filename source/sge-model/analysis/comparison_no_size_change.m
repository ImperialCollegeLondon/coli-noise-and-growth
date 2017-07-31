%% parse data
demands = { 'Q' , 'R' , 'P' } ;
for i_d = 1:length(demands)
    data.(demands{i_d}).km.no_size_change = readtable ( [ '../output-data/sim-result_demand-' demands{i_d} '_response-km_no_size_change.csv' ] ) ;
    data.(demands{i_d}).kp.no_size_change = readtable ( [ '../output-data/sim-result_demand-' demands{i_d} '_response-kp_no_size_change.csv' ] ) ;
    data.(demands{i_d}).km.size_change = readtable ( [ '../output-data/sim-result_demand-' demands{i_d} '_response-km.csv' ] ) ;
    data.(demands{i_d}).kp.size_change = readtable ( [ '../output-data/sim-result_demand-' demands{i_d} '_response-kp.csv' ] ) ;
end

%% some plotting parameters
marker_size = 20 ;
line_width = 4 ;

%% plot CV prot conc vs div rate
for i_d = 1:length(demands)
    plot ( data.(demands{i_d}).km.size_change.mu .* 60 , data.(demands{i_d}).km.size_change.prot_birth_conc_CV , ...
        'b' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
    plot ( data.(demands{i_d}).kp.size_change.mu .* 60 , data.(demands{i_d}).kp.size_change.prot_birth_conc_CV , ...
        'r' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
    plot ( data.(demands{i_d}).km.no_size_change.mu .* 60 , data.(demands{i_d}).km.no_size_change.prot_birth_conc_CV , ...
        '--b' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
    plot ( data.(demands{i_d}).kp.no_size_change.mu .* 60 , data.(demands{i_d}).kp.no_size_change.prot_birth_conc_CV , ...
        '--r' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
    xlabel ( 'Division rate (dblgs/hr)' ) ; ylabel ( 'CV birth protein concentration' ) ;
    legend ( { 'Transcriptional response' , 'Translational response' ,  'No size change, transcriptional response' , 'No size change, translational response' } , 'Location' , 'NorthWest' ) ;
    set ( gcf , 'Color' , 'w' ) ; export_fig ( [ '../plots/comparison_size-change_no-size-change/CV-birth-prot-conc_' demands{i_d} '.png' ] ) ; close ;
end