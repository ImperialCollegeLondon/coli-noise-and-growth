

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
size_no_size = { 'no_size_change' , 'size_change' } ;
for i_s = 1:length(size_no_size)
    for i_d = 1:length(demands)
        plot ( data.(demands{i_d}).km.(size_no_size{i_s}).mu .* 60 , data.(demands{i_d}).km.(size_no_size{i_s}).prot_birth_conc_CV , ...
            '-bo' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
        plot ( data.(demands{i_d}).kp.(size_no_size{i_s}).mu .* 60 , data.(demands{i_d}).kp.(size_no_size{i_s}).prot_birth_conc_CV , ...
            '-ro' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
        xlabel ( 'Division rate (dblgs/hr)' ) ; ylabel ( 'CV birth protein concentration' ) ;
        legend ( { 'Transcriptional response' , 'Translational response' } , 'Location' , 'NorthWest' ) ;
        set ( gcf , 'Color' , 'w' ) ; export_fig ( [ '../plots/CV-birth-prot-conc_' demands{i_d} '_' size_no_size{i_s} '.csv' ] ) ; close ;
    end
end

%% plot mRNA_avg vs div rate
size_no_size = { 'no_size_change' , 'size_change' } ;
for i_s = 1:length(size_no_size)
    for i_d = 1:length(demands)
        plot ( data.(demands{i_d}).km.(size_no_size{i_s}).mu .* 60 , data.(demands{i_d}).km.(size_no_size{i_s}).mRNA_birth_avg , ...
            '-bo' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
        plot ( data.(demands{i_d}).kp.(size_no_size{i_s}).mu .* 60 , data.(demands{i_d}).kp.(size_no_size{i_s}).mRNA_birth_avg , ...
            '-ro' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
        xlabel ( 'Division rate (dblgs/hr)' ) ; ylabel ( 'Avg birth mRNA number' ) ;
        legend ( { 'Transcriptional response' , 'Translational response' } , 'Location' , 'NorthWest' ) ;
        set ( gcf , 'Color' , 'w' ) ; export_fig ( [ '../plots/mRNA_avg_' demands{i_d} '_' size_no_size{i_s} '.png' ] ) ; close ;
    end
end