%% parse data
demands = { 'Q' , 'R' , 'P' } ;
for i_d = 1:length(demands)
    data_no_lnm.(demands{i_d}).km.size_change = readtable ( [ '../../sge-model/output-data/sim-result_demand-' demands{i_d} '_response-km.csv' ] ) ;
    data_no_lnm.(demands{i_d}).kp.size_change = readtable ( [ '../../sge-model/output-data/sim-result_demand-' demands{i_d} '_response-kp.csv' ] ) ;
    data.(demands{i_d}).km.size_change = readtable ( [ '../output-data/sim-result_demand-' demands{i_d} '_response-km.csv' ] ) ;
    data.(demands{i_d}).kp.size_change = readtable ( [ '../output-data/sim-result_demand-' demands{i_d} '_response-kp.csv' ] ) ;
end

%% verify same reference parameters for the normal / lnm data
i_mu_ref = round ( size(data_no_lnm.Q.km.size_change,1) / 2 ) ;
if ( data_no_lnm.Q.km.size_change{i_mu_ref,'mu'} ~= data.Q.km.size_change{i_mu_ref,'mu'} )
    error ( 'different mu_refs btw lnm / no lnm data !' ) ;
end
if ( data_no_lnm.Q.km.size_change{i_mu_ref,'km'} ~= data.Q.km.size_change{i_mu_ref,'km'} )
    error ( 'different km_refs btw lnm / no lnm data !' ) ;
end
if ( data_no_lnm.Q.km.size_change{i_mu_ref,'kp'} ~= data.Q.km.size_change{i_mu_ref,'kp'} )
    error ( 'different kp_refs btw lnm / no lnm data !' ) ;
end
if ( data_no_lnm.Q.km.size_change{i_mu_ref,'rm'} ~= data.Q.km.size_change{i_mu_ref,'rm'} )
    error ( 'different rm_refs btw lnm / no lnm data !' ) ;
end
if ( data_no_lnm.Q.km.size_change{i_mu_ref,'V_birth'} ~= data.Q.km.size_change{i_mu_ref,'V_birth'} )
    error ( 'different V_birth_refs btw lnm / no lnm data !' ) ;
end

%% some plotting parameters
marker_size = 20 ;
line_width = 4 ;

%% plot CV prot conc vs div rate
for i_d = 1:length(demands)
    plot ( data.(demands{i_d}).km.size_change.mu .* 60 , data.(demands{i_d}).km.size_change.prot_birth_conc_CV , ...
        '-b' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
    plot ( data.(demands{i_d}).kp.size_change.mu .* 60 , data.(demands{i_d}).kp.size_change.prot_birth_conc_CV , ...
        '-r' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
    plot ( data_no_lnm.(demands{i_d}).km.size_change.mu .* 60 , data_no_lnm.(demands{i_d}).km.size_change.prot_birth_conc_CV , ...
        '--b' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
    plot ( data_no_lnm.(demands{i_d}).kp.size_change.mu .* 60 , data_no_lnm.(demands{i_d}).kp.size_change.prot_birth_conc_CV , ...
        '--r' , 'MarkerSize' , marker_size , 'LineWidth' , line_width ) ; hold on ;
    xlabel ( 'Division rate (dblgs/hr)' ) ; ylabel ( 'CV birth protein concentration' ) ;
    legend ( { 'LNM, transcriptional response' , 'LNM, translational response' , ...
        'no LNM, transcriptional response' , 'no LNM, translational response' } , 'Location' , 'NorthWest' ) ;
    set ( gcf , 'Color' , 'w' ) ; export_fig ( [ '../plots/comparison_lnm_no-lnm/CV-birth-prot-conc_' demands{i_d} '.png' ] ) ; close ;
end

