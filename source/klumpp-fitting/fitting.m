
% load data
ksm_data = load ( 'klumpp_ksm.dat' ) ;
gene_dosage_data = load ( 'klumpp_gene_dosage.dat' ) ;
vol_data = load ( 'klumpp_volume.dat' ) ;
ksp_data = load ( 'klumpp_ksp.dat' ) ;
rm_data = load ( 'klumpp_mrna_deg.dat' ) ;

% fitting range
gr_vec = linspace ( 0.5 , 3.5 , 100 ) ;

% do fits
[ksm_vec,pars_fit_ksm] = simple_fit ( ksm_data(:,1) , ksm_data(:,2) , gr_vec , 'sig' ) ;
[gd_vec,pars_fit_gd] = simple_fit ( gene_dosage_data(:,1) , gene_dosage_data(:,2) , gr_vec , 'exp' ) ;
[vol_vec,pars_fit_vol] = simple_fit ( vol_data(:,1) , vol_data(:,2) , gr_vec , 'exp' ) ;
[ksp_vec,pars_fit_ksp] = simple_fit ( ksp_data(:,1) , ksp_data(:,2) , gr_vec , 'lin' ) ;
[rm_vec,pars_fit_rm] = simple_fit ( rm_data(:,1) , rm_data(:,2) , gr_vec , 'lin' ) ;

% do plots
subplot ( 3 , 2 , 1 ) ;
plot ( gr_vec , ksm_vec ) ; hold on ; plot ( ksm_data(:,1) , ksm_data(:,2) , 'ro' ) ;
xlim ( [ 0.5 , 3.2 ] ) ;
ylabel ( 'transcription rate' ) ;
subplot ( 3 , 2 , 2 ) ;
plot ( gr_vec , gd_vec ) ; hold on ; plot ( gene_dosage_data(:,1) , gene_dosage_data(:,2) , 'ro' ) ;
xlim ( [ 0.5 , 3.2 ] ) ;
ylabel ( 'gene dosage' ) ;
subplot ( 3 , 2 , 3 ) ;
plot ( gr_vec , vol_vec ) ; hold on ; plot ( vol_data(:,1) , vol_data(:,2) , 'ro' ) ;
xlim ( [ 0.5 , 3.2 ] ) ;
ylabel ( 'cell volume' ) ;
subplot ( 3 , 2 , 4 ) ;
plot ( gr_vec , rm_vec ) ; hold on ; plot ( rm_data(:,1) , rm_data(:,2) , 'ro' ) ;
xlim ( [ 0.5 , 3.2 ] ) ;
ylabel ( 'mRNA degradation' ) ;
subplot ( 3 , 2 , 5 ) ;
plot ( gr_vec , ksp_vec ) ; hold on ; plot ( ksp_data(:,1) , ksp_data(:,2) , 'ro' ) ;
xlim ( [ 0.5 , 3.2 ] ) ;
ylabel ( 'translation rate' ) ;
xlabel ( 'doublings per hour' ) ;

% save plots
set ( gcf , 'Color' , 'None' ) ;
set ( gcf , 'Position' , [ 0 0 1200 900 ] ) ;
export_fig ( gcf , 'epsc' , 'klumpp_fit_plots.pdf' ) ;
close ;

%% also plot the case where km = km_per_gene * gene_dosage

i_ref_vec = find ( gr_vec > 2 , 1 , 'first' ) ;
i_ref_data_km_kp_vol = 4 ;

% km * gd
subplot ( 1 , 4 , 1 ) ; % subplot ( 2 , 2 , 1 ) ;
plot ( gr_vec , ksm_vec .* gd_vec ./ (ksm_vec(i_ref_vec) * gd_vec(i_ref_vec)) ) ; hold on ; 
plot ( ksm_data(:,1) , ksm_data(:,2) .* gene_dosage_data(:,2) ./ (ksm_data(i_ref_data_km_kp_vol,2) .* gene_dosage_data(i_ref_data_km_kp_vol,2)) , 'ro' ) ; hold on ;
plot ( [0 4] , [1 1] , 'k' , 'LineWidth' , 1.0 ) ;
xlim ( [ 0. , 4 ] ) ; ylim ( [ 0 2.5 ] ) ;
ylabel ( 'Rel. transcription rate' ) ;
xlabel ( 'Division rate (dblg/hr)' ) ;

% vol
subplot (  1, 4 , 2 ) ; % subplot ( 2 , 2 , 2 ) ;
plot ( gr_vec , vol_vec ./ vol_vec(i_ref_vec) ) ; hold on ;
plot ( vol_data(:,1) , vol_data(:,2) ./ vol_data(i_ref_data_km_kp_vol,2) , 'ro' ) ; hold on ;
plot ( [0 4] , [1 1] , 'k' , 'LineWidth' , 1.0 ) ;
xlim ( [ 0. , 4 ] ) ; ylim ( [ 0 2.5 ] ) ;
ylabel ( 'Rel. cell volume' ) ;
xlabel ( 'Division rate (dblg/hr)' ) ;

% kp
subplot (  1, 4 , 3 ) ; % subplot ( 2 , 2 , 3 ) ;
plot ( gr_vec , ksp_vec ./ ksp_vec(i_ref_vec) ) ; hold on ;
plot ( ksp_data(:,1) , ksp_data(:,2) ./ ksp_data(i_ref_data_km_kp_vol,2) , 'ro' ) ; hold on ;
plot ( [0 4] , [1 1] , 'k' , 'LineWidth' , 1.0 ) ;
xlim ( [ 0. , 4 ] ) ; ylim ( [ 0 1.2 ] ) ;
ylabel ( 'Rel. translation rate' ) ;
xlabel ( 'Division rate (dblg/hr)' ) ;

subplot (  1, 4 , 4 ) ; % subplot ( 2 , 2 , 4 ) ;
plot ( gr_vec , rm_vec ./ rm_vec(i_ref_vec) ) ; hold on ; plot ( rm_data(:,1) , rm_data(:,2) ./ rm_vec(i_ref_vec) , 'ro' ) ; hold on ;
plot ( [0 4] , [1 1] , 'k' , 'LineWidth' , 1.0 ) ;
xlim ( [ 0.4 , 4 ] ) ; ylim ( [ 0 1.4] ) ;
ylabel ( 'Rel. mRNA degrad. rate' ) ;
xlabel ( 'Division rate (dblg/hr)' ) ;
% save plots
set ( gcf , 'Color' , 'None' ) ;
set ( gcf , 'Position' , [ 0 0 2000 400 ] ) ;
export_fig ( gcf , 'epsc' , 'klumpp_fit_plots_no_gene_dosage.pdf' ) ;
close ;

