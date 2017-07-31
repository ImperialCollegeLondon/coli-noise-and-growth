
%% parameters
% all rates relative to alpha_max
pars.alpha_max = 1 ;
% fix P_toxin (the scale for prot concentration)
pars.P_toxin = 1 ;
% for now fix hill exponent
pars.hill = 2 ;
% only sigma and rp remain !
sigma_log10_vec = -1.2:0.015:0.2 ;
rp_log10_vec = -2.3:0.015:-0.9 ;


%% explore all combinations
num_roots_mat = zeros (  length(rp_log10_vec) , length(sigma_log10_vec) ) ;
dist_extreme_roots_mat = - ones (  length(rp_log10_vec) , length(sigma_log10_vec) ) ;
for i_sig=1:length(sigma_log10_vec)
    for i_rp=1:length(rp_log10_vec)
        pars.sigma = 10 ^ sigma_log10_vec(i_sig) ;
        pars.rp = 10 ^ rp_log10_vec(i_rp) ;
        roots = count_steady_states (pars) ;
        num_roots_mat(i_rp,i_sig) = roots.num_roots ;
        if ( roots.num_roots == 3 )
            dist_extreme_roots_mat(i_rp,i_sig) = roots.values(3)-roots.values(1) ;
        end
    end
end

%% define two example points, plot their nullclines intersections
give_ref_pars ;
draw_nullclines ( pars_1 , 2 , 1 ) ; 
draw_nullclines ( pars_2 , 2 , 2 ) ;
for i=1:4
    subplot ( 2 , 2 , i ) ; xlim( [0.2 60] ) ;
end
set ( gcf , 'Position' , [0 0 1600 1000] ) ;
export_fig ( gcf , '../plots/bistability_det_model_region_example_points' ) ;
close ;

%% plot the exploration result
size_points = 50 ;
% regions
subplot ( 1 , 2 , 1 ) ;
[SIG,RP] = meshgrid (sigma_log10_vec,rp_log10_vec) ;
scatter ( SIG(:) , RP(:) , size_points , num_roots_mat(:) , 's' , 'filled' ) ; hold on ;
plot ( log10(pars_1.sigma) , log10(pars_1.rp) , 'ro' , 'MarkerFaceColor' , 'r' ) ; hold on ;
plot ( log10(pars_2.sigma) , log10(pars_2.rp) , 'go' , 'MarkerFaceColor' , 'g' ) ;
xlim ( minmax(sigma_log10_vec) ) ; ylim ( minmax(rp_log10_vec) ) ; 
xlabel ( 'log10(\sigma)' ) ; ylabel ( 'log10(\gamma_p)' ) ;
colorbar ;
% distance btw stable states
subplot ( 1 , 2 , 2 ) ;
scatter ( SIG(:) , RP(:) , size_points , dist_extreme_roots_mat(:) , 's' , 'filled' ) ; hold on ;
plot ( log10(pars_1.sigma) , log10(pars_1.rp) , 'ro' , 'MarkerFaceColor' , 'r' ) ; hold on ;
plot ( log10(pars_2.sigma) , log10(pars_2.rp) , 'go' , 'MarkerFaceColor' , 'g' ) ;
xlim ( minmax(sigma_log10_vec) ) ; ylim ( minmax(rp_log10_vec) ) ; 
xlabel ( 'log10(\sigma)' ) ; ylabel ( 'log10(\gamma_p)' ) ;
colorbar ;
set ( gcf , 'Position' , [0 0 1500 500] ) ;
export_fig ( gcf , '../plots/bistability_exploration_det_model' ) ;
close ;
