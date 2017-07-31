
% load param sets
input_pars = readtable ( [ '../input-data/km_scaling_with_mu/set_1.csv' ] ) ;

% prepare plots
plot_folder = [ '../plots/km_scaling_with_mu/set_1/' ] ;
mkdir ( plot_folder ) ;

% make a subplot of alpha trajs
for i_sim = 1:size(input_pars,1)
    lineage_trajs = readtable ( [ '../output-data/km_scaling_with_mu/set_1/parset-' num2str(i_sim) '_lineage-trajs.csv' ] ) ;
    time_vec = linspace ( 0 , size(lineage_trajs,1) * input_pars.update_period(1) * input_pars.num_updates_per_output(1) , size(lineage_trajs,1) ) ;
    subplot ( 3 , 4 , i_sim ) ;
    plot ( time_vec , lineage_trajs.alpha_fast ./ input_pars.alpha_max(i_sim) , 'b' ) ;
    xlim ( minmax(time_vec) ) ; ylim ( [0 1] ) ;
    xlabel ( 'Time' , 'FontSize' , 10 ) ; ylabel ( '\alpha / \alpha_{max}' , 'FontSize' , 10 ) ;
    set ( gca , 'FontSize' , 15 ) ;
    title ( [ '\mu = ' num2str(input_pars.alpha_max(i_sim)) ] , 'FontSize' , 15 ) ;
end
set ( gcf , 'Position' , [0 0 2000 2000] ) ;
export_fig ( [ plot_folder '/all-alpha-trajs' ] ) ; close ;

% % make a subplot of alpha distributions V_birth vs EM_birth
% alpha_bins = 0:0.15:1 ;
% V_birth_vec = unique ( input_pars.V_birth ) ;
% EM_birth_vec = unique ( input_pars.EM_birth ) ;
% for i_sim = 1:size(input_pars,1)
%     lineage_trajs = readtable ( [ '../output-data/prot_mRNA_noise_exploration/sig-' num2str(sig) '_rp-' num2str(rp) '/parset-' num2str(i_sim) '_lineage-trajs.csv' ] ) ;
%     V_birth = input_pars{i_sim,'V_birth'} ; EM_birth = input_pars{i_sim,'EM_birth'} ;
%     i_V = find ( V_birth_vec == V_birth ) ; i_EM = find ( EM_birth_vec == EM_birth ) ;
%     subplot ( length(EM_birth_vec) , length(V_birth_vec) , sub2ind([length(V_birth_vec),length(EM_birth_vec)],i_V,i_EM) ) ;
%     freqs_fast = histc ( lineage_trajs.alpha_fast , alpha_bins ) ./ length(lineage_trajs.alpha_fast) ;
%     freqs_slow = histc ( lineage_trajs.alpha_slow , alpha_bins ) ./ length(lineage_trajs.alpha_slow) ;
%     bar ( alpha_bins , [freqs_fast,freqs_slow] ) ;
%     xlim ( [-0.1 1.1] ) ; ylim ( [0 1] ) ;
%     xlabel ( '\alpha' , 'FontSize' , 10 ) ;
%     set ( gca , 'FontSize' , 15 ) ;
%     title ( [ 'V_{birth} = ' num2str(V_birth) ' , EM_{birth} = ' num2str(EM_birth) ] , 'FontSize' , 15 ) ;
% end
% set ( gcf , 'Position' , [0 0 2000 2000] ) ;
% export_fig ( [ plot_folder '/all-alpha-distribs' ] ) ; close ;


