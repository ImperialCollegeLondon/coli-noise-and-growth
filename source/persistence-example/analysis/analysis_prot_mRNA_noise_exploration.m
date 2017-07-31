
% vectors of sigma and rp
sigmas = [ 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 ] ;
rps = [ 1 1.2 1.4 1.6 1.8 2.0 2.2 ] ;

% subplots per sigma x rp couple
for i_sig = 1:length(sigmas)
    for i_rp = 1:length(rps)
        
        % get sig and rp
        sig = sigmas(i_sig) ;
        rp = rps(i_rp) ;
        
        % load param sets
        input_pars = readtable ( [ '../input-data/prot_mRNA_noise_exploration/sig-' num2str(sig) '_rp-' num2str(rp) '.csv' ] ) ;
        
        % prepare plots
        plot_folder = [ '../plots/prot_mRNA_noise_exploration/sig-' num2str(sig) '_rp-' num2str(rp) ] ;
        mkdir ( plot_folder ) ;
        
        % make a subplot of V_birth vs EM_birth
        V_birth_vec = unique ( input_pars.V_birth ) ;
        EM_birth_vec = unique ( input_pars.EM_birth ) ;
        for i_sim = 1:size(input_pars,1)
            lineage_trajs = readtable ( [ '../output-data/prot_mRNA_noise_exploration/sig-' num2str(sig) '_rp-' num2str(rp) '/parset-' num2str(i_sim) '_lineage-trajs.csv' ] ) ;
            time_vec = linspace ( 0 , size(lineage_trajs,1) * input_pars.update_period(1) * input_pars.num_updates_per_output(1) , size(lineage_trajs,1) ) ;
            V_birth = input_pars{i_sim,'V_birth'} ; EM_birth = input_pars{i_sim,'EM_birth'} ;
            i_V = find ( V_birth_vec == V_birth ) ; i_EM = find ( EM_birth_vec == EM_birth ) ;
            subplot ( length(EM_birth_vec) , length(V_birth_vec) , sub2ind([length(V_birth_vec),length(EM_birth_vec)],i_V,i_EM) ) ;
            plot ( time_vec , lineage_trajs.alpha_fast , 'b' ) ;
            xlim ( minmax(time_vec) ) ; ylim ( [0 1] ) ;
            xlabel ( 'Time' , 'FontSize' , 10 ) ; ylabel ( '\alpha' , 'FontSize' , 10 ) ;
            set ( gca , 'FontSize' , 15 ) ;
            title ( [ 'V_{birth} = ' num2str(V_birth) ' , EM_{birth} = ' num2str(EM_birth) ] , 'FontSize' , 15 ) ;
        end
        set ( gcf , 'Position' , [0 0 2000 2000] ) ;
        export_fig ( [ plot_folder '/all-alpha-trajs' ] ) ; close ;
        
        % make a subplot of alpha distributions V_birth vs EM_birth
        alpha_bins = 0:0.15:1 ;
        V_birth_vec = unique ( input_pars.V_birth ) ;
        EM_birth_vec = unique ( input_pars.EM_birth ) ;
        for i_sim = 1:size(input_pars,1)
            lineage_trajs = readtable ( [ '../output-data/prot_mRNA_noise_exploration/sig-' num2str(sig) '_rp-' num2str(rp) '/parset-' num2str(i_sim) '_lineage-trajs.csv' ] ) ;
            V_birth = input_pars{i_sim,'V_birth'} ; EM_birth = input_pars{i_sim,'EM_birth'} ;
            i_V = find ( V_birth_vec == V_birth ) ; i_EM = find ( EM_birth_vec == EM_birth ) ;
            subplot ( length(EM_birth_vec) , length(V_birth_vec) , sub2ind([length(V_birth_vec),length(EM_birth_vec)],i_V,i_EM) ) ;
            freqs_fast = histc ( lineage_trajs.alpha_fast , alpha_bins ) ./ length(lineage_trajs.alpha_fast) ;
            freqs_slow = histc ( lineage_trajs.alpha_slow , alpha_bins ) ./ length(lineage_trajs.alpha_slow) ;
            bar ( alpha_bins , [freqs_fast,freqs_slow] ) ;
            xlim ( [-0.1 1.1] ) ; ylim ( [0 1] ) ;
            xlabel ( '\alpha' , 'FontSize' , 10 ) ;
            set ( gca , 'FontSize' , 15 ) ;
            title ( [ 'V_{birth} = ' num2str(V_birth) ' , EM_{birth} = ' num2str(EM_birth) ] , 'FontSize' , 15 ) ;
        end
        set ( gcf , 'Position' , [0 0 2000 2000] ) ;
        export_fig ( [ plot_folder '/all-alpha-distribs' ] ) ; close ;
        
    end
end

% subplot per EM_birth * V_birth couple
for i_V = 1:length(V_birth_vec)
    for i_EM = 1:length(EM_birth_vec)
        
        % prepare plots
        plot_folder = [ '../plots/prot_mRNA_noise_exploration/EM_birth-' num2str(EM_birth_vec(i_EM)) '_V_birth-' num2str(V_birth_vec(i_V)) ] ;
        mkdir ( plot_folder ) ;
        
        % iterate on couple
        for i_sig = 1:length(sigmas)
            for i_rp = 1:length(rps)
                
                subplot ( length(sigmas) , length(rps) , sub2ind([length(rps),length(sigmas)],i_rp,i_sig) ) ;
                
                % get sig and rp and pars for this couple
                sig = sigmas(i_sig) ;
                rp = rps(i_rp) ;
                input_pars = readtable ( [ '../input-data/prot_mRNA_noise_exploration/sig-' num2str(sig) '_rp-' num2str(rp) '.csv' ] ) ;
                
                % get the good traj (corresponding to that i_V and i_EM)
                for i_sim = 1:size(input_pars,1)
                    if ( (input_pars{i_sim,'V_birth'} == V_birth_vec(i_V)) && (input_pars{i_sim,'EM_birth'} == EM_birth_vec(i_EM)) )
                        lineage_trajs = readtable ( [ '../output-data/prot_mRNA_noise_exploration/sig-' num2str(sig) '_rp-' num2str(rp) '/parset-' num2str(i_sim) '_lineage-trajs.csv' ] ) ;
                        time_vec = linspace ( 0 , size(lineage_trajs,1) * input_pars.update_period(1) * input_pars.num_updates_per_output(1) , size(lineage_trajs,1) ) ;
                    end
                end
                
                % plot it
                plot ( time_vec , lineage_trajs.alpha_fast , 'b' ) ;
                xlim ( minmax(time_vec) ) ; ylim ( [0 1] ) ;
                xlabel ( 'Time' , 'FontSize' , 10 ) ; ylabel ( '\alpha' , 'FontSize' , 10 ) ;
                set ( gca , 'FontSize' , 15 ) ;
                title ( [ '\sigma = 10^{-' num2str(sig) '} , rp = 10^{-' num2str(rp) '}' ] , 'FontSize' , 15 ) ;

            end
        end
        
        % save the plot
        set ( gcf , 'Position' , [0 0 2000 2000] ) ;
        export_fig ( [ plot_folder '/EMVbirth_all-alpha-trajs' ] ) ; close ;

    end
end
