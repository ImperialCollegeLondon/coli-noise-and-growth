
function sim_data = basic_no_feedback_study()

% sim_data = vary_nlm_a();
sim_data = vary_nlm_sigma_2();

end


function pars = give_ref_model()
%%% model type definition
pars.GE_model = 'dependent_rates';
pars.alpha_model = 'norm_distrib';
pars.partitioning_type = 'normal';

%%% model parameters that don't change
% Dependent GE parameters
pars.km_0 = 0.;
pars.km_per_size = 0.278 * 0.8;
pars.km_per_alpha = 0;
pars.rm = 0.139;
pars.kp = 0.939;
pars.rp = 0.0;
% lnm parameters
pars.lnm_a = 0.9;
pars.V_birth = 0.88;
pars.lnm_b = (2-pars.lnm_a) * pars.V_birth;
pars.lnm_sigma_1 = 0.18;
pars.lnm_sigma_2 = 0.035;
% growth media and toxin parameters
pars.mu_max = 2; % in doublings / hours
pars.alpha_CV = 0.08;
% simulation parameters
pars.random_seed = 0;
pars.num_lineages = 1;
sim_duration_ref = 100000 * 60;
pars.sim_duration = sim_duration_ref;
pars.update_period = 0.1;
pars.num_updates_per_output = 10;
end



function sim_data = vary_alpha_CV()

alpha_CV_vec = 0.5 ; %[ linspace( 0. , 0.2 , 5 ) , linspace(0.2 , 0.6 , 20) ] ;
pars = give_ref_model(); 

for i_a=1:length(alpha_CV_vec)
    
    pars.alpha_CV = alpha_CV_vec(i_a);
    
    sim_data = do_single_sim(pars);
    
    I_birth = find(diff(sim_data.traj_size)<0)+1;
    cell_size.birth.points = sim_data.traj_size(I_birth);
    prot.birth.points = sim_data.traj_prot(I_birth);
    prot_conc.birth.points = prot.birth.points ./ cell_size.birth.points;
    cc_time.points = diff(sim_data.traj_time(I_birth));
    
    cc_time.avg(i_a) = mean(cc_time.points);
    cc_time.compute_CV(i_a) = mean(cc_time.points);
    prot.birth.avg(i_a) = mean(prot.birth.points);
    prot.birth.CV(i_a) = compute_CV(prot.birth.points);
    prot_conc.birth.avg(i_a) = mean(prot_conc.birth.points);
    prot_conc.birth.CV(i_a) = compute_CV(prot_conc.birth.points);
    cell_size.birth.avg(i_a) = mean(cell_size.birth.points);
    cell_size.birth.CV(i_a) = compute_CV(cell_size.birth.points);
    
end

% plotting
% marker_size = 4;
% subplot(2,3,1); plot( alpha_CV_vec,prot_conc.birth.avg,'-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(prot_conc.birth.avg)]); ylabel('[prot] avg'); xlabel('alpha CV'); grid();
% subplot(2,3,[4 5]); plot( alpha_CV_vec , prot_conc.birth.CV , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); hold on;
% subplot(2,3,2); plot( alpha_CV_vec , prot.birth.avg , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(prot.birth.avg)]); ylabel('prot # avg'); xlabel('alpha CV'); grid();
% subplot(2,3,[4 5]); plot( alpha_CV_vec , prot.birth.CV , '-bs','MarkerSize',marker_size,'MarkerFaceColor','b'); ylim([0 1.5*max(prot.birth.CV)]); xlim(minmax(alpha_CV_vec)); ylabel('CV'); xlabel('alpha CV'); grid(); legend({'prot conc.','prot #'});
% subplot(2,3,3); plot( alpha_CV_vec , cell_size.birth.avg , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(cell_size.birth.avg)]); ylabel('birth size avg'); xlabel('alpha CV'); grid();
% subplot(2,3,6); plot( alpha_CV_vec , cell_size.birth.CV , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(cell_size.birth.CV)]); ylabel('birth size CV'); xlabel('alpha CV'); grid();
% for i=[1 2 3 6]; subplot(2,3,i); set(gca,'LineWidth',2,'FontSize',25); end
% subplot(2,3,[4 5]); set(gca,'LineWidth',2,'FontSize',25);
% set(gcf,'Color','w','Position',[200 163 1630 837]);
% export_fig(gcf,['lnm-study_alpha-CV-variation_sigma1-' num2str(pars.lnm_sigma_1) '_sigma2-' num2str(pars.lnm_sigma_2) '_lnm-a-' num2str(pars.lnm_a) '_lnm-b-' num2str(pars.lnm_b) '.pdf']); close;
end



function sim_data = vary_nlm_a()

a_vec = linspace( 0 , 1.5 , 10 );
pars = give_ref_model(); 

for i_a=1:length(a_vec)
    
    pars.lnm_a = a_vec(i_a);
    pars.lnm_b = (2-pars.lnm_a) * pars.V_birth; % to get the same average birth size
    
    sim_data = do_single_sim(pars);
    
    I_birth = find(diff(sim_data.traj_size)<0)+1;
    cell_size.birth.points = sim_data.traj_size(I_birth);
    prot.birth.points = sim_data.traj_prot(I_birth);
    prot_conc.birth.points = prot.birth.points ./ cell_size.birth.points;
    cc_time.points = diff(sim_data.traj_time(I_birth));
    
    cc_time.avg(i_a) = mean(cc_time.points);
    cc_time.compute_CV(i_a) = mean(cc_time.points);
    prot.birth.avg(i_a) = mean(prot.birth.points);
    prot.birth.CV(i_a) = compute_CV(prot.birth.points);
    prot_conc.birth.avg(i_a) = mean(prot_conc.birth.points);
    prot_conc.birth.CV(i_a) = compute_CV(prot_conc.birth.points);
    cell_size.birth.avg(i_a) = mean(cell_size.birth.points);
    cell_size.birth.CV(i_a) = compute_CV(cell_size.birth.points);
    
end

% plotting
marker_size = 4;
subplot(2,3,1); plot( a_vec,prot_conc.birth.avg,'-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(prot_conc.birth.avg)]); ylabel('[prot] avg'); xlabel('nlm a'); grid();
subplot(2,3,[4 5]); plot( a_vec , prot_conc.birth.CV , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); hold on;
subplot(2,3,2); plot( a_vec , prot.birth.avg , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(prot.birth.avg)]); ylabel('prot # avg'); xlabel('nlm a'); grid();
subplot(2,3,[4 5]); plot( a_vec , prot.birth.CV , '-bs','MarkerSize',marker_size,'MarkerFaceColor','b'); ylim([0 1.5*max(prot.birth.CV)]); xlim(minmax(a_vec)); ylabel('CV'); xlabel('nlm a'); grid(); legend({'prot conc.','prot #'});
subplot(2,3,3); plot( a_vec , cell_size.birth.avg , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(cell_size.birth.avg)]); ylabel('birth size avg'); xlabel('nlm a'); grid();
subplot(2,3,6); plot( a_vec , cell_size.birth.CV , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(cell_size.birth.CV)]); ylabel('birth size CV'); xlabel('nlm a'); grid();
for i=[1 2 3 6]; subplot(2,3,i); set(gca,'LineWidth',2,'FontSize',25); end
subplot(2,3,[4 5]); set(gca,'LineWidth',2,'FontSize',25);
set(gcf,'Color','w','Position',[200 163 1630 837]);
export_fig(gcf,['lnm-study_a-variation_km-prop-size_sigma1-' num2str(pars.lnm_sigma_1) '_sigma2-' num2str(pars.lnm_sigma_2) '_alpha-CV-' num2str(pars.alpha_CV) '_lnm-b-' num2str(pars.lnm_b) '.pdf']); close;
end



function sim_data = vary_nlm_sigma_1()

s1_vec = linspace( 0 , 0.8 , 10 );
pars = give_ref_model(); 

for i_s=1:length(s1_vec)
    
    pars.lnm_sigma_1 = s1_vec(i_s);
    
    
    sim_data = do_single_sim(pars);
    
    I_birth = find(diff(sim_data.traj_size)<0)+1;
    cell_size.birth.points = sim_data.traj_size(I_birth);
    prot.birth.points = sim_data.traj_prot(I_birth);
    prot_conc.birth.points = prot.birth.points ./ cell_size.birth.points;
    cc_time.points = diff(sim_data.traj_time(I_birth));
    
    cc_time.avg(i_s) = mean(cc_time.points);
    cc_time.compute_CV(i_s) = mean(cc_time.points);
    prot.birth.avg(i_s) = mean(prot.birth.points);
    prot.birth.CV(i_s) = compute_CV(prot.birth.points);
    prot_conc.birth.avg(i_s) = mean(prot_conc.birth.points);
    prot_conc.birth.CV(i_s) = compute_CV(prot_conc.birth.points);
    cell_size.birth.avg(i_s) = mean(cell_size.birth.points);
    cell_size.birth.CV(i_s) = compute_CV(cell_size.birth.points);
    
end

% plotting
marker_size = 4;
subplot(2,3,1); plot( s1_vec,prot_conc.birth.avg,'-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(prot_conc.birth.avg)]); ylabel('[prot] avg'); xlabel('sigma 1'); grid();
subplot(2,3,[4 5]); plot( s1_vec , prot_conc.birth.CV , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); hold on;
subplot(2,3,2); plot( s1_vec , prot.birth.avg , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(prot.birth.avg)]); ylabel('prot # avg'); xlabel('sigma 1'); grid();
subplot(2,3,[4 5]); plot( s1_vec , prot.birth.CV , '-bs','MarkerSize',marker_size,'MarkerFaceColor','b'); ylim([0 1.5*max(prot.birth.CV)]); xlim(minmax(s1_vec)); ylabel('CV'); xlabel('sigma 1'); grid(); legend({'prot conc.','prot #'});
subplot(2,3,3); plot( s1_vec , cell_size.birth.avg , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(cell_size.birth.avg)]); ylabel('birth size avg'); xlabel('sigma 1'); grid();
subplot(2,3,6); plot( s1_vec , cell_size.birth.CV , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(cell_size.birth.CV)]); ylabel('birth size CV'); xlabel('sigma 1'); grid();
for i=[1 2 3 6]; subplot(2,3,i); set(gca,'LineWidth',2,'FontSize',25); end
subplot(2,3,[4 5]); set(gca,'LineWidth',2,'FontSize',25);
set(gcf,'Color','w','Position',[200 163 1630 837]);
export_fig(gcf,['lnm-study_sigma-1-variation_km-prop-size_lnm-a-' num2str(pars.lnm_a) '_sigma2-' num2str(pars.lnm_sigma_2) '_alpha-CV-' num2str(pars.alpha_CV) '_lnm-b-' num2str(pars.lnm_b) '.pdf']); close;
end


function sim_data = vary_nlm_sigma_2()

s2_vec = linspace( 0 , 0.15 , 10 );
pars = give_ref_model(); 

for i_s=1:length(s2_vec)
    
    pars.lnm_sigma_2 = s2_vec(i_s);
    
    
    sim_data = do_single_sim(pars);
    
    I_birth = find(diff(sim_data.traj_size)<0)+1;
    cell_size.birth.points = sim_data.traj_size(I_birth);
    prot.birth.points = sim_data.traj_prot(I_birth);
    prot_conc.birth.points = prot.birth.points ./ cell_size.birth.points;
    cc_time.points = diff(sim_data.traj_time(I_birth));
    
    cc_time.avg(i_s) = mean(cc_time.points);
    cc_time.compute_CV(i_s) = mean(cc_time.points);
    prot.birth.avg(i_s) = mean(prot.birth.points);
    prot.birth.CV(i_s) = compute_CV(prot.birth.points);
    prot_conc.birth.avg(i_s) = mean(prot_conc.birth.points);
    prot_conc.birth.CV(i_s) = compute_CV(prot_conc.birth.points);
    cell_size.birth.avg(i_s) = mean(cell_size.birth.points);
    cell_size.birth.CV(i_s) = compute_CV(cell_size.birth.points);
    
end

% plotting
marker_size = 4;
subplot(2,3,1); plot( s2_vec,prot_conc.birth.avg,'-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(prot_conc.birth.avg)]); ylabel('[prot] avg'); xlabel('sigma 2'); grid();
subplot(2,3,[4 5]); plot( s2_vec , prot_conc.birth.CV , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); hold on;
subplot(2,3,2); plot( s2_vec , prot.birth.avg , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(prot.birth.avg)]); ylabel('prot # avg'); xlabel('sigma 2'); grid();
subplot(2,3,[4 5]); plot( s2_vec , prot.birth.CV , '-bs','MarkerSize',marker_size,'MarkerFaceColor','b'); ylim([0 1.5*max(prot.birth.CV)]); xlim(minmax(s2_vec)); ylabel('CV'); xlabel('sigma 2'); grid(); legend({'prot conc.','prot #'});
subplot(2,3,3); plot( s2_vec , cell_size.birth.avg , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(cell_size.birth.avg)]); ylabel('birth size avg'); xlabel('sigma 2'); grid();
subplot(2,3,6); plot( s2_vec , cell_size.birth.CV , '-ko','MarkerSize',marker_size,'MarkerFaceColor','k'); ylim([0 1.2 * max(cell_size.birth.CV)]); ylabel('birth size CV'); xlabel('sigma 2'); grid();
for i=[1 2 3 6]; subplot(2,3,i); set(gca,'LineWidth',2,'FontSize',25); end
subplot(2,3,[4 5]); set(gca,'LineWidth',2,'FontSize',25);
set(gcf,'Color','w','Position',[200 163 1630 837]);
export_fig(gcf,['lnm-study_sigma-2-variation_km-prop-size_lnm-a-' num2str(pars.lnm_a) '_sigma1-' num2str(pars.lnm_sigma_1) '_alpha-CV-' num2str(pars.alpha_CV) '_lnm-b-' num2str(pars.lnm_b) '.pdf']); close;
end
