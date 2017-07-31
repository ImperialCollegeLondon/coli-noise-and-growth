function result_table = do_set_of_sims ( par_file )
%DO_SET_OF_SIMS Simulate all param set in a csv file
%   par_file in input-data

par_data = readtable ( par_file ) ;
sim_pars = table2struct ( readtable ( '../input-data/simulation_params.csv' ) ) ;
lnm_trends = readtable ( '../input-data/linear_trends_lnm_params.csv' , 'ReadRowNames' , true ) ;

for i_sim = 1:size(par_data,1)

    % build the struct to pass to the simulation function
    par_struct = table2struct ( par_data(i_sim,:) ) ;
    names = [ fieldnames(par_struct) ; fieldnames(sim_pars) ] ;
    pars = cell2struct ( [struct2cell(par_struct); struct2cell(sim_pars)] , names , 1 ) ;
    
    % compute the lnm parameters
    pars.lnm_a = lnm_trends{'a_LNM','intercept'} + pars.mu * 60 * lnm_trends{'a_LNM','slope'} ;
    pars.lnm_b = pars.V_birth * ( 2 - pars.lnm_a ) ; % b such that average size at birth respected
    lnm_b_orig_unit = lnm_trends{'b_LNM','intercept'} + pars.mu * 60 * lnm_trends{'b_LNM','slope'} ;
    lnm_sigma_1_orig_unit = lnm_trends{'sig1_LNM','intercept'} + pars.mu * 60 * lnm_trends{'sig1_LNM','slope'} ;
    pars.lnm_sigma_1 = lnm_sigma_1_orig_unit * pars.lnm_b / lnm_b_orig_unit ; 
    pars.lnm_sigma_2 = lnm_trends{'sig2_LNM','intercept'} + pars.mu * 60 * lnm_trends{'sig2_LNM','slope'} ;
    pars.elong_rate_CV = lnm_trends{'alpha_CV','intercept'} + pars.mu * 60 * lnm_trends{'alpha_CV','slope'} ;

    % do the simulation
    sim_data = do_single_sim(pars) ;
    mRNA.birth.avg(i_sim) = sim_data.mRNA_avg(1) ;
    mRNA.birth.CV(i_sim) = sim_data.mRNA_CV(1) ;
    prot.birth.avg(i_sim) = sim_data.prot_avg(1) ;
    prot.birth.CV(i_sim) = sim_data.prot_CV(1) ;
    prot_conc.birth.avg(i_sim) = sim_data.prot_conc_avg(1) ;
    prot_conc.birth.CV(i_sim) = sim_data.prot_conc_CV(1) ;
    cell_size.birth.avg(i_sim) = sim_data.size_avg(1) ;
    cell_size.birth.CV(i_sim) = sim_data.size_CV(1) ;
end

% give result as a table with pars,sim_pars,mRNA_ ...
mu = par_data.mu ;
km = par_data.km ;
kp = par_data.kp ;
rm = par_data.rm ;
V_birth = par_data.V_birth ;
if isfield(par_data,'prot_conc_demand')
    prot_conc_demand = par_data.prot_conc_demand ;
end
num_lineages = repmat(sim_pars.num_lineages,size(par_data,1),1) ;
num_generations = repmat(sim_pars.num_generations,size(par_data,1),1) ;
mRNA_birth_avg = mRNA.birth.avg' ;
mRNA_birth_CV = mRNA.birth.CV' ;
prot_birth_avg = prot.birth.avg' ;
prot_birth_CV = prot.birth.CV' ;
prot_birth_conc_avg = prot_conc.birth.avg' ;
prot_birth_conc_CV = prot_conc.birth.CV' ;
% size_birth_avg = cell_size.birth.avg' ;
% size_birth_CV = cell_size.birth.CV' ;

if isfield (par_data,'prot_conc_demand')
    result_table = table ( mu , km , rm , kp , V_birth , prot_conc_demand , ...
        num_lineages , num_generations , ...
        mRNA_birth_avg , mRNA_birth_CV , prot_birth_avg , prot_birth_CV , ...
        prot_birth_conc_avg , prot_birth_conc_CV ) ;
else
    result_table = table ( mu , km , rm , kp , V_birth , ...
        num_lineages , num_generations , ...
        mRNA_birth_avg , mRNA_birth_CV , prot_birth_avg , prot_birth_CV , ...
        prot_birth_conc_avg , prot_birth_conc_CV ) ;
end

end
