function result_table = do_set_of_sims ( par_file )
%DO_SET_OF_SIMS Simulate all param set in a csv file
%   par_file in input-data

par_data = readtable ( par_file ) ;
sim_pars = table2struct ( readtable ( '../input-data/simulation_params.csv' ) ) ;

for i_sim = 1:size(par_data,1)
    
    % build the struct to pass to the simulation function
    par_struct = table2struct ( par_data(i_sim,:) ) ;
    names = [ fieldnames(par_struct) ; fieldnames(sim_pars) ] ;
    pars = cell2struct ( [struct2cell(par_struct); struct2cell(sim_pars)] , names , 1 ) ;
    
    % do the simulation
    pars.partitioning_type = 'normal' ;
    sim_data = do_single_sim(pars) ;
    mRNA.birth.avg(i_sim) = sim_data.mRNA_avg(1) ;
    mRNA.birth.CV(i_sim) = sim_data.mRNA_CV(1) ;
    prot.birth.avg(i_sim) = sim_data.prot_avg(1) ;
    prot.birth.CV(i_sim) = sim_data.prot_CV(1) ;

    % do the simulation with equal partitioning
    pars.partitioning_type = 'equal' ;
    sim_data = do_single_sim(pars) ;
    equal_part.mRNA.birth.avg(i_sim) = sim_data.mRNA_avg(1) ;
    equal_part.mRNA.birth.CV(i_sim) = sim_data.mRNA_CV(1) ;
    equal_part.prot.birth.avg(i_sim) = sim_data.prot_avg(1) ;
    equal_part.prot.birth.CV(i_sim) = sim_data.prot_CV(1) ;
end

% give result as a table with pars,sim_pars,mRNA_
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
prot_birth_conc_avg = prot_birth_avg ./ V_birth ;
prot_birth_conc_CV = prot_birth_CV ;
equal_part_mRNA_birth_avg = equal_part.mRNA.birth.avg';
equal_part_mRNA_birth_CV = equal_part.mRNA.birth.CV';
equal_part_prot_birth_avg = equal_part.prot.birth.avg' ;
equal_part_prot_birth_CV = equal_part.prot.birth.CV' ;
equal_part_prot_birth_conc_avg = equal_part_prot_birth_avg ./ V_birth ;
equal_part_prot_birth_conc_CV = equal_part_prot_birth_CV ;
if isfield (par_data,'prot_conc_demand')
    result_table = table ( mu , km , rm , kp , V_birth , prot_conc_demand , ...
        num_lineages , num_generations , ...
        mRNA_birth_avg , mRNA_birth_CV , prot_birth_avg , prot_birth_CV , ...
        prot_birth_conc_avg , prot_birth_conc_CV , ...
        equal_part_prot_birth_conc_avg , equal_part_prot_birth_conc_CV , equal_part_mRNA_birth_avg , equal_part_mRNA_birth_CV ) ;
else
    result_table = table ( mu , km , rm , kp , V_birth , ...
        num_lineages , num_generations , ...
        mRNA_birth_avg , mRNA_birth_CV , prot_birth_avg , prot_birth_CV , ...
        prot_birth_conc_avg , prot_birth_conc_CV , ...
        equal_part_prot_birth_conc_avg , equal_part_prot_birth_conc_CV , equal_part_mRNA_birth_avg , equal_part_mRNA_birth_CV ) ;
end
