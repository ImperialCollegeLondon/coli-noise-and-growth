

%% ref params
mu_ref = 0.0333333333333333 ;
km_ref = 0.279459339463851 ;
rm_ref = 0.138629436111989 ;
kp_ref = 0.938742978928143 ;
Vb_ref = 0.885197525716027 ;

%% params of the explo
fold_change = 2.0 ;

%% write the param vectors
mu = [ mu_ref ; fold_change * mu_ref ; fold_change * mu_ref ; fold_change * mu_ref ; fold_change * mu_ref ; fold_change * mu_ref ; fold_change * mu_ref ] ;
km = [ km_ref ; km_ref ; fold_change * km_ref ; km_ref ; km_ref ; fold_change * km_ref ; km_ref ] ;
rm = [ rm_ref ; rm_ref ; rm_ref ; rm_ref / fold_change ; rm_ref ; rm_ref ; rm_ref ] ;
kp = [ kp_ref ; kp_ref ; kp_ref ; kp_ref ; fold_change * kp_ref ; kp_ref ; fold_change * kp_ref ] ;
V_birth = [ Vb_ref ; Vb_ref ; Vb_ref ; Vb_ref ; Vb_ref ; fold_change * Vb_ref ; fold_change * Vb_ref ] ;

%% write the table
result_table = table ( mu , km , rm , kp , V_birth ) ;
writetable ( result_table , 'single_change_explo.csv' ) ;
