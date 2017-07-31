%% 
data = readtable ( '../output-data/sim-result_single_change_explo.csv' ) ;

%% 
plot ( data{1,'prot_birth_conc_avg'} , data{1,'prot_birth_conc_CV'} , 'ks' , 'MarkerSize' , 20 , 'MarkerFaceColor' , 'k' ) ; hold on ;
plot ( data{2,'prot_birth_conc_avg'} , data{2,'prot_birth_conc_CV'} , 'ko' , 'MarkerSize' , 20 , 'MarkerFaceColor' , 'k' ) ; hold on ;
plot ( data{3,'prot_birth_conc_avg'} , data{3,'prot_birth_conc_CV'} , 'bo' , 'MarkerSize' , 20 , 'MarkerFaceColor' , 'b' ) ; hold on ;
plot ( data{4,'prot_birth_conc_avg'} , data{4,'prot_birth_conc_CV'} , 'go' , 'MarkerSize' , 20 , 'MarkerFaceColor' , 'g' ) ; hold on ;
plot ( data{5,'prot_birth_conc_avg'} , data{5,'prot_birth_conc_CV'} , 'ro' , 'MarkerSize' , 20 , 'MarkerFaceColor' , 'r' ) ; hold on ;
plot ( data{6,'prot_birth_conc_avg'} , data{6,'prot_birth_conc_CV'} , 'mo' , 'MarkerSize' , 20 , 'MarkerFaceColor' , 'm' ) ; hold on ;
plot ( data{7,'prot_birth_conc_avg'} , data{7,'prot_birth_conc_CV'} , 'co' , 'MarkerSize' , 20 , 'MarkerFaceColor' , 'c' ) ; hold on ;

%%
ylim ( [0.15 0.5] ) ; xlim ( [0 70] ) ;
xlabel ( 'average protein conc. (at birth)' ) ; ylabel ( 'CV protein conc. (at birth)') ;
title ( 'effect of growth rate X2' ) ;
legend ( { 'origin' , 'no change' , 'km change (X2)' , 'rm change (/2)' , 'kp change (X2)' , 'V_{birth} and km change (X2)' , 'V_{birth} and kp change (X2)' } , ...
    'Location' , 'SouthWest' ) ;
export_fig ( gcf , '../plots/single_change_explo' ) ; 
close ;



