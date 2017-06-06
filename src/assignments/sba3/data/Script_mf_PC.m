
name = 'CADH1_HUMAN_e3_n2_m40'
protein = 'CADH1_HUMAN'
number = '3'
exper = '2O72.pdb'

% uncomment the lines below for the second protein
% name = 'LACI_ECOLI_e40_n2_m40'
% protein = 'LACI_ECOLI'
% number = '40'
% exper = '1JWL.pdb'

covariance_analysis2(horzcat(name, '.fas'), protein, horzcat(name, '_MI_DIs.txt' ), 0.5,'N','Y');
get_constraints_PC(name,'M','FROB')
test_constraints_PC(name, protein, number, exper, 'M')