
name = 'DYR_ECOLI_e3_n2_m40'
protein = 'DYR_ECOLI'
number = '3'
exper = '1RX2.pdb'

covariance_analysis2( horzcat(name, '.fas'), protein, horzcat(name, '_MI_DIs.txt' ), 0.5,'N','Y');

get_constraints_PC(name,'M','FROB')

test_constraints_PC(name, protein, number, exper, 'M')
