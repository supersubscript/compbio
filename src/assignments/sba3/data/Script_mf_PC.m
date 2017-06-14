% uncomment the lines below for the first protein
% name = 'CADH1_HUMAN_e3_n2_m40';
% protein = 'CADH1_HUMAN';
% number = '3';
% exper = '2O72.pdb';

% uncomment the lines below for the second protein
% name = 'LACI_ECOLI_e40_n2_m40';
% protein = 'LACI_ECOLI';
% number = '40';
% exper = '1JWL.pdb';

% uncomment the lines below for the third protein
% name = 'DYR_ECOLI_e3_n2_m40';
% protein = 'DYR_ECOLI';
% number = '3';
% exper = '1RX2.pdb';

% last param is theta
tt = [0:0.1:1];
parfor indxone = 1:length(tt)
    for indxtwo = 1:length(tt)
        theta = tt(indxone);
        pc_weight = tt(indxtwo);
        covariance_analysis2(horzcat(name, '.fas'), protein, horzcat(name, '_MI_DIs_', protein, '_theta_', num2str(theta), '_pc_weight_', num2str(pc_weight), '.txt' ), pc_weight,'N','Y', theta, pc_weight, protein);
        get_constraints_PC(name,'M','FROB')
        test_constraints_PC(name, protein, number, exper, 'M', num2str(theta), num2str(pc_weight))
    end
end