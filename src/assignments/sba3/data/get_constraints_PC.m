% this code reads in the results of the covariance_analysis.m programme and
% extracts the highest scoring residue pairs for which the two residues
% involved are not in close proximity in sequence. 

% this code does not consider SS. 	

function get_constraints_PC(path, type,SCORING_METHOD_FOR_RANKING_CONSTRAINTS)
% get_constraints_PSEUDO('DYR_ECOLI_e3_n2_m40')
set(0,'DefaultFigureColor','w')

% input file name
if type == 'G'
	AllPairConstraintScoreFile = horzcat(path,'_GREMLIN_PROCESSED');
else
    if type == 'P'
	AllPairConstraintScoreFile = horzcat(path,'_L1_PSEUDO_PROCESSED');
    else
        AllPairConstraintScoreFile =   horzcat(path,'_MI_DIs.txt');
    end
end
%IndexMappingFile = horzcat(path,'.indextable');
IndexMappingFile = horzcat(path,'.indextable');

% output file name
DIScoreFilename = horzcat(path,'_DIScores.csv');
if type == 'G'
	DIScoreFilename = horzcat(path,'_GREMLIN_DIScores.csv');
end
if type == 'P'
	DIScoreFilename = horzcat(path,'_L1_PSEUDO_DIScores.csv');
end
     
Predicted_Constraints_name = horzcat(path,'_Predicted_Constraints');


% parameter values
%SCORING_METHOD_FOR_RANKING_CONSTRAINTS = 'FROB';
PAIR_MIN_RESIDUE_OFFSET_DISTANCE = 5;
TOTAL_NUMBER_CONSTRAINTS = 1000;
Nconstraints_runs = 200;


%read in the index mapping file
fid = fopen(IndexMappingFile);
fgetl(fid); %discard header
index_map_file_contents = textscan(fid,'%*s %c %c %s %s %s %c %c %c'); %aa ss_code ss_conf msa_ind msa_cons% msa_sig in_const
fclose(fid);
uniprot_residue = index_map_file_contents{1};
ssCode = index_map_file_contents{2};
uniprot_msa_cons = zeros(1,size(index_map_file_contents{5},1));
for i = 1:size(index_map_file_contents{5})
	if ~strcmp(index_map_file_contents{5}{i},'-')
		uniprot_msa_cons(i) = str2num(index_map_file_contents{5}{i});
	end
end
SeqLength = size(uniprot_residue,1);

%read in the DI score text file.
AllPairConstraintScoreFile
fid = fopen(AllPairConstraintScoreFile);
FileRead_temp = fscanf(fid,'%d %*s %d %*s %g %g %g',[5,Inf]);
fclose(fid);
DIScores_allpairs = FileRead_temp';
temp = size(DIScores_allpairs)
if (strcmp(SCORING_METHOD_FOR_RANKING_CONSTRAINTS,'DI'))
	scoreColumnIndex = 4;
elseif(strcmp(SCORING_METHOD_FOR_RANKING_CONSTRAINTS,'FROB'))
	scoreColumnIndex = 5;
else
	error('Unknown scoring method selected');
end
SCORING_METHOD_FOR_RANKING_CONSTRAINTS
for i = 1:temp(1)
	DIMatrix(DIScores_allpairs(i,1),DIScores_allpairs(i,2)) = DIScores_allpairs(i,scoreColumnIndex);
end
size(DIMatrix)
%TODO: matrix is now non-packed (there are uniprot offsets which are not used) FIX or map indicies to a reduced matrix -- or limit range to [min, max] coordinates
DIStart = DIScores_allpairs(1,1); %may not start at 1
DIEnd = size(DIMatrix,2);
DIMatrix(DIEnd,:) = 0; %adds a row of zeros at end

% remove those pairs that contain two residues in close sequence proximity,
% e.g. (i,j) where i = 10 and j = 12. 
for i = DIStart:DIEnd
	for j = i:min(i+PAIR_MIN_RESIDUE_OFFSET_DISTANCE,DIEnd)
		DIMatrix(i,j) = 0;
	end
end

[HiScores,IDX] = sort(DIMatrix(:),'descend');
for i = 1:TOTAL_NUMBER_CONSTRAINTS
	DIScores(i,1) = ceil(IDX(i)/DIEnd);
	DIScores(i,2) = IDX(i) - floor(IDX(i)/DIEnd)*DIEnd;
	DIScores(i,3) = HiScores(i);
end
%csvwrite(DIScoreFilename,DIScores);

% %if (strcmp(BATCH_MODE,'yes') ~= 1)
% 	figurehandle1 = figure();
% 	subplot(1,2,1), hist(DIMatrix(:),50)
% 	title(horzcat('DI Score histogram for ', PROTEIN_NAME),'Interpreter','none')
% 	subplot(1,2,2), hist(HiScores(1:500),100), line([HiScores(100) HiScores(100)],[0 50])
% 	title(horzcat('Top 500 DI Scores for ', PROTEIN_NAME),'Interpreter','none')
% 	eval(['print -dpdf -f' num2str(figurehandle1) ' ' protein_name '_DIScoreHist']);
%end

DIScores(DIScores(:,2) == 0, : ) = []; %removes 0 values
DIScores(DIScores(:,1) >=  length(uniprot_msa_cons), :  ) = []; %removes 0 values
uniprot_msa_cons 
max(DIScores(:,2))

DIScores(:,5) = uniprot_msa_cons(DIScores(:,1));
DIScores(:,6) = uniprot_msa_cons(DIScores(:,2));

csvwrite(DIScoreFilename, DIScores);

% make predicted constraint lots
predicted_constraint = zeros(SeqLength,SeqLength);
%if (strcmp(BATCH_MODE,'yes') ~= 1)
	figurehandle4 = figure();
	rankPosition = 0;
	i = 1;
	while i <= TOTAL_NUMBER_CONSTRAINTS && rankPosition < Nconstraints_runs(length(Nconstraints_runs))
		%if ~ssClashFlag(i) && ~highConservationFlag(i) && ~cys_filter_flag(i)
			t1 = DIScores(i,1);
			t2 = DIScores(i,2);
			predicted_constraint(t1,t2) = 1;
			predicted_constraint(t2,t1) = 1;
			rankPosition = rankPosition + 1;
		%end
		i = i + 1;
	end
	spy(predicted_constraint(1:SeqLength,1:SeqLength),'r.',8)
	rectangle('position',[DIStart,DIStart,DIEnd-DIStart,DIEnd-DIStart]);
	eval(['print -dpdf -f'  ' ' Predicted_Constraints_name]);
%end
end
