% this code reads in the output of get_constraints.m and compares it to the
% experimentally solved crystal structure of a protein sequence from the
% alignment. 

% command: test_constraints('DYR_ECOLI_e3_n2_m40','DYR_ECOLI','3','1RX2.pdb')
% here evalue needs to be a string, like '3'. 

function test_constraints_PC(path,protein_name,evalue,pdbFilename, type) %type should be 'P', 'G' or 'M'
set(0,'DefaultFigureColor','w')
 
% input file names
indexMappingFilename = horzcat(path,'.indextableplus')

if type == 'P'
    DIScoreFilename = horzcat(path,'_L2_DIScores.csv');  
else
    if type == 'G'
    DIScoreFilename = horzcat(path,'_GREMLIN_DIScores.csv');  
    else
    DIScoreFilename = horzcat(path,'_DIScores.csv');  
    end
end

% output file names

if type == 'P'
    DIoutputfile = horzcat(path,'_L1_PSEUDO_DIScoresCompared.csv');  
else
    if type == 'G'
    DIoutputfile = horzcat(path,'_GREMLIN_DIScoresCompared.csv');  
    else
    DIoutputfile = horzcat(path,'_DIScoresCompared.csv');  
    end
end

FPoutputfile = horzcat(path,'_FPplot');
FPoutputfile2 = horzcat(path,'_SeqDists');
Cmapoutputfile = horzcat(path,'_Cmap');

% parameters
BATCH_MODE = 'no'; %use 'yes' to inhibit window dependant functions
CONTACT_MAP_EIC_COUNT = [200];
FALSE_POSITIVE_PLOT_EIC_COUNT = 200;
RESIDUE_PROXIMITY_THRESHOLD = 5;
RESIDUE_PROXIMITY_METHOD = 'nearest_atom';
%RESIDUE_PROXIMITY_METHOD = 'carbon_alpha';

% plot colour specification
NO_DATA_COLOR = [0.9,0.9,1.0];
CRYSTAL_CONTACT_COLOR = [.82,.82,.82];
EIC_COLOR = [1.0,.0,.0];
NON_EIC_COLOR = [.0,.0,.0];
PREDICTION_CONTACT_COLOR = [.0,1.0,.0];

%read in the index mapping file
fid = fopen(indexMappingFilename)
fgetl(fid); %discard header
uniprotIndexTable = textscan(fid,'%*s %c %c %s %s %s %c %c %s %s %s %s %s %s %s'); %aa ss_code ss_conf msa_ind msa_cons% msa_cons in_const atom# chain prot_coord aa xpos ypos zpos
fclose(fid);
proteinLength = size(uniprotIndexTable{1},1);
alignmentStartOffset = find(uniprotIndexTable{6} == '*',1);
alignmentEndOffset = find(uniprotIndexTable{6} == '*',1,'last');
constructStartOffset = find(uniprotIndexTable{7} == '*',1);
constructEndOffset = find(uniprotIndexTable{7} == '*',1,'last');
crystalStartOffset = find(~strcmp('-',uniprotIndexTable{10}),1);
crystalEndOffset = find(~strcmp('-',uniprotIndexTable{10}),1,'last');
chain_in_alignment = uniprotIndexTable{9}{find(cell2mat(uniprotIndexTable{9})' ~= '-',1)};
crystalGapOffset = zeros(1,proteinLength);
crystalGapOffset(cell2mat(uniprotIndexTable{9})' == '-') = 1;
%crystalGapOffset(find(strcmp('-',uniprotIndexTable{8})')) = 1;

DIScoreFilename
%load DI list
fid = fopen(DIScoreFilename);
DIScores = textscan(fid,'%f %f %f %f %f %f','Delimiter',','); %j i score del del del 
fclose(fid);
filtered_DI_index = ~(DIScores{4} == 0);
constraintTableLength = size(DIScores{1},1);

%load pdb
pdbStruct = pdbread([pdbFilename]);
atomStruct = pdbStruct.Model(1).Atom(:);

%TODO we may be ignoring needed positions in the HETATM array

% plots min all atom distance against DI ranking (x axis).
% Create Nx3 matrix of 3D coords
pdb_coords = [atomStruct.X;
	atomStruct.Y;
	atomStruct.Z;]';
AtomDistance = pdist(pdb_coords);
AtomDistanceMatrix = squareform(AtomDistance);
residue_is_CA = [atomStruct.resSeq]' == 1e+99;
for atom_i = 1:size(atomStruct,1)
	if (~strcmp(chain_in_alignment,atomStruct(atom_i).chainID))
		continue %atoms on different chains are not relevant, but can have duplicate resSeq offsets
	end
	residue_is_CA(atom_i,1) = strcmp(atomStruct(atom_i).AtomName,'CA');
	residue_combined_string = strcat(num2str(atomStruct(atom_i).resSeq),atomStruct(atom_i).iCode);
	for uniprot_i = 1:size(uniprotIndexTable{1},1)
		if (strcmp(residue_combined_string,uniprotIndexTable{10}{uniprot_i}))
			atom_to_uniprot_offset(atom_i) = uniprot_i;
		end
	end
end
MISSING_POSITION_INFO_MAGIC_NUMBER = 1e+6;
FALSE_POSITIVE_MISSING_POINT_Y = 0;
residueDistance = MISSING_POSITION_INFO_MAGIC_NUMBER * ones(proteinLength,proteinLength);
if (strcmp(RESIDUE_PROXIMITY_METHOD,'nearest_atom'))
	for uniprot_i = crystalStartOffset:crystalEndOffset
		Res_i_Atoms = atom_to_uniprot_offset == uniprot_i;
		for uniprot_j = uniprot_i + 1:crystalEndOffset
			Res_j_Atoms = atom_to_uniprot_offset == uniprot_j;
			DistHolder = AtomDistanceMatrix(Res_i_Atoms,Res_j_Atoms);
			minDistance = min(DistHolder(:));
			if ~isempty(minDistance)
				residueDistance(uniprot_i,uniprot_j) = min(DistHolder(:));
			end
			residueDistance(uniprot_j,uniprot_i) = residueDistance(uniprot_i,uniprot_j);
			clear DistHolder
		end
	end
elseif (strcmp(RESIDUE_PROXIMITY_METHOD,'carbon_alpha'))
	for uniprot_i = crystalStartOffset:crystalEndOffset
		Res_i_Atoms = atom_to_uniprot_offset == uniprot_i & residue_is_CA;
		for uniprot_j = uniprot_i + 1:crystalEndOffset
			Res_j_Atoms = atom_to_uniprot_offset == uniprot_j & residue_is_CA;
			DistHolder = AtomDistanceMatrix(Res_i_Atoms,Res_j_Atoms);
			minDistance = min(DistHolder(:));
			if ~isempty(minDistance)
				residueDistance(uniprot_i,uniprot_j) = min(DistHolder(:));
			end
			residueDistance(uniprot_j,uniprot_i) = residueDistance(uniprot_i,uniprot_j);
			clear DistHolder
		end
	end
else
	error('Error: unkonwn RESIDUE_PROXIMITY_METHOD')
end

% puts the structure distances into the DIScores table
for constraint_i = 1:constraintTableLength
	smaller_index = min(DIScores{1}(constraint_i),DIScores{2}(constraint_i));
	bigger_index = max(DIScores{1}(constraint_i),DIScores{2}(constraint_i));
	DIScores{4}(constraint_i) = residueDistance(smaller_index,bigger_index);
    SeqDist(constraint_i) = bigger_index-smaller_index;
end

%save('DIScoresCompared','DIScores')
%csvwrite(horzcat(proteinName,'_DIScoresCompared.csv'),DIScores);
csvwrite(DIoutputfile,DIScores);

distanceDescriptor = 'all atom';
if (strcmp(RESIDUE_PROXIMITY_METHOD,'carbon_alpha') == 1)
	distanceDescriptor = 'CA';
end

% plots the distance in structure for each pair vs the rank of the pair
if (strcmp(BATCH_MODE,'yes') ~= 1)
	ats = num2cell(1:1:FALSE_POSITIVE_PLOT_EIC_COUNT);
	EIC_indices = find(filtered_DI_index == 0);
	EIC_pointCount = min(size(EIC_indices,1), FALSE_POSITIVE_PLOT_EIC_COUNT);
	scatterEIC_X = EIC_indices(1:EIC_pointCount);
	scatterEIC_Y = DIScores{4}(scatterEIC_X);
	scatterEIC_Y(scatterEIC_Y == MISSING_POSITION_INFO_MAGIC_NUMBER) = FALSE_POSITIVE_MISSING_POINT_Y;
	nonEIC_indices = find(filtered_DI_index);
	scatterNonEIC_X = nonEIC_indices(nonEIC_indices < scatterEIC_X(EIC_pointCount));
	scatterNonEIC_Y = DIScores{4}(scatterNonEIC_X);
	scatterNonEIC_Y(scatterNonEIC_Y == MISSING_POSITION_INFO_MAGIC_NUMBER) = FALSE_POSITIVE_MISSING_POINT_Y;
	figurehandle2 = figure();
	axeshandle = axes('Parent',figurehandle2,'FontWeight','bold','FontSize',12);
	hold(axeshandle,'on');
	scatter(scatterEIC_X,scatterEIC_Y,'o','MarkerFaceColor',EIC_COLOR,'MarkerEdgeColor',EIC_COLOR)
	scatter(scatterNonEIC_X,scatterNonEIC_Y,'o','MarkerFaceColor',NON_EIC_COLOR,'MarkerEdgeColor',NON_EIC_COLOR)
	legend({'EIC pair','non-EIC DI pair'},'Location','Best','FontWeight','bold','FontSize',14);
	xlabel('DI rank','FontWeight','bold','FontSize',14);
	ylabel('Minimum distance in crystal structure','FontWeight','bold','FontSize',14);
	title(horzcat(protein_name,' / ',' Minimum ', distanceDescriptor, ' distance for top ', num2str(FALSE_POSITIVE_PLOT_EIC_COUNT), ' EIC pairs'),'Interpreter','none')
    eval(['print -dpdf -f' num2str(figurehandle2) ' ' FPoutputfile]);
	%saveas(figurehandle2,[protein_name '_FP_Plot.fig'],'fig');
end

% make contact map overlay plots
if (strcmp(BATCH_MODE,'yes') ~= 1)
	closeResidues = residueDistance < RESIDUE_PROXIMITY_THRESHOLD;
	[r,c] = find(closeResidues); % row and column indices of those atoms that are close
	crystalInteractions = sparse(r,c,1);
	crystalInteractions(proteinLength+1,proteinLength+1) = 0;
	if (exist('predictedPDBFilename'))
		closePredictionResidues = predictedResidueDistance < RESIDUE_PROXIMITY_THRESHOLD;
		[predicted_r,predicted_c] = find(closePredictionResidues);
		predictionInteractions = sparse(predicted_r,predicted_c,1);
		predictionInteractions(proteinLength+1,proteinLength+1) = 0;
	end
	figurehandle4 = figure()
	hold on
	residueOffsetMin = min(constructStartOffset,crystalStartOffset);
	residueOffsetMax = max(constructEndOffset,crystalEndOffset);
	noCrystalDataEnd = residueOffsetMin;
	noCrystalDataBegin = residueOffsetMin - 1;
	while (noCrystalDataEnd <= residueOffsetMax)
		if (crystalGapOffset(noCrystalDataEnd))
			if (noCrystalDataBegin == residueOffsetMin - 1)
				noCrystalDataBegin = noCrystalDataEnd; %set nodata start if not set yet
			end
			if ((noCrystalDataEnd == residueOffsetMax) && ~(noCrystalDataBegin > alignmentEndOffset) && ~(noCrystalDataEnd < alignmentStartOffset))
				barBegin = max(noCrystalDataBegin,alignmentStartOffset);
				barEnd = min(noCrystalDataEnd,alignmentEndOffset);
				rectangle('Position',[barBegin - 0.5,alignmentStartOffset - 0.5,(barEnd - barBegin) + 1,(alignmentEndOffset - alignmentStartOffset) + 1],'FaceColor',NO_DATA_COLOR,'LineStyle','None');
				rectangle('Position',[alignmentStartOffset - 0.5,barBegin - 0.5,(alignmentEndOffset - alignmentStartOffset) + 1,(barEnd - barBegin) + 1],'FaceColor',NO_DATA_COLOR,'LineStyle','None');
			end
		else
			if ((noCrystalDataBegin ~= residueOffsetMin - 1) && ~(noCrystalDataBegin > alignmentEndOffset) && ~(noCrystalDataEnd < alignmentStartOffset))
				barBegin = max(noCrystalDataBegin,alignmentStartOffset);
				barEnd = min(noCrystalDataEnd,alignmentEndOffset);
				if (barBegin < barEnd)
					rectangle('Position',[barBegin - 0.5,alignmentStartOffset - 0.5,(barEnd - barBegin),(alignmentEndOffset - alignmentStartOffset) + 1],'FaceColor',NO_DATA_COLOR,'LineStyle','None');
					rectangle('Position',[alignmentStartOffset - 0.5,barBegin - 0.5,(alignmentEndOffset - alignmentStartOffset) + 1,(barEnd  - barBegin)],'FaceColor',NO_DATA_COLOR,'LineStyle','None');
				end
				noCrystalDataBegin = residueOffsetMin - 1;
			end
		end
		noCrystalDataEnd = noCrystalDataEnd + 1;
	end
	EIC_Pair = zeros(size(crystalInteractions));
	foundConstraints = 0;
	index = 1;
	while (index <= size(DIScores{1},1) && foundConstraints < CONTACT_MAP_EIC_COUNT)
		if ~filtered_DI_index(index)
			t1 = DIScores{1}(index);
			t2 = DIScores{2}(index);
			EIC_Pair(t1,t2) = 1; % this should be fine as it just
			%reads in the numbers from DIScores, so should avoid missing numbers.
			EIC_Pair(t2,t1) = 1;
			foundConstraints = foundConstraints + 1;
		end
		index = index + 1;
	end
	axis([residueOffsetMin - 0.5,residueOffsetMax + 0.5,residueOffsetMin - 0.5,residueOffsetMax + 0.5]);
	set(gca,'Box','on');
	set(gca,'YDir','reverse');
	set(gca,'DataAspectRatio',[1,1,1]);
	[plotRows,plotCols] = find(crystalInteractions);
	plot(plotRows,plotCols,'o','MarkerSize',6,'MarkerFaceColor',CRYSTAL_CONTACT_COLOR,'MarkerEdgeColor',CRYSTAL_CONTACT_COLOR);
	if (exist('predictedPDBFilename'))
		[plotRows,plotCols] = find(predictionInteractions);
		plot(plotRows,plotCols,'s','MarkerSize',8,'MarkerEdgeColor',PREDICTION_CONTACT_COLOR);
	end
	[plotRows,plotCols] = find(EIC_Pair);
    % plotRows = plotRows(33:end) - 9;
    % plotCols = plotCols(33:end) - 9;
	plot(plotRows,plotCols,'*','MarkerSize',4,'MarkerFaceColor',EIC_COLOR,'MarkerEdgeColor',EIC_COLOR);
	xlabel(['number of constraints = ' num2str(size(plotRows,1)/2)]);
	rectangle('position',[alignmentStartOffset - 0.5,alignmentStartOffset - 0.5,alignmentEndOffset - alignmentStartOffset + 1,alignmentEndOffset - alignmentStartOffset + 1]);
	title(sprintf('%s, red e%s alignment, grey %s residues < %.2f A apart',protein_name,evalue,pdbFilename,RESIDUE_PROXIMITY_THRESHOLD),'Interpreter','none')
	eval(['print -dpdf -f' num2str(figurehandle4) ' ' horzcat(Cmapoutputfile,'_',num2str(CONTACT_MAP_EIC_COUNT))]);
    %eval(['print -dpdf -f' num2str(figurehandle2) ' ' FPoutputfile]);
	%saveas(figurehandle4,[protein_name '_ContactMap.fig'],'fig');
end
%exit
end %function

