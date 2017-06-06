% this is the first file that you run. 
% It takes in a protein sequence alignmnet, which will always be a .fas
% file, and analyzes the covariance structure of the sequences. 

function Frob = covariance_analysis2(msa_fasta_filename, seqid_of_interest, outputfile, pseudocount_weight, l2, apc)
% example of how to call the function: covariance_analysis('DYR_ECOLI_e3_n2_m40.fas','DYR_ECOLI','DYR_ECOLI_e3_n2_m40_MI_DIs.txt')

% parameters
theta = 0.3; % this is for the sequence weighting - if theta = 0.3 then more than 70% identical. 

% function calls
[Pij_true, Pi_true, alignment_width, q, encoded_seq_of_interest, focus_to_uniprot_offset_map, W] = read_alignment(msa_fasta_filename, seqid_of_interest, theta);
[Pij, Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q);
number2letter_map = create_number2letter_map();
C = Compute_C_gap(Pij, Pi, alignment_width, q);  

lambda = 0.001;
if l2 == 'Y'
	penalty =  lambda * eye(length(C),length(C));%l2 reg
	C = (penalty + 0.25 *C^2)^0.5 + 0.5 *C ;
end

invC = inv(C);
for i=1:(alignment_width)
    for j=(i):alignment_width
        W_frob = ReturnW_frob_gap(invC, i, j, q);
		% this just pulls out Cij 20 by 20 block. 
        [Frob(i,j), ~] = frob_link(W_frob);       
        [Frob(j,i), ~] = frob_link(W_frob);       
    end
end

%APC
if apc == 'Y'
    row_mean = mean(Frob,2);
    all_mean = abs(mean(abs(row_mean)));
    Frob = Frob - (row_mean * (row_mean'))/all_mean;
    Frob = Frob - min(min(Frob)); % A constant is added so that there are no values included in the actual plot that are less than zero
end

% print results file
fp = fopen(outputfile, 'w');
for i=1:(alignment_width)
    for j=(i):alignment_width
        % print mutual information
        [MI_true, ~, ~] = calculate_mi(i, j, Pij_true, Pi_true, q);
        fprintf(fp, '%d %s %d %s %g ', focus_to_uniprot_offset_map(i), number2letter_map(encoded_seq_of_interest(i)), focus_to_uniprot_offset_map(j), number2letter_map(encoded_seq_of_interest(j)), MI_true);
        fprintf(fp, '%g ', Frob(i,j));
        fprintf(fp, '%g ', Frob(i,j));
        fprintf(fp, '\n');
    end
end
fclose(fp);
end

% these first two functions (read_alignment and read_alignment_fasta) simply
% count the frequency with which each amino acid occurs in each column, and
% with which each pair of amino acids occurs in each pair of columns. 

function [Pij_true, Pi_true, alignment_width, q, encoded_seq_of_interest, focus_to_uniprot_offset_map, W] = read_alignment(msa_fasta_filename, seqid_of_interest, theta)
% this is the call to the data cleaning function read_alignmnet_fasta
[encoded_focus_alignment, focus_index_of_interest, focus_to_uniprot_offset_map] = read_alignment_fasta(msa_fasta_filename, seqid_of_interest);
encoded_seq_of_interest = encoded_focus_alignment(focus_index_of_interest,:);
[alignment_height,alignment_width] = size(encoded_focus_alignment);
W = ones(1, alignment_height);
if(theta > 0.0)   
    W = (1./(1+sum(squareform(pdist(encoded_focus_alignment, 'hamm') < theta))));    
end

Meff=sum(W); % effective number of sequence in alignment counted using sequence weights

hist(W,200)
q = max(max(encoded_focus_alignment));
Pij_true = zeros(alignment_width, alignment_width, q, q);
Pi_true = zeros(alignment_width, q);

% single columns counts
for j=1:alignment_height
	for i=1:alignment_width
		Pi_true(i, encoded_focus_alignment(j, i)) = Pi_true(i, encoded_focus_alignment(j, i)) + W(j);
    end
end
Pi_true = Pi_true/Meff;

% pair column counts
for l=1:alignment_height
	for i=1:alignment_width-1
		for j=i+1:alignment_width
			Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j)) = Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j)) + W(l);
			Pij_true(j, i, encoded_focus_alignment(l, j), encoded_focus_alignment(l, i)) = Pij_true(i, j, encoded_focus_alignment(l, i), encoded_focus_alignment(l, j));
		end
	end
end
Pij_true = Pij_true/Meff;

scra = eye(q, q);
for i=1:alignment_width
	for alpha=1:q
		for beta=1:q
			Pij_true(i, i, alpha, beta) = Pi_true(i, alpha) * scra(alpha, beta);
		end
	end
end
end

% (this just does some preprocessing of the data) . 
function [encoded_focus_alignment, focus_index_of_interest, focus_to_uniprot_offset_map] = read_alignment_fasta(msa_fasta_filename, seqid_of_interest)
% sorts out ambiguous resiudes and how to handle them. 
METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = getenv('DI_METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES'); % 1 = change them to gaps .. 2 = mask entire sequence
if (size(METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES,2) == 0)
	METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = 2;
else
	METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = str2num(METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES);
end
% read in alignment. 
full_alignment = fastaread(msa_fasta_filename);
alignment_width = size(full_alignment(1).Sequence, 2);
alignment_height = size(full_alignment, 1);
letter2number_map = create_letter2number_map();
[full_index_of_interest, range_of_interest_start, range_of_interest_end] = find_seq_of_interest(full_alignment, seqid_of_interest);
encoded_focus_alignment = [];
skipped_sequence_counter = 0;
[focuscolumnlist, focus_to_uniprot_offset_map] = scan_sequence_of_interest_for_focus_columns(full_alignment(full_index_of_interest).Sequence, range_of_interest_start, letter2number_map);
for full_alignment_index=1:alignment_height
	focus_alignment_row = full_alignment(full_alignment_index).Sequence(focuscolumnlist);
	encoded_focus_alignment_row = letter2number_map(focus_alignment_row);
	if (size(find(encoded_focus_alignment_row == 0),2) > 0)
		error(['Error: sequence in alignment has illegal characters: ' full_alignment(full_alignment_index).Sequence]);
	end
	if (size(find(encoded_focus_alignment_row <= -2),2) > 0)
		error(['Error: sequence in alignment has dot or lowercase in conserved position: ' full_alignment(full_alignment_index).Sequence]);
	end
	if (size(find(encoded_focus_alignment_row == -1),2) > 0)
		if (METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES == 1)
			encoded_focus_alignment_row(find(encoded_focus_alignment_row == -1)) = 1;
		else
			if (METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES == 2)
				continue %skip sequences with ambiguous residues
			else
				error('Internal Error');
			end
		end
	end
	encoded_focus_alignment(size(encoded_focus_alignment,1) + 1,:) = encoded_focus_alignment_row;
	if (full_alignment_index == full_index_of_interest)
		focus_index_of_interest = size(encoded_focus_alignment,1);
	end
end
end

% this adds pseudocounts to the frequency counts from the data. 
% because the input data, the set of aligned sequences, is undersampled, we
% add random sequences to the input data set - these are called
% pseudocounts. 
% with_pc = with_pseudocounts. 
function [Pij, Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q)
Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(alignment_width, alignment_width, q, q);
Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(alignment_width, q);
scra = eye(q);
for i=1:alignment_width
	for alpha = 1:q
		for beta = 1:q
			Pij(i, i, alpha, beta) = (1.-scra(alpha, beta)*pseudocount_weight)*Pij_true(i, i, alpha, beta) + pseudocount_weight/q*scra(alpha, beta);
		end
	end
end
end

function C = Compute_C_gap(Pij, Pi, alignment_width, q)
C=zeros(alignment_width*(q-1), alignment_width*(q-1));
for i=1:alignment_width
	for j=1:alignment_width
		for alpha=1:q-1
			for beta=1:q-1
				 C(mapkey(i, alpha, q), mapkey(j, beta, q)) = Pij(i, j, alpha+1, beta+1) - Pi(i, alpha+1)*Pi(j, beta+1);
			end
		end
	end
end
end

function A=mapkey(i, alpha, q)
A = (q-1)*(i-1)+alpha;
end

function W=ReturnW_gap(C, i, j, q)
W = ones(q, q);
W(2:q, 2:q) = exp(-C(mapkey(i, 1:q-1, q), mapkey(j, 1:q-1, q)));
end

function W=ReturnW_frob_gap(C, i, j, q)
W = zeros(q, q);
W(2:q, 2:q) = C(mapkey(i, 1:q-1, q), mapkey(j, 1:q-1, q));
end

function [Frob, W] = frob_link(W)
q = 21;


for i = 1:q
    W(:,i) = W(:,i) - mean(W(:,i));
end
for i = 1:q
    W(i,:) = W(i,:) - mean(W(i,:));
end

Frob = norm(W,'fro'); %This calculates the 2-norm of the column vector, W(:).
% The 2-norm is equal to the Euclidean length of the vector.
return;
end

function [M, s1, s2] = calculate_mi(i, j, P2, P1, q)
M = 0.;
for alpha=1:q
	for beta = 1:q
		 if(P2(i, j, alpha, beta)>0)
			M = M + P2(i, j, alpha, beta)*log(P2(i, j, alpha, beta) / P1(i, alpha)/P1(j, beta));
		end
	end
end

s1=0.;
s2=0.;
for alpha=1:q
	if(P1(i, alpha)>0)
		s1 = s1 - P1(i, alpha) * log(P1(i, alpha));
	end
	if(P1(j, alpha)>0)
		s2 = s2 - P1(j, alpha) * log(P1(j, alpha));
	end
end
end

% this is all data processing. 

function [index_of_interest, range_start, range_end] = find_seq_of_interest(full_alignment, seqid_of_interest)
index_of_interest = -1;
for scan_index = 1:size(full_alignment,1)
	[seqid, range_start, range_end] = split_uniprot_id(full_alignment(scan_index).Header);
	if (strcmp(seqid,seqid_of_interest) == 1)
		index_of_interest = scan_index;
		break
	end
end
if (index_of_interest == -1)
	error(['Error: could not find sequence of interest (' seqid_of_interest ') in multiple sequence alignment']);
end
end

function [seqid, range_start, range_end] = split_uniprot_id(pfam_uniprot_range_line)
slashposition = findstr('/', pfam_uniprot_range_line);
if (size(slashposition,2) ~= 1 || slashposition == 1 || slashposition == size(pfam_uniprot_range_line,2))
	error(['Error: could not parse (slash error) uniprot range line in pfam alignment : ' pfam_uniprot_range_line]);
end
seqid = pfam_uniprot_range_line(1:slashposition - 1);
rangestring = pfam_uniprot_range_line(slashposition + 1:size(pfam_uniprot_range_line,2));
hyphenposition = findstr('-', rangestring);
if (size(hyphenposition,2) ~= 1 || hyphenposition == 1 || hyphenposition == size(rangestring,2))
	error(['Error: could not parse (hyphen error) uniprot range line in pfam alignment : ' pfam_uniprot_range_line]);
end
range_start = str2num(rangestring(1:hyphenposition - 1));
range_end = str2num(rangestring(hyphenposition + 1 : size(rangestring,2)));
if (isempty(range_start) || isempty(range_end))
	error(['Error: could not parse (range start/end) uniprot range line in pfam alignment : ' pfam_uniprot_range_line]);
end
end

function [focuscolumnlist, uniprotoffsetlist] = scan_sequence_of_interest_for_focus_columns(sequence_of_interest, range_of_interest_start, letter2number_map)
focuscolumnlist = [];
uniprotoffsetlist = [];
next_uniprotoffset = range_of_interest_start;
for pos=1:size(sequence_of_interest,2)
	residuecode = letter2number_map(sequence_of_interest(pos));
	if (residuecode == 0)
		error(['Error: sequence of interest contains undefined residues:' sequence_of_interest]);
	end
	if (residuecode == -1)
		error(['Error: sequence of interest contains ambiguous residues:' sequence_of_interest]);
	end
	if (residuecode > 1)
		focuscolumnlist = [focuscolumnlist pos];
		uniprotoffsetlist = [uniprotoffsetlist next_uniprotoffset];
	end
	if (residuecode == -2 || residuecode > 1)
		next_uniprotoffset = next_uniprotoffset + 1;
	end
end
end

function letter2number_map = create_letter2number_map()
letter2number_map(256) = 0; %initiallize all bytes to 0
letter2number_map('-') = 1;
letter2number_map('A') = 2;
letter2number_map('C') = 3;
letter2number_map('D') = 4;
letter2number_map('E') = 5;
letter2number_map('F') = 6;
letter2number_map('G') = 7;
letter2number_map('H') = 8;
letter2number_map('I') = 9;
letter2number_map('K') = 10;
letter2number_map('L') = 11;
letter2number_map('M') = 12;
letter2number_map('N') = 13;
letter2number_map('P') = 14;
letter2number_map('Q') = 15;
letter2number_map('R') = 16;
letter2number_map('S') = 17;
letter2number_map('T') = 18;
letter2number_map('V') = 19;
letter2number_map('W') = 20;
letter2number_map('Y') = 21;
letter2number_map('B') = -1; %ambiguous : skip sequences containing these
letter2number_map('Z') = -1; %ambiguous : skip sequences containing these
letter2number_map('J') = -1; %ambiguous : skip sequences containing these
letter2number_map('X') = -1; %ambiguous : skip sequences containing these
letter2number_map('U') = -1; %non-standard : skip sequences containing these
letter2number_map('O') = -1; %non-standard : skip sequences containing these
letter2number_map('a') = -2; %non-conserved: skip in seq of interest
letter2number_map('c') = -2; %non-conserved: skip in seq of interest
letter2number_map('d') = -2; %non-conserved: skip in seq of interest
letter2number_map('e') = -2; %non-conserved: skip in seq of interest
letter2number_map('f') = -2; %non-conserved: skip in seq of interest
letter2number_map('g') = -2; %non-conserved: skip in seq of interest
letter2number_map('h') = -2; %non-conserved: skip in seq of interest
letter2number_map('i') = -2; %non-conserved: skip in seq of interest
letter2number_map('k') = -2; %non-conserved: skip in seq of interest
letter2number_map('l') = -2; %non-conserved: skip in seq of interest
letter2number_map('m') = -2; %non-conserved: skip in seq of interest
letter2number_map('n') = -2; %non-conserved: skip in seq of interest
letter2number_map('p') = -2; %non-conserved: skip in seq of interest
letter2number_map('q') = -2; %non-conserved: skip in seq of interest
letter2number_map('r') = -2; %non-conserved: skip in seq of interest
letter2number_map('s') = -2; %non-conserved: skip in seq of interest
letter2number_map('t') = -2; %non-conserved: skip in seq of interest
letter2number_map('v') = -2; %non-conserved: skip in seq of interest
letter2number_map('w') = -2; %non-conserved: skip in seq of interest
letter2number_map('y') = -2; %non-conserved: skip in seq of interest
letter2number_map('b') = -2; %non-conserved: skip in seq of interest
letter2number_map('z') = -2; %non-conserved: skip in seq of interest
letter2number_map('j') = -2; %non-conserved: skip in seq of interest
letter2number_map('x') = -2; %non-conserved: skip in seq of interest
letter2number_map('u') = -2; %non-conserved: skip in seq of interest
letter2number_map('o') = -2; %non-conserved: skip in seq of interest
letter2number_map('.') = -3; %non-conserved: skip in seq of interest, do not advance position
end

function number2letter_map = create_number2letter_map()
number2letter_map(1) = '-';
number2letter_map(2) = 'A';
number2letter_map(3) = 'C';
number2letter_map(4) = 'D';
number2letter_map(5) = 'E';
number2letter_map(6) = 'F';
number2letter_map(7) = 'G';
number2letter_map(8) = 'H';
number2letter_map(9) = 'I';
number2letter_map(10) = 'K';
number2letter_map(11) = 'L';
number2letter_map(12) = 'M';
number2letter_map(13) = 'N';
number2letter_map(14) = 'P';
number2letter_map(15) = 'Q';
number2letter_map(16) = 'R';
number2letter_map(17) = 'S';
number2letter_map(18) = 'T';
number2letter_map(19) = 'V';
number2letter_map(20) = 'W';
number2letter_map(21) = 'Y';
end
