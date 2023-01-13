function [A,B,WT_vector, KO_vector] = computingReaction_WT_KO_vectorversion(WT_vector, KO_vector,t)
%This function computes the ranked list of reactions 
%
%Input:
%
%WT_vector = a vector containing the mean of the wild-type fluxes accross the different metabolic objective
%
%KO_vector = a vector containing the mean of the disease fluxes accross the different metabolic objective 
%
%t = string for final txt output
%
%Example = if t = 'PH1' then the output will be ReactionsRank_PH1.txt
%
%Output: 
%
%A txt file contains the metabolites rank according to DFA



%Loading model for final print
Hepa2 = load('HepaforRankedLists');

HepaModel2 = Hepa2.Hepa2;


%Computing vector difference
A = WT_vector - KO_vector;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computing the dimension of matrix A
[n,m] = size(A)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:n
%     
%     j = 1;
%         
%         if ((WT_vector(i,j)) + (KO_vector(i,j))) == 0
%             
%             B(i,j) = 0;
%             
%         else
%             
%             B(i,j) = ((WT_vector(i,j)) - (KO_vector(i,j)))/((WT_vector(i,j)) + (KO_vector(i,j)));
%             
%         end
%         
% end

%Computing the dimension of stoichiometrix matrix of HepaModel
[ns,ms] = size(HepaModel2.S);

%Sorting 
[Sort_Mean, IX_Mean] = sort(A, 'descend');

fid = fopen(strcat('ReactionsRank_',t,'.txt'), 'w');

fprintf( fid, '%s\t%s\t%s\t%s\t%s\t%s\n', 'Reactions', 'Enzyme/Transport', 'Delta Flux','Abs Delta Flux','WT Reaction Sign','PH1 Reaction Sign');
    

for i = 1:n
    
    x = printRxnFormula(HepaModel2, HepaModel2.rxns{IX_Mean(i)});

    fprintf( fid, '%s\t%s\t%u\t%u\t%u\t%u\n', x{1}, HepaModel2.rxnGeneMat{IX_Mean(i)}, Sort_Mean(i),abs(Sort_Mean(i)),sign(WT_vector(IX_Mean(i))),sign(KO_vector(IX_Mean(i))));
    
end

%Closing text files
fclose(fid);

