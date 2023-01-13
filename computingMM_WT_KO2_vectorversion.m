function [SumAbsDelta] = computingMM_WT_KO2_vectorversion(WT_vector, KO_vector,t)
%This function computes the ranked list of metabolites according to DFA
%
%Input:
%
%WT_vector = a vector containing the mean of the wild-type fluxes accross the different metabolic objective
%
%KO_vector = a vector containing the mean of the disease fluxes accross the different metabolic objective 
%
%t = string for final txt output
%
%Example = if t = 'PH1' then the output will be MetabolitesRank_PH1.txt
%
%Output: 
%
%A txt file contains the metabolites rank according to DFA

%Loading model for final print
Hepa2 = load('HepaforRankedLists');

HepaModel2 = Hepa2.Hepa2;


%Computing differences using mean 2
A = WT_vector - KO_vector;

%Computing the dimension of matrix A
[n,m] = size(A);

%Computing the dimension of stoichiometrix matrix of HepaModel
[ns,ms] = size(HepaModel2.S);

Mean = A;

AbsMean = (abs(A));

for i = 1:ns
    
    Pipp = find(HepaModel2.S(i,:)~=0);
    
    SumDelta(i,1) = 0;
    
    SumAbsDelta(i,1) = 0;
    
    %SumAbsDelta_WT(i,1) = 0;
    
    %SumAbsDelta_KO(i,1) = 0;
    
    for j = 1:length(Pipp)
        
            %SumDelta(i,1) = SumDelta(i,1) + Mean(Pipp(j),1);
            
            SumAbsDelta(i,1) = SumAbsDelta(i,1) + AbsMean(Pipp(j),1);
            
            %SumAbsDelta_WT(i,1) = SumAbsDelta_WT(i,1) + abs(WT_vector(Pipp(j),1));
            
            %SumAbsDelta_KO(i,1) = SumAbsDelta_KO(i,1) + abs(KO_vector(Pipp(j),1));
            
            
    end
        
            
    
end


    %%%Sorting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %[Sort_SumDelta,IX_Sum_Delta] = sort(SumDelta,'descend');

    [Sort_SumAbsDelta,IX_Sum_AbsDelta] = sort(SumAbsDelta,'descend');

    %fid_SortSumDelta = fopen(strcat('MetSortSumDelta2_',t,'.txt'),'w');

    fid_SortAbsSumDelta = fopen(strcat('MetabolitesRank_',t,'.txt'),'w');

    fprintf(fid_SortAbsSumDelta, '%s\t%s\n','Metabolote','DFA Value');
    
    %%%Writing on txt%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(Sort_SumAbsDelta)
        
        fprintf(fid_SortAbsSumDelta, '%s\t%u\n', HepaModel2.metNames{IX_Sum_AbsDelta(i)},Sort_SumAbsDelta(i));
    
        
    end

    
    fclose(fid_SortAbsSumDelta);

    

