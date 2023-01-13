function [Results ,ResultsRed] = AnalysingHepaNew(EnzymeDisorderString,EnzymeDisorderVal,EnzymeKOOverString,EnzymeKOOverValLower,EnzymeKOOverValUpper,DoniniLowerVal,DoniniUpperVal,CreatingModel)
%This function computes the fluxes for the 442 metabolic objectives.
%Input:
%
%EnzymeDisorderString: string representing the enzyme related to a disorder.
%Example: a) for WT simulations, EnzymeDisorderString = 'WT'
%         b) for PH1 simulations, EnzymeDisorderString = 'EC:2.6.1.44'
%
%EnzymeDisorderVal: value for the UpperBound of reaction/reactions related 
%to disorder. This value is considered only if EnzymeDisorderString is not equal to 'WT'.
%Example: for PH1 simulations, if you set EnzymeDisorderString = 'EC:2.6.1.44'
%         and EnzymeDisorderVal = 0, then you'll simulate the KO of 'EC:2.6.1.44'
%
%EnzymeKOOverString = string representing enzyme for overexpression and
%downregulation.
%
%EnzymeKOOverValLower and EnzymeKOOverValUpper = values for simulating
%overexpression and downregulation.
%
%Example: if you set EnzymeDisorderString = 'EC:2.6.1.44'
%                    EnzymeDisorderVal = 0
%                    EnzymeKOOverString = 'EC2.6.1.2'
%                    EnzymeKOOverValLower = 500 
%                    EnzymeKOOverValUpper = 1500
%youl'll simulate the overexpression of GPT in PH1, because the lover value
%for a rection is set to 0 while the upper value is set to 1000.
%
%DoniniLowerVal and DoniniUpperVal = values to set the upper bounds and
%lower bounds for the non canonical AGT-like reaction.
%
%Ouput: 
%
%Results: a matrix storing the results for the extend flux minimization problems for each metabolic objective
%
%ResultsRed =  a matrix storing a value for each of the reaction of the model according
%with the principle of flux minimization of H. G. Holzhutter.



%%%Selecting optimization problem solver%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%changeCobraSolver('glpk','LP');

%%%Output variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results = [];

ResultsRed = [];

%%%Variables for simulations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InsertedReaction = 52;

n_Model = 442;

%%%Loading Equilibrium Constants
EquilibriumConstants = load('EquilibriumConstants');

Keq = EquilibriumConstants.Keq;

dGo = EquilibriumConstants.dGo;

if strcmp(CreatingModel,'Y')
    
    StrModelPath = 'CreatedMetabolicObjectivesMatlab/Hepa1MetabolicObjSlim_';
    
else
    
    StrModelPath = 'MetabolicObjectivesMatLab/Hepa1MetabolicObjSlim_';
    
end

%%%Computing solution for each of the 442 objective functions%%%%%%%%%%%%%%
for i = 1:n_Model
    
    %%%Loading model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Model = load(strcat(StrModelPath,num2str(i)),'Model');
    
    %Model = load(strcat('Hepa1MetabolicObjSlim_',num2str(i)),'Model');
    %rxnGeneMat = load('MetabolicObjectivesPhayton//Enzymes.mat');
    
    %%%Setting AGXT lower bound for WT simulations%%%%%%%%%%%%%%%%%%%%%%%%%
    SimModel = Model.Model;
    
    %SimModel.rxnGeneMat = rxnGeneMat.rxnGeneMat;
    
    SimModel.lb(5080) = 500;
    
    SimModel.lb(5079) = 0;
    
    SimModel.ub(5079) = 1000;
    
    
    %%%Enzyme Disorder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(EnzymeDisorderString,'WT') == 0
        
        %Enzymes Disorder
        [IndexesEnzymes,IndexesEnzymes2] = linkingEnzymeReactions(SimModel,EnzymeDisorderString,dGo,InsertedReaction);
        
        for xyz = 1:length(IndexesEnzymes)
            
            SimModel = changeRxnBounds(SimModel,SimModel.rxns{IndexesEnzymes(xyz),1},EnzymeDisorderVal,'b');
            
        end
        
        for xyz = 1:length(IndexesEnzymes2)
            
            SimModel = changeRxnBounds(SimModel,SimModel.rxns{IndexesEnzymes2(xyz),1},EnzymeDisorderVal,'b');
            
        end
        
    end
        
    %%%%%%%%%%%%DoniniLower,UpperBound%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SimModel.lb(2539*2+25) = DoniniLowerVal;
    
    SimModel.ub(2539*2+25) = DoniniUpperVal;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Enzyme KO Over EnzymeKOOverVal == 0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(EnzymeKOOverString,'WT') == 0
        
        %Enzymes Overexpression
        [IndexesEnzymes1, IndexesEnzymes2] = linkingEnzymeReactions(SimModel,EnzymeKOOverString, dGo,InsertedReaction);
        
        if EnzymeKOOverValUpper == 0
    
            for xyz = 1:length(IndexesEnzymes1)
            
                SimModel = changeRxnBounds(SimModel,SimModel.rxns{IndexesEnzymes1(xyz),1},EnzymeKOOverValUpper,'b');
            
            end
            
            for xyz = 1:length(IndexesEnzymes2)
            
                SimModel = changeRxnBounds(SimModel,SimModel.rxns{IndexesEnzymes2(xyz),1},EnzymeKOOverValUpper,'b');
            
            end
            
        else
            
            for xyz = 1:length(IndexesEnzymes1)
            
                SimModel = changeRxnBounds(SimModel,SimModel.rxns{IndexesEnzymes1(xyz),1},EnzymeKOOverValUpper,'u');
                
                SimModel = changeRxnBounds(SimModel,SimModel.rxns{IndexesEnzymes1(xyz),1},EnzymeKOOverValLower,'l');
            
            end
            
            for xyz = 1:length(IndexesEnzymes2)
            
                SimModel = changeRxnBounds(SimModel,SimModel.rxns{IndexesEnzymes2(xyz),1},0,'b');
            
            end
            
        end
        
    end
    
    %%%Computing solution and saving solutions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nS,mS] = size(SimModel.S);
    
    Solution = optimizeCbModel(SimModel,'min');

    if isempty(Solution.x) ~= 1
        
        for j = 1:length(Solution.x)
        
            Results(j,i) = Solution.x(j);
            
        end
        
        for z = 1:2539
                
                if Solution.x(z) > 0
                
                    %Result contains the intere vectors of fluxes
                    %the first mExt/2 reactions represent the Hepatonet3
                    %reactions
                    ResultsRed(z,i) = Solution.x(z);
                    
                else
                    
                    ResultsRed(z,i) = -Solution.x(z+2539);
                    
                end
                
                if abs(ResultsRed(z,i)) <= 1e-10
                    
                    ResultsRed(z,i) = 0;
                    
                end
                
        end
        
        t = 1;
        
        for z = 2539*2+1:length(Solution.x)
            
            if abs(Solution.x(z)) > 1e-10
                
                    ResultsRed(2539+t,i) = Solution.x(z);
                    
            else
                
                ResultsRed(2539+t,i) = Solution.x(z);
                
            end
            
            t = t  + 1; 
            
        end
                
    
    end
    
end
