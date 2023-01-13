%This script creates the 442 metabolic objectives for the simulations

mkdir('CreatedMetabolicObjectivesMatlab')

Hepatonet2 = load('HepatonetModels.mat');

Hepatonet2 = Hepatonet2.Hepatonet2;

EvalPerc = 1;

ReactionNumber = 5131;

BoundsChange = 0;

%reading the text file
[NumTestBase,TxtTestBase,RawTestBase] = xlsread('HepaTestsBase.xls');

NumList = PosNegObjNumerical('HepaPosNegObj.xls');

for i = 3:length(RawTestBase)
    
        [Model] = creatingHepaObjective(Hepatonet2, ReactionNumber, EvalPerc, RawTestBase,i, BoundsChange,NumList);
        
        save(strcat('Hepa1MetabolicObjSlim_',num2str(i-2)),'Model');
        
        movefile(strcat('Hepa1MetabolicObjSlim_',num2str(i-2),'.mat'),'CreatedMetabolicObjectivesMatlab');
end
            
    
    
