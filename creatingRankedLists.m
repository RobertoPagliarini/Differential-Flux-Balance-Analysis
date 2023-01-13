function [MeanFluxesWT, MeanFluxesPH1] = creatingRankedLists(CreatingModel,OutString)
%Input: 
%CreatingModel = string that allows to create objective functions or to
%upload the existent models. 
%If CreatingModel = 'Y' then the function creates models, otherwise the
%function uploads existent models.
%
%OutString = string for final txt output
%
%Example = if OutString = 'PH1' then the output will be
%ReactionsRank_PH1.txt and MetabolitesRank_PH1.txt
%
%Output =
%
%Two vectors conyaining mean fluxes in WT and PH1


%%%paths to Cobra toolbox and solver%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath(genpath('/Applications/MATLAB_R2012a.app/toolbox/glpkmex'));

%addpath(genpath('/Applications/MATLAB_R2012a.app/toolbox/cobra'));


%%%Model creation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(CreatingModel,'Y')
    
    creatingHepaMetabolicObjectivesStrcut;
    
end

%%%Computing WT and PH1 flux distributions accros the 442 objective functions
[ResultsWT, ResultsRedWT] = AnalysingHepaNew('WT',0,'WT',0,0,0,1000,CreatingModel);

[ResultsPH1, ResultsRedPH1] = AnalysingHepaNew('EC:2.6.1.44',0, 'WT',0,0,0,1000,CreatingModel);

%%%Computing mean of fluxdistributions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MeanFluxesWT(1:2539+53,1) = mean([ResultsRedWT(1:2539+53,:)]')';

MeanFluxesPH1(1:2539+53,1) = mean([ResultsRedPH1(1:2539+53,:)]')';

%%% Computing ranked lists%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
computingMM_WT_KO2_vectorversion(MeanFluxesWT,MeanFluxesPH1,OutString);

computingReaction_WT_KO_vectorversion(MeanFluxesWT, MeanFluxesPH1,OutString);