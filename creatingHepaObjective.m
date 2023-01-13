function [Prova1, foundOut,ProvaNoEval,ConReactionNumber,number_inserted_reactions,prodreac,ReactionNumber] = creatingHepaObjective(Hepatonet1,ReactionNumber,EvalPerc,RawTestBase,i,BoundsChange ,NumList)

        %copying Hepatonet model
        Prova1 = Hepatonet1;
        ProvaNoEval = Hepatonet1;
        
        %changing the links to the metabolites
        Prova1.mets = Prova1.metNames;
        ProvaNoEval.mets = ProvaNoEval.metNames;
        
        %weight for metabolic objective
        %[f] = createHepaWeights(Prova1.rxns);

        %Prova1.c = f;
        
        %ProvaNoEval.c = f;
    
        %our objective is constraints in xls file
        ReadingObj = RawTestBase(i,3);
    
        ReadingObj = ReadingObj{1};
    
        if isnan(ReadingObj) ~= 1
       
            MetabolicObj = textscan(ReadingObj,'%s');
            
        else
            
            MetabolicObj = {};
            
            MetabolicObj{1,1}{1,1} = NaN(1,1);
            
        end
    
        %reading evaluator. This is evaluator also in xls file
        ReadingEval = RawTestBase(i,4);
        
        ReadingEval = ReadingEval{1};
        
        ReadingEval = '';
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%da qui in su abbiamo commentato%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
    
        %we now read all the reactions in the constraints
        for j = 1:length(MetabolicObj{1,1})
        
            %each cell in the MetabolicObj is a new reaction linking the
            %cell with its environment
            ReactionNumber = ReactionNumber + 1;
            
            %the j-esimo metabolite in the cell
            p = MetabolicObj{1,1}{j,1};
        
            %if the metabolite has a - then we insert a reaction which
            %introduces these metabolites
            if (p(1) == '-')
            
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
            
                ConReaction = strcat('->', p);
            
                ConReaction(3) = ' ';
            
                Prova1 = addReaction(Prova1, ConReactionNumber, ConReaction);
                
                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber, ConReaction);
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                 
%                 if strcmp(ReadingEval,p(2:end))
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                     
%                 end
%             
            
            
            end
        
            if (p(1) == '=')
            
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
            
                ConReaction = strcat(p(2:end), ' <=>');
            
                Prova1 = addReaction(Prova1, ConReactionNumber, ConReaction);
                
                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber, ConReaction);
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%                 if strcmp(ReadingEval,p(2:end))
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
%                 
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
%                 
%                 end
%             
            end
        
        
            if (p(1) == '+')
            
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
            
                ConReaction = strcat(p(2:end), ' ->');
            
                Prova1 = addReaction(Prova1, ConReactionNumber, ConReaction);
                
                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber, ConReaction);
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%                 if strcmp(ReadingEval,p(2:end))
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
%             

            
            end
            
            %%%%%%controllare questa parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (isnan(p(1)) & strcmp(ReadingEval,'exclude(s)') == 0)
            
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
            
                ConReaction = strcat(ReadingEval, ' <=>');
            
                Prova1 = addReaction(Prova1, ConReactionNumber, ConReaction);
            
                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber, ConReaction);
                
%                 Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
%                 
%                 Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                
            end
            %%%%%%controllare questa parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%
            %%%%%%%%%%%MES
            if (strcmp(p,'MES'))
            
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2O(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' H2O(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval,'H2O(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
%                 
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
%                 
%                 end
%                 
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> O2(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> O2(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval,'O2(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Pi(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Pi(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval,'Pi(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
%                 
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
%                 
%                 end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' CO2(s) ->');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' CO2(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval,'CO2(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> NH3(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> NH3(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'NH3(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Glucose(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Glucose(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Glucose(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Lysine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Lysine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Lysine(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Sulfate(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Sulfate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Sulfate(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
%                     
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
%                 
%                 end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Methionine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Methionine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Methionine(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Tryptophan(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Tryptophan(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Tryptophan(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Phenylalanine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Phenylalanine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Phenylalanine(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urea(s) ->');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Urea(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Urea(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Choline(c)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Choline(c)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Choline(c)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 
%                 end
%                 
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Leucine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Leucine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Leucine(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Histidine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Histidine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Histidine(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                     
%                 end

                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Nicotinamide(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Nicotinamide(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
%                 if strcmp(ReadingEval, 'Nicotinamide(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                     
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Valine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Valine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Valine(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Threonine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Threonine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Threonine(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                     
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Ethanolamine(c)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Ethanolamine(c)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Ethanolamine(c)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 
%                 end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Arachidonate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Arachidonate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Arachidonate(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Riboflavin(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Riboflavin(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Riboflavin(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                                     
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pyridoxine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pyridoxine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pyridoxine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urate(s) ->');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Urate(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Urate(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Ubiquinone(m)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Ubiquinone(m)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Ubiquinone(m)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Isoleucine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Isoleucine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Isoleucine(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Folate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Folate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%                 if strcmp(ReadingEval, 'Folate(s)')
%             
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
%                 
%                     
%                 
%                 end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pantothenate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pantothenate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pantothenate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                 %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Linoleate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Linoleate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Linoleate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                
                end
                
                 %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Palmitolate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Palmitolate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Palmitolate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                
            end
            %%%%%%%%%%%%%%%%%% MES
            %%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%
            %%%%%%%%%%%HES
            if (strcmp(p,'HES'))
            
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2O(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' H2O(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'H2O(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Pi(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Pi(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'Pi(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' CO2(s) ->');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' CO2(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'CO2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> NH3(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> NH3(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'NH3(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Glucose(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Glucose(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Glucose(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Lysine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Lysine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Lysine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Sulfate(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Sulfate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Sulfate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Methionine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Methionine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Methionine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Tryptophan(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Tryptophan(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Tryptophan(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Phenylalanine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Phenylalanine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Phenylalanine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urea(s) ->');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Urea(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urea(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Choline(c)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Choline(c)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Choline(c)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Leucine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Leucine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Leucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Histidine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Histidine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Histidine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                end

                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Nicotinamide(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Nicotinamide(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                if strcmp(ReadingEval, 'Nicotinamide(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Valine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Valine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Valine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Threonine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Threonine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Threonine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Ethanolamine(c)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Ethanolamine(c)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Ethanolamine(c)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Arachidonate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Arachidonate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Arachidonate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Riboflavin(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Riboflavin(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Riboflavin(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                                    
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pyridoxine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pyridoxine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pyridoxine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urate(s) ->');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Urate(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Ubiquinone(m)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Ubiquinone(m)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Ubiquinone(m)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Isoleucine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Isoleucine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Isoleucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Folate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Folate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Folate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pantothenate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pantothenate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pantothenate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                 %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Linoleate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Linoleate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Linoleate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                
                end
                
                 %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Palmitolate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Palmitolate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Palmitolate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                
            end
            %%%%%%%%%%%%%%%%%% HES
            %%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%
            %%%%%%%%%%%DES
            if (strcmp(p,'DES'))
            
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2O(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' H2O(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'H2O(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' O2(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' O2(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'O2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
               
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Pi(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Pi(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'Pi(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' CO2(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' CO2(s) <=>');
                
                if strcmp(ReadingEval,'CO2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' NH3(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' NH3(s) <=>');

                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'NH3(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Glucose(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Glucose(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Glucose(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
               
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,'Lysine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,'Lysine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Lysine(s)')
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Sulfate(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Sulfate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Sulfate(s)')
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
               
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Methionine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Methionine(s) <=>');

                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Methionine(s)')
            
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Tryptophan(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Tryptophan(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Tryptophan(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Phenylalanine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Phenylalanine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Phenylalanine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urea(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Urea(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urea(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Choline(c) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Choline(c) <=>');

                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Choline(c)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,'Leucine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,'Leucine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Leucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Histidine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Histidine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Histidine(s)')
            
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end

                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,'Nicotinamide(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,'Nicotinamide(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Nicotinamide(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                    
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,'Valine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,'Valine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Valine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Cholesterol(b) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Cholesterol(b) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Cholesterol(b)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,'Threonine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,'Threonine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Threonine(s)')
            
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Ethanolamine(c) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Ethanolamine(c) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Ethanolamine(c)')
            
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Arachidonate(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Arachidonate(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Arachidonate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,'Riboflavin(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,'Riboflavin(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Riboflavin(s)')
            
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,'Pyridoxine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,'Pyridoxine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pyridoxine(s)')
            
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urate(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Urate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urate(s)')
            
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Ubiquinone(m) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Ubiquinone(m) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Ubiquinone(m)')
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Isoleucine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Isoleucine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Isoleucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                
               %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Folate(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Folate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Folate(s)')
            
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Oleate(c) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Oleate(c) ->');
                
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Oleate(c)')
            
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Pantothenate(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Pantothenate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pantothenate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                 %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Linoleate(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Linoleate(s) ->');
                
                
                if strcmp(ReadingEval, 'Linoleate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                    
                end
                
                 %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Palmitolate(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Palmitolate(s) ->');
                
                
                if strcmp(ReadingEval, 'Palmitolate(s)')
            
                   Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                   Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' ATP-energy(c) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' ATP-energy(c) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'ATP-energy(c)')
            
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' ATP-energy(m) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' ATP-energy(m) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'ATP-energy(m)')
            
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                
            end
            %%%%%%%%%%%%%%%%%% DES
            %%%%%%%%%%%%%%%%%%
        
        
            %%%%%%%%%%%%%%
            %%%%%%%%%%%WES
            if (strcmp(p,'WES'))
            
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2O(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' H2O(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               
                if strcmp(ReadingEval,'H2O(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Pi(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Pi(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval,'Pi(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' CO2(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' CO2(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval,'CO2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' NH3(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' NH3(s) ->');

                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval, 'NH3(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Sulfate(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Sulfate(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval, 'Sulfate(s)')
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urea(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Urea(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urea(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %%%%%%%%%%%%WES
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            end
            
            %%%%%%%%%%%%%%
            %%%%%%%%%%%MIMES
            if (strcmp(p,'MIMES'))
            
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2O(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' H2O(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval,'H2O(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> O2(s)');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' -> O2(s)');
                
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval,'O2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
               
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Pi(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Pi(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber, 1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber, -1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber, 1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber, 1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if strcmp(ReadingEval,'Pi(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' CO2(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' CO2(s) ->');
               
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval,'CO2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> NH3(s)');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' -> NH3(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                
                if strcmp(ReadingEval, 'NH3(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Glucose(s)');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' -> Glucose(s)');

                
                 %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval, 'Glucose(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Lysine(s)');
                
               ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' -> Lysine(s)');

               %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
               
                if strcmp(ReadingEval, 'Lysine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
               
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Sulfate(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Sulfate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Sulfate(s)')
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
               
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urea(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Urea(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urea(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Methionine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Methionine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Methionine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Tryptophan(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Tryptophan(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                if strcmp(ReadingEval, 'Tryptophan(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Phenylalanine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Phenylalanine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Phenylalanine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Choline(c) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Choline(c) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Choline(c)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                    
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Leucine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Leucine(s)');
                
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Leucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Histidine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Histidine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Histidine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                end

                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Nicotinamide(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Nicotinamide(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Nicotinamide(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Valine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Valine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Valine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Threonine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Threonine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Threonine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Ethanolamine(c) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Ethanolamine(c) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Ethanolamine(c)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Arachidonate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Arachidonate(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Arachidonate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Riboflavin(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Riboflavin(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Riboflavin(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                                    
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pyridoxine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pyridoxine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pyridoxine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urate(s) ->');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Urate(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Isoleucine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Isoleucine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Isoleucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Folate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Folate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Folate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pantothenate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pantothenate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pantothenate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                 %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Linoleate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Linoleate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                if strcmp(ReadingEval, 'Linoleate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Palmitolate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Palmitolate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Palmitolate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                 
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Cholesterol(b) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Cholesterol(b) -> ');
                
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Cholesterol(b)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2S(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' H2S(s) -> ');
                
                if strcmp(ReadingEval, 'H2S(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Fe2+(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Fe2+(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Fe2+(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                 
                
                end
                
                %%%%%%%%%%%%MIMES
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            end
            
            
            %%%%%%%%%%%%%%
            %%%%%%%%%%%MIPES
            if (strcmp(p,'MIPES'))
            
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2O(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' H2O(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'H2O(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> O2(s)');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' -> O2(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'O2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
               
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Pi(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Pi(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'Pi(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' CO2(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' CO2(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'CO2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> NH3(s)');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' -> NH3(s)');

                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'NH3(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Glucose(s)');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' -> Glucose(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval, 'Glucose(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Lysine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Lysine(s) <=>');

                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Lysine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
               
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Sulfate(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Sulfate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Sulfate(s)')
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
               
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urea(s) ->');
                
                 ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Urea(s) ->');
                 
                 %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urea(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,'Methionine(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Methionine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Methionine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Tryptophan(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Tryptophan(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Tryptophan(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Phenylalanine(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Phenylalanine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Phenylalanine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Choline(c) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Choline(c) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Choline(c)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                    
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Leucine(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Leucine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Leucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Histidine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Histidine(s) <=> ');
                
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Histidine(s)')
            
                   Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                   Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                    
                end

                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Nicotinamide(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Nicotinamide(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Nicotinamide(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Valine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Valine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Valine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Threonine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Threonine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Threonine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Ethanolamine(c) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Ethanolamine(c) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Ethanolamine(c)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Arachidonate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,'Arachidonate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Arachidonate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Riboflavin(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Riboflavin(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Riboflavin(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                                    
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pyridoxine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pyridoxine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pyridoxine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urate(s) ->');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Urate(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Isoleucine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Isoleucine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Isoleucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Folate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Folate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Folate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pantothenate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pantothenate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pantothenate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Linoleate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Linoleate(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Linoleate(s)')
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Palmitolate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Palmitolate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Palmitolate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Cholesterol(b) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Cholesterol(b) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Cholesterol(b)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2S(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' H2S(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'H2S(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Fe2+(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Fe2+(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Fe2+(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Glutamate(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Glutamate(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Glutamate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Glycine(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Glycine(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Glycine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Alanine(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Alanine(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Alanine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Aspartate(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Aspartate(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Aspartate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Arginine(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Arginine(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Arginine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Glutamine(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Glutamine(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Glutamine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Serine(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Serine(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Serine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Tyrosine(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Tyrosine(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Tyrosine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Cysteine(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Cysteine(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Cysteine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Proline(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Proline(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Proline(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Asparagine(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Asparagine(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Asparagine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' L-Lactate(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' L-Lactate(s) -> ');
                
                if strcmp(ReadingEval, 'L-Lactate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Cystine(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Cystine(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Cystine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %%%%%%%%%%%%MIPES
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            end
            
            %%%%%%%%%%%%%%
            %%%%%%%%%%%PIPES
            if (strcmp(p,'PIPES'))
            
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2O(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' H2O(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'H2O(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> O2(s)');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' -> O2(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'O2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
               
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Pi(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Pi(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'Pi(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' CO2(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' CO2(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval,'CO2(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> NH3(s)');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' -> NH3(s)');
                 
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
                if strcmp(ReadingEval, 'NH3(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Glucose(s) <=> ');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Glucose(s) <=> ');

                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Glucose(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Lysine(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Lysine(s) <=> ');

                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Lysine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
               
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Sulfate(s) <=>');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Sulfate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Sulfate(s)')
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
               
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urea(s) ->');
                
                ProvaNoEval = addReaction(ProvaNoEval,ConReactionNumber,' Urea(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urea(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,'Methionine(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Methionine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Methionine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Tryptophan(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Tryptophan(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Tryptophan(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Phenylalanine(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Phenylalanine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Phenylalanine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Choline(c) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Choline(c) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                if strcmp(ReadingEval, 'Choline(c)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                    
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Leucine(s) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Leucine(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Leucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Histidine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Histidine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Histidine(s)')
            
                   Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                   Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                    
                end

                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Nicotinamide(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Nicotinamide(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Nicotinamide(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Valine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Valine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                if strcmp(ReadingEval, 'Valine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Threonine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Threonine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Threonine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber , 1000*EvalPerc,'u');
                
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Ethanolamine(c) <=>');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Ethanolamine(c) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Ethanolamine(c)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                
                end
                
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Arachidonate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,'Arachidonate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Arachidonate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Riboflavin(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Riboflavin(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Riboflavin(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                                    
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pyridoxine(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pyridoxine(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pyridoxine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Urate(s) ->');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Urate(s) ->');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Urate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Isoleucine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Isoleucine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Isoleucine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Folate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Folate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Folate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                    
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' -> Pantothenate(s)');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' -> Pantothenate(s)');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Pantothenate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                
                end
                
                 %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Linoleate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Linoleate(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Linoleate(s)')
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Palmitolate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Palmitolate(s) <=>');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Palmitolate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                               
                end
                
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' H2S(s) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' H2S(s) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'H2S(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Fe2+(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Fe2+(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Fe2+(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                               
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Glutamate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Glutamate(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Glutamate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                    
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Glycine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Glycine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Glycine(s)')
            
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Alanine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Alanine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Alanine(s)')
            
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Aspartate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Aspartate(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Aspartate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Arginine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Arginine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Arginine(s)')
            
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Glutamine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Glutamine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Glutamine(s)')
            
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Serine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Serine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Serine(s)')
            
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Tyrosine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Tyrosine(s) <=> ');
                
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Tyrosine(s)')
            
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Cysteine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Cysteine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Cysteine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Proline(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Proline(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                if strcmp(ReadingEval, 'Proline(s)')
            
                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Asparagine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Asparagine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Asparagine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' L-Lactate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' L-Lactate(s) <=> ');
                
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'L-Lactate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');

                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Cystine(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Cystine(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Cystine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Cholesterol(b) -> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' CHolesterol(b) -> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Cystine(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'b');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Palmitate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Palmitate(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Palmitate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Oleate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Oleate(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if strcmp(ReadingEval, 'Oleate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                %
                ReactionNumber = ReactionNumber+1;
                
                StrReactionNumber = num2str(ReactionNumber);
            
                ConReactionNumber = strcat('r',StrReactionNumber);
                
                Prova1 = addReaction(Prova1,ConReactionNumber,' Stearate(s) <=> ');

                ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Stearate(s) <=> ');
                
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
                
                ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
                %%%%%%%%%%%%%%%%%nuona parte%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                if strcmp(ReadingEval, 'Stearate(s)')
            
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
                    
                    Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
                
                end
                
                
                %%%%%%%%%%%%PIPES
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            end
           
           
           
       end
   
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%da qui in su abbiamo cancellato%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%inserimento reazioni di trasporto dei metaboliti
        %%%%%%%%%%%%%%%%%nello spazio sinusoidale%%%%%%%%%%%%%%%%%%%%%%%%%
%         for i = 1:length(Hepatonet1.metNames)
%             
%             Metab = Hepatonet1.metNames{i};
%             
%             if strcmp(Metab(end-1) , 's');
%                 
%                 ReactionNumber = ReactionNumber+1;
%                 
%                 StrReactionNumber = num2str(ReactionNumber);
%             
%                 ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%                 ConReaction = strcat(Metab, ' <=>');
%                 
%                 [Prova1,rxnIDEexists] = addReaction(Prova1,ConReactionNumber, ConReaction);
% 
%                 if isempty(rxnIDEexists) == 0
%                     
%                     ReactionNumber = ReactionNumber - 1;
%         
%                 else
%                 
%                 
%                     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber, ConReaction);
%                 
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
%                 
%                     Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
%                 
%                     ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
%                 
%                     ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
%                     
%                 end
%                 
%                 
%             end
%             
%         end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%inserimento reazioni di trasporto dei metaboliti
            %%%%%%%%%%%%%%%%%nello spazio sinusoidale%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%inserimento reazioni di trasporto conosciute%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         %reading the text file
%         [NumTestBase2,TxtTestBase2,RawTestBase2] = xlsread('HepaTestsProd.xls');
% 
%         
%         for i = 1:length(RawTestBase2)
%             
%             %reading evaluator. This is evaluator also in xls file
%             Metab = RawTestBase2{i,2};
%             
%             if (strcmp(Metab(end-1),'s') == 0)
%         
%             ReactionNumber = ReactionNumber+1;
%                 
%             StrReactionNumber = num2str(ReactionNumber);
%             
%             ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%             ConReaction = strcat(Metab, ' <=>');
%                 
%             [Prova1,rxnIDEexists] = addReaction(Prova1,ConReactionNumber, ConReaction);
% 
%             if isempty(rxnIDEexists) == 0
%                     
%                 ReactionNumber = ReactionNumber - 1;
%         
%             else
%                 
%                 
%                 ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber, ConReaction);
%                 
%                 Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,1000*EvalPerc,'u');
%                 
%                 Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,-1000*EvalPerc,'l');
%                 
%                 ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,1000*EvalPerc,'u');
%                 
%                 ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,-1000*EvalPerc,'l');
%                     
%             end
%             
%             end
%             
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%inserimento reazioni di trasporto dei metaboliti
%         %%%%%%%%%%%%%%%%%nello spazio sinusoidale%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
                    
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%prove inserimento reazioni%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
    %massimizzazione della reazione: Glycine(p) + H2O(p) + O2(p) -> Glyoxylate(p) + H202(p) + NH4+(p)
    
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' Glyoxylate(p) + O2(p) ->  H2O2(p) + Oxalate(p) ');
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber, 'Glyoxylate(p) + O2(p) ->  H2O2(p) + Oxalate(p)');
%     
%     if isempty(rxnIDEexists)
%         
%         %Prova1.c(end) = -1;
%         
%         %ProvaNoEval.c(end) = -1;
%         
%         prodreac = length(ProvaNoEval.c);
%         
%         
%     else
%         
%         ConReactionNumber = Prova1.rxns{rxnIDEexists};
%     
%         %Prova1.c(rxnIDEexists) = -1;
%     
%         %ProvaNoEval.c(rxnIDEexists) = -1;
%         
%         prodreac = rxnIDEexists;
%         
%         ReactionNumber = ReactionNumber-1;
%     
%     end

prodreac = 5084;

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end nuova prova
     
     
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' NH3(p) <=> NH3(c)');
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' NH3(p) <=> NH3(c) ');
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' H+(PG)(p) <=> H+(PG)(c) '); 
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' H+(PG)(p) <=> H+(PG)(c) '); 
%     
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end nuova prova
%      
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' H+(PG)(c) -> '); 
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' H+(PG)(c) ->  '); 
%     
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end nuova prova
%      
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' Alanine(c) <=> Alanine(p) ');
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Alanine(c) <=> Alanine(p) ');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' Pyruvate(p) <=> Pyruvate(c) ');
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Pyruvate(p) <=> Pyruvate(c) ');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end nuova prova
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' Serine(p) <=> Serine(c) ');
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,' Serine(p) <=> Serine(c) ');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end nuova prova
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' NADH(c) <=> ');
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,'  NADH(c) <=> ');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end nuova prova

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' NADP+(c) <=> ');
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,'  NADP+(c) <=> ');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end nuova prova

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber,' NAD+(c) <=> ');
%     
%     ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber,'  NAD+(c) <=> ');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end nuova prova


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber ,' -> Dgcholcoa(p) ');
%     
%     [ProvaNoEval] = addReaction(ProvaNoEval,ConReactionNumber ,' -> Dgcholcoa(p) ');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nuova prova
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber ,' -> Choloyl-CoA(p) ');
%     
%     [ProvaNoEval] = addReaction(ProvaNoEval,ConReactionNumber ,' -> Choloyl-CoA(p) ');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber ,' Glycochenodeoxycholate(p) <=> Glycochenodeoxycholate(c)');
%     
%     [ProvaNoEval] = addReaction(ProvaNoEval,ConReactionNumber ,'  Glycochenodeoxycholate(p) <=> Glycochenodeoxycholate(c) ');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber ,' Glycocholate(p) <=> Glycocholate(c)');
%     
%     [ProvaNoEval] = addReaction(ProvaNoEval,ConReactionNumber ,'  Glycocholate(p) <=> Glycocholate(c) ');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber ,' Gdchola(p) <=> Gdchola(c)');
%     
%     [ProvaNoEval] = addReaction(ProvaNoEval,ConReactionNumber ,'  Gdchola(p) <=> Gdchola(c) ');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber ,' CoA(p) <=> CoA(c)');
%     
%     [ProvaNoEval] = addReaction(ProvaNoEval,ConReactionNumber ,'  CoA(p) <=> CoA(c) ');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber ,' Sarcosine(p) <=> Sarcosine(c)');
%     
%     [ProvaNoEval] = addReaction(ProvaNoEval,ConReactionNumber ,' Sarcosine(p) <=> Sarcosine(c) ');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ReactionNumber = ReactionNumber+1;
%                 
%     StrReactionNumber = num2str(ReactionNumber);
%             
%     ConReactionNumber = strcat('r',StrReactionNumber);
%                 
%     [Prova1, rxnIDEexists] = addReaction(Prova1,ConReactionNumber ,' Formaldehyde(p) <=> Formaldehyde(c)');
%     
%     [ProvaNoEval] = addReaction(ProvaNoEval,ConReactionNumber ,' Formaldehyde(p) <=> Formaldehyde(c) ');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     
    %Bounds do not change   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if BoundsChange == 0
        
        for z = 1:length(Prova1.lb)
    
            Prova1.lb(z) = Prova1.lb(z)*EvalPerc/EvalPerc;
    
            Prova1.ub(z) = Prova1.ub(z)*EvalPerc/EvalPerc;
            
            ProvaNoEval.lb(z) = ProvaNoEval.lb(z)*EvalPerc/EvalPerc;
    
            ProvaNoEval.ub(z) = ProvaNoEval.ub(z)*EvalPerc/EvalPerc;
            
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %constraints target reaction flux. In xls file this is the column
    %objective
    ReadingConstraint = RawTestBase(i,2);
    
    %ReadingConstraint = ReadingConstraint{1};
    
    StrReactionNumber = num2str(ReactionNumber+1);
   
    ConReactionNumber = strcat('r',StrReactionNumber);
    
    if NumList(i-2) == 0 
    
        %ConReaction = strcat(ReadingConstraint, ' -> ');
        
        
        [Prova1,rxnIDExists] = addReaction(Prova1, ConReactionNumber,ReadingConstraint, [-1],'false');
        
        Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,500*EvalPerc,'b');
        
        %%%%%%%%%%%%%%%%%%%%%%%new 26/11/2015%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MetCol = Prova1.S(:,end);
        
        MetRow = find(abs(MetCol) == 1);
        
        [nP,mP] = size(Prova1.S);
        
        ii = 1;
        
        found = 0;
        
        foundOut = [];
        
        while found == 0 & ii < (mP-1)
            
            
            if sum(abs(Prova1.S(:,ii))) == 1 && abs(Prova1.S(MetRow,ii)) == 1
                
                Prova1 = removeRxns(Prova1,Prova1.rxnNames(ii));
                
                found = 1;
                
                foundOut = ii;
                
            end
            
            ii = ii + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            
        
    else
        
        %ConReaction = strcat(' -> ', ReadingConstraint);
        
        [Prova1,rxnIDExists] = addReaction(Prova1, ConReactionNumber,ReadingConstraint, [1],'false');
        
        Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,500*EvalPerc,'b');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MetCol = Prova1.S(:,end)
        
        MetRow = find(abs(MetCol) == 1)
        
        [nP,mP] = size(Prova1.S);
        
        ii = 1;
        
        found = 0;
        
        foundOut = [];
        
        while found == 0 & ii < (mP-1)
            
            
            if sum(abs(Prova1.S(:,ii))) == 1 && abs(Prova1.S(MetRow,ii)) == 1
                
                Prova1 = removeRxns(Prova1,Prova1.rxnNames(ii));
                
                found = 1;
                
            end
            
            ii = ii + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
            
    %[Prova1,rxnIDExists] = addReaction(Prova1, ConReactionNumber, ConReaction);
    
%     if isempty(rxnIDExists) == 0
%         
%         ConReactionNumber = strcat('r',rxnIDExists);
%         
%         Prova1 = changeRxnBounds(Prova1, Prova1.rxnNames{rxnIDExists} ,500*EvalPerc,'b');
%     
%         ProvaNoEval = changeRxnBounds(ProvaNoEval, Prova1.rxnNames{rxnIDExists} ,500*EvalPerc,'b');
%         
%     else
%         
%         
%         ProvaNoEval = addReaction(ProvaNoEval, ConReactionNumber, ConReaction);
%     
%         Prova1 = changeRxnBounds(Prova1, ConReactionNumber ,500*EvalPerc,'b');
%     
%         ProvaNoEval = changeRxnBounds(ProvaNoEval, ConReactionNumber ,500*EvalPerc,'b');
%         
%     end
%     
    %Counting the number of reactions inserted by this function
    
    ReactionNumber = ReactionNumber+1;
    
    [nHepa1,mHepa1] = size(Hepatonet1.S);
    
    [nProva1,mProva1] = size(Prova1.S);
    
    number_inserted_reactions = mProva1 - mHepa1;
    
    %Chesi bounds%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prova1.lb(16) = 0;
     %Prova1.ub(16) = 0;
%     
     %Prova1.lb(2555) = 0;
     %Prova1.ub(2555) = 0;
%     
%     Prova1.lb(17) = 1500;
%     Prova1.ub(17) = 2000;
%     
%     Prova1.lb(2556) = 1500;
%     Prova1.ub(2556) = 2000;
%     
%     Prova1.lb(1452) = 1500;
%     Prova1.ub(1452) = 2000;
%     
%     Prova1.lb(3991) = 1500;
%     Prova1.ub(3991) = 2000;
  
    %Chesi bounds%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %Gaucher disease type I bounds%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Prova1 = changeRxnBounds(Prova1,'r1406',111,'b');
%     Prova1 = changeRxnBounds(Prova1,'r3945',111,'b');
    
    
%     Prova1.lb(1406) = 111;
%     Prova1.ub(1406) = 111;
%      
%     Prova1.lb(3945) = 111;
%     Prova1.ub(3945) = 111;
    %Gaucher disease type I bounds%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %GSD1A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %Prova1.lb(396) = 0;
     %Prova1.ub(396) = 0;
     
     %Prova1.lb(2935) = 0;
     %Prova1.ub(2935) = 0;
    %Gaucher disease type I bounds%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%