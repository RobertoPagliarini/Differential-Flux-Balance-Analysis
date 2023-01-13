function [Indexes1, Indexes2] = linkingEnzymeReactions(model,Enzyme,dGo,InsertedReaction);
%This function links enzymes to reactions for simulations by considering
%the most probable according with Gibb's free energy in the standard
%conditions

Indexes1 = [];

Indexes2 = [];

[n,m] = size(model.rxnGeneMat);

for i = 1:n
    
    if (i >= 1 & i <= 2539) 
    
        if strcmp(model.rxnGeneMat{i,1},Enzyme)
        
            if dGo(i) <= 0
        
                Indexes1(end+1,1) = i;
                
                Indexes2(end+1,1) = 2539 + i;
                
            else
                
                Indexes1(end+1,1) = 2539 + i;
                
                Indexes2(end+1,1) =  i;
                
            end
            
        end
        
    end
       
    if (i > 2539*2 & i <= 2539*2+InsertedReaction)
        
        if strcmp(model.rxnGeneMat{i,1},Enzyme)
            
            Indexes1(end+1,1) = i;
            
        end
        
    end
    
end
