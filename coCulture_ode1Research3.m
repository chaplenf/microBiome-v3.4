function dy = coCulture_ode1Research3(t, y, ~)
global community bucket HenrysLawCoefficients kLa Temp
%ODE framework
    
dy = zeros(size(y));
    numSpecies = length(community.species);
    numLiquidSpecies = length(community.Lfeed);
    numLiquidNMSpecies = length(community.LNMfeed);
    numGasSpecies = length(community.Gfeed);
    numBioSpecies= length(community.mets);
    numBBuckets=length(bucket.numbbucket);
    numReac= length(bucket.flux);

V = y(1);                     %Bioreactor liquid volume [L]
Vg = y(2);                    %Bioreactor gas volume [L]
for j = 1:numLiquidSpecies
    L(j) = y(2+numSpecies+j);  %Substrates [mmol]
end
for k = 1:numLiquidNMSpecies
    LNM(k)=y(2+numSpecies+numLiquidSpecies+k);
end
% assigning growth rates and metabolic production/consumption rates
    % here, the rates are calculated using FBA through linear programming models 
    % run in the GAMS environment

    
vs =zeros(numSpecies,((numReac+numBBuckets)/numSpecies));
mu = zeros(1, numSpecies);

reactionIdentifiers = [community.mets bucket.reac];
subFlux = [L LNM];

uptakeRates131;  % updating FBA organism uptake constraints

% calculating the growth rates and production/consumption rates of the organisms

for i = 1:numSpecies
try
    switch community.species{1,i}.description
        case 'Nitrosomonas europaea'
               rxnID={'EX_cpd00013_e' 'EX_cpd00013y_e'};
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00418x_c',0,'u');
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00418_e',0,'u');                
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00418_e',0,'l');                 
               community.species{1,i} = changeObjective(community.species{1,i},rxnID,-1);                 
    [vNeu]=optimizeCbModel(community.species{1,i},'max',1e-5);
            if (vNeu.stat ~= 1)    
                mu(i) = 0;
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416_e',0,'l');  
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416_e',0,'u');             
                disp('deathphase');
            else
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416_e',0.0,'l');  
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416_e',0.031,'u');   
                [~,reacIDModel]=ismember(reactionIdentifiers,community.species{1,i}.rxns);
                [~,loc]=ismember('EX_cpd11416_e',community.species{1,i}.rxns); 
                vs(i,:)=vNeu.v(reacIDModel);
                mu(i)=vNeu.v(loc);
            end
       case 'Nitrobacter winogradskyi'            
           [vNwi]=optimizeCbModel(community.species{1,i},'max',1e-5);
            if (vNwi.stat ~= 1)               
                mu(i) = 0;
              community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416_e',0,'l');  
              community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416_e',0,'u');             
                disp('deathphase');
            else
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416_e',0.0,'l');  
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416_e',0.031,'u'); 
                [~,reacIDModel]=ismember(reactionIdentifiers,community.species{1,i}.rxns);
                [~,loc]=ismember('EX_cpd11416_e',community.species{1,i}.rxns);
                vs(i,:)=vNwi.v(reacIDModel);
                mu(i)=vNwi.v(loc);
            end
    end
catch
    if (community.species{1,i}.ub < community.species{1,i}.lb)
    for j=1:length(community.species{1,j}.ub)
        if (community.species{1,i}.ub(j) < community.species{1,i}.lb(j))
        community.species{1,i} = changeRxnBounds(community.species{1,i},community.species.rxns(j),0,'b');    
        end
    end    
    end
   disp('catch');    
end    
HONOModelODEResearch3;  %updating abiotic and biotic reaction sets
t, vs(:,:)
end
