function [t, V, Vg, X, L, LNM, G, AB, BB, MT, R] = Co3Research3(tspan, initialConditions)
global community bucket HenrysLawCoefficients kLa Temp
%   tspan: simulation time
%   intialConditions: the initial reactor volume, biomass and metabolite concentrations

    numSpecies = length(community.species);
    numLiquidSpecies = length(bucket.liquidreactormassbalance);
    numLNM = length(bucket.NMliquidreactormassbalance);
    numGasSpecies = length(bucket.gasreactormassbalance);
    numABucket = length(bucket.numabucket);
    numBBucket = length(bucket.numbbucket);
    numMTbucket = length(bucket.masstransfer);
    numRbucket = length(bucket.flux);
    
sparseJacobian = ones(151);

nonNeg=linspace(1,151,151);

        DO=@coCulture_ode3Research3;
        options=odeset('JPattern',sparseJacobian,'NonNegative',nonNeg,'RelTol',1e-7,'AbsTol',1e-5);
%        options=odeset('Mass',M,'MassSingular','yes','RelTol',1e-10,'AbsTol',1e-8);
        [t,y] = ode15s(DO,tspan,initialConditions,options);

V = y(:,1);
Vg = y(:,2);
X = y(:,3:2+numSpecies);
L = y(:,2+numSpecies+1: 2+numSpecies+numLiquidSpecies);
LNM = y(:,2+numSpecies+numLiquidSpecies+1: 2+numSpecies+numLiquidSpecies+numLNM);
G = y(:,2+numSpecies+numLiquidSpecies+numLNM+1:2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies);
AB = y(:,2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+1:2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket);
BB= y(:,2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+1:2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+numBBucket);
MT= y(:,2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+numBBucket+1:2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+numBBucket+numMTbucket);
R= y(:,2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+numBBucket+numMTbucket+1:2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+numBBucket+numMTbucket+numRbucket);
end
