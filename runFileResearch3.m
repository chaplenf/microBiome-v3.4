function     [ttot,Vtot,Vgtot,Xtot,Ltot,LNMtot,Gtot, ABtot, BBtot, MTtot, Rtot]= runFileResearch3(pH, NO2, NH3, CO3, PO4, O2, OD)
global community bucket HenrysLawCoefficients kLa Temp
% set the variable "model" to a COBRA Toolbox model      

Neuro=readCbModel('Neuro6_e1.xml');
     model1 = Neuro;
     model1.description = 'Nitrosomonas europaea';
     Nwino=readCbModel('Nwino3G.xml');  
     model2 = Nwino;
     model2.description = 'Nitrobacter winogradskyi';

%Advise not modifying anything beyond here for the moment.

% Initiatize variable community
%   community:  description of the microbial community and its environment
%       community.species:      species in the community
%       community.mets:         identifiers for substrates and products in bioreactor
%       community.Xfeed:        biomass in the feed (gDCW/L)
%       community.Lfeed:        metabolites in the feed (mmol/L)
%       community.LNMfeed:      metabolites in the feed (mmol/L)
%       community.Gfeed         partial pressure of species in the gas feed (atm) 
%       community.Fin:          liquid flowrate in (L/hr)
%       community.Fout:         liquid flowrate out (L/hr)
%       community.Vin:          gas flowrate in (L/hr)
%       community.Vout:         gas flowrate out (L/hr)

    community.species = {model1 model2};
    community.mets = {'EX_cpd00075_e','EX_cpd00209_e','EX_cpd00418_e','EX_cpd00007_e','EX_cpd00013y_e','EX_cpd00047_e','EX_cpd00659_e','EX_cpd00528_e','EX_cpd00011_e','EX_cpd00242_e','EX_cpd00165_e','EX_cpd00418x_c','EX_cpd00067_e','EX_cpd00009_e','EX_cpd00013_e'};
    community.Xfeed= [0 0];                                 
    community.Lfeed = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];            
    community.LNMfeed = [0 0 0 0 0 0 0 0 0 0];                                  
    community.Gfeed = [0 0 0 0.21 0 0 0.78 0.0004 ];         
    community.Fin = 0;                                      
    community.Fout = 0;                                     
    community.Vgin = 0;                                     
    community.Vgout = 0;
    community.B = [0,0,0,0];

    %Initialize variable bucket    
%   bucket: cumulative fluxes from selected species and internal cell reactions 
%       bucket.reac                         identifiers for internal reaction fluxes being monitored (mmol/L)
%       bucket.liquidreactormassbalance     metabolite species in liquid phase (mmol/L)
%       bucket.NMliquidreactormassbalance   non-metabolite species in liquid phase (mmol/L)
%       bucket.gasreactormassbalance        species in gas phase (mmol/L)
%       bucket.masstransfer                 species being transferred between phases (mmol/L)
%       bucket.numbbucket                   biotic buckets (cumulative flux) for selected species (mmol/L)
%       bucket.numabucket                   abiotic buckets for selected species (mmol/L)
%       bucket.flux                         internal reaction flux buckets for selected reactions (mmol/L)

    bucket.reac = {'rxn14202a','rxn14202b','rxn14202c','rxn14202d','rxn14202e','rxn14202f','rxn14202g','rxn01806mod','rxn05795mod','rxnP460a','rxnP460b','rxn10113','rxn10125mod','rxn08975','rxn10122mod','rxn10043mod','rxn00058','rxn10042','rxn12750','rxn14173a','rxn14173b','rxn00567a','rxn00567b','rxn00567c','rxn09008','rxn12750mod','rxn14010','rxn10811','Cytochrome_Exchange1','Cytochrome_Exchange2','Cytochrome_Exchange3'};
    bucket.liquidreactormassbalance = zeros(1,15);
    bucket.NMliquidreactormassbalance =zeros(1,10);
    bucket.gasreactormassbalance = zeros(1,8);
    bucket.masstransfer = zeros(1,8);
    bucket.numbbucket = zeros(1,30);
    bucket.numabucket = zeros(1,14);
    bucket.flux=zeros(1,62);
    
%    pH=7.8;
%Henderson-Hasselbalch Equation for impact of initial pH on split between
%acid and base in buffers.

%PO4 ... pKa=7.21  pH=pK' + log {[Proton Acceptor]/[Proton Donor]} ... solve
%for proton acceptor and donor with [Proton Acceptor]+[Proton Donor] = PO4

PO4Ratio = 10^(pH-12.32);
HPO4Ratio = 10^(pH-6.86);
H2PO4Ratio = 10^(pH-2.21);

%CO3 ... pKa=6.33  pH=pK' + log {[Proton Acceptor]/[Proton Donor]} ... solve
%for proton acceptor and donor with [Proton Acceptor]+[Proton Donor] = CO3

H2CO3Ratio = 10^(pH-3.33);
HCO3Ratio = 10^(pH-10.32);

%Henry's Law Coefficient Temperature Variation

%CO3 ... pKa=10.33  pH=pK' + log {[Proton Acceptor]/[Proton Donor]} ... solve
%for proton acceptor and donor with [Proton Acceptor]+[Proton Donor] = CO3

NH3Ratio = 10^(pH-9.24);

%Henry's Law Coefficient Temperature Variation

HenrysLawParameters=[298.15 298.15 298.15 298.15 298.15 298.15 298.15 298.15;                                                                   %Tref
                    (273.15+Temp) (273.15+Temp) (273.15+Temp) (273.15+Temp) (273.15+Temp) (273.15+Temp) (273.15+Temp) (273.15+Temp);            %Texpt
                    0.0012667 53.25 0.0014252 0.036 0.02425 0.00063 0.023 49;                                                                   %Constant 1
                    1650 3708.3 2100 2393 2675 1300 2150 4820];                                                                                 %Constant 2

HenrysLawCoefficients=exp(log(HenrysLawParameters(3,:))-((1./HenrysLawParameters(1,:)-1./HenrysLawParameters(2,:)).*HenrysLawParameters(4,:)));

%Matrix of parameters for adjusting mass transfer coefficients with
%temperature

massTransferParameters=[0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05;                                                        %Vessel Diameter (d; m)
                        3.333 3.333 3.333 3.333 3.333 3.333 3.333 3.333;                                                %Rotational Rate (n; s^-1)
                        0.03 0.03 0.03 0.03 0.03 0.03 0.03 0.03;                                                        %Diameter of Rotation (do; m)
                        0.000005 0.000005 0.000005 0.000005 0.000005 0.000005 0.000005 0.000005;                        %Liquid Volume (Vl; m^3)
                        2.1e-9 2.1e-9 2.1e-9 2.1e-9 2.1e-9 2.1e-9 2.1e-9 2.1e-9;                                        %Liquid Diffusivity (D; m^2s^-1)
                        8e-7 8e-7 8e-7 8e-7 8e-7 8e-7 8e-7 8e-7;                                                        %Kinematic Viscosity (v; m^2s^-1)
                        9.81 9.81 9.81 9.81 9.81 9.81 9.81 9.81;                                                        %Acceleration due to Gravity (g; m^1s^-2)
                        0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01;                                                        %Osmolarity (Osm)
                        3.72e-7 3.72e-7 3.72e-7 3.72e-7 3.72e-7 3.72e-7 3.72e-7 3.72e-7];                               %Constant

kLa=3600*(massTransferParameters(9,:)./HenrysLawCoefficients).*massTransferParameters(8,:).^0.05.*massTransferParameters(2,:).^(1.18-(massTransferParameters(8,:)/10.1)).* ...
                                  massTransferParameters(4,:).^-0.74.*massTransferParameters(3,:).^0.33.*massTransferParameters(1,:).^1.88;                           
                              
%Set simulation time and call solver
%Segment 1

%New Model Switch
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567a',1.3,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567a',-1.3,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',0,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',0,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',0,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',0,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',8,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',8,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Transporter',10000,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'EX_cpd11416_e',0.031,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'EX_cpd11416_e',0.031,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn05795mod',47.88,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn05795mod',-47.88,'l');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn01806mod',0.6,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460a',1.2,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460b',0.6,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202g',0,'l');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202g',0.0,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202b',0,'l');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202b',0.04788,'u');    

    %Initial conditions for system variables

    initV = 0.005;                                                                                                                                                          %L
    initVg = 0.155;                                                                                                                                                         %L
    initX = [OD(1,1)*0.16*initV OD(1,2)*0.251*initV];                                                                                                                       %gDCW%
    initL = [NO2 0 0 O2 0.023 0 0 0 0 CO3 0 0 (initV*1000*10^-pH) 0.065 0.002];                                                                      %mmol
    initLNM = [0 0 0 0.00375 0 (initV*1000*10^-(14-pH)) 0 0 initV*0.1 0 ]; 
    initG=[0 0 0 1.308 0 0 4.86 0.0025];%mmol
    initAB=zeros(1,14);
    initBB=zeros(1,30);
    initMT=zeros(1,8);
    initR=zeros(1,62);

initialConditions = [initV initVg initX initL initLNM initG initAB initBB initMT initR];   
    tspan = [0:0.01:0.5];
  
    [t, V, Vg, X, L, LNM, G, AB, BB, MT, R] = Co1Research3(tspan,initialConditions);    

ttot=[t];
Vtot=[V];
Vgtot=[Vg];
Xtot=[X];
Ltot=[L];
LNMtot = [LNM];
Gtot=[G];
ABtot=[AB];
BBtot = [BB];
MTtot=[MT];
Rtot=[R];

lastV = size(V);
lastVg = size(Vg);
lastX = size(X);
lastL = size(L);
lastLNM = size(LNM);
lastG = size(G);
lastAB = size(AB);
lastBB = size(BB);
lastMT = size(MT);
lastR = size(R);

initV1 = V(lastV(1),1); %L
initVg1 = Vg(lastVg(1),1); %L
initX1 = X(lastX(1),:); %gDCW
initL1 = L(lastL(1),:);   %mmol
initLNM1 = LNM(lastLNM(1),:); %mmol
initG1 = G(lastG(1),:);%mmol
initAB1 = AB(lastAB(1),:);   %mmol
initBB1 = BB(lastBB(1),:); %mmol
initMT1 = MT(lastMT(1),:);%mmol
initR1 = R(lastR(1),:);

initialConditions = [initV1 initVg1 initX1 initL1 initLNM1 initG1 initAB1 initBB1 initMT1 initR1];

%Set simulation time and call solver
%Segment 2

%New Model Switch

    community.species{1,1} = changeRxnBounds(community.species{1,1},'EX_cpd11416_e',0.031,'u');    
    community.species{1,2} = changeRxnBounds(community.species{1,2},'EX_cpd11416_e',0.031,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn05795mod',47.88,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn01806mod',0.6,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460a',1.2,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460b',0.6,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202g',0,'l');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202g',0.0,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202b',0,'l');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202b',0.04788,'u');    

    tspan = [0.5:0.01:1.];
  
    [t, V, Vg, X, L, LNM, G, AB, BB, MT, R] = Co2Research3(tspan,initialConditions);    

ttot=[ttot; t];
Vtot=[Vtot; V];
Vgtot=[Vgtot; Vg];
Xtot=[Xtot; X];
Ltot=[Ltot; L];
LNMtot = [LNMtot;LNM];
Gtot=[Gtot; G];
ABtot=[ABtot; AB];
BBtot = [BBtot;BB];
MTtot=[MTtot; MT];
Rtot=[Rtot; R];

lastV = size(V);
lastVg = size(Vg);
lastX = size(X);
lastL = size(L);
lastLNM = size(LNM);
lastG = size(G);
lastAB = size(AB);
lastBB = size(BB);
lastMT = size(MT);
lastR = size(R);

initV2 = V(lastV(1),1); %L
initVg2 = Vg(lastVg(1),1); %L
initX2 = X(lastX(1),:); %gDCW
initL2 = L(lastL(1),:);   %mmol
initLNM2 = LNM(lastLNM(1),:); %mmol
initG2 = G(lastG(1),:);%mmol
initAB2 = AB(lastAB(1),:);   %mmol
initBB2 = BB(lastBB(1),:); %mmol
initMT2 = MT(lastMT(1),:);%mmol
initR2 = R(lastR(1),:);

initialConditions = [initV2 initVg2 initX2 initL2 initLNM2 initG2 initAB2 initBB2 initMT2 initR2];

%Segment 2
%Poughon Model Switch
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567a',1.3,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567a',-1.3,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',0,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',0,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',0,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',0,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',8,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',8,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Transporter',10000,'u');
tspan = [1:0.01:3];       

   [t, V, Vg, X, L, LNM, G, AB, BB, MT, R] = Co3Research3(tspan,initialConditions);    

ttot=[ttot; t];
Vtot=[Vtot; V];
Vgtot=[Vgtot; Vg];
Xtot=[Xtot; X];
Ltot=[Ltot; L];
LNMtot = [LNMtot;LNM];
Gtot=[Gtot; G];
ABtot=[ABtot; AB];
BBtot = [BBtot;BB];
MTtot=[MTtot; MT];
Rtot = [Rtot; R];

lastV = size(V);
lastVg = size(Vg);
lastX = size(X);
lastL = size(L);
lastLNM = size(LNM);
lastG = size(G);
lastAB = size(AB);
lastBB = size(BB);
lastMT = size(MT);
lastR = size(R);

initV3 = V(lastV(1),1); %L
initVg3 = Vg(lastVg(1),1); %L
initX3 = X(lastX(1),:); %gDCW
initL3 = L(lastL(1),:);   %mmol
initLNM3 = LNM(lastLNM(1),:); %mmol
initG3 = G(lastG(1),:);%mmol
initAB3 = AB(lastAB(1),:);   %mmol
initBB3 = BB(lastBB(1),:); %mmol
initMT3 = MT(lastMT(1),:);%mmol
initR3 = R(lastR(1),:);

initialConditions = [initV3 initVg3 initX3 initL3 initLNM3 initG3 initAB3 initBB3 initMT3 initR3];

%Segment 3
   
tspan = [3:0.01:10];
    
    [t, V, Vg, X, L, LNM, G, AB, BB, MT, R] = Co4Research3(tspan,initialConditions);       

ttot=[ttot; t];
Vtot=[Vtot; V];
Vgtot=[Vgtot; Vg];
Xtot=[Xtot; X];
Ltot=[Ltot; L];
LNMtot = [LNMtot;LNM];
Gtot=[Gtot; G];
ABtot=[ABtot; AB];
BBtot = [BBtot;BB];
MTtot=[MTtot; MT];
Rtot=[Rtot;R];

end
