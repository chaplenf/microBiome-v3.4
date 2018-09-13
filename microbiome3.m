clear;
global community HenrysLawCoefficients kLa Temp
%sets initial pH, NH3 concentration, O2 concentration, and carbonate
%reactions

pH=4;%suggested range for initial pH: 6.8-7.8
NH3=0.025; %suggested range for initial NH3: 0.0 to 0.125 mmol
NO2=0.0; %suggested range for initial NO2: 0.0 to 0.125 mmol
O2=0.001172; %suggested range for initial O2: 0 to 0.001172 mmol 
CO3=0.0115; %suggested range for initial carbonate: 0 to 0.0115 mmol
PO4=0.06875; %suggested range for initial phosphate buffer: 0.025 to 0.125 mmol
Temp=30; %suggested range for temperature; 10-30 degrees C
OD(1,1)=0.2; %OD of Nitrosomonas europaea. Suggested range 0 to 0.25.
OD(1,2)=0.05; %OD of Nitrobacter winogradskyi. Suggested range 0 to 0.10.

%Decisions with respect to scavenging NO reactions and nitrite reductase, an enzyme involved in nitric oxide gas production
% [NO scavenger;oxidase reactions;nitrite reductase mutant] ... 1 switches on reactions ... 0 switches off reactions

%Could easily add additional inhibitors/mutants; May be able to allow partial downregulation (intermediate states)

    [ttot,Vtot,Vgtot,Xtot,Ltot,LNMtot,Gtot, ABtot, BBtot, MTtot, Rtot] = runFileResearch3(pH, NO2, NH3, CO3, PO4, O2, OD); 
    
%ttot,Vtot,Vgtot,Xtot,Ltot,LNMtot,Gtot,ABtot,BBtot,Rtot represent model variable outputs.  
    
%Final N2O concentration from experimental run and column vectors
    
    N2O=1000000.*(Gtot(:,6)).*(0.00008206*(Temp+273.15)./0.155);

%Figures of key variables to check for integration consistency.  Nitrate, nitrite and
%Ammonia should not go above 5 mM (0.025 mmol for a 0.005 liter liquid volume).

finalOD(:,1)=Xtot(:,1)./(0.16*0.005);
finalOD(:,2)=Xtot(:,2)./(0.251*0.005);
X=Xtot./0.005;
L=Ltot./0.005;
AB=ABtot/0.005;
BB=BBtot/0.005;
MT=MTtot/0.005;
LNM=LNMtot./0.005;
NOxppmv=1000000.*(Gtot(:,1)+Gtot(:,2)).*(0.00008206*(Temp+273.15)./0.155);
N2Oppmv=1000000.*(Gtot(:,6)).*(0.00008206*(Temp+273.15)./0.155);
pH=-(log10(0.001*L(:,13)));

    figure;
subplot(5,1,1);plot(ttot,finalOD);ylabel('Optical Density (OD)'); 
axis ([0 10 0 0.25]);
subplot(5,1,2);plot(ttot,X);ylabel('Biomass (gDCW/L)');xlabel('Time (Hours)');
axis ([0 10 0 0.05]);
    figure;
subplot(5,1,1);plot(ttot,L(:,11));ylabel('NH2OH (mMol)');
axis ([0 10 0 0.25]);
subplot(5,1,2);plot(ttot,NOxppmv);ylabel('NOx (ppmV)');
axis ([0 10 0 200]);
subplot(5,1,3);plot(ttot,N2Oppmv);ylabel('N2O (ppmV)');
axis ([0 10 0 100]);
subplot(5,1,4);plot(ttot,pH);ylabel('pH');xlabel('Time (Hours)');
axis ([0 10 6 9]);
    figure;
subplot(5,1,1);plot(ttot,L(:,5));ylabel('Ammonia (mMol)');
axis ([0 10 0 10]);
subplot(5,1,2);plot(ttot,L(:,1));ylabel('Nitrite (mMol)');
axis ([0 10 0 10]);
subplot(5,1,3);plot(ttot,L(:,2));ylabel('Nitrate (mMol)');
axis ([0 10 0 10]);
subplot(5,1,4);plot(ttot,L(:,4));ylabel('Oxygen (mMol)');xlabel('Time (Hours)');
axis ([0 10 0 0.3]);