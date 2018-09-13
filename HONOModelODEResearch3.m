 %Rate constants and other parameters
kga1 = 19.764; 
kga2=  7416; 
KsNO = 1.3*10^-2;
KsNO2 = 1.3*10^-2; 
kminus1 = 10^11; %Assume equilibrium is achieved for each time step. 3600 is conversion from s-1 to h-1
k1=(10^14)*10^-3.15; %Assume equilibrium is achieved for each time step. 1 is concentration of water in mM 3.3 is pKa of HONO.
k2 = 8640;  % 8640 %liquid oxidation of 4 NO to 4 NO2- AB(2)
k3 = 24.84; 
k4 = 102.96; 
kminus4 = 6.012*10^8; 
k5 = 17280000;
kminus6=10^11; %Assume equilibrium is achieved for each time step
k6=(10^14)*10^-2.148; %Assume equilibrium is achieved for each time step. 1 is concentration of water in mM 2.12 is pKa of H2PO4-.
kminus7=10^11; %Assume diffusion limited reaction rate
k7=(10^14)*10^-6.84; %Assume equilibrium is achieved for each time step. 1 is concentration of water in mM 7.21 is pKa of HPO42-.
kminus8=10^11; %Assume diffusion limited reaction rate
k8=(10^14)*10^-12.375; %Assume equilibrium is achieved for each time step. 1 is concentration of water in mM 12.32 is pKa of PO43-.
kminus9=10^11; %Assume diffusion limited reaction rate
k9=(10^14)*10^-14; %Assume diffusion limited reaction rate
k10=3600*0.06;
kminus10=3600*20;
kminus11=10^11;%Assume diffusion limited reaction rate
k11=(10^14)*10^-3.33;%Assume diffusion limited reaction rate
kminus12=10^11;%Assume diffusion limited reaction rate
k12=(10^14)*10^-10.29;%Assume diffusion limited reaction rate
k13=3600*2*10^-3;
kminus13=3600*0.001; 
kminus14=10^11; %Assume diffusion limited reaction rate
k14=(10^14)*10^-9.25; %Assume diffusion limited reaction rate

%Liquid Volume (L)
Vdum(1)=y(1);

%Gas volume (L)
Vdum(2)=y(2);

%Biomass Concentration Transformation (gDCW/L)
Xdum=y(3:2+numSpecies)/Vdum(1);

%Liquid Phase Species Concentration Transformation (mmol/L)
Ldum=y(2+numSpecies+1:2+numSpecies+numLiquidSpecies)/Vdum(1);

%Liquid Phase Non-Metabolite Species Concentration Transformation (mmol/L)
LNMdum=y(2+numSpecies+numLiquidSpecies+1:2+numSpecies+numLiquidSpecies+numLiquidNMSpecies)/Vdum(1);

%Gas Phase Species Concentration Transformation (atm)
Gdum=y(2+numSpecies+numLiquidSpecies+numLiquidNMSpecies+1:2+numSpecies+numLiquidSpecies+numLiquidNMSpecies+numGasSpecies)*(0.00008206*(Temp+273.15))/Vdum(2);

%Liquid species reactor balance (mmol)
LBdum=((community.Fin/Vdum(1))*(community.Lfeed' - Ldum)*Vdum(1));

%Liquid Non-Metabolite Species reactor balance (mmol)
LNMBdum=(community.Fin/Vdum(1)*(community.LNMfeed' - LNMdum)*Vdum(1));

%Gas species reactor balance (mmol)
GBdum=(community.Vgin/Vdum(2)*(community.Gfeed' - Gdum)*(Vdum(2)/0.00008206*(Temp+273.15)));

%Biotic Reaction Buckets (mmol/(L))
BBdum(1)= Xdum(1)*vs(1,1); %NO2- B(1)
BBdum(2)= Xdum(2)*vs(2,1); %NO2- B(2)
BBdum(3)= Xdum(1)*vs(1,2); %nitrate B(3)
BBdum(4)= Xdum(2)*vs(2,2); %nitrate B(4)
BBdum(5)= Xdum(1)*vs(1,3); %NO B(5)
BBdum(6)= Xdum(2)*vs(2,3); %NO B(6)
BBdum(7)= Xdum(1)*vs(1,4); %O2 B(7)
BBdum(8)= Xdum(2)*vs(2,4); %O2 B(8)
BBdum(9)= Xdum(1)*vs(1,5); %NH3 B(9)
BBdum(10)= Xdum(2)*vs(2,5); %NH3 B(10)
BBdum(11)= Xdum(1)*vs(1,6); %CHO2 B(11)
BBdum(12)= Xdum(2)*vs(2,6); %CHO2 B(12)
BBdum(13)= Xdum(1)*vs(1,7); %N2O B(13)
BBdum(14)= Xdum(2)*vs(2,7); %N2O B(14)
BBdum(15)= Xdum(1)*vs(1,8); %N2 B(15)
BBdum(16)= Xdum(2)*vs(2,8); %N2 B(16)
BBdum(17)= Xdum(1)*vs(1,9); %CO2 B(17)
BBdum(18)= Xdum(2)*vs(2,9); %CO2 B(18)
BBdum(19)= Xdum(1)*vs(1,10); %HCO3 B(19)
BBdum(20)= Xdum(2)*vs(2,10); %HCO3 B(20)
BBdum(21)= Xdum(1)*vs(1,11); %NH2OH B(21)
BBdum(22)= Xdum(2)*vs(2,11); %NH2OH B(22)
BBdum(23)= Xdum(1)*vs(1,12); %NOhao B(23)
BBdum(24)= Xdum(2)*vs(2,12); %NOhao B(24)
BBdum(25)= Xdum(1)*vs(1,13); %H+ B(25)
BBdum(26)= Xdum(2)*vs(2,13); %H+ B(26)
BBdum(27)= Xdum(1)*vs(1,14); %HPO42- B(27)
BBdum(28)= Xdum(2)*vs(2,14); %HPO42- B(28)
BBdum(29)= Xdum(1)*vs(1,15); %NH4 B(29)
BBdum(30)= Xdum(2)*vs(2,15); %NH4 B(30)

%Mass transfer Buckets(mmol/(L))
MBdum(1)= kLa(3)*(1000*HenrysLawCoefficients(3)*Gdum(1)-Ldum(3)); %NO
MBdum(2)= -kga1*(Gdum(2)-(1/(1000*HenrysLawCoefficients(7)))*LNMdum(1))*Vdum(2)/(0.00008206*(Temp+273.15)); %NO2
MBdum(3)= -kga2*(Gdum(3)-(1/(1000*HenrysLawCoefficients(8)))*LNMdum(2))*Vdum(2)/(0.00008206*(Temp+273.15)); %HONO
MBdum(4)= kLa(1)*(1000*HenrysLawCoefficients(1)*Gdum(4)-Ldum(4)); %O2
MBdum(5)= kLa(2)*(1000*HenrysLawCoefficients(2)*Gdum(5)-Ldum(15)); %NH3
MBdum(6)= kLa(5)*(1000*HenrysLawCoefficients(5)*Gdum(6)-Ldum(7)); %N2O
MBdum(7)= kLa(6)*(1000*HenrysLawCoefficients(6)*Gdum(7)-Ldum(8)); %N2
MBdum(8)= kLa(4)*(1000*HenrysLawCoefficients(4)*Gdum(8)-Ldum(9)); %CO2

%Abiotic reaction buckets (mmol/(L))
ABdum(1)= k1*LNMdum(2)*1-(kminus1*Ldum(1)*Ldum(13)); %liquid formation of NO2- from HONO AB(1)
ABdum(2)= 4*k2*Ldum(3)*Ldum(3)*Ldum(4); %liquid oxidation of 4 NO to 4 NO2- AB(2)
ABdum(3)= 2*k3*Gdum(1)*Gdum(1)*Gdum(4)/((0.00008206*(Temp+273.15))^3); %gas oxidation of 2 NO2 to 2 NO AB(3)
ABdum(4)= k4*LNMdum(2)*LNMdum(2)-kminus4*((Ldum(3)))*LNMdum(1); %liquid side oxidation of 2 HONO to NO and NO2 AB(4)
ABdum(5)= k5*LNMdum(1)*LNMdum(1); %liquid side oxidation of 2NO2 to NO3- and HNO2 AB(5)
ABdum(6)= k6*LNMdum(3)*1-kminus6*LNMdum(4)*Ldum(13); %formation of H2PO4- from H3PO4 AB(6)
ABdum(7)= k7*LNMdum(4)*1-kminus7*Ldum(14)*Ldum(13); %formation of HPO42- from H2PO4- AB(7)
ABdum(8)= k8*Ldum(14)*1-kminus8*LNMdum(5)*Ldum(13); %formation of PO43- from HPO42- AB(8)
ABdum(9)= k9*1*1-kminus9*LNMdum(6)*Ldum(13); %consumption of water to form OH- and H+ AB(9)
ABdum(10)= k10*Ldum(9)-kminus10*LNMdum(7); %formation of H2CO3 from aqueous CO2 AB(10)
ABdum(11)= k11*LNMdum(7)-kminus11*Ldum(13)*Ldum(10); %liquid consumption of H2CO3 to form HCO3- and H+ AB(11)
ABdum(12)= k12*Ldum(10)-kminus12*LNMdum(8)*Ldum(13); %liquid consumption of HCO3- to form CO32- and H+ AB(12)
ABdum(13) = k13*LNMdum(8)*LNMdum(9)-kminus13*LNMdum(10); %liquid formation of CaCO3 from CO32- and Ca2+ AB(13)
ABdum(14) = k14*Ldum(5)*1-kminus14*Ldum(13)*Ldum(15); %liquid formation of NH3 from NH4+ AB(14)

pH=-log10(0.001*Ldum(13))

%Bioreactor volume changes (L)
dy(1) = community.Fin - community.Fout;
dy(2) = community.Vgin - community.Vgout;

%Biomass changes (gDCW)
dy(3)= Vdum(1)*(mu(1)*Xdum(1)+community.Fin/Vdum(1)*(community.Xfeed(1) - Xdum(1))); %Neuro Biomass X(1)
dy(4)= Vdum(1)*(mu(2)*Xdum(2)+community.Fin/Vdum(1)*(community.Xfeed(2) - Xdum(2))); %Nwino Biomass X(2)
   
%Liquid phase metabolite reaction network (mmol)
dy(5)= Vdum(1)*(ABdum(1)+ABdum(2)+BBdum(1)+BBdum(2))+LBdum(1); %NO2- L(1)
dy(6)= Vdum(1)*(ABdum(5)+BBdum(3)+BBdum(4))+LBdum(2); %nitrate L(2)
dy(7)= Vdum(1)*MBdum(1)+Vdum(1)*(ABdum(4)-ABdum(2)+BBdum(5)+BBdum(23)+BBdum(6)+BBdum(24))+LBdum(3); %NO L(3)
dy(8)= Vdum(1)*MBdum(4)+Vdum(1)*(-0.25*ABdum(2)+BBdum(7)+BBdum(8))+LBdum(4); %O2 L(4)
dy(9)= Vdum(1)*(-ABdum(14)+BBdum(9)+BBdum(10))+LBdum(5); %NH4+ L(5) 
dy(10)= Vdum(1)*(BBdum(10)+ BBdum(12))+LBdum(6); %CHO2 L(6)
dy(11)= Vdum(1)*MBdum(6)+Vdum(1)*(BBdum(13)+ BBdum(14))+LBdum(7); %N2O L(7)%
dy(12)= Vdum(1)*MBdum(7)+Vdum(1)*(BBdum(15)+ BBdum(16))+LBdum(8); %N2 L(8)
dy(13)= Vdum(1)*MBdum(8)+Vdum(1)*(-ABdum(10)+BBdum(17)+BBdum(18))+LBdum(9); %CO2 L(9)
dy(14)= Vdum(1)*(ABdum(11)-ABdum(12)+BBdum(19)+ BBdum(20))+LBdum(10); %HCO3 L(10)
dy(15)= Vdum(1)*(BBdum(21)+BBdum(22))+LBdum(11); %L(11) NH2OH
dy(16)= Vdum(1)*(BBdum(23)+BBdum(24))+LBdum(12); %L(12) NOHAO
dy(17)= Vdum(1)*(ABdum(1)+ABdum(2)+ABdum(5)+ABdum(6)+ABdum(7)+ABdum(8)+ABdum(9)+ABdum(11)+ABdum(12)+ABdum(14)+BBdum(25)+BBdum(26))+LBdum(13); %H+ L(13)
dy(18)= Vdum(1)*(ABdum(7)-ABdum(8)+BBdum(27)+BBdum(28))+LBdum(14); %HPO42- L(14)
dy(19)=Vdum(1)*(MBdum(5)+ABdum(14)+BBdum(29)+BBdum(30))+LBdum(15); %NH3 L(15)

%Liquid phase non-metabolite reaction network (mmol)
dy(20)= Vdum(1)*(ABdum(4)-2*ABdum(5))-Vdum(2)*MBdum(2)+LNMBdum(1); %NO2 LNM(1)
dy(21)= Vdum(1)*(-ABdum(1)-2*ABdum(4)+ABdum(5))-Vdum(2)*MBdum(3)+LNMBdum(2);%HONO LNM(2)
dy(22)= Vdum(1)*(-ABdum(6))+LNMBdum(3); %H3PO4 LNM(3)
dy(23)= Vdum(1)*(ABdum(6)-ABdum(7))+LNMBdum(4); %H2PO4- LNM(4)
dy(24)= Vdum(1)*(ABdum(8))+LNMBdum(5); %PO43- LNM(5)
dy(25)= Vdum(1)*(ABdum(9))+LNMBdum(6); %OH- LNM(6)
dy(26)= Vdum(1)*(ABdum(10)-ABdum(11))+LNMBdum(7); %H2CO3 LNM(7)
dy(27)= Vdum(1)*(ABdum(12)-ABdum(13))+LNMBdum(8); %CO32- LNM(8)
dy(28)=Vdum(1)*(-ABdum(13))+LNMBdum(9); %Ca2+ LNM(9) 
dy(29)=Vdum(1)*(ABdum(13))+LNMBdum(10); %CaCO3(s)LNM(10)

%Gas phase reaction/mass transfer network (mmol)
dy(30)= -Vdum(1)*MBdum(1)-Vdum(2)*ABdum(3)+GBdum(1)*Vdum(2)/(0.00008206*(Temp+273.15)); %NO(g)G(1)
dy(31)= Vdum(2)*MBdum(2)+Vdum(2)*ABdum(3)+GBdum(2)*Vdum(2)/(0.00008206*(Temp+273.15));%NO2(g) G(2)
dy(32)= Vdum(2)*MBdum(3)+ GBdum(3)*Vdum(2)/(0.00008206*(Temp+273.15)); %HONO (g) G(3)
dy(33)= -Vdum(1)*MBdum(4)-0.5*Vdum(2)*ABdum(3)+GBdum(4)*Vdum(2)/(0.00008206*(Temp+273.15)); %O2 (g) G(4)
dy(34)= -Vdum(1)*MBdum(5)+GBdum(5)*Vdum(2)/(0.00008206*(Temp+273.15)); %NH3 (g) G(5)
dy(35)= -Vdum(1)*MBdum(6)+GBdum(6)*Vdum(2)/(0.00008206*(Temp+273.15)); %N2O (g) G(6)
dy(36)= -Vdum(1)*MBdum(7)+GBdum(7)*Vdum(2)/(0.00008206*(Temp+273.15)); %N2 (g) G(7)
dy(37)= -Vdum(1)*MBdum(8)+GBdum(8)*Vdum(2)/(0.00008206*(Temp+273.15)); %CO2 (g) G(8)

%Abiotic Buckets Retrieval (mmol/(L))
dy(38)= 0;
dy(39)= 0;
dy(40)= 0;
dy(41)= 0;
dy(42)= 0;
dy(43) = 0;
dy(44)= 0;
dy(45)= 0;
dy(46)= 0;
dy(47)= 0;
dy(48)= 0;
dy(49)= 0;
dy(50)= 0;
dy(51)= 0;

%Biotic Buckets Retrieval (mmol/(L))
dy(52)= Vdum(1)*Xdum(1)*vs(1,1); %NO2- B(1)
dy(53)= Vdum(1)*Xdum(2)*vs(2,1); %NO2- B(2)
dy(54)= Vdum(1)*Xdum(1)*vs(1,2); %nitrate B(3)
dy(55)= Vdum(1)*Xdum(2)*vs(2,2); %nitrate B(4)
dy(56)= Vdum(1)*Xdum(1)*vs(1,3); %NO B(5)
dy(57)= Vdum(1)*Xdum(2)*vs(2,3); %NO B(6)
dy(58)= Vdum(1)*Xdum(1)*vs(1,4); %O2 B(7)
dy(59)= Vdum(1)*Xdum(2)*vs(2,4); %O2 B(8)
dy(60)= Vdum(1)*Xdum(1)*vs(1,5); %NH4+ B(9)
dy(61)= Vdum(1)*Xdum(2)*vs(2,5); %NH4+ B(10)
dy(62)= Vdum(1)*Xdum(1)*vs(1,6); %CHO2 B(11)
dy(63)= Vdum(1)*Xdum(2)*vs(2,6); %CHO2 B(12)
dy(64)= Vdum(1)*Xdum(1)*vs(1,7); %N2O B(13)
dy(65)= Vdum(1)*Xdum(2)*vs(2,7); %N2O B(14)
dy(66)= Vdum(1)*Xdum(1)*vs(1,8); %N2 B(15)
dy(67)= Vdum(1)*Xdum(2)*vs(2,8); %N2 B(16)
dy(68)= Vdum(1)*Xdum(1)*vs(1,9); %CO2 B(17)
dy(69)= Vdum(1)*Xdum(2)*vs(2,9); %CO2 B(18)
dy(70)= Vdum(1)*Xdum(1)*vs(1,10); %HCO3 B(19)
dy(71)= Vdum(1)*Xdum(2)*vs(2,10); %HCO3 B(20)
dy(72)= Vdum(1)*Xdum(1)*vs(1,11); %NH2OH B(21)
dy(73)= Vdum(1)*Xdum(2)*vs(2,11); %NH2OH B(22)
dy(74)= Vdum(1)*Xdum(1)*vs(1,12); %NOhao B(23)
dy(75)= Vdum(1)*Xdum(2)*vs(2,12); %NOhao B(24)
dy(76)= Vdum(1)*Xdum(1)*vs(1,13); %H+ B(25)
dy(77)= Vdum(1)*Xdum(2)*vs(2,13); %H+ B(26)
dy(78)= Vdum(1)*Xdum(1)*vs(1,14); %PO43- B(27)
dy(79)= Vdum(1)*Xdum(2)*vs(2,14); %PO43- B(28)
dy(80)= Vdum(1)*Xdum(1)*vs(1,15); %NH3 B(29)
dy(81)= Vdum(1)*Xdum(2)*vs(2,15); %NH3 B(30)

%Mass Transfer Buckets (mmol/L)
dy(82)=0;
dy(83)=0;
dy(84)=0;
dy(85)=0;
dy(86)=0;
dy(87)=0; 
dy(88)=0;
dy(89)=0;

%Internal Flux Retrieval (mmol/(L))
dy(90)=Xdum(1)*vs(1,16); %R(1) rxn14202a
dy(91)=Xdum(2)*vs(2,16); %R(1) 
dy(92)=Xdum(1)*vs(1,17); %R(2) rxn14202b
dy(93)=Xdum(2)*vs(2,17); %R(2) 
dy(94)=Xdum(1)*vs(1,18); %R(3) rxn14202c
dy(95)=Xdum(2)*vs(2,18); %R(3) 
dy(96)=Xdum(1)*vs(1,19); %R(4) rxn14202d
dy(97)=Xdum(2)*vs(2,19); %R(4) 
dy(98)=Xdum(1)*vs(1,20); %R(5) rxn014202e
dy(99)=Xdum(2)*vs(2,20); %R(5) 
dy(100)=Xdum(1)*vs(1,21); %R(6) rxn14202f
dy(101)=Xdum(2)*vs(2,21); %R(6)
dy(102)=Xdum(1)*vs(1,22); %R(7) rxn14202g
dy(103)=Xdum(2)*vs(2,22); %R(7) 
dy(104)=Xdum(1)*vs(1,23); %R(8) rxn1806mod
dy(105)=Xdum(2)*vs(2,23); %R(8) 
dy(106)=Xdum(1)*vs(1,24); %R(9) rxn05795mod
dy(107)=Xdum(2)*vs(2,24); %R(9) 
dy(108)=Xdum(1)*vs(1,25); %R(10) rxnP460a
dy(109)=Xdum(2)*vs(2,25); %R(10) 
dy(110)=Xdum(1)*vs(1,26); %R(11) rxnP460b
dy(111)=Xdum(2)*vs(2,26); %R(11) 
dy(112)=Xdum(1)*vs(1,27); %R(12) rxn10113
dy(113)=Xdum(2)*vs(2,27); %R(12) 
dy(114)=Xdum(1)*vs(1,28); %R(13) rxn10125mod
dy(115)=Xdum(2)*vs(2,28); %R(13)
dy(116)=Xdum(1)*vs(1,29); %R(14) rxn08975
dy(117)=Xdum(2)*vs(2,29); %R(14) 
dy(118)=Xdum(1)*vs(1,30); %R(15) rxn10122mod
dy(119)=Xdum(2)*vs(2,30); %R(15) 
dy(120)=Xdum(1)*vs(1,31); %R(16) rxn10043mod
dy(121)=Xdum(2)*vs(2,31); %R(16) 
dy(122)=Xdum(1)*vs(1,32); %R(17) rxn00058
dy(123)=Xdum(2)*vs(2,32); %R(17) 
dy(124)=Xdum(1)*vs(1,33); %R(18) rxn10042
dy(125)=Xdum(2)*vs(2,33); %R(18) 
dy(126)=Xdum(1)*vs(1,34); %R(19) rxn12750
dy(127)=Xdum(2)*vs(2,34); %R(19) 
dy(128)=Xdum(1)*vs(1,35); %R(20) rxn14173a
dy(129)=Xdum(2)*vs(2,35); %R(20)
dy(130)=Xdum(1)*vs(1,36); %R(21) rxn14173b
dy(131)=Xdum(2)*vs(2,36); %R(21) 
dy(132)=Xdum(1)*vs(1,37); %R(22) rxn00567a
dy(133)=Xdum(2)*vs(2,37); %R(22) 
dy(134)=Xdum(1)*vs(1,38); %R(23) rxn00567b
dy(135)=Xdum(2)*vs(2,38); %R(23) 
dy(136)=Xdum(1)*vs(1,39); %R(24) rxn00567c
dy(137)=Xdum(2)*vs(2,39); %R(24) 
dy(138)=Xdum(1)*vs(1,40); %R(25) rxn09008
dy(139)=Xdum(2)*vs(2,40); %R(25) 
dy(140)=Xdum(1)*vs(1,41); %R(26) rxn12750mod
dy(141)=Xdum(2)*vs(2,41); %R(26) 
dy(142)=Xdum(1)*vs(1,42); %R(27) rxn14010
dy(143)=Xdum(2)*vs(2,42); %R(27) 
dy(144)=Xdum(1)*vs(1,43); %R(28) rxn10811
dy(145)=Xdum(2)*vs(2,43); %R(28) 
dy(146)=Xdum(1)*vs(1,44); %R(29) Cytochrome Exchange1
dy(147)=Xdum(2)*vs(2,44); %R(29) 
dy(148)=Xdum(1)*vs(1,45); %R(30) Cytochrome Exchange2
dy(149)=Xdum(2)*vs(2,45); %R(30) 
dy(150)=Xdum(1)*vs(1,46); %R(31) Cytochrome Exchange3
dy(151)=Xdum(2)*vs(2,46); %R(31) 




