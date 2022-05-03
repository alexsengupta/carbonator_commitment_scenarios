function [A,AF,AM,Alk, Cat,Clo,Co,Cs,Cup,Ksp,L,N,OM,P,So,VF,a2,alb0,alphaRF_ch4,alpha_ch4,d,dr0,e,eps,g,k1,k2,ka,kd,kh,m,tau_ch4_pi,vtau,tau_N20,alphaRF_N2O,tau_CFC12,alphaRF_CFC12]=load_constants(parameter_set,fn);
% parameters based on an 1850 start date

load(fn)
clearvars -except OPT* parameter_set

% Constants for temperature model
Cs=9;       % *60*60*24*365; %W yr/m2/K (Geoffroy uses 7.3 - 9 gives a better volcanic response
Co=106;     % *60*60*24*365; %W yr/m2/K

L=1;        % 1.13 from GEOFFROY, but 1 matches more closely to CMIP5 multi-model mean; %Feedback parameter W/m2/K %L ~ 3.5 for no feedback
g=0.73;     % heat exchange parameter W/m2/K
eps=0.6;    % Internal variavility magnitude
vtau=1.2;   % yr Volcanic aerosol decay timescale
alb0=0.31;   % initial albedo
AF=-0.0052; % conversion from SO2 emissions to RF W/m2/TgSO2/yr
VF=-20;     % conversion from optical thickness of volcanic aerosols to RF

% Constants for ocean carbon model (glotter)
d=50;       % ratio of upper to deep ocean volume
ka=1/5;     % Glotter paper suggests 1/5
kd=1/20;    % Glotter paper suggests 1/20
OM=7.8e22;  % moles of water in ocean
AM=1.77e20; % moles of air in atmosphere
Alk=767; %GtC   % Alkalinity; assumed constant here but buffered on long timescales O[10ka] from the dissolution of CaCO3
kh=1.0548e+03;  % at 15oC; ratio of the molar concentrations of CO2 in atmosphere and ocean (Henry's Law)
k1=8.7184e-07;  % at 15oC; dissociation constant
k2=5.4426e-10;  % at 15oC; dissociation constant
A=kh*AM/(OM/(d+1)); % A is the ratio of mass of CO2 in atmospheric to upper ocean dissolved CO2,% i.e. A is inversely proportional to
% CO2 solubility. A is temperature dependent, however the % effect is small and has been neglected here.
A=132.216074; % Slightly different to above in order to have equilibrium at 1850

% Constants for Terrestrial model (Svirezhev)
m=8.7e-2;   % /yr  1/residence time of carbon in vegetation
a2=4.7e-4;  % /GtC %strength of CO2 fertilisation
dr0=0.024625; %0.025;  % /yr decomposition rate of soil
e=0.5;      % proportion of vegetation that forms soil (remainder goes to atmosphere)

% constant for Aragonite calculation
Ksp=10^(-6.19); %=6.4565e-7 (aragonite dissolution constant at 25degC)

%N2O model
tau_N20=109;
alphaRF_N2O=0.12;

%Halogen Model
tau_CFC12=102;%trs
alphaRF_CFC12=0.32; % The radiative forcing due to CFC-12 of 0.32 Wm−2 ppbv−1 used in WMO (1999) is retained [Ch 6 AR4]



% OPTIMISED Constants for methane model
% tau_ch4_pi=8; % pre industrial methane decay
% alpha_ch4=-0.12;
alphaRF_ch4=0.0316;
tau_ch4_pi= 7.8865;
alpha_ch4= -0.1538;


% Carbon cycle parameters - OPTIMISED
A=OPT_A(parameter_set);                     
a2=OPT_a2(parameter_set);                    
d=OPT_d(parameter_set);                     
dr0=OPT_dr0(parameter_set);                   
e=OPT_e(parameter_set);                     
ka=OPT_ka(parameter_set);                    
kd=OPT_kd(parameter_set);                    
m=OPT_m(parameter_set);                     

%Initial conditions - OPTIMISED
Clo=OPT_Clo(parameter_set);                   
Cup=OPT_Cup(parameter_set); 
Cat= 596; %OPTIMISED

N=OPT_N(parameter_set);                     
P=OPT_P(parameter_set);                     
So=OPT_So(parameter_set);  


figure(55);clf
plot(0:1/12:100,OPT_CO2_100yr_ts','linewidth',2)
legend('1','2','3','4','5','6','7')
set(gca,'fontsize',15)
% print -f55 -depsc 'figures/Cparam_decay'

