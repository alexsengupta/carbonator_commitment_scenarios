function [years,meinhausenRF_OZO,meinhausenRF_AER,TSI,alb,emissionCH4,emissionCO2,emissionN2O,emissionSO2,emission_volc,radf_other,radf_halo,radf_N2O,radf_CO2,radf_CH4,conc_CO2,conc_CH4,conc_N2O,mTSI,conc_Hal] = load_RCP(myrcp,startyear,endyear);

switch myrcp
    case 'RCP85'
        eCO2=59;
        eCH4=60;
        eN2O=61;
        eSO2=67;
        cCO2=43;
        cCH4=44;
        cN2O=45;
        rCO2=21;
        rCH4=22;
        rN2O=23;
        rhal=24;
        roth=25;
    case 'RCP6'
        eCO2=50;
        eCH4=51;
        eN2O=52;
        eSO2=64;
        cCO2=28;
        cCH4=29;
        cN2O=30;
        rCO2=3;
        rCH4=4;
        rN2O=5;
        rhal=6;
        roth=7;
    case 'RCP3'
        eCO2=56;
        eCH4=57;
        eN2O=58;
        eSO2=66;
        cCO2=38;
        cCH4=39;
        cN2O=40;
        rCO2=15;
        rCH4=16;
        rN2O=17;
        rhal=18;
        roth=19;
    case 'RCP45'
        eCO2=53;
        eCH4=54;
        eN2O=55;
        eSO2=65;
        cCO2=33;
        cCH4=34;
        cN2O=35;
        rCO2=9;
        rCH4=10;
        rN2O=11;
        rhal=12;
        roth=13;
    case 'COM20'
        eCO2=59;
        eCH4=60;
        eN2O=61;
        eSO2=67;
        cCO2=43;
        cCH4=44;
        cN2O=45;
        rCO2=21;
        rCH4=22;
        rN2O=23;
        rhal=24;
        roth=25;
    case 'COM30'
        eCO2=59;
        eCH4=60;
        eSO2=67;
        eN2O=61;
        cCO2=43;
        cCH4=44;
        cN2O=45;
        rCO2=21;
        rCH4=22;
        rN2O=23;
        rhal=24;
        roth=25;
    case 'COM40'
        eCO2=59;
        eCH4=60;
        eN2O=61;
        eSO2=67;
        cCO2=43;
        cCH4=44;
        cN2O=45;
        rCO2=21;
        rCH4=22;
        rN2O=23;
        rhal=24;
        roth=25;
end

% Load RCP85 emissions and RF
disp(myrcp)

DT=1/12;
default_constants_optimised


% load Greenhouse gas and aerosol emission, concentration, forcing from
% PCMDI and interpolate in time
rcp_hist_RF_CO2e=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/rcp_hist_RF_CO2e.xls');
years=startyear:DT:endyear;
ind1870=find(years>1861 & years<=1880);
ind2012=find(years>2006 & years<=2018);

% Volcanic
% load volcanic optical thickness (http://data.giss.nasa.gov/modelforce/strataer/)
clear OT* E_*
load ('/Users/z3045790/Dropbox/Simple Climate Model/EV.txt')
emission_volc=interp1(EV(1,:),EV(2,:),years);
emission_volc(isnan(emission_volc))=0;
emission_volc(emission_volc<0)=0;

% Solar
% load TSI - ANNUAL MEAN TSI: Lean (GRL 2000) with Wang Lean Sheeley (ApJ 2005)
TSI=xlsread('/Users/z3045790/Dropbox/Simple Climate Model/Code/Code in Matlab/TSI_WLS_ann_1610_2008.xls');
TSI=interp1(TSI(:,1),TSI(:,3),years);
mTSI=nanmean(TSI);
TSI(isnan(TSI))=mTSI;

alb(1:length(years))=alb0;

% interpolate data to DT timestep
for n=1:length(rcp_hist_RF_CO2e(:,1))
    if ~all(isnan(rcp_hist_RF_CO2e(n,:)))
        rcp_hist_RF_CO2eI(n,:)=interp1(rcp_hist_RF_CO2e(1,:),rcp_hist_RF_CO2e(n,:),years);
    end
end

emissionN2O=rcp_hist_RF_CO2eI(eN2O,:);
emissionCH4=rcp_hist_RF_CO2eI(eCH4,:);
emissionCO2=rcp_hist_RF_CO2eI(eCO2,:);
emissionSO2=rcp_hist_RF_CO2eI(eSO2,:);

conc_N2O=rcp_hist_RF_CO2eI(cN2O,:);
conc_CH4=rcp_hist_RF_CO2eI(cCH4,:);
conc_CO2=rcp_hist_RF_CO2eI(cCO2,:);

radf_CH4=rcp_hist_RF_CO2eI(rCH4,:);
radf_CO2=rcp_hist_RF_CO2eI(rCO2,:);
radf_N2O=rcp_hist_RF_CO2eI(rN2O,:);
radf_halo=rcp_hist_RF_CO2eI(rhal,:);
radf_other=rcp_hist_RF_CO2eI(roth,:);


% create emission data from 1750
ind1850_60 =find(years>1850 & years<=1860);
ind1855=find(years==1855);
emissionN2O(1:ind1855)=(1:ind1855)/(ind1855)*mean(emissionN2O(ind1850_60));
emissionCH4(1:ind1855)=(1:ind1855)/(ind1855)*mean(emissionCH4(ind1850_60));
emissionCO2(1:ind1855)=(1:ind1855)/(ind1855)*mean(emissionCO2(ind1850_60));
emissionSO2(1:ind1855)=(1:ind1855)/(ind1855)*mean(emissionSO2(ind1850_60));
emission_volc(1:ind1855)=0;
TSI(1:ind1855)=mTSI;



ind1850=find(years==1850);
conc_N2O(1:ind1850)=linspace(273,conc_N2O(ind1850),ind1850);
conc_CO2(1:ind1850)=linspace(280,conc_CO2(ind1850),ind1850);
conc_CH4(1:ind1850)=linspace(722,conc_CH4(ind1850),ind1850);

% RF data from Meinshausen et al 2011

myrcp_=myrcp; 
if strcmp(myrcp,'RCP3');myrcp_='RCP3PD';end
if strcmp(myrcp,'COM20') | strcmp(myrcp,'COM30')| strcmp(myrcp,'COM40');myrcp_='RCP85';end
data=xlsread(['/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/',myrcp_,'_MIDYEAR_RADFORCING.xls'],[myrcp_,'_MIDYEAR_RADFORCING']);
myyear=data(:,1);
ind=find(myyear<=2100);
myyear=myyear(ind);
RF_OZO=squeeze(data(ind,50)+data(ind,51));
RF_AER=squeeze(data(ind,42))+squeeze(data(ind,49));
mRF_OZO=interp1(myyear,RF_OZO,years);
mRF_AER=interp1(myyear,RF_AER,years);
mRF_OZO(isnan(mRF_OZO))=0;
mRF_AER(isnan(mRF_AER))=0;

data=xlsread(['/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/',myrcp_,'_MIDYR_CONC.xlsx']);
data(1:36,:)=[];
data(:,1)=[];
conc_Hal=interp1(myyear,data(ind,8),years);
conc_Hal(isnan(conc_Hal))=0;


% scale forcing to 0.35W/m2
% scale forcing to 0.4W/m2 (updated)
meinhausenRF_OZO=mRF_OZO; rfd=mean(meinhausenRF_OZO(ind2012))-mean(meinhausenRF_OZO(ind1870));
meinhausenRF_OZO=(meinhausenRF_OZO-mean(meinhausenRF_OZO(ind1870)))*0.4/rfd+mean(meinhausenRF_OZO(ind1870)) ;
meinhausenRF_AER=mRF_AER;

if strcmp(myrcp,'COM20')
    ind=find(years==2020);
    emissionN2O(ind:end)=0;
    emissionCH4(ind:end)=0;
    emissionCO2(ind:end)=0;
    emissionSO2(ind:end)=0;
    meinhausenRF_OZO(ind:end)=0;
    meinhausenRF_AER(ind:end)=0;
end
if strcmp(myrcp,'COM30')
    ind=find(years==2020);
    ind2=find(years==2030);
    emissionN2O(ind:end)=emissionN2O(ind);
    emissionCH4(ind:end)=emissionCH4(ind);
    emissionCO2(ind:end)=emissionCO2(ind);
    emissionSO2(ind:end)=emissionSO2(ind);
    meinhausenRF_OZO(ind:end)=meinhausenRF_OZO(ind);
    meinhausenRF_AER(ind:end)=meinhausenRF_AER(ind);
    emissionN2O(ind2:end)=0;
    emissionCH4(ind2:end)=0;
    emissionCO2(ind2:end)=0;
    emissionSO2(ind2:end)=0;
    meinhausenRF_OZO(ind2:end)=0;
    meinhausenRF_AER(ind2:end)=0;
end
if strcmp(myrcp,'COM40')
    ind=find(years==2020);
    ind2=find(years==2040);
    emissionN2O(ind:end)=emissionN2O(ind);
    emissionCH4(ind:end)=emissionCH4(ind);
    emissionCO2(ind:end)=emissionCO2(ind);
    emissionSO2(ind:end)=emissionSO2(ind);
    meinhausenRF_OZO(ind:end)=meinhausenRF_OZO(ind);
    meinhausenRF_AER(ind:end)=meinhausenRF_AER(ind);
    emissionCH4(ind2:end)=0;
    emissionN2O(ind2:end)=0;
    emissionCO2(ind2:end)=0;
    emissionSO2(ind2:end)=0;
    meinhausenRF_OZO(ind2:end)=0;
    meinhausenRF_AER(ind2:end)=0;
end





% %% test code
% clear;close all
% startyear=1750;
% endyear=2100;
% 
% for s=5
%     scen={'RCP85','RCP45','RCP3','RCP6','COM20'}
%     myrcp=scen{s};
%     
%     
%     
%     [years,meinhausenRF_OZO,TSI,alb,emissionCH4,emissionCO2,emissionSO2,emission_volc,radf_other,radf_halo,radf_N2O,radf_CO2,radf_CH4,conc_CO2,conc_CH4] = load_RCP(myrcp,startyear,endyear);
%     figure(1)
%     subplot(2,2,1)
%     plot(years,emissionCO2);hold on
%     subplot(2,2,2)
%     plot(years,emissionCH4);hold on
%     subplot(2,2,3)
%     plot(years,emissionSO2);hold on
%     subplot(2,2,4)
%     plot(years,meinhausenRF_OZO);hold on
%     
%     figure(2)
%     subplot(2,2,1)
%     plot(years,conc_CO2);hold on
%     subplot(2,2,2)
%     plot(years,conc_CH4);hold on
%     subplot(2,2,3)
%     plot(years,radf_CO2);hold on
%     subplot(2,2,4)
%     plot(years,radf_CH4);hold on
%     pause
%     
% end
