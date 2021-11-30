
load(fnC);
figure(1)
subplot(1,2,1)
plot(0:1/12:100,OPT_CO2_100yr_ts);title('CO2 decay following emission termination after 1000Gt emission')
subplot(1,2,2)
plot(0:1/12:100,OPT_CO2fraction_100yr_ts);title('Associated CO2 fraction')

csvwrite('CO2_change.csv',OPT_CO2_100yr_ts');
csvwrite('CO2_fraction.csv',OPT_CO2fraction_100yr_ts');

% Large parameter sweep. Simulations where the 1850-2020 CO2 RMSE was
% <6ppm and the 1861-1880 ro 2006-2018 CO2 change was within +/-1ppm of

clear OPT_CO2* OPT_rmsE OPT_best*
var=who('OPT_*')

disp('Carbon model parameters')
for n=1:length(var)
    V=var{n};
    eval(['disp([''',V(5:end),','',num2str(',V,')])'])
end

disp('OM, 7.8e22  % moles of water in ocean')
disp('AM, 1.77e20 % moles of air in atmosphere')
disp('Alk, 767 %GtC   % Alkalinity; assumed constant here but buffered on long timescales O[10ka] from the dissolution of CaCO3')
disp('kh, 1.0548e+03  % at 15oC; ratio of the molar concentrations of CO2 in atmosphere and ocean (Henrys Law)')
disp('k1, 8.7184e-07  % at 15oC; dissociation constant')
disp('k2, 5.4426e-10  % at 15oC; dissociation constant')


disp('Volcanic aerosol parameters')
disp('vtau, 1.2   % yr Volcanic aerosol decay timescale')

disp('Albedo parameters')
disp('alb0, 0.31   % initial albedo')

disp('Aerosol parameters')
disp('VF, -20     % conversion from optical thickness of volcanic aerosols to RF')
disp('AF, -0.0052 % conversion from SO2 emissions to RF W/m2/TgSO2/yr')

disp('N2O parameters')
disp('tau_N20, 109 %yrs')
disp('alphaRF_N2O, 0.12;')

disp('Holaogen parameters')
disp('tau_CFC12, 102 %yrs')
disp('alphaRF_CFC12, 0.32')

disp('CH4 parameters')
disp('alphaRF_ch4, 0.0316')
disp('tau_ch4_pi,  7.8865 %yrs')
disp('alpha_ch4,  -0.1538')

figure(101);clf
for n=1:4
    plot(yearall,squeeze(ALL_bestC_CO2(ALLsubset(n),:,:))','r')
    hold on
end
for s=1:4
    switch s
        case 1
            myrcp='RCP3';
        case 2
            myrcp='RCP45';
        case 3
            myrcp='RCP6';
        case 4
            myrcp='RCP85';
    end
    [years,meinhausenRF_OZO,meinhausenRF_AER,TSI,alb,emissionCH4,emissionCO2,emissionN2O,emissionSO2,emission_volc,radf_other,radf_halo,radf_N2O,radf_CO2,radf_CH4,conc_CO2,conc_CH4,conc_N2O,mTSI,conc_Hal] = load_RCP(myrcp,1850,2100);
    plot(years,conc_CO2,'k--')
    hold on
end
set(gca,'fontsize',16)
title('CO2 concentratio, for 4 subset models (red) and RCP data (black)')


clear ALL* OPT*

print -f1 -depsc 'images/subset_carbon_models'
print -f101 -depsc 'images/subset_carbon_models_RCP'

%% Radiative forcing compared to AR6
data=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/AR6_ERF_1750-2019.xls');
AR6_year=data(:,1);
AR6_ERF_co2(1,:)=data(:,2);
AR6_ERF_ch4(1,:)=data(:,3);
AR6_ERF_n2o(1,:)=data(:,4);
AR6_ERF_other_wmghg(1,:)=data(:,5);
AR6_ERF_o3(1,:)=data(:,6);
AR6_ERF_h2o_stratospheric(1,:)=data(:,7);
AR6_ERF_contrails(1,:)=data(:,8);
AR6_ERF_aerosol_radiation_interactions(1,:)=data(:,9);
AR6_ERF_aerosol_cloud_interactions(1,:)=data(:,10);
AR6_ERF_bc_on_snow(1,:)=data(:,11);
AR6_ERF_land_use(1,:)=data(:,12);
AR6_ERF_volcanic(1,:)=data(:,13);
AR6_ERF_solar(1,:)=data(:,14);
AR6_ERF_nonco2_wmghg(1,:)=data(:,15);
AR6_ERF_aerosol(1,:)=data(:,16);
AR6_ERF_chapter2_other_anthro(1,:)=data(:,17);
AR6_ERF_total_anthropogenic(1,:)=data(:,18);
AR6_ERF_total_natural(1,:)=data(:,19);
AR6_ERF_total(1,:)=data(:,20);

clear data
data(1,:,:)=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/RCP3PD_MIDYEAR_RADFORCING.xls','RCP3PD_MIDYEAR_RADFORCING');
data(2,:,:)=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/RCP45_MIDYEAR_RADFORCING.xls','RCP45_MIDYEAR_RADFORCING');
data(3,:,:)=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/RCP6_MIDYEAR_RADFORCING.xls','RCP6_MIDYEAR_RADFORCING');
data(4,:,:)=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/RCP85_MIDYEAR_RADFORCING.xls','RCP85_MIDYEAR_RADFORCING');
myyear=data(1, :,1);
ind=find(myyear>=1850 & myyear<=2100);
RCP_year=myyear(ind);
RCP_RF_TOT=squeeze(data(:,ind,2));
RCP_RF_SOL=squeeze(data(:,ind,4));
RCP_RF_VOLC=squeeze(data(:,ind,3));
RCP_RF_CO2=squeeze(data(:,ind,9));
RCP_RF_CH4=squeeze(data(:,ind,10));
RCP_RF_AER=squeeze(data(:,ind,42))+squeeze(data(:,ind,49)); %49 is cloud albedo effect
RCP_RF_OZO=squeeze(data(:,ind,50)+data(:,ind,51));
RCP_RF_ALB=squeeze(data(:,ind,53)+data(:,ind,49));
RCP_RF_TOT_minAER=RCP_RF_TOT-RCP_RF_AER;




figure(2);clf
ind1870=find(yearall>=1861 & yearall<1880);
ind2012=find(yearall>=2006 & yearall<2018);

for myscen =1:4
    carbonator_tot_min_aerosol=squeeze(projdat_RCO2(myscen,:,:))+projdat_RCH4+projdat_RN2O+projdat_RCOZO+projdat_RHalo+projdat_RSOL+projdat_RVOLC;
    plot(yearall,carbonator_tot_min_aerosol(1:4,:)-mean(carbonator_tot_min_aerosol(1,ind1870)),'r')
    hold on
end
plot(yearall,RCP_RF_TOT_minAER-mean(RCP_RF_TOT_minAER(1,ind1870)),'k')

ind1870ar6=find(AR6_year>=1861 & AR6_year<1880);
ind2012ar6=find(AR6_year>=2006 & AR6_year<2018);
plot(AR6_year,AR6_ERF_total-AR6_ERF_aerosol-mean(AR6_ERF_total(ind1870ar6)-AR6_ERF_aerosol(ind1870ar6)),'b','linewidth',2)
xlim([1850 2100])
set(gca,'fontsize',16)
title('Total radiative forcing, excluding aerosols')
legend('Carbonator 1','Carbonator 2','Carbonator 3','Carbonator 4','RCP Meinshausen2011','AR6 estimate')

print -f2 -depsc 'images/RF_noAer'

% based on RCP85 additional years
tmp=squeeze(mean(projdat_RCO2(:,4,ind2012),3))-squeeze(mean(projdat_RCO2(:,4,ind1870),3));tmp2=mean(AR6_ERF_co2(ind2012ar6))-mean(AR6_ERF_co2(ind1870ar6));  disp(['CO2: ',num2str(range(tmp(:))),' AR6: ',num2str(tmp2)])
tmp=squeeze(mean(projdat_RCH4(4,ind2012),2))-squeeze(mean(projdat_RCH4(4,ind1870),2));    tmp2=mean(AR6_ERF_ch4(ind2012ar6))-mean(AR6_ERF_ch4(ind1870ar6));  disp(['CH4: ',num2str(tmp),' AR6: ',num2str(tmp2)])
tmp=squeeze(mean(projdat_RN2O(4,ind2012),2))-squeeze(mean(projdat_RN2O(4,ind1870),2));    tmp2=mean(AR6_ERF_n2o(ind2012ar6))-mean(AR6_ERF_n2o(ind1870ar6));  disp(['N2O: ',num2str(tmp),' AR6: ',num2str(tmp2)])
tmp=squeeze(mean(projdat_RCOZO(4,ind2012),2))-squeeze(mean(projdat_RCOZO(4,ind1870),2));    tmp2=mean(AR6_ERF_o3(ind2012ar6))-mean(AR6_ERF_o3(ind1870ar6));  disp(['O3: ',num2str(tmp),' AR6: ',num2str(tmp2)])
tmp=squeeze(mean(projdat_RHalo(4,ind2012),2))-squeeze(mean(projdat_RHalo(4,ind1870),2));    tmp2=mean(AR6_ERF_other_wmghg(ind2012ar6))-mean(AR6_ERF_other_wmghg(ind1870ar6));  disp(['Halogen: ',num2str(tmp),' AR6: ',num2str(tmp2)])
tmp=squeeze(mean(projdat_RSOL(4,ind2012),2))-squeeze(mean(projdat_RSOL(4,ind1870),2));    tmp2=mean(AR6_ERF_solar(ind2012ar6))-mean(AR6_ERF_solar(ind1870ar6));  disp(['Solar: ',num2str(tmp),' AR6: ',num2str(tmp2)])
tmp=squeeze(mean(projdat_RVOLC(4,ind2012),2))-squeeze(mean(projdat_RVOLC(4,ind1870),2));    tmp2=mean(AR6_ERF_volcanic(ind2012ar6))-mean(AR6_ERF_volcanic(ind1870ar6));  disp(['Volcanic: ',num2str(tmp),' AR6: ',num2str(tmp2)])

figure(201);clf
cols={'g','b','k','r'}
for n=1:4
    plot(yearall,squeeze(prctile(projdat_RCSO2(:,n,:),[50])),cols{n},'linewidth',2)
    hold on
    plot(yearall,RCP_RF_AER(n,:),cols{n})
end
set(gca,'fontsize',16)

figure(202);clf
cols={'g','b','k','r'}
for n=1:4
    plot(yearall,squeeze(mean(projdat_RCSO2(:,n,:))),cols{n},'linewidth',2)
    hold on
    plot(yearall,RCP_RF_AER(n,:),cols{n})
end
set(gca,'fontsize',16)






%% Distribution of critical parameters
figure(3);clf
for v=1:6
    switch v
        case 1;
            var='g';
        case 2;
            var='Co';
        case 3;
            var='Cs';
        case 4;
            var='climate_sensitivity';
        case 5;
            var='deltaT';
        case 6;
            var='RFA';
    end
    subplot(3,2,v)
    eval(['tmp=alldat_',var,';']) % all simulations
    [N,X]=hist(tmp,100);
    plot(X,N/max(N),'r')
    hold on
    eval(['tmp=dat_',var,';']) % simulations that meet dT criteria
    [N,X]=hist(tmp,100);
    hold on; plot(X,N/max(N),'k')
    
    eval(['tmp=datJD_',var,';']) % JD simulations
    [N,X]=hist(tmp,100);
    hold on; plot(X,N/max(N),'b')
    title(var)
end

%%
[sw_ecs,sw_arf]=sherwood_ecs_arf(20000);
load_2d_histogram;
dhistogram=dhistogram';
ecs=0.1:.2:7.9;
arf=-2.9:.2:.1;
dhistogram=dhistogram/max(dhistogram(:));
dhistogram=reshape(dhistogram,42,16);
dhistogram(41:42,:)=[];

figure(4);clf
sfh2=subplot(2,2,2)
% Red All simulations that meet dT
% Blak simulations that follow JD
% Blue sherwood JD
scatter(dat_RFA, dat_climate_sensitivity,3,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
hold on
scatter(datJD_RFA,datJD_climate_sensitivity,3,'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
contour(arf,ecs,dhistogram,[0:.1:1],'b','linewidth',1.5)

xlabel('Aerosol RF')
ylabel('ECS')
xlim([-2 0.25]);ylim([1 8])
sfh1=subplot(2,2,1)
cla
[N,X]=hist(datJD_climate_sensitivity,1:.1:8); plot(N/sum(N),X,'k')
hold on
[N,X]=hist(sw_ecs,1:.1:8); plot(N/sum(N),X,'b--')
[N,X]=hist(dat_climate_sensitivity,1:.1:8); plot(N/sum(N),X,'r')

sfh4=subplot(2,2,4)
cla
[N,X]=hist(datJD_RFA,-2:.1:0.5); plot(X,N/sum(N),'k')
hold on
[N,X]=hist(sw_arf,-2:.1:0.5); plot(X,N/sum(N),'b--')
[N,X]=hist(dat_RFA,-2:.1:0.5); plot(X,N/sum(N),'r')

pos1 = get(sfh1, 'Position')
pos2 = get(sfh2, 'Position')
pos3 = get(sfh4, 'Position')
pos1_=[0.05    0.3    0.2    0.6]
pos2_=[0.3    0.3    0.6    0.6]
pos4_=[.3    0.05    0.6    0.2]
set(sfh1, 'Position',pos1_ )
set(sfh2, 'Position',pos2_ )
set(sfh4, 'Position',pos4_ )

print(['images/distributions_',num2str(test)], '-f4', '-djpeg100', '-r250' )


%% JD colored by parameter
for CP=1:length(unique(datJD_Cparam))
    ind=find(datJD_Cparam==CP);
    figure(40+CP);clf
    subplot(2,2,1)
    scatter(datJD_RFA(ind), datJD_climate_sensitivity(ind),20,datJD_g(ind),'filled')
    hold on
    contour(arf,ecs,dhistogram,[0:.1:1],'k')
    xlabel('Aerosol RF')
    ylabel('ECS')
    xlim([-1.5 0.1]);ylim([1.5 5])
    colorbar
    title('g');set(gca,'fontsize',15)
    subplot(2,2,2)
    scatter(datJD_RFA(ind), datJD_climate_sensitivity(ind),20,datJD_Co(ind),'filled')
    hold on
    contour(arf,ecs,dhistogram,[0:.1:1],'k')
    xlabel('Aerosol RF')
    ylabel('ECS')
    xlim([-1.5 0.1]);ylim([1.5 5])
    colorbar
    title('Co');set(gca,'fontsize',15)
    
    subplot(2,2,3)
    scatter(datJD_RFA(ind), datJD_climate_sensitivity(ind),20,datJD_Cs(ind),'filled')
    hold on
    contour(arf,ecs,dhistogram,[0:.1:1],'k')
    xlabel('Aerosol RF')
    ylabel('ECS')
    xlim([-1.5 0.1]);ylim([1.5 5])
    colorbar
    title('Cs');set(gca,'fontsize',15)
    
    subplot(2,2,4)
    plot([0.5 0.5],[0 5])
    hold on
    text(1,1,['ECS: ',num2str(min(alldat_climate_sensitivity)),' - ',num2str(max(alldat_climate_sensitivity))])
    
    text(1,2,['g: ',num2str(min(alldat_g)),' - ',num2str(max(alldat_g))])
    text(1,3,['Co: ',num2str(min(alldat_Co)),' - ',num2str(max(alldat_Co))])
    text(1,4,['Cs: ',num2str(min(alldat_Cs)),' - ',num2str(max(alldat_Cs))])
    xlim([0 3]);ylim([0 5]);set(gca,'fontsize',15)
end

print(['images/col_distributions_',num2str(test)], '-f41', '-djpeg100', '-r250' )

%% Parameter values for JD simulations
figure(4);clf
subplot(2,2,1)
[N,X]=hist(datJD_g,0:.01:3);
bar(X,N/sum(N))
hold on
plot([range_g(1) range_g(2)],[1 1]*max(N/sum(N))*1.05,'r','linewidth',2)
subplot(2,2,2)
[N,X]=hist(datJD_Co,0:300);
bar(X,N/sum(N))
hold on
plot([range_Co(1) range_Co(2)],[1 1]*max(N/sum(N))*1.05,'r','linewidth',2)
subplot(2,2,3)
[N,X]=hist(datJD_Cs,0:.1:20);
bar(X,N/sum(N))
hold on
plot([range_Cs(1) range_Cs(2)],[1 1]*max(N/sum(N))*1.05,'r','linewidth',2)
hold on
print(['images/PDF_g_Co_Cs_',num2str(test)], '-f4', '-djpeg100', '-r250' )
subplot(2,2,4)
[N,X]=hist(datJD_Cparam,0.5:7.5);
bar(X,N)
title('Carbon model parameter set')



figure(41);clf
cols={'g','b','k','r'}
ind1870=find(yearall>=1861 & yearall<1880);
for myscen=1:4
    for CM=1:4
        plot(yearall,squeeze(projdat_RCO2(CM,myscen,:))-mean(projdat_RCO2(myscen,ind1870)),cols{myscen})
        hold on
    end
end
ind1870=find(AR6_year>=1861 & AR6_year<1880);
plot(AR6_year,AR6_ERF_co2-mean(AR6_ERF_co2(ind1870)))




%% RCP temperature distributions
figure(30);clf
col={'b','k','g','r'} ;
subplot(2,1,1)
for s=1:4
    plot(yearall,prctile(squeeze(projdat_Ts(:,s,:)),[10 50 90],1),col{s})
    hold on
end
xlim([1850 2100])

subplot(2,1,2)
ind0=find(yearall>=1861 & yearall<=1880);
ind1=find(yearall>=2006 & yearall<=2018);
ind2=find(yearall>=2080 & yearall<=2100);

tmp0=squeeze(mean(projdat_Ts(:,:,ind0),3));
tmp1=squeeze(mean(projdat_Ts(:,:,ind1),3));
tmp2=squeeze(mean(projdat_Ts(:,:,ind2),3));

for s=1:4
    tmp=tmp2(:,s)-tmp0(:,s)
    [N2,X]=hist(tmp,0:.1:7);
    plot(0:.1:7,N2,col{s},'linewidth',2)
    hold on
    plot([1 1]*median(tmp),[0 800],col{s},'linewidth',1)
    
end
print(['images/RCP_SAT_1860to2090_',num2str(test)], '-f30', '-djpeg100', '-r250' )

%% RCP temperature distributions relative to 1861-1880



projdat_Ts_1870=projdat_Ts*NaN;
for n=1:length(projdat_Ts)
    n
    for s=1:7
        projdat_Ts_1870(n,s,:)=projdat_Ts(n,s,:)-squeeze(mean(projdat_Ts(n,s,ind0),3));
    end
end
figure(5);clf
col={'b','k','g','r'} ;
subplot(2,1,1)
for s=1:4
    plot(yearall,prctile(squeeze(projdat_Ts_1870(:,s,:)),[10 50 90],1),col{s})
    hold on
end
set(gca,'fontsize',15)
xlim([1850 2100])
title('Global SAT anomaly relatibe to 1861-1880')
pbaspect([1.3 1 1])
subplot(2,1,2)
ind2=find(yearall>=2080 & yearall<=2100);
tmp2=squeeze(mean(projdat_Ts_1870(:,:,ind2),3));
for s=1:4
    tmp=tmp2(:,s);
    [N2,X]=hist(tmp,0:.1:7);
    plot(0:.1:7,N2,col{s},'linewidth',2)
    hold on
    plot([1 1]*median(tmp),[0 8000],col{s},'linewidth',1)
    
end
xlim([1 6])
set(gca,'fontsize',15)
pbaspect([1.3 1 1])
title('Global SAT anomaly 2080-2100 minus 1861-1880')

print -depsc -f5 'RCP_SAT_pdfs'

%% % sanity check of recen SST trend

figure(33);clf

for s=1
    plot(yearall,prctile(squeeze(projdat_Ts_1870(:,s,:)),[5 50 95],1),col{s})
    hold on
    
end
load('/Users/z3045790/CCRC/DATASETS/global_SATA_SSTA/global_TAS.mat')

ind=find(dyr_HC5>=1861 & dyr_HC5<=1880);plot(dyr_HC5,TAS_HC5-mean(TAS_HC5(ind)),'k')
% ind=find(dyr_GT4>=1861 & dyr_GT4<=1880);plot(dyr_GT4,TAS_GT4,'g')
ind=find(dyr_BE>=1861 & dyr_BE<=1880);plot(dyr_BE,TAS_BE-mean(TAS_BE(ind)),'r')
xlim([1850 2020]);ylim([-1 1.5]);pbaspect([1.5 1 1])
print -f33 -depsc 'images/SAT_carbonator_HC_BE'

figure(34);clf
ind=find(dyr_BE>=2010 & dyr_BE<=2020)
plot(dyr_BE(ind),TAS_BE(ind));hold on
plot_polyfit(dyr_BE(ind),TAS_BE(ind),'r')
trend(TAS_BE(ind),dyr_BE(ind))*10

ind=find(dyr_BE>=2000 & dyr_BE<=2020)
plot(dyr_BE(ind),TAS_BE(ind));hold on
plot_polyfit(dyr_BE(ind),TAS_BE(ind),'r--')
trend(TAS_BE(ind),dyr_BE(ind))*10


%%
clear data
data(1,:,:)=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/RCP3PD_MIDYEAR_RADFORCING.xls','RCP3PD_MIDYEAR_RADFORCING');
data(2,:,:)=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/RCP45_MIDYEAR_RADFORCING.xls','RCP45_MIDYEAR_RADFORCING');
data(3,:,:)=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/RCP6_MIDYEAR_RADFORCING.xls','RCP6_MIDYEAR_RADFORCING');
data(4,:,:)=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/RCP85_MIDYEAR_RADFORCING.xls','RCP85_MIDYEAR_RADFORCING');
myyear=data(1, :,1);
ind=find(myyear>=1850 & myyear<=2100);
myyear=myyear(ind);
RF_TOT=squeeze(data(:,ind,2));
RF_SOL=squeeze(data(:,ind,4));
RF_VOLC=squeeze(data(:,ind,3));
RF_CO2=squeeze(data(:,ind,9));
RF_CH4=squeeze(data(:,ind,10));
RF_AER=squeeze(data(:,ind,42))+squeeze(data(:,ind,49));
RF_OZO=squeeze(data(:,ind,50)+data(:,ind,51));
RF_ALB=squeeze(data(:,ind,53)+data(:,ind,49));
RF_OTHER=RF_TOT-(RF_SOL+RF_VOLC+RF_CO2+RF_CH4+RF_AER+RF_OZO);

ind1870=find(myyear>=1861 & myyear<1880);
clear R*s
for s=1:4
    RF_TOTi(s,:)=interp1(myyear,RF_TOT(s,:),yearall)-mean(RF_TOT(s,ind1870));
    RF_SOLi(s,:)=interp1(myyear,RF_SOL(s,:),yearall)-mean(RF_SOL(s,ind1870));
    RF_VOLCi(s,:)=interp1(myyear,RF_VOLC(s,:),yearall)-mean(RF_VOLC(s,ind1870));
    RF_CO2i(s,:)=interp1(myyear,RF_CO2(s,:),yearall)-mean(RF_CO2(s,ind1870));
    RF_CH4i(s,:)=interp1(myyear,RF_CH4(s,:),yearall)-mean(RF_CH4(s,ind1870));
    RF_AERi(s,:)=interp1(myyear,RF_AER(s,:),yearall)-mean(RF_AER(s,ind1870));
    RF_OZOi(s,:)=interp1(myyear,RF_OZO(s,:),yearall)-mean(RF_OZO(s,ind1870));
    RF_ALBi(s,:)=interp1(myyear,RF_ALB(s,:),yearall)-mean(RF_ALB(s,ind1870));
    RF_OTHERi(s,:)=interp1(myyear,RF_OTHER(s,:),yearall)-mean(RF_OTHER(s,ind1870));
end
figure(40);clf
col={'b','k','y','r'};
cold={'b--','k--','g--','r--'} ;

for s=1:4
    plot(yearall,prctile(squeeze(projdat_RF(:,s,:)),[10 50 90],1),col{s})
    hold on
    plot(yearall,RF_TOTi(s,:),cold{s})
end
xlim([1850 2100])



% %%
% C_CH4pi=719;
% figure(41);clf
% figure(42);clf
% figure(43);clf
% figure(44);clf
% figure(45);clf
% for s=1:4
%     switch s
%         case 1
%             myrcp='RCP3';
%         case 2
%             myrcp='RCP45';
%         case 3
%             myrcp='RCP6';
%         case 4
%             myrcp='RCP85';
%     end
%     [years,meinhausenRF_OZO,meinhausenRF_AER,TSI,alb,emissionCH4,emissionCO2,emissionSO2,emission_volc,radf_other,radf_halo,radf_N2O,radf_CO2,radf_CH4,conc_CO2,conc_CH4,mTSI] = load_RCP(myrcp,startyear,endyear);
%
%     figure(41)
%     plot(years,radf_N2O,'k')
%     hold on
%     plot(years,radf_halo,'r')
%     plot(years,radf_other,'b')
%     plot(yearall,RF_AERi(s,:)+RF_OZOi(s,:),'b--')
%
%     set(gca,'fontsize',14)
%
%     figure(42);
%     subplot(2,2,1)
%     plot(yearall,projdat_CCO2(s,:),'k')
%     hold on
%     plot(years,conc_CO2,'k--');title('concentration CO2')
%     set(gca,'fontsize',14);xlim([1850 2100])
%     subplot(2,2,3)
%     plot(yearall,projdat_RCO2(s,:),'k')
%     hold on
%     plot(years,radf_CO2,'k--');title('RF CO2')
%     plot(yearall,RF_CO2i(s,:),'r--')
%     set(gca,'fontsize',14);xlim([1850 2100])
%     subplot(2,2,2)
%     plot(yearall,C_CH4pi+projdat_CCH4(s,:),'k')
%     hold on
%     plot(years,conc_CH4,'k--');title('concentration CH4')
%     set(gca,'fontsize',14);xlim([1850 2100])
%     subplot(2,2,4)
%     plot(yearall,projdat_RCH4(s,:),'k')
%     hold on
%     plot(years,radf_CH4,'k--');title('RF CO2')
%     plot(yearall,RF_CH4i(s,:),'r--')
%     set(gca,'fontsize',14);xlim([1850 2100])
%
%     figure(43);
%     plot(yearall,projdat_RCOZO(s,:),'k')
%     hold on
%     plot(yearall,RF_OZOi(s,:),'k--')
%     set(gca,'fontsize',14);title('all ozone')
%
%     figure(44)
%     subplot(2,1,1)
%     plot(yearall,prctile(squeeze(projdat_RCSO2(:,s,:)),[50],1),'k')
%     hold on
%     plot(yearall,RF_AERi(s,:),'k--');grid on;xlim([1850 2100])
%     subplot(2,1,2)
%     plot(years,emissionSO2);
%     hold on
%     grid on;xlim([1850 2100])
%
%     %     figure(45)
%     %     plot(yearall,projdat_RCOZO(s,:)+prctile(squeeze(projdat_RCSO2(:,s,:)),[50],1),'k')
%     %     hold on
%     %     plot(years,radf_other,'k--')
%
% end
% xlim([1750 2100])
%
% figure(41);legend('NO2','HALO','Aerosol/Ozone etc')
%
% % print -f42 -depsc 'images/all_co2_RF_modelVSmeinhausen'
cols={'g','b','k','r','b--','k--','r--'}
figure(8);clf
for myscen=1:7
    subplot(2,3,1)
    plot(yearall,squeeze(projdat_RCO2(1,myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('CO2')
    subplot(2,3,2)
    plot(yearall,squeeze(projdat_RCH4(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('CH4')
    subplot(2,3,3)
    plot(yearall,squeeze(projdat_RN2O(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('N2O')
    subplot(2,3,4)
    plot(yearall,squeeze(projdat_RCOZO(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('O3')
    subplot(2,3,5)
    plot(yearall,squeeze(projdat_RHalo(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('Halogens')
    subplot(2,3,6)
    plot(yearall,squeeze(median(projdat_RCSO2(:,myscen,:))),cols{myscen});hold on;set(gca,'fontsize',15);title('Aerosol (median)')
end

print -f8 -depsc 'all_radiative_forcing'

%% Committment
ind2020=find(yearall==2020)
figure(100);clf
SATcommit=projdat_Ts(:,1:3,:)*NaN;
SATmax=projdat_Ts(:,1:3)*NaN;
SATmax_year=projdat_Ts(:,1:3)*NaN;
for n=1:length(projdat_Ts)
    n
    for s=5:7
        SATcommit(n,s-4,:)= (projdat_Ts(n,s,:)-projdat_Ts(n,s,ind2020));
        SATmax(n,s-4)= max(projdat_Ts(n,s,:));
        ind=find(projdat_Ts(n,s,:)==SATmax(n,s-4));SATmax_year(n,s-4)=yearall(ind);
        SATmax2020(n,s-4)= max(projdat_Ts(n,s,:))-projdat_Ts(n,s,ind2020);
    end
end

n=1;
scen=5;
figure(100);clf
subplot(2,1,1)
plot(yearall(ind2020:end),prctile(squeeze(SATcommit(:,scen-4,ind2020:end)),[10 50 90],1))
subplot(2,1,2)
plot(yearall(ind2020:end),prctile(squeeze(projdat_Ts(:,scen,ind2020:end)),[10 50 90],1))


figure(105);clf
plot(yearall(ind2020-10*12:end),prctile(squeeze(SATcommit(:,scen-4,ind2020-10*12:end)),[5 50 95],1),'r','linewidth',2)
xlim([2010 2100])
set(gca,'ytick',[-.5:.25:1],'fontsize',16)
ylim([-.6 .75])
grid on
print -f105 -depsc 'commitSAT_overlay'

%%
SATtmp=squeeze(SATcommit(:,1,ind2020:end));
yeartmp=yearall(ind2020:end);
[M,I] =max(SATtmp,[],2);

figure(110);clf
subplot(2,1,1)
[N,X]=hist(M,-.1:.025:3);hold on
plot(X,N,'k');hold on
xlabel('Max commited warming relative to 2020')
set(gca,'fontsize',16)
subplot(2,1,2)
[N,X]=hist(yeartmp(I),2019.5:2100.5);
plot(X,N,'k');hold on
set(gca,'fontsize',16)
xlabel('Year of Max commited warming ')

SATtmp=squeeze(SATcommit(:,2,ind2020:end));
yeartmp=yearall(ind2020:end);
[M,I] =max(SATtmp,[],2);
subplot(2,1,1)
[N,X]=hist(M,-.1:.025:3);
plot(X,N,'b')
subplot(2,1,2)
[N,X]=hist(yeartmp(I),2019.5:2100.5)
plot(X,N,'b');

SATtmp=squeeze(SATcommit(:,3,ind2020:end));
yeartmp=yearall(ind2020:end);
[M,I] =max(SATtmp,[],2);
xlim([0 1.5])
subplot(2,1,1)
[N,X]=hist(M,-.1:.025:3);
plot(X,N,'r')
subplot(2,1,2)
[N,X]=hist(yeartmp(I),2019.5:2100.5)
plot(X,N,'r');
xlim([2020 2060])


figure(120);clf
[N,X]=hist(SATmax(:,1),0.9:.05:4);
plot(X,N/sum(N),'k');
hold on
[N,X]=hist(SATmax(:,2),0.9:.05:4); plot(X,N/sum(N),'b')
[N,X]=hist(SATmax(:,3),0.9:.05:4); plot(X,N/sum(N),'r')
xlim([0.9 2.5])


figure(1200);clf
for CP=1:4
    ind=find(datJD_Cparam==CP);
    [N,X]=hist(SATmax(ind,1),0.9:.05:4);
    plot(X,N/sum(N),'k');
    hold on
    [N,X]=hist(SATmax(ind,2),0.9:.05:4); plot(X,N/sum(N),'b')
    [N,X]=hist(SATmax(ind,3),0.9:.05:4); plot(X,N/sum(N),'r')
    xlim([0.9 2.5])
end
pbaspect([2 1 1])




figure(121);clf
for CP=1:4
    ind=find(datJD_Cparam==CP);
    [N,X]=hist(SATmax(ind,1),0.9:.05:4); plot(X,cumsum(N)/sum(N),'k')
    hold on
    [N,X]=hist(SATmax(ind,2),0.9:.05:4); plot(X,cumsum(N)/sum(N),'b')
    [N,X]=hist(SATmax(ind,3),0.9:.05:4); plot(X,cumsum(N)/sum(N),'r')
end
xlim([0.9 2.5])
grid on

print -f120 -depsc 'commitSAT_pdf'
print -f1200 -depsc 'commitSAT_pdf_4CP'
print -f121 -depsc 'commitSAT_cumpdf_4CP'



%%


figure(111);clf
coldc={'b','k','g','r','b--','k--','r--'} ;
for dr=1:7
    plot(yearall(ind2020-10*12:end),prctile(squeeze(SATcommit(find(datJD_Cparam==dr),scen-4,ind2020-10*12:end)),[5 50 95],1),coldc{dr},'linewidth',1)%fast decay
    hold on
end
xlim([2010 2100])
set(gca,'ytick',[-.5:.25:1],'fontsize',16)
ylim([-.6 1])
grid on
title('SAT relative to 2020')

figure(112);clf
for dr=1:7
    SATtmp=squeeze(SATcommit(find(datJD_Cparam==dr),1,ind2020:end));
    yeartmp=yearall(ind2020:end);
    [M,I] =max(SATtmp,[],2);
    [N,X]=hist(M,0:.1:2);
    plot(X,N,coldc{dr})
    hold on
end

figure(113);clf
vars={'g','Co','Cs','climate_sensitivity','RFA'};
for var=1:length(vars)
    subplot(2,3,var)
    eval(['tmp=datJD_',vars{var},';'])
    xx=linspace(min(tmp),max(tmp),20);
    for dr=1:7
        ind=find(datJD_Cparam==dr);
        disp([vars{var},': ',num2str(median(tmp(ind)))]);
        [N,X]=hist(tmp(ind),xx);
        plot(X,N/sum(N),coldc{dr})
        hold on
        plot([1 1]*median(tmp(ind)),[0 1]*max(N)/sum(N),coldc{dr})
    end
    title(vars{var})
end

print -f111 -depsc 'comitted_warming_4CMparameters';





%% RF
ind1870=find(yearall>=1861 & yearall<1880);
ind2012=find(yearall>=2006 & yearall<2018);
projdat_RCSO2_bs=projdat_RCSO2*NaN;
projdat_RF_bs=projdat_RCSO2*NaN;
RF_AERi_bs=RF_AERi*NaN;
for n=1:length(projdat_RCSO2)
    n
    for s=1:7
        projdat_RCSO2_bs(n,s,:)=projdat_RCSO2(n,s,:)-squeeze(mean(projdat_RCSO2(n,s,ind1870),3));
        projdat_RF_bs(n,s,:)=projdat_RF(n,s,:)-squeeze(mean(projdat_RF(n,s,ind1870),3));
    end
end
for s=1:4
    RF_AERi_bs(s,:)=RF_AERi(s,:)-mean(RF_AERi(s,ind1870),2);
    RF_TOTi_bs(s,:)=RF_TOTi(s,:)-mean(RF_TOTi(s,ind1870),2);
    RF_OZOi_bs(s,:)=RF_OZOi(s,:)-mean(RF_OZOi(s,ind1870),2);
end


figure(60);clf
plot(yearall,squeeze(prctile(projdat_RCSO2_bs(:,1,:),[5 50 95])),'r')
hold on
plot(yearall,squeeze(prctile(projdat_RCSO2_bs(:,2,:),[5 50 95])),'b')
plot(yearall,squeeze(prctile(projdat_RCSO2_bs(:,3,:),[5 50 95])),'c')
plot(yearall,squeeze(prctile(projdat_RCSO2_bs(:,4,:),[5 50 95])),'m')
plot(yearall,RF_AERi_bs,'k--')
xlim([1850 2100])

figure(61);clf
plot(yearall,squeeze(prctile(projdat_RF_bs(:,1,:),[5 50 95])),'r')
hold on
plot(yearall,squeeze(prctile(projdat_RF_bs(:,2,:),[5 50 95])),'b')
plot(yearall,squeeze(prctile(projdat_RF_bs(:,3,:),[5 50 95])),'c')
plot(yearall,squeeze(prctile(projdat_RF_bs(:,4,:),[5 50 95])),'m')
plot(yearall,RF_TOTi_bs,'k--')
xlim([1850 2100])


figure(62);clf
plot(yearall,squeeze(prctile(projdat_RF_bs(:,1,:)-projdat_RCSO2_bs(:,1,:),[5 50 95])),'r')
hold on
plot(yearall,squeeze(prctile(projdat_RF_bs(:,2,:)-projdat_RCSO2_bs(:,2,:),[5 50 95])),'b')
plot(yearall,squeeze(prctile(projdat_RF_bs(:,3,:)-projdat_RCSO2_bs(:,3,:),[5 50 95])),'c')
plot(yearall,squeeze(prctile(projdat_RF_bs(:,4,:)-projdat_RCSO2_bs(:,4,:),[5 50 95])),'m')
plot(yearall,RF_TOTi_bs-RF_AERi_bs,'k--')
xlim([1850 2100])

squeeze(mean(RF_TOTi_bs(:,ind2012)-RF_AERi_bs(:,ind2012),2))

mean(squeeze(mean(projdat_RF_bs(:,1:4,ind2012)-projdat_RCSO2_bs(:,1:4,ind2012),3)))



%% Paper Figure
figure(1000);clf
ind=find(datJD_Cparam==1);
scatter(datJD_RFA(ind), datJD_climate_sensitivity(ind),7,datJD_g(ind),'filled')
hold on
contour(arf,ecs,dhistogram,[.01 .05 .1 .25 .5 .75],'k')
%contour(arf,ecs,smooth2a(dhistogram,1,1),[.01 .05 .1 .25 .5 .75],'k')
xlabel('Aerosol Radiative Forcing [W/m^2]')
ylabel('Equilibrium Climate Sensitivity [^oC]')
xlim([-2 0.2]);ylim([1.2 7])
colorbar;colormap(lbmap(21,'BlueRed'))
title('g [W/m^2/K]');set(gca,'fontsize',15)

figure(1001);clf
scatter(datJD_RFA(ind), datJD_climate_sensitivity(ind),7,datJD_Co(ind),'filled')
hold on
contour(arf,ecs,dhistogram,[.01 .05 .1 .25 .5 .75],'k')
%contour(arf,ecs,smooth2a(dhistogram,1,1),[.01 .05 .1 .25 .5 .75],'k')
xlabel('Aerosol Radiative Forcing [W/m^2]')
ylabel('Equilibrium Climate Sensitivity [^oC]')
xlim([-2 0.2]);ylim([1.2 7])
colorbar;colormap(lbmap(21,'BlueRed'))
title('Co [Wyr/m2/K]');set(gca,'fontsize',15)

figure(1002);clf
scatter(datJD_RFA(ind), datJD_climate_sensitivity(ind),7,datJD_Cs(ind),'filled')
hold on
contour(arf,ecs,dhistogram,[.01 .05 .1 .25 .5 .75],'k')
%contour(arf,ecs,smooth2a(dhistogram,1,1),[.01 .05 .1 .25 .5 .75],'k')
xlabel('Aerosol Radiative Forcing [W/m^2]')
ylabel('Equilibrium Climate Sensitivity [^oC]')
xlim([-2 0.2]);ylim([1.2 7])
colorbar;colormap(lbmap(21,'BlueRed'))
title('Cs [Wyr/m2/K]');set(gca,'fontsize',15)

figure(1003);clf
scatter(datJD_RFA(ind), datJD_climate_sensitivity(ind),7,SATmax(ind,1),'filled')
hold on
contour(arf,ecs,dhistogram,[.01 .05 .1 .25 .5 .75],'k')
%contour(arf,ecs,smooth2a(dhistogram,1,1),[.01 .05 .1 .25 .5 .75],'k')
xlabel('Aerosol Radiative Forcing [W/m^2]')
ylabel('Equilibrium Climate Sensitivity [^oC]')
xlim([-2 0.2]);ylim([1.2 7])
colorbar;caxis([1 2]);colormap(lbmap(21,'BlueRed'))
title('Maximim Temperature');set(gca,'fontsize',15)

print -f1000 -depsc 'figures/Fig1a_JD_g'
print -f1001 -depsc 'figures/Fig1b_JD_Co'
print -f1002 -depsc 'figures/FigRef_JD_Cs'
print -f1003 -depsc 'figures/FigRef_JD_Tmax'

SATtmp=squeeze(SATcommit(:,1,ind2020:end));
yeartmp=yearall(ind2020:end);
[M,I] =max(SATtmp,[],2);

figure(1003);clf
subplot(1,2,1)
plot(datJD_RFA(ind),SATmax(ind,1),'.');hold on
xlabel('Aerosol Radiative Forcing [W/m^2]')
ylabel('Committed Warming above pre-industrial [^oC]')
set(gca,'fontsize',15) ;pbaspect([1 1 1])
plot([-2.5 0.5],[1.5 1.5],'k--');plot([-2.5 0.5],[2 2],'k--')
subplot(1,2,2)
plot(datJD_climate_sensitivity(ind),SATmax(ind,1),'.');hold on
xlabel('Equilibrium Climate Sensitivity [^oC]')
ylabel('Committed Warming above pre-industrial [^oC]')
set(gca,'fontsize',15) ;pbaspect([1 1 1])
plot([1 8],[1.5 1.5],'k--');plot([1 8],[2 2],'k--')

print -f1003 -depsc 'figures/FigRef_warming_vs_ARF'

ind=find(datJD_Cparam==1);
tmp=SATmax(ind,1);
tmpRFA=datJD_RFA(ind);
tmpECS=datJD_climate_sensitivity(ind);
cc=1;
for n=0.5:-.1:-2.5
    ind=find(tmpRFA<n);
    p15(cc)=length(find(tmp(ind)>1.5))/length(ind);
    p2(cc)=length(find(tmp(ind)>2))/length(ind);
    cc=cc+1;
end
cc=1;
for n=1:.25:8
    ind=find(tmpECS>n);
    p15ecs(cc)=length(find(tmp(ind)>1.5))/length(ind);
    p2ecs(cc)=length(find(tmp(ind)>2))/length(ind);
    cc=cc+1;
end
figure(1004);clf
subplot(1,2,1)
plot(0.5:-.1:-2.5,p15);hold on
plot(0.5:-.1:-2.5,p2)
xlabel('Aerosol Radiative Forcing [W/m^2]')
ylabel('Proportion')
set(gca,'fontsize',15) ;pbaspect([1 1 1])
xlim([-2.5 .5])
subplot(1,2,2)
plot(1:.25:8,p15ecs);hold on
plot(1:.25:8,p2ecs)
xlabel('Equilibrium Climate Sensitivity [^oC]')
ylabel('Proportion')
set(gca,'fontsize',15) ;pbaspect([1 1 1])
xlim([1 8])

print -f1004 -depsc 'figures/FigRef_CDFwarming_vs_ARF'



figure(1005);clf
coldc={'b','k','y','r'} ;
for dr=1:4
    plot(yearall(ind2020-10*12:end),prctile(squeeze(SATcommit(find(datJD_Cparam==dr),scen-4,ind2020-10*12:end)),[1 5 50 95 99],1),coldc{dr},'linewidth',1)%fast decay
    hold on
    pause
end
xlim([2010 2100])
set(gca,'ytick',[-.5:.25:1],'fontsize',16)
ylim([-.6 1])
grid on
ylabel('SAT relative to 2020 [^oC]')

print -f1005 -depsc 'figures/committedSAT_rel2020'



figure(1006);clf
cols={'g','b','k','r','b--','k--','r--'}
for myscen=4:7
    figure(1006);
    plot(yearall,squeeze(projdat_RCO2(1,myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('CO2');pbaspect([1.3 1 1])
    plot(yearall,squeeze(projdat_RCO2(2,myscen,:)),cols{myscen})
    plot(yearall,squeeze(projdat_RCO2(3,myscen,:)),cols{myscen})
    plot(yearall,squeeze(projdat_RCO2(4,myscen,:)),cols{myscen})
    figure(1007);
    plot(yearall,squeeze(projdat_RCH4(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('CH4');pbaspect([1.3 1 1]);ylim([0 1])
    figure(1008);
    plot(yearall,squeeze(projdat_RN2O(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('N2O');pbaspect([1.3 1 1]);ylim([0 1])
    figure(1009);
    plot(yearall,squeeze(projdat_RCOZO(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('O3');pbaspect([1.3 1 1]);ylim([0 1])
    figure(1010);
    plot(yearall,squeeze(projdat_RHalo(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('Halogens');pbaspect([1.3 1 1]);ylim([0 1])
    figure(1011);
    plot(yearall,squeeze(prctile(projdat_RCSO2(:,myscen,:),[ 50 ])),cols{myscen});hold on;set(gca,'fontsize',15);title('Aerosol (median)');pbaspect([1.3 1 1]);ylim([-1 0])
end

print -f1006 -depsc 'figures/all_radiative_forcingCO2'
print -f1007 -depsc 'figures/all_radiative_forcingCH4'
print -f1008 -depsc 'figures/all_radiative_forcingN2O'
print -f1009 -depsc 'figures/all_radiative_forcingO3'
print -f1010 -depsc 'figures/all_radiative_forcingHal'
print -f1011 -depsc 'figures/all_radiative_forcingAer'

% Radiative forcing relative to 2020

ind2020=find(yearall==2020);
s=5;
projdat_RF_2020=projdat_RF*NaN;
for n=1:length(projdat_RF)
    projdat_RF_2020(n,s,:)=projdat_RF(n,s,:)-projdat_RF(n,s,ind2020);
end

figure(300);clf
coldc={'b','k','y','r'} ;
for CP=1:4
    plot(yearall(ind2020:end),squeeze(projdat_RCO2(CP,s,ind2020:end))-squeeze(projdat_RCO2(CP,s,ind2020)),coldc{CP})
    hold on
end
plot(yearall(ind2020:end),squeeze(projdat_RCH4(s,ind2020:end))-squeeze(projdat_RCH4(s,ind2020)),'k')
plot(yearall(ind2020:end),squeeze(projdat_RN2O(s,ind2020:end))-squeeze(projdat_RN2O(s,ind2020)),'b')
plot(yearall(ind2020:end),squeeze(projdat_RHalo(s,ind2020:end))-squeeze(projdat_RHalo(s,ind2020)),'g')

%plot(yearall(ind2020:end),prctile(squeeze(projdat_RF_2020(:,s,ind2020:end)),50),'m')

dRF_OZO=squeeze(projdat_RCOZO(s,ind2020)-projdat_RCOZO(s,ind2020-1));
dRF_AER=prctile(squeeze(projdat_RCSO2(:,s,ind2020)-projdat_RCSO2(:,s,ind2020-1)),[5 50 95]);
plot([2020 2100],[1 1]*(dRF_OZO+dRF_AER(1)),'k--')
plot([2020 2100],[1 1]*(dRF_OZO+dRF_AER(2)),'k--')
plot([2020 2100],[1 1]*(dRF_OZO+dRF_AER(3)),'k--')
xlim([2020 2050])

set(gca,'fontsize',16)
legend('CO2','CO2','CO2','CO2','CH4','N2O','Halogens','Aer+Ozone')
pbaspect([1 1 1])

figure(301);clf
[N,X]=hist(SATmax_year(:,1),2020:2100);
plot(X,N)
pbaspect([1.5 1 1])
set(gca,'fontsize',16)
xlim([2020 2100])

print -depsc -f300 'images/RF_rel_to_2020'
print -depsc -f301 'images/year_of_max_T'


figure(300);clf
for CP=1:4
    plot(yearall(ind2020:end),squeeze(projdat_RCO2(CP,s,ind2020:end))-squeeze(projdat_RCO2(CP,s,ind2020)),'r')
    hold on
    pause
end
plot(yearall(ind2020:end),squeeze(projdat_RCH4(s,ind2020:end))-squeeze(projdat_RCH4(s,ind2020)),'k')
plot(yearall(ind2020:end),squeeze(projdat_RN2O(s,ind2020:end))-squeeze(projdat_RN2O(s,ind2020)),'b')
plot(yearall(ind2020:end),squeeze(projdat_RHalo(s,ind2020:end))-squeeze(projdat_RHalo(s,ind2020)),'g')
xlim([2020 2100])
ylim([-.9 0])
set(gca,'fontsize',16)
legend('CO2','CO2','CO2','CO2','CH4','N2O','Halogens','Aer+Ozone')
pbaspect([2 1 1])
print -depsc -f300 'images/RF_rel_to_2020_simple'

% figure(302)
% subplot(3,1,1)
% plot(yearall,squeeze(projdat_RCSO2(1,5,:)));xlim([1850 2050])
% subplot(3,1,2)
% plot(yearall,squeeze(projdat_RCOZO(5,:)));xlim([1850 2050])
% subplot(3,1,3)
% plot(yearall,squeeze(projdat_RF(1,5,:)));xlim([1850 2050])

%% look at concentration
cdat=xlsread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/RadiativeForcing/CO2 mixing ratio.xlsx')
cyr1=cdat(:,1);
cco21=cdat(:,2);ind=find(isnan(cyr1));cyr1(ind)=[];cco21(ind)=[];
cyr2=cdat(:,4);
cco22=cdat(:,5);;ind=find(isnan(cyr2));cyr2(ind)=[];cco22(ind)=[];
cyr3=cdat(:,7);
cco23=cdat(:,8);;ind=find(isnan(cyr3));cyr3(ind)=[];cco23(ind)=[];


colCP={'b','k','y','r'} ;
figure(500);clf
ind=find(yearall>=1861 & yearall<1880);
indc=find(cyr3>=1861 & cyr3<1880);
for CP=1:4
    plot(yearall,squeeze(projdat_CCO2(CP,4:7,:))-squeeze(mean(projdat_CCO2(CP,4,ind),3)),colCP{CP})
    %plot(yearall,squeeze(projdat_CCO2(CP,4:7,:)),colCP{CP})
    hold on
end
xlim([1850 2100]);set(gca,'fontsize',16)
ind=find(years>=1861 & years<1880);
plot(years,conc_CO2-mean(conc_CO2(ind)),'k--')
% plot(years,conc_CO2,'k--')
pbaspect([1.7 1 1])
plot(cyr1,cco21-mean(cco23(indc)))
plot(cyr2,cco22-mean(cco23(indc)))
plot(cyr3,cco23-mean(cco23(indc)))
;ylim([-15 230])

figure(501);clf
ind=find(yearall>=1861 & yearall<1880);
plot(yearall,squeeze(projdat_CCH4(4:7,:))-squeeze(mean(projdat_CCH4(4,ind),2)),'k')
hold on
xlim([1850 2100]);ylim([-300 1800]);set(gca,'fontsize',16)
ind=find(years>=1861 & years<1880);
plot(years,conc_CH4-mean(conc_CH4(ind)),'k--')
pbaspect([1.7 1 1])

figure(502);clf
ind=find(yearall>=1861 & yearall<1880);
plot(yearall,squeeze(projdat_CN2O(4:7,:))-squeeze(mean(projdat_CN2O(4,ind),2)),'k')
hold on
xlim([1850 2100]);ylim([-20 100]);set(gca,'fontsize',16)
ind=find(years>=1861 & years<1880);
plot(years,conc_N2O-mean(conc_N2O(ind)),'k--')
pbaspect([1.7 1 1])

figure(503);clf
ind=find(yearall>=1861 & yearall<1880);
plot(yearall,squeeze(projdat_CHalo(4:7,:))-squeeze(mean(projdat_CHalo(4,ind),2)),'k')
hold on
xlim([1850 2100]);ylim([-200 1200]);set(gca,'fontsize',16)
ind=find(years>=1861 & years<1880);
plot(years,conc_Hal-mean(conc_Hal(ind)),'k--')
pbaspect([1.7 1 1])

print -f500 -depsc 'images/CO2conc'
print -f501 -depsc 'images/CH4conc'

% figure(1003);clf
% subplot(1,2,1)
% plot(datJD_RFA(ind),M(ind),'.')
% xlabel('Aerosol Radiative Forcing [W/m^2]')
% ylabel('Committed Warming above 2020 [^oC]')
% set(gca,'fontsize',15) ;pbaspect([1 1 1])
% subplot(1,2,2)
% plot(datJD_climate_sensitivity(ind),M(ind),'.')
% xlabel('Equilibrium Climate Sensitivity [^oC]')
% ylabel('Committed Warming above 2020 [^oC]')
% set(gca,'fontsize',15) ;pbaspect([1 1 1])

%% Fig 3
figure(600);clf
[N,X]=hist(SATmax(:,1),0.9:.05:4);
plot(X,N/sum(N),'k');
hold on
[N,X]=hist(SATmax(:,2),0.9:.05:4); plot(X,N/sum(N),'b')
[N,X]=hist(SATmax(:,3),0.9:.05:4); plot(X,N/sum(N),'r')
xlim([1 2.5])
pbaspect([1.6 1 1])




figure(601);clf

[N,X]=hist(SATmax(:,1),0.9:.05:4); plot(X,cumsum(N)/sum(N),'k')
hold on
[N,X]=hist(SATmax(:,2),0.9:.05:4); plot(X,cumsum(N)/sum(N),'b')
[N,X]=hist(SATmax(:,3),0.9:.05:4); plot(X,cumsum(N)/sum(N),'r')
pbaspect([1.6 1 1])
xlim([1 2.5])
grid on


print -f600 -depsc 'commitSAT_pdf_corrected'
print -f601 -depsc 'commitSAT_cumpdf_corrected'