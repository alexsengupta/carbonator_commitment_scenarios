% based on round 3 of optimiseCO2_historicalANDdecay3_updated3.m


clear;close all
% Try 1
%     MY_Cup = 713
%     MY_P=60
%     MY_dr0=0.024625
%     MY_a2=4.7e-4
%     uncertaintyCup=0.05;

%     % Try 2
%     MY_Cup = 691
%     MY_P=60
%     MY_dr0=0.024625
%     MY_a2=4.7e-4
%     uncertaintyCup=0.025;
%
%     % Try 3
%     MY_Cup = 691
%     MY_P=75
%     MY_dr0=0.03
%     MY_a2=3.5e-4
%     uncertaintyCup=0.025;


RERUN_OPTIMISATION=1
if RERUN_OPTIMISATION
    % Try 4
    MY_Cup = 687
    MY_P=80
    MY_dr0=0.035
    MY_a2=3.5e-4
    uncertaintyCup=0.025;
    
    
    cases=500000
    usePre2020toOptimise=1;
    uncertainty=0.2; % 25% & +/- unvertanty on selected variables
    opt_yr_end=2020;
    opt_yr_start=1850;
    
    
    DT=1/12;
    
    years=1750:DT:2100;
    yearsX=years;
    
    
    % % Constants for ocean carbon model (glotter)
    Cs=9;       % *60*60*24*365; %W yr/m2/K (Geoffroy uses 7.3 - 9 gives a better volcanic response
    Co=106;     % *60*60*24*365; %W yr/m2/K
    %
    L=1;        % 1.13 from GEOFFROY, but 1 matches more closely to CMIP5 multi-model mean; %Feedback parameter W/m2/K %L ~ 3.5 for no feedback
    g=0.73;     % heat exchange parameter W/m2/K
    %
    d_=50;       % ratio of upper to deep ocean volume
    % ka=1/5;     % Glotter paper suggests 1/5
    % kd=1/20;    % Glotter paper suggests 1/20
    % OM=7.8e22;  % moles of water in ocean
    % AM=1.77e20; % moles of air in atmosphere
    Alk=767; %GtC   % Alkalinity; assumed constant here but buffered on long timescales O[10ka] from the dissolution of CaCO3
    % kh=1.0548e+03;  % at 15oC; ratio of the molar concentrations of CO2 in atmosphere and ocean (Henry's Law)
    k1=8.7184e-07;  % at 15oC; dissociation constant
    k2=5.4426e-10;  % at 15oC; dissociation constant
    % A=kh*AM/(OM/(d_+1)); % A is the ratio of mass of CO2 in atmospheric to upper ocean dissolved CO2,% i.e. A is inversely proportional to
    % % CO2 solubility. A is temperature dependent, however the % effect is small and has been neglected here.
    A_=132.216074 % Slightly different to above in order to have equilibrium at 1850
    %
    % % Constants for Terrestrial model (Svirezhev)
    % m=8.7e-2;   % /yr  1/residence time of carbon in vegetation
    % a2=4.7e-4;  % /GtC %strength of CO2 fertilisation
    % dr0=0.024625; %0.025;  % /yr decomposition rate of soil
    % e=0.5;      % proportion of vegetation that forms soil (remainder goes to atmosphere)
    
    
    allC_CO2=zeros(4,cases,round(length(years)/12)+1)*NaN;
    % allR_CO2=zeros(4,cases,length(years))*NaN;
    
    for scenario=1:4
        switch scenario
            case 1
                col='m';scen='RCP6';
            case 2
                col='k';scen='RCP45';
            case 3
                col='y';scen='RCP3';
            case 4
                col='r';scen='RCP85';
        end
%         [years,meinhausenRF_OZO,meinhausenRF_AER,TSI,alb,emissionCH4(scenario,:),emissionCO2(scenario,:),emissionSO2(scenario,:),emission_volc(scenario,:),radf_other,radf_halo,radf_N2O,radf_CO2(scenario,:),radf_CH4(scenario,:),concCO2(scenario,:),concCH4(scenario,:),mTSI] = load_RCP(scen,1750,2100);
        [years,meinhausenRF_OZO,meinhausenRF_AER,TSI,alb,emissionCH4(scenario,:),emissionCO2(scenario,:),emissionN2O(scenario,:),emissionSO2(scenario,:),emission_volc(scenario,:),radf_other,radf_halo,radf_N2O(scenario,:),radf_CO2(scenario,:),radf_CH4(scenario,:),concCO2(scenario,:),concCH4(scenario,:),concN2O(scenario,:),mTSI,conc_Hal] = load_RCP(scen,1750,2100);;
    end
    for NN=1:cases
        NN
        if NN==1
            rnd=1+randn(12,1)*0;
        else
            rnd=1+randn(12,1)*uncertainty;
            rnd(8)=1+randn*uncertaintyCup;
        end
        
        all_rnd(NN,:)=rnd;
        
        d=rnd(1)*d_;       % ratio of upper to deep ocean volume
        ka=rnd(2)*1/5;     % Glotter paper suggests 1/5
        kd=rnd(3)*1/20;    % Glotter paper suggests 1/20
        A=A_*(d+1)/(50+1);
        
        m=rnd(4)*8.7e-2;   % /yr  1/residence time of carbon in vegetation
        %%%a2=rnd(5)*4.7e-4;  % /GtC %strength of CO2 fertilisation
        a2=rnd(5)*MY_a2;
        %%%dr0=rnd(6)*0.024625; %0.025;  % /yr decomposition rate of soil
        dr0=rnd(6)*MY_dr0;
        e=rnd(7)*0.5;      % proportion of vegetation that forms soil (remainder goes to atmosphere)
        
        Cat=596;
        Cup=rnd(8)*MY_Cup;
        %Clo(1)=rnd(9)*Cup(1)*d;
        Clo=Cup(1)*d;
        
        %%%P(1)=rnd(10)*60; %Gt/yr; NPP P=P(CO2,T,N); here we assume that P only changes due to CO2 fertilisation
        P=rnd(10)*MY_P;
        N=rnd(11)*689.6552; %original value 700; %Gt; Carbon in vegetation
        So=rnd(12)*1218.3; %original value 1200; %Gt  Carbon in Soil
        
        
        for scenario=1:4
            switch scenario
                case 1
                    col='m';scen='RCP6';
                case 2
                    col='k';scen='RCP45';
                case 3
                    col='y';scen='RCP3';
                case 4
                    col='r';scen='RCP85';
            end
            C_CO2pi=Cat(1)/2.13; % % PI concentration; NB glotter used 285 (but need this value for equilibrium)
            C_CO2(1)=C_CO2pi;
            R_CO2(1)=0;
            
            Ts(1)=0;
            To(1)=0;
            
            for y=2:length(years)
                % ocean from Glotter
                
                a=Cup(y-1)/Alk;
                H=(-k1*(1-a) + sqrt( k1^2*(1-a)^2 - 4*k1*k2*(1-2*a) ) )/2; % hydrogen concentration
                pH(y)=-log10(H); % ocean pH
                B=1/( 1 + k1/H + k1*k2/H^2 ); % ratio of dissolved CO2 to total oceanic carbon (B=B(pH))
                Btmp(y)=B;
                
                % terrestrial from Svirezhev
                P(y)=P(1)*(1+a2*(Cat(y-1)-Cat(1))); %NPP with ferilisation effects
                So(y) = So(y-1) + DT*(e*m*N(y-1) - dr0*So(y-1)); %Soil carbon
                N(y) = N(y-1) + DT*(P(y-1) - m*N(y-1)); %Vegetation carbon
                
                
                % Carbon budget
                Cat(y) = Cat(y-1) + DT*emissionCO2(scenario,y-1) + DT*( - ka* (Cat(y-1) - A*B*Cup(y-1) ) ) ... %glotter
                    + DT*( -P(y-1) +(1-e)*m*N(y-1) + dr0*So(y-1) );        %svirezhev
                Cup(y) = Cup(y-1)                       + DT*( + ka* (Cat(y-1) - A*B*Cup(y-1) ) - kd*( Cup(y-1) - Clo(y-1)/d ) );
                Clo(y) = Clo(y-1) +                       DT*(                                    kd*( Cup(y-1) - Clo(y-1)/d ) );
                
                C_CO2(y)=Cat(y)/2.13;
                %R_CO2(y)=5.35*log((C_CO2(y))/C_CO2pi);
                
                %Ts(y)=Ts(y-1)+ DT*( R_CO2(y) - L*Ts(y-1) - g*(Ts(y-1)-To(y-1)) )/Cs ;
                %To(y)=To(y-1)+ DT*( g*(Ts(y-1)-To(y-1)) )/Co;
                
            end
            
            allC_CO2(scenario,NN,:)=C_CO2(1:12:end);
            
            %allR_CO2(scenario,NN,:)=R_CO2;
        end
        ind=find(years==1850);
        ALL1850_Cup(NN)=Cup(ind);
        ALL1850_Clo(NN)=Clo(ind);
        ALL1850_Cat(NN)=Cat(ind);
        ALL1850_P(NN)=P(ind);
        ALL1850_N(NN)=N(ind);
        ALL1850_So(NN)=So(ind);
        ALL1850_B(NN)=Btmp(ind);
    end
    
    concCO2=concCO2(:,1:12:end);    
    myears=years(1:12:end);
    figure(1);
    plot(myears,squeeze(allC_CO2(:,1,:)),'k')
    hold on
    if cases>1;plot(myears,squeeze(allC_CO2(:,2,:)),'r');end
    plot(myears,squeeze(concCO2),'--')
    
    ind_errorEnd=find(myears==opt_yr_end);
    ind_errorStart=find(myears==opt_yr_start);
    clear rmsE
    for n=1:cases
        if usePre2020toOptimise
            error=squeeze(allC_CO2(:,n,ind_errorStart:ind_errorEnd))-concCO2(:,ind_errorStart:ind_errorEnd);
        else
            error=squeeze(allC_CO2(:,n,ind_errorStart:end))-concCO2(:,ind_errorStart:end);
        end
        rmsE(n)=mean(sqrt(mean(error.^2,2)));
    end
    
    figure(2);clf
    plot(rmsE)
    min(rmsE)
    
    ind1850=find(myears==1850);
    ind1870=find(myears>=1861 & myears<1880);
    ind2012=find(myears>=2006 & myears<2018);
    
    CO21870=mean(mean(concCO2(4,ind1870),2)); % Use historical + RCP85
    CO22012=mean(mean(concCO2(4,ind2012),2));
    dCO2_target=CO22012-CO21870;
    
    myCO2_1870=mean(mean(allC_CO2(4,:,ind1870),3),1);
    myCO2_2012=mean(mean(allC_CO2(4,:,ind2012),3),1);
    dCO2=myCO2_2012-myCO2_1870;
    
    %     ind_optim=find(abs(dCO2-dCO2_target)<2 & abs(myCO2_2012-CO22012)<2)
    
    %%
    ind1860_2018=find(myears>1860 & myears<2018)
    my_meanCO2=mean(mean(allC_CO2(4,:,ind1860_2018),3),1);
    target_meanCO2=mean(mean(concCO2(4,ind1860_2018),2));
    %ind_optim=find(abs(dCO2-dCO2_target)<2 & abs(my_meanCO2-target_meanCO2)<2);
    ind_optim=find(abs(dCO2-dCO2_target)<1 & rmsE<6);
    
    %%
    
    figure(22);clf
    for n=1:4
        plot(myears,squeeze(allC_CO2(n,ind_optim,:)),'k')
        hold on
    end
    plot(myears,squeeze(concCO2),'r--')
    xlim([1850 2100])
    
    alld=all_rnd(:,1)*d_;       % ratio of upper to deep ocean volume
    allka=all_rnd(:,2)*1/5;     % Glotter paper suggests 1/5
    allkd=all_rnd(:,3)*1/20;    % Glotter paper suggests 1/20
    allA=A_*(alld+1)/(50+1);
    allm=all_rnd(:,4)*8.7e-2;   % /yr  1/residence time of carbon in vegetation
    alla2=all_rnd(:,5)*MY_a2;  % /GtC %strength of CO2 fertilisation
    alldr0=all_rnd(:,6)*MY_dr0; %0.025;  % /yr decomposition rate of soil
    alle=all_rnd(:,7)*0.5;      % proportion of vegetation that forms soil (remainder goes to atmosphere)
    allP=all_rnd(:,10)*MY_P; %Gt/yr; NPP P=P(CO2,T,N); here we assume that P only changes due to CO2 fertilisation
    allN=all_rnd(:,11)*689.6552; %original value 700; %Gt; Carbon in vegetation
    allSo=all_rnd(:,12)*1218.3; %original value 1200; %Gt  Carbon in Soil
    allCup=all_rnd(:,8)*MY_Cup;
    allClo=allCup.*alld;
    
    save 'optimised_parameters_v4_dump.mat'
else
    load 'optimised_parameters_v4_dump.mat'
end


%% RUN FROM 1850
clear emission* radf_* conc*
for scenario=1:4
        switch scenario
            case 1
                col='m';scen='RCP6';
            case 2
                col='k';scen='RCP45';
            case 3
                col='y';scen='RCP3';
            case 4
                col='r';scen='RCP85';
        end
%         [years,meinhausenRF_OZO,meinhausenRF_AER,TSI,alb,emissionCH4(scenario,:),emissionCO2(scenario,:),emissionSO2(scenario,:),emission_volc(scenario,:),radf_other,radf_halo,radf_N2O,radf_CO2(scenario,:),radf_CH4(scenario,:),concCO2(scenario,:),concCH4(scenario,:),mTSI] = load_RCP(scen,1750,2100);
        [years,meinhausenRF_OZO,meinhausenRF_AER,TSI,alb,emissionCH4(scenario,:),emissionCO2(scenario,:),emissionN2O(scenario,:),emissionSO2(scenario,:),emission_volc(scenario,:),radf_other,radf_halo,radf_N2O(scenario,:),radf_CO2(scenario,:),radf_CH4(scenario,:),concCO2(scenario,:),concCH4(scenario,:),concN2O(scenario,:),mTSI,conc_Hal] = load_RCP(scen,1850,2100);;
end
    

ALLrmsE = rmsE(ind_optim);
clear ALL_*
for SS=1:length(ind_optim)
    
    d=alld(ind_optim(SS));       % ratio of upper to deep ocean volume
    ka=allka(ind_optim(SS));     % Glotter paper suggests 1/5
    kd=allkd(ind_optim(SS));    % Glotter paper suggests 1/20
    A=allA(ind_optim(SS));
    m=allm(ind_optim(SS));   % /yr  1/residence time of carbon in vegetation
    %%%a2=rnd(5)*4.7e-4;  % /GtC %strength of CO2 fertilisation
    a2=alla2(ind_optim(SS));
    %%%dr0=rnd(6)*0.024625; %0.025;  % /yr decomposition rate of soil
    dr0=alldr0(ind_optim(SS));
    e=alle(ind_optim(SS));      % proportion of vegetation that forms soil (remainder goes to atmosphere)
    
    clear Cat Clo Cup P N So C_CO2 C_CH4 To Ts
    
    Cat=ALL1850_Cat(ind_optim(SS));
    %%%Cup(1)=rnd(8)*713;
    Cup=ALL1850_Cup(ind_optim(SS));
    Clo=ALL1850_Clo (ind_optim(SS));
    %%%P(1)=rnd(10)*60; %Gt/yr; NPP P=P(CO2,T,N); here we assume that P only changes due to CO2 fertilisation
    P= ALL1850_P(ind_optim(SS));
    N=ALL1850_N(ind_optim(SS)); %original value 700; %Gt; Carbon in vegetation
    So=ALL1850_So(ind_optim(SS)); %original value 1200; %Gt  Carbon in Soil
    
    
    To(1)=0;
    Ts(1)=0;
    
    
    ALL_d(SS)=d;
    ALL_ka(SS)=ka;
    ALL_kd(SS)=kd;
    ALL_m(SS)=m;
    ALL_a2(SS)=a2;
    ALL_dr0(SS)=dr0;
    ALL_e(SS)=e;
    ALL_Cup(SS)=Cup(1);
    ALL_Clo(SS)=Clo(1);
    ALL_P(SS)=P(1);
    ALL_N(SS)=N(1);
    ALL_So(SS)=So(1);
    ALL_A(SS)=A(1);
    
    for scenario=1:4
        clear C_CO2 R_CO2
        C_CO2pi=Cat(1)/2.13; % % PI concentration; NB glotter used 285 (but need this value for equilibrium)
        C_CO2(1)=C_CO2pi;
        R_CO2(1)=0;
        
        for y=2:length(years)
            % ocean from Glotter
            a=Cup(y-1)/Alk;
            H=(-k1*(1-a) + sqrt( k1^2*(1-a)^2 - 4*k1*k2*(1-2*a) ) )/2; % hydrogen concentration
            pH(y)=-log10(H); % ocean pH
            B=1/( 1 + k1/H + k1*k2/H^2 ); % ratio of dissolved CO2 to total oceanic carbon (B=B(pH))
            
            % terrestrial from Svirezhev
            P(y)=P(1)*(1+a2*(Cat(y-1)-Cat(1))); %NPP with ferilisation effects
            So(y) = So(y-1) + DT*(e*m*N(y-1) - dr0*So(y-1)); %Soil carbon
            N(y) = N(y-1) + DT*(P(y-1) - m*N(y-1)); %Vegetation carbon
            
            
            % Carbon budget
            Cat(y) = Cat(y-1) + DT*emissionCO2(scenario,y-1) + DT*( - ka* (Cat(y-1) - A*B*Cup(y-1) ) ) ... %glotter
                + DT*( -P(y-1) +(1-e)*m*N(y-1) + dr0*So(y-1) );        %svirezhev
            Cup(y) = Cup(y-1)                       + DT*( + ka* (Cat(y-1) - A*B*Cup(y-1) ) - kd*( Cup(y-1) - Clo(y-1)/d ) );
            Clo(y) = Clo(y-1) +                       DT*(                                    kd*( Cup(y-1) - Clo(y-1)/d ) );
            
            C_CO2(y)=Cat(y)/2.13;
            R_CO2(y)=5.35*log((C_CO2(y))/C_CO2pi);
            
            Ts(y)=Ts(y-1)+ DT*( R_CO2(y) - L*Ts(y-1) - g*(Ts(y-1)-To(y-1)) )/Cs ;
            To(y)=To(y-1)+ DT*( g*(Ts(y-1)-To(y-1)) )/Co;
        end
        bestC_CO2(scenario,:)=C_CO2(1:12:end);
        bestR_CO2(scenario,:)=R_CO2(1:12:end);
        bestTs(scenario,:)=Ts(1:12:end);
    end
    ALL_bestC_CO2(SS,:,:)=bestC_CO2;
    ALL_bestR_CO2(SS,:,:)=bestR_CO2;
    ALL_bestTs(SS,:,:)=bestTs;
    clear bestC_CO2 bestR_CO2 bestTs
    
    % 1pc CO2 simulation, until 100Ht released and then stop emission
    F1=1;
    F2=0;
    F3=0;
    F4=0;
    F5=0;
    F6=0;
    randomise_parameters=0;
    
    Ts(1)=0;
    To(1)=0;
    
    C_CO2pi=Cat(1)/2.13; % % PI concentration; NB glotter used 285 (but need this value for equilibrium)
    C_CO2(1)=C_CO2pi;
    R_CO2(1)=0;
    C_CH4pi=791; % PI concentration
    C_CH4(1)=0;
    R_CH4(1)=0;
    R_SO2(1)=0;
    OpTkhkV(1)=0;
    SL(1)=0;
    
    CC_CO2=C_CO2(1);
    
    standard_emission=0;
    DC=[];
    forced_emision(1)=0;
    co2switch=0;
    r=nthroot(1.01,12);
    % Loop through timesteps
    for y=2:length(yearsX)
        
        r=nthroot(1.01,12);
        
        % Carbon Cycle
        % ocean from Glotter
        a=Cup(y-1)/Alk;
        H=(-k1*(1-a) + sqrt( k1^2*(1-a)^2 - 4*k1*k2*(1-2*a) ) )/2; % hydrogen concentration
        pH(y)=-log10(H); % ocean pH
        B=1/( 1 + k1/H + k1*k2/H^2 ); % ratio of dissolved CO2 to total oceanic carbon (B=B(pH))
        
        % terrestrial from Svirezhev
        P(y)=P(1)*(1+a2*(Cat(y-1)-Cat(1))); %NPP with ferilisation effects
        So(y) = So(y-1) + F1*DT*(e*m*N(y-1) - dr0*So(y-1)); %Soil carbon
        N(y) = N(y-1) + F1*DT*(P(y-1) - m*N(y-1)); %Vegetation carbon
        
        CC_CO2(y)=CC_CO2(y-1)*r;
        if sum(DC)<=1000 & co2switch==0
            % calculate emission required for a 1% CO2 increase
            forced_emision(y) = (CC_CO2(y)*2.13 - Cat(y-1))/DT - (( - ka* (Cat(y-1) - A*B*Cup(y-1) ) )   +    ( -P(y-1) +(1-e)*m*N(y-1) + dr0*So(y-1) ));
        else
            forced_emision(y)=0;
            co2switch=1;
        end
        
        % Carbon budget
        if standard_emission
            Cat(y) = Cat(y-1) + F1*DT*emissionCO2(scenario,y-1) + F1*DT*( - ka* (Cat(y-1) - A*B*Cup(y-1) ) ) ... %glotter
                + F1*DT*( -P(y-1) +(1-e)*m*N(y-1) + dr0*So(y-1) );        %svirezhev
            Cup(y) = Cup(y-1)                       + F1*DT*( + ka* (Cat(y-1) - A*B*Cup(y-1) ) - kd*( Cup(y-1) - Clo(y-1)/d ) );
            Clo(y) = Clo(y-1) +                       F1*DT*(                                    kd*( Cup(y-1) - Clo(y-1)/d ) );
        else
            Cat(y) = Cat(y-1) + F1*DT*forced_emision(y) + F1*DT*( - ka* (Cat(y-1) - A*B*Cup(y-1) ) ) ... %glotter
                + F1*DT*( -P(y-1) +(1-e)*m*N(y-1) + dr0*So(y-1) );        %svirezhev
            Cup(y) = Cup(y-1)                       + F1*DT*( + ka* (Cat(y-1) - A*B*Cup(y-1) ) - kd*( Cup(y-1) - Clo(y-1)/d ) );
            Clo(y) = Clo(y-1) +                       F1*DT*(                                    kd*( Cup(y-1) - Clo(y-1)/d ) );
        end
        
        % convert CO2 emissions to RF
        C_CO2(y)=Cat(y)/2.13;
        R_CO2(y)=5.35*log((C_CO2(y))/C_CO2pi);
        
        DC(y)=(Cat(y)+Cup(y)+Clo(y)+P(y)+So(y)+N(y))-(Cat(y-1)+Cup(y-1)+Clo(y-1)+P(y-1)+So(y-1)+N(y-1));
        
        Ts(y)=Ts(y-1)+ DT*( R_CO2(y) - L*Ts(y-1) - g*(Ts(y-1)-To(y-1)) )/Cs ;
        To(y)=To(y-1)+ DT*( g*(Ts(y-1)-To(y-1)) )/Co;
        
        
    end
    
    
    ind=find(C_CO2==max(C_CO2));
    ALL_decay_CO2(SS,:)=C_CO2;
    figure(3);
    plot(yearsX-yearsX(ind),C_CO2-C_CO2(ind),'r-')
    hold on
    xlim([0 100]);ylim([-120 10])
    
    ALL_CO2_100yr(SS)=C_CO2(ind+100*12)-C_CO2(ind);
    ALL_CO2_100yr_ts(SS,:)=C_CO2(ind:ind+100*12)-C_CO2(ind);
    ALL_CO2fraction_100yr_ts(SS,:)=(C_CO2(ind:ind+100*12)-C_CO2(1))/(C_CO2(ind)-C_CO2(1));
    
    
    clear C_CO2 R_CO2
    
end
concCO2=concCO2(:,1:12:end);
myears=years(1:12:end);
figure(499);clf
plot(myears,squeeze(ALL_bestC_CO2(:,1,:)));hold on
plot(myears,squeeze(ALL_bestC_CO2(:,2,:)))
plot(myears,squeeze(ALL_bestC_CO2(:,3,:)))
plot(myears,squeeze(ALL_bestC_CO2(:,4,:)))
plot(myears,concCO2,'r--','linewidth',2)
xlim([1850 2100])


i1870=find(myears>=1861 & myears<1880);
i2012=find(myears>=2006 & myears<2018);
figure(29);clf
hist(squeeze(mean(ALL_bestC_CO2(:,1,i2012),3)-mean(ALL_bestC_CO2(:,1,i1870),3)))
mean(concCO2(1,i2012))-mean(concCO2(1,i1870))

figure(3);clf
subplot(3,4,1)
tmp=allC_CO2(4,:,end); % CO2 at 2100
plot(tmp,allP,'.');hold on
plot(tmp(ind_optim),allP(ind_optim),'ro');title(['allP ',num2str(corr(tmp(ind_optim)',allP(ind_optim)))])
subplot(3,4,2)
plot(tmp,allN,'.');hold on
plot(tmp(ind_optim),allN(ind_optim),'ro');title(['allN ',num2str(corr(tmp(ind_optim)',allN(ind_optim)))])
subplot(3,4,3)
plot(tmp,allSo,'.');hold on
plot(tmp(ind_optim),allSo(ind_optim),'ro');title(['allSo ',num2str(corr(tmp(ind_optim)',allSo(ind_optim)))])
subplot(3,4,4)
plot(tmp,allCup,'.');hold on
plot(tmp(ind_optim),allCup(ind_optim),'ro');title(['allCup ',num2str(corr(tmp(ind_optim)',allCup(ind_optim)))])
subplot(3,4,5)
plot(tmp,allClo,'.');hold on
plot(tmp(ind_optim),allClo(ind_optim),'ro');title(['allClo ',num2str(corr(tmp(ind_optim)',allClo(ind_optim)))])
subplot(3,4,6)
plot(tmp,alle,'.');hold on
plot(tmp(ind_optim),alle(ind_optim),'ro');title(['alle ',num2str(corr(tmp(ind_optim)',alle(ind_optim)))])
subplot(3,4,7)
plot(tmp,alldr0,'.');hold on
plot(tmp(ind_optim),alldr0(ind_optim),'ro');title(['alldr0 ',num2str(corr(tmp(ind_optim)',alldr0(ind_optim)))])
subplot(3,4,8)
plot(tmp,alla2,'.');hold on
plot(tmp(ind_optim),alla2(ind_optim),'ro');title(['alla2 ',num2str(corr(tmp(ind_optim)',alla2(ind_optim)))])
subplot(3,4,9)
plot(tmp,allm,'.');hold on
plot(tmp(ind_optim),allm(ind_optim),'ro');title(['allm ',num2str(corr(tmp(ind_optim)',allm(ind_optim)))])
subplot(3,4,10)
plot(tmp,allkd,'.');hold on
plot(tmp(ind_optim),allkd(ind_optim),'ro');title(['allkd ',num2str(corr(tmp(ind_optim)',allkd(ind_optim)))])
subplot(3,4,11)
plot(tmp,allka,'.');hold on
plot(tmp(ind_optim),allka(ind_optim),'ro');title(['allka ',num2str(corr(tmp(ind_optim)',allka(ind_optim)))])
subplot(3,4,12)
plot(tmp,alld,'.');hold on
plot(tmp(ind_optim),alld(ind_optim),'ro');title(['alld ',num2str(corr(tmp(ind_optim)',alld(ind_optim)))])


figure(4);clf
subplot(3,4,1)
plot(myCO2_2012,allP,'.');hold on
plot(myCO2_2012(ind_optim),allP(ind_optim),'ro');plot(3000,allP(1),'kx','linewidth',4)
subplot(3,4,2)
plot(myCO2_2012,allN,'.');hold on
plot(myCO2_2012(ind_optim),allN(ind_optim),'ro');plot(3000,allN(1),'kx','linewidth',4)
subplot(3,4,3)
plot(myCO2_2012,allSo,'.');hold on
plot(myCO2_2012(ind_optim),allSo(ind_optim),'ro');plot(3000,allSo(1),'kx','linewidth',4)
subplot(3,4,4)
plot(myCO2_2012,allCup,'.');hold on
plot(myCO2_2012(ind_optim),allCup(ind_optim),'ro');plot(3000,allCup(1),'kx','linewidth',4);%plot(ind_optim./ind_optim*2800,Cup1850(ind_optim),'k.')
subplot(3,4,5)
plot(myCO2_2012,allClo,'.');hold on
plot(myCO2_2012(ind_optim),allClo(ind_optim),'ro');plot(3000,allClo(1),'kx','linewidth',4);%plot(ind_optim./ind_optim*2800,Clo1850(ind_optim),'k.')
subplot(3,4,6)
plot(myCO2_2012,alle,'.');hold on
plot(myCO2_2012(ind_optim),alle(ind_optim),'ro');plot(3000,alle(1),'kx','linewidth',4)
subplot(3,4,7)
plot(myCO2_2012,alldr0,'.');hold on
plot(myCO2_2012(ind_optim),alldr0(ind_optim),'ro');plot(3000,alldr0(1),'kx','linewidth',4)
subplot(3,4,8)
plot(myCO2_2012,alla2,'.');hold on
plot(myCO2_2012(ind_optim),alla2(ind_optim),'ro');plot(3000,alla2(1),'kx','linewidth',4)
subplot(3,4,9)
plot(myCO2_2012,allm,'.');hold on
plot(myCO2_2012(ind_optim),allm(ind_optim),'ro');plot(3000,allm(1),'kx','linewidth',4)
subplot(3,4,10)
plot(myCO2_2012,allkd,'.');hold on
plot(myCO2_2012(ind_optim),allkd(ind_optim),'ro');plot(3000,allkd(1),'kx','linewidth',4)
subplot(3,4,11)
plot(myCO2_2012,allka,'.');hold on
plot(myCO2_2012(ind_optim),allka(ind_optim),'ro');plot(3000,allka(1),'kx','linewidth',4)
subplot(3,4,12)
plot(myCO2_2012,alld,'.');hold on
plot(myCO2_2012(ind_optim),alld(ind_optim),'ro');plot(3000,alld(1),'kx','linewidth',4)


figure(5);clf
subplot(3,4,1);plot(ALL_CO2_100yr,ALL_d,'.');title(['ALL_d  ',num2str(corr(ALL_CO2_100yr',ALL_d'))])
subplot(3,4,2);plot(ALL_CO2_100yr,ALL_ka,'.');title(['ALL_ka  ',num2str(corr(ALL_CO2_100yr',ALL_ka'))])
subplot(3,4,3);plot(ALL_CO2_100yr,ALL_kd,'.');title(['ALL_kd  ',num2str(corr(ALL_CO2_100yr',ALL_kd'))])
subplot(3,4,4);plot(ALL_CO2_100yr,ALL_m,'.');title(['ALL_m  ',num2str(corr(ALL_CO2_100yr',ALL_m'))])
subplot(3,4,5);plot(ALL_CO2_100yr,ALL_a2,'.');title(['ALL_a2  ',num2str(corr(ALL_CO2_100yr',ALL_a2'))])
subplot(3,4,6);plot(ALL_CO2_100yr,ALL_dr0,'.');title(['ALL_dr0  ',num2str(corr(ALL_CO2_100yr',ALL_dr0'))])
subplot(3,4,7);plot(ALL_CO2_100yr,ALL_e,'.');title(['ALL_e  ',num2str(corr(ALL_CO2_100yr',ALL_e'))])
subplot(3,4,8);plot(ALL_CO2_100yr,ALL_Cup,'.');title(['ALL_Cup  ',num2str(corr(ALL_CO2_100yr',ALL_Cup'))])
subplot(3,4,9);plot(ALL_CO2_100yr,ALL_Clo,'.');title(['ALL_Clo  ',num2str(corr(ALL_CO2_100yr',ALL_Clo'))])
subplot(3,4,10);plot(ALL_CO2_100yr,ALL_P,'.');title(['ALL_P  ',num2str(corr(ALL_CO2_100yr',ALL_P'))])
subplot(3,4,11);plot(ALL_CO2_100yr,ALL_N,'.');title(['ALL_N  ',num2str(corr(ALL_CO2_100yr',ALL_N'))])
subplot(3,4,12);plot(ALL_CO2_100yr,ALL_A,'.');title(['ALL_A  ',num2str(corr(ALL_CO2_100yr',ALL_A'))])


figure(51);clf
subplot(3,4,1);plot(ALL_bestC_CO2(:,4,end),ALL_d,'.');title(['ALL_d  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_d'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_d) min(ALL_d)])
subplot(3,4,2);plot(ALL_bestC_CO2(:,4,end),ALL_ka,'.');title(['ALL_ka  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_ka'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_ka) min(ALL_ka)])
subplot(3,4,3);plot(ALL_bestC_CO2(:,4,end),ALL_kd,'.');title(['ALL_kd  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_kd'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_kd) min(ALL_kd)])
subplot(3,4,4);plot(ALL_bestC_CO2(:,4,end),ALL_m,'.');title(['ALL_m  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_m'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_m) min(ALL_m)])
subplot(3,4,5);plot(ALL_bestC_CO2(:,4,end),ALL_a2,'.');title(['ALL_a2  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_a2'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_a2) min(ALL_a2)])
subplot(3,4,6);plot(ALL_bestC_CO2(:,4,end),ALL_dr0,'.');title(['ALL_dr0  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_dr0'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_dr0) min(ALL_dr0)])
subplot(3,4,7);plot(ALL_bestC_CO2(:,4,end),ALL_e,'.');title(['ALL_e  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_e'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_e) min(ALL_e)])
subplot(3,4,8);plot(ALL_bestC_CO2(:,4,end),ALL_Cup,'.');title(['ALL_Cup  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_Cup'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_Cup) min(ALL_Cup)])
subplot(3,4,9);plot(ALL_bestC_CO2(:,4,end),ALL_Clo,'.');title(['ALL_Clo  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_Clo'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_Clo) min(ALL_Clo)])
subplot(3,4,10);plot(ALL_bestC_CO2(:,4,end),ALL_P,'.');title(['ALL_P  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_P'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_P) min(ALL_P)])
subplot(3,4,11);plot(ALL_bestC_CO2(:,4,end),ALL_N,'.');title(['ALL_N  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_N'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_N) min(ALL_N)])
subplot(3,4,12);plot(ALL_bestC_CO2(:,4,end),ALL_A,'.');title(['ALL_A  ',num2str(corr(ALL_bestC_CO2(:,4,end),ALL_A'))]);hold on;plot([1 1]*concCO2(4,end),[max(ALL_A) min(ALL_A)])


figure(6);clf
[N,X]=hist(allCup,[600:10:800]);plot(X,N/max(N),'b');hold on
[N,X]=hist(allCup(ind_optim),[600:10:800]);plot(X,N/max(N),'r');


figure(59);clf
I = imread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/decayCO2.jpg');
yr=0:1/12:100
imshow(I)
hold on
plot(180+7.85*yr,50-4.76*ALL_CO2_100yr_ts(:,:)','r','linewidth',2)



figure(60);clf
for SS=1:length(ALL_A)
    figure(60);cla
    I = imread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/decayCO2.jpg');
    imshow(I)
    hold on
    
    %if ALL_CO2_100yr_ts(SS,end)>-130 & ALL_CO2_100yr_ts(SS,ind40)>-90 & ALLrmsE(SS)<8
    plot(180+7.85*yr,50-4.76*ALL_CO2_100yr_ts(SS,:)','r','linewidth',2)
    figure(61);clf
    plot(myears,squeeze(ALL_bestC_CO2(SS,1,:)));hold on
    plot(myears,squeeze(ALL_bestC_CO2(SS,2,:)))
    plot(myears,squeeze(ALL_bestC_CO2(SS,3,:)))
    plot(myears,squeeze(ALL_bestC_CO2(SS,4,:)))
    plot(myears,concCO2,'r--','linewidth',2)
    SS
    pause
    %end
end


ALLsubset=[5,16,21,41];

figure(60);clf
I = imread('/Users/z3045790/OneDrive - UNSW/My Projects/sherwood_committement/decayCO2.jpg');
yr=0:1/12:100
imshow(I)
hold on
for SS=ALLsubset
    figure(60);hold on
    plot(180+7.85*yr,50-4.76*ALL_CO2_100yr_ts(SS,:)','r','linewidth',2);hold on
    
    figure(61)
    plot(myears,squeeze(ALL_bestC_CO2(SS,1,:)));hold on
    plot(myears,squeeze(ALL_bestC_CO2(SS,2,:)))
    plot(myears,squeeze(ALL_bestC_CO2(SS,3,:)))
    plot(myears,squeeze(ALL_bestC_CO2(SS,4,:)))
    plot(myears,concCO2,'r--','linewidth',2)
    SS
    pause
end


OPT_CO2_100yr=ALL_CO2_100yr(ALLsubset);
OPT_CO2_100yr_ts=ALL_CO2_100yr_ts(ALLsubset,:);
OPT_CO2fraction_100yr_ts=ALL_CO2fraction_100yr_ts(ALLsubset,:)
OPT_decay_CO2=ALL_decay_CO2(ALLsubset);
OPT_rmsE=ALLrmsE(ALLsubset);
OPT_bestC_CO2=ALL_bestC_CO2(ALLsubset,:,:);
OPT_bestR_CO2=ALL_bestR_CO2(ALLsubset,:,:);
OPT_bestTs=ALL_bestTs(ALLsubset,:,:);

OPT_A=ALL_A(ALLsubset);
OPT_d=ALL_d(ALLsubset);
OPT_dr0=ALL_dr0(ALLsubset);
OPT_e=ALL_e(ALLsubset);
OPT_ka=ALL_ka(ALLsubset);
OPT_Clo=ALL_Clo(ALLsubset);
OPT_kd=ALL_kd(ALLsubset);
OPT_Cup=ALL_Cup(ALLsubset);
OPT_m=ALL_m(ALLsubset);
OPT_N=ALL_N(ALLsubset);
OPT_P=ALL_P(ALLsubset);
OPT_So=ALL_So(ALLsubset);
OPT_a2=ALL_a2(ALLsubset);


save 'optimised_parameters_v4.mat' ALL* myears OPT*
