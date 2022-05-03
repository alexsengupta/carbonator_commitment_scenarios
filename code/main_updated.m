% updated vesrio ensure aerosol and ozone forcing starts at 0
clear; close all
fnC='optimised_parameters_v4.mat';
load(fnC)
clearvars -except ALLsubset fnC

randomise_parameters=1; %Randomise g,L,ARF,Co,Cs
% HadCRU4,Berkely Earth: DL/DO=1.32, 1.43.
% ((0.3)DL+(0.7)DO)/DO = 1.096,1.129
% Use average value 1.11
DTscale=1.11; 
Trange=[0.89 1.17]; % change between 1861–1880 and 2006–2018

no_of_good_simulations=300000;
CparameterSet=length(ALLsubset); % coresponds to the number of carbon model paramter sets
 

for test=[4]; % Use range from glotter
    
    % tasks
    find_simulations=0; % find parameter set for models that meet teperature criteria
    plot_simulations=0;
    subset_based_onJointPDF=0;
    jointPDF_projections=0;
    jointPDF_projections_const_conc=0;
    plot_jointPDF_projections=1;
    plot_jointPDF_projections_const_conc=0;
    
    check_RF=0;
    
    if test==0
        load_glotter % load g, C from glotter
    end
    if test==1
        load_glotter % load g, C from glotter
        range_g(1)=range_g(1)*.5;
        range_g(2)=range_g(2)*2;
        range_Co(1)=range_Co(1)*.5;
        range_Co(2)=range_Co(2)*2;
        range_Cs(1)=range_Cs(1)*.5;
        range_Cs(2)=range_Cs(2)*2;
    end 
    if test==2
        load_glotter % load g, C from glotter
        range_g(1)=range_g(1)*.1;
        range_g(2)=range_g(2)*2;
        range_Co(1)=range_Co(1)*.5;
        range_Co(2)=range_Co(2)*2;
        range_Cs(1)=range_Cs(1)*.1;
        range_Cs(2)=range_Cs(2)*2;
    end
    if test==3
        load_glotter % load g, C from glotter
        range_g(1)=range_g(1)*.05;
        range_g(2)=range_g(2)*1.7;
        range_Co(1)=range_Co(1)*.5;
        range_Co(2)=range_Co(2)*2;
        range_Cs(1)=range_Cs(1)*.05;
        range_Cs(2)=range_Cs(2)*2;
    end
    if test==4
        load_glotter % load g, C from glotter
        range_g(1)=range_g(1)*.05;
        range_g(2)=range_g(2)*2;
        range_Co(1)=range_Co(1)*.5;
        range_Co(2)=range_Co(2)*2;
        range_Cs(1)=range_Cs(1)*.05;
        range_Cs(2)=range_Cs(2)*2;
    end
    
    datafolder=['data_',num2str(test)];
    mkdir(datafolder)
    
    
    if check_RF
         uniform_dist=1; % pick RF and climate sensitivity from a uniform distribution based on ranges from the sherwood joint distribution
        
        startyear=1850;
        endyear=2100;
        RF_SO2_2006_2018=-0.50; % default rediative forcing change between 1861–1880 and 2006–2018 (will be scaled by sf)
        [years,meinhausenRF_OZO,meinhausenRF_AER,TSI,alb,emissionCH4,emissionCO2,emissionN2O,emissionSO2,emission_volc,radf_other,radf_halo,radf_N2O,radf_CO2,radf_CH4,conc_CO2,conc_CH4,conc_N2O,mTSI,conc_Hal]=load_RCP('RCP85',startyear,endyear) % Load RCP85 emission & radiative forcing
        
        plotRF=0;
        
        
    end
    if find_simulations
        
        uniform_dist=1; % pick RF and climate sensitivity from a uniform distribution based on ranges from the sherwood joint distribution
        
        startyear=1850; %<<UPDATE 2>>
        endyear=2018;
        RF_SO2_2006_2018=-0.50; % default rediative forcing change between 1861–1880 and 2006–2018 (will be scaled by sf)
        [years,meinhausenRF_OZO,meinhausenRF_AER,TSI,alb,emissionCH4,emissionCO2,emissionN2O,emissionSO2,emission_volc,radf_other,radf_halo,radf_N2O,radf_CO2,radf_CH4,conc_CO2,conc_CH4,conc_N2O,mTSI,conc_Hal]=load_RCP('RCP85',startyear,endyear) % Load RCP85 emission & radiative forcing
        
        plotRF=1;
        % NB faster to process in chunks
        % total number of simulations match temperature criteria = no_of_good_simulations x no_iterations

        %select_simulations_in_Trange
        select_simulations_in_Trange_updated
        combine_select_simulations_in_Trange
        
        figure(1);clf
        subplot(2,2,1)
        scatter(alldat_climate_sensitivity,alldat_RF)
        subplot(2,2,2)
        scatter(alldat_Co,alldat_Cs)
        subplot(2,2,3)
        scatter(alldat_g,alldat_sf)
        subplot(2,2,4)
        hist(alldat_deltaT,100)
        
        figure(2);clf
        subplot(2,2,1)
        scatter(dat_climate_sensitivity,dat_RFA)
        subplot(2,2,2)
        scatter(dat_Co,dat_Cs)
        subplot(2,2,3)
        scatter(dat_g,dat_sf)
        subplot(2,2,4)
        hist(dat_deltaT,100)
        
        figure(3);plot(alldat_RF)
        
        
    end
    
    if plot_simulations
        
        load(['data_',num2str(test),'/alldat_select_simulations_in_Trange'])
        
        figure(1+test);clf
        for v=1:6
            switch v
                case 1;
                    var='g';
                case 2;
                    var='Co';
                case 3;
                    var='Cs';
                case 4;
                    var='L';
                case 5;
                    var='deltaT';
                case 6;
                    var='RFA';
            end
            subplot(3,2,v)
            eval(['tmp=alldat_',var,';'])
            [N,X]=hist(tmp,100)
            plot(X,N/max(N),'r')
            hold on
            eval(['tmp=dat_',var,';'])
            [N,X]=hist(tmp,100)
            hold on; plot(X,N/max(N),'k')
            title(var)
        end
        
        figure(10+test);clf
        load_2d_histogram;
        dhistogram=dhistogram';
        ecs=0.1:.2:7.9;
        arf=-2.9:.2:.1;
        dhistogram=dhistogram/max(dhistogram(:));
        dhistogram=reshape(dhistogram,42,16);
        dhistogram(41:42,:)=[];
        scatter(dat_RFA, dat_climate_sensitivity,3,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.9)
        hold on
        contour(arf,ecs,dhistogram,[0:.1:1])
        xlabel('Aerosol RF')
        ylabel('ECS')
        xlim([-2 0.5])
        
    end
    %%
    if subset_based_onJointPDF
        
        load(['data_',num2str(test),'/alldat_select_simulations_in_Trange'])
        
        no_of_jointDist_members=10000;
        
        subset_based_on_jointPDF_byCparam
        
        save(['data_',num2str(test),'/JD_pdf'],'datJD_*')
    end
    %%
    if jointPDF_projections
        'jointPDF_projections'
        
        load(['data_',num2str(test),'/JD_pdf'])
        startyear=1850; %<<UPDATE 2>>
        endyear=2100;
        all_scenarios=[1:7]; % RCP+COMMIT
        
        %JD_projections
        JD_projections_updated
        
        save(['data_',num2str(test),'/JD_projections_pdf'],'projdat_*','yearall','datJD_*')
    end
    
    if jointPDF_projections_const_conc
        'jointPDF_projections_const_conc'
        load(['data_',num2str(test),'/JD_pdf'])
        startyear=1850; %<<UPDATE 2>>
        endyear=2100;
        %JD_projections_constant_conc
        JD_projections_constant_conc_updated
        save(['data_',num2str(test),'/JD_projections_const_conc_pdf'],'projdat_CC_*','yearall')
    end
    
    if plot_jointPDF_projections
        'plot_jointPDF_projections'
        load(['data_',num2str(test),'/alldat_select_simulations_in_Trange'])
        load(['data_',num2str(test),'/JD_pdf'])
        load(['data_',num2str(test),'/JD_projections_pdf'])
        
        startyear=1750;
        endyear=2100;
        plot_final
    end
    
    if plot_jointPDF_projections_const_conc
        'plot_jointPDF_projectionsconst_conc'
        load(['data_',num2str(test),'/JD_projections_pdf'])
        load(['data_',num2str(test),'/JD_projections_const_conc_pdf'])
        plot_final_const_conc
    end
end

%%

