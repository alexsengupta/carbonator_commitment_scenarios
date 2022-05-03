% Set Input parameters

% time step (years)
DT=1/12; % needs to be small otherwise ocean chemistry calculations become unstable

ind1870=find(years>1861 & years<=1880);
ind2012=find(years>2006 & years<=2018);
mRF_AER_2006_2018=-0.89;


for NI=1:CparameterSet
    
    disp(['ITERATION....',num2str(NI),'/',num2str(CparameterSet)])
    
    %switches for forcing
    F1=1; %CO2
    F2=1; %CH4
    F3=1; %SO2
    F4=1; %Volcanics
    F5=1; %Solar
    F6=1; %ozone
    
    nn=no_of_good_simulations*10;
    % Load in Sherwood joint distribution pairs
    [sw_ecs,sw_arf]=sherwood_ecs_arf(nn); % choose ARF and climate sensitivity from sherwoord jPDF
    if uniform_dist==1 % choose ARF and climate sensitivity from uniform distibution with ranges from sherwoord jPDF
        sw_ecs=rand(nn,1)*(max(sw_ecs)-min(sw_ecs))+min(sw_ecs);
        sw_arf=rand(nn,1)*(max(sw_arf)-min(sw_arf))+min(sw_arf);
    end
    
    % load model parameters and initial conditions
    [A,AF,AM,Alk, Cat,Clo,Co,Cs,Cup,Ksp,L,N,OM,P,So,VF,a2,alb0,alphaRF_ch4,alpha_ch4,d,dr0,e,eps,g,k1,k2,ka,kd,kh,m,tau_ch4_pi,vtau,tau_N20,alphaRF_N2O,tau_CFC12,alphaRF_CFC12]=load_constants(NI,fnC);%<<UPDATE 2>>
    
    % Model iterations
    counter=1;
    counter2=1;
    
    while counter<=no_of_good_simulations
        
        if mod(counter2,100)==0;
            disp([num2str(counter),' -- ',num2str(counter2)])
        end
        
        R_OZO=meinhausenRF_OZO;
        R_AER=meinhausenRF_AER;
        if randomise_parameters
            g=rand*diff(range_g)+range_g(1);
            Cs=rand*diff(range_Cs)+range_Cs(1);
            Co=rand*diff(range_Co)+range_Co(1);
            L =5.35*log(2)/sw_ecs(counter2);
            sf=sw_arf(counter2); %<<UPDATE>>
            R_AER=(R_AER-mean(R_AER(ind1870)))*sf/mRF_AER_2006_2018 + mean(R_AER(ind1870));
            R_AER=R_AER-R_AER(1);  %<<UPDATE 2>>
            R_OZO=R_OZO-R_OZO(1);  %<<UPDATE 3>>
        end
      
       
        % Temperature (Geoffroy)
        Ts(1)=0;
        Ts_(1)=0;
        To(1)=0;
        deltaT(1)=0;
        
        C_CO2pi=Cat(1)/2.13; % % PI concentration; NB glotter used 285 (but need this value for equilibrium)
        C_CO2(1)=C_CO2pi;
        R_CO2(1)=0;
        C_CH4pi=791; % PI concentration
        C_CH4(1)=0;
        R_CH4(1)=0;
        R_SO2(1)=0;
        OpTkhkV(1)=0;
        
        C_N2Opi=275.4251; %UPDATE5 ;% in 1750 275.4 in 1850; % PI concentration
        C_N2O(1)=0;
        R_N2O(1)=0;
        
        R_Halo(1)=0;
        C_Halo(1)=0;
        
        stop_year=999999;
        %% Loop through timesteps
        
        for y=2:length(years)
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
            
            % Carbon budget
            Cat(y) = Cat(y-1) + F1*DT*emissionCO2(y-1) + F1*DT*( - ka* (Cat(y-1) - A*B*Cup(y-1) ) ) ... %glotter
                + F1*DT*( -P(y-1) +(1-e)*m*N(y-1) + dr0*So(y-1) );        %svirezhev
            Cup(y) = Cup(y-1)                       + F1*DT*( + ka* (Cat(y-1) - A*B*Cup(y-1) ) - kd*( Cup(y-1) - Clo(y-1)/d ) );
            Clo(y) = Clo(y-1) +                       F1*DT*(                                    kd*( Cup(y-1) - Clo(y-1)/d ) );
            %           if years(y)>2000;pause;end
            % convert CO2 emissions to RF
            C_CO2(y)=Cat(y)/2.13;
            R_CO2(y)=5.35*log((C_CO2(y))/C_CO2pi);
            
            % convert CH4 emissions to RF
            T=tau_ch4_pi*( C_CH4pi/(C_CH4(y-1)+C_CH4pi) )^alpha_ch4;
            C_CH4(y)=C_CH4(y-1)  + DT*( F2*emissionCH4(y-1)/2.78 - (1/T)*C_CH4(y-1) );
            %R_CH4(y)=0.66*log((C_CH4pi+C_CH4(y))/C_CH4pi);
            R_CH4(y)=alphaRF_ch4*(sqrt(C_CH4pi+C_CH4(y))-sqrt(C_CH4pi)); % based on RF proportional sqrt(CH4)-sqrt(CH4pi) see IPCC TAR chp6
            
            % anthropogenic aerosols
            R_SO2(y)=R_AER(y);%<<UPDATED>>
            
            % N2O %<<UPDATE 2>>
            if years(y)<stop_year
                C_N2O(y)=conc_N2O(y)-C_N2Opi;
            else
                C_N2O(y)=C_N2O(y-1)  + DT*( - (1/tau_N20)*C_N2O(y-1) );
            end
            R_N2O(y)=alphaRF_N2O*(sqrt(C_N2Opi+C_N2O(y))-sqrt(C_N2Opi));
            
            % Halogen %<<UPDATE 2>>
            C_Halo(y)=conc_Hal(y);
%             if years(y)<stop_year
%                 C_Halo(y)=conc_Hal(y);
%             else
%                 C_Halo(y)=C_Halo(y-1)  + DT*( - (1/tau_CFC12)*C_Halo(y-1) );
%             end
            R_Halo(y)=alphaRF_CFC12/1000*C_Halo(y); %(convert to ppb)
            
            % volcanic aerosols
            OpTkhkV(y)=OpTkhkV(y-1) + DT*( emission_volc(y) - OpTkhkV(y-1)/vtau );
            R_volc(y)= VF*OpTkhkV(y);
            
            % Ozone
            %R_OZO
            
            % Solar
            R_sol(y)=((TSI(y)-mTSI)/4)*(1-alb(y)) - TSI(y)/4*(alb(y)-alb0);
            
            %temperature model
            RF_GG = F1*R_CO2(y-1) + F2*R_CH4(y-1) + F6*R_OZO(y-1) ...
                +R_N2O(y-1) + R_Halo(y-1); % %<<UPDATE 2>>
            
            RF_(y) = RF_GG + F3*R_SO2(y-1) + F4*R_volc(y-1) + F5*R_sol(y);
            Ts_(y)=Ts_(y-1)+ DT*( RF_(y-1) - L*Ts_(y-1) - g*(Ts_(y-1)-To(y-1)) )/Cs ;
            To(y)=To(y-1)+ DT*( g*(Ts_(y-1)-To(y-1)) )/Co;
            Ts(y)=Ts_(y)*DTscale;
            
        end
        
        % 1861–1880 to 2006–2018
        %RF_SO2_2006_2018=mean(R_SO2(ind2))-mean(R_SO2(ind1))
        
        %calculate RF_SO2_2006_2018, need to set sf=1 prior to main loop
        T1870=mean(Ts(ind1870));
        T2012=mean(Ts(ind2012));
        deltaT=T2012-T1870;
        
        
        if randomise_parameters
            alldat_g(counter2)=g;
            alldat_L(counter2)=L;
            alldat_Co(counter2)=Co;
            alldat_Cs(counter2)=Cs;
            alldat_RF(counter2)=mean(RF_(ind2012))-mean(RF_(ind1870));
            alldat_RFA(counter2)=mean(R_SO2(ind2012))-mean(R_SO2(ind1870));
            alldat_sf(counter2)=sf;
            alldat_deltaT(counter2)=deltaT;
            alldat_Cparam(counter2)=NI;
            counter2=counter2+1;
            
            if deltaT>Trange(1) & deltaT<Trange(2)
                dat_g(counter)=g;
                dat_L(counter)=L;
                dat_Co(counter)=Co;
                dat_Cs(counter)=Cs;
                dat_RF(counter)=mean(RF_(ind2012))-mean(RF_(ind1870));
                dat_RFA(counter)=mean(R_SO2(ind2012))-mean(R_SO2(ind1870));
                dat_sf(counter)=sf;
                dat_deltaT(counter)=deltaT;
                dat_Cparam(counter)=NI;
                
                
                if plotRF
                    if counter==1
                        figure(100);clf
                        mean(R_CO2(ind2012))-mean(R_CO2(ind1870))
                        mean(R_CH4(ind2012))-mean(R_CH4(ind1870))
                        mean(R_N2O(ind2012))-mean(R_N2O(ind1870))
                        mean(R_Halo(ind2012))-mean(R_Halo(ind1870))
                        mean(R_sol(ind2012))-mean(R_sol(ind1870))
                        mean(R_volc(ind2012))-mean(R_volc(ind1870))
                        mean(R_OZO(ind2012))-mean(R_OZO(ind1870))
                        
                        
                        plot(years,R_CO2,'r')
                        hold on
                        plot(years,R_CH4,'b')
                        plot(years,R_volc,'g')
                        plot(years,R_sol,'k')
                        pbaspect([2.5 1 1])
                        title('Radiative Forcing')
                        set(gca,'fontsize',12)
                        xlim([1850 2018])
                        plot(years,R_SO2,'color',[.7 .7 .7])
                    end
                end
                counter=counter+1;
            end
        end
        Cup(2:end)=[];
        Cat(2:end)=[];
        Clo(2:end)=[];
        P(2:end)=[];
        So(2:end)=[];
        N(2:end)=[];
        clear R_* C_* OpTkhkV To Ts_ Ts deltaT RF_ aragonite_saturation
    end
    
    
%     CsdT/dt=RF-LT-g(T-To)
%     At equlibrium, subject to RF associated with a doubleing of CO2 (dRF)
%     dT/dt=0, T=To
%     So, equilibriam temperature change = dFR/L
    
    % Need to multiply DTscale to obtain SAT climate sensitivity
    dat_climate_sensitivity =DTscale*( 5.35*log(2)./dat_L );
    alldat_climate_sensitivity=DTscale*( 5.35*log(2)./alldat_L );
    
    % compute TCR
    F2=0;
    F3=0;
    F4=0;
    F5=0;
    
    NN=length(dat_g);
    for counter=1:length(dat_g)
        default_constants_optimised
        g=dat_g(counter);
        L =dat_L(counter);
        Co=dat_Co(counter);
        Cs=dat_Cs(counter);
        sf=dat_sf(counter);
        R_OZO=meinhausenRF_OZO;
        R_AER=meinhausenRF_AER;
        R_AER=(R_AER-mean(R_AER(ind1870)))*sf/mRF_AER_2006_2018 + mean(R_AER(ind1870)); %<<UPDATE>>
        
        % Temperature (Geoffroy)
        Ts(1)=0;
        Ts_(1)=0;
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
        
        % Loop through timesteps
        y=1;
        while C_CO2(end)<=C_CO2(1)*2
            y=y+1;
            
            r=nthroot(1.01,12);
            C_CO2(y)=C_CO2(y-1)*r;
            R_CO2(y)=5.35*log((C_CO2(y))/C_CO2pi);
            
% %             % convert CH4 emissions to RF
% %             T=tau_ch4_pi*( C_CH4pi/(C_CH4(y-1)+C_CH4pi) )^alpha_ch4;
% %             C_CH4(y)=C_CH4(y-1)  + DT*( F2*emissionCH4(y-1)/2.78 - (1/T)*C_CH4(y-1) );
% %             R_CH4(y)=0.66*log((C_CH4pi+C_CH4(y))/C_CH4pi);
% %             
% %             % anthropogenic aerosols
% %             R_SO2(y)=R_AER(y);%<<UPDATED>>
% %             %             R_SO2(y)=AF*(emissionSO2(y-1)-emissionSO2(1));
% %             %             R_SO2(y)=R_SO2(y)*sf; % scale anthropogenic aerosols <<<<<<
% %             
% %             % volcanic aerosols
% %             OpTkhkV(y)=OpTkhkV(y-1) + DT*( emission_volc(y) - OpTkhkV(y-1)/vtau );
% %             R_volc(y)= VF*OpTkhkV(y);
% %             
% %             
% %             % Solar
% %             R_sol(y)=((TSI(y)-mTSI)/4)*(1-alb(y)) - TSI(y)/4*(alb(y)-alb0);
% %             
% %             %temperature model
% %             RF_GG = (F1*R_CO2(y-1) + F2*R_CH4(y-1))  ;
% %             RF_(y) = RF_GG + F3*R_SO2(y-1) + F4*R_volc(y-1) + F5*R_sol(y);
            
            RF_GG = R_CO2(y-1);
            RF_(y) = RF_GG;
            Ts_(y)=Ts_(y-1)+ DT*( RF_(y-1) - L*Ts_(y-1) - g*(Ts_(y-1)-To(y-1)) )/Cs ;
            To(y)=To(y-1)+ DT*( g*(Ts_(y-1)-To(y-1)) )/Co;
            Ts(y)=Ts_(y)*DTscale;
        end
        
        dat_TCR(counter)=Ts(end);
        
        if mod(counter,100)==0;
            disp(['TCR loop ',num2str(counter),' / ',num2str(NN)])
        end
        
        Cup(2:end)=[];
        Cat(2:end)=[];
        Clo(2:end)=[];
        P(2:end)=[];
        So(2:end)=[];
        N(2:end)=[];
        clear R_* C_* OpTkhkV To Ts Ts deltaT RF_ aragonite_saturation
        
    end

    save(['data_',num2str(test),'/dat_select_simulations_in_Trange_',num2str(NI)],'dat_*','alldat_*');
    clear dat_* alldat_*
end


