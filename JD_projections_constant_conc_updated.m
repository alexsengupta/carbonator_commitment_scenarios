% Set Input parameters

% time step (years)
DT=1/12; % needs to be small otherwise ocean chemistry calculations become unstable
years=startyear:DT:endyear;
ind1870=find(years>1861 & years<=1880);
ind2012=find(years>2006 & years<=2018);
mRF_AER_2006_2018=-0.89;

no_iterations=length(datJD_g);

projdat_CC_CCO2=zeros(CparameterSet,3,length(years(1:12:end)))*NaN;
projdat_CC_RCO2=zeros(CparameterSet,3,length(years(1:12:end)))*NaN;
projdat_CC_CCH4=zeros(3,length(years(1:12:end)))*NaN;
projdat_CC_RCH4=zeros(3,length(years(1:12:end)))*NaN;
projdat_CC_RCOZO=zeros(3,length(years(1:12:end)))*NaN;
projdat_CC_RCSO2=zeros(no_iterations,3,length(years(1:12:end)))*NaN;
projdat_CC_RHalo=zeros(3,length(years(1:12:end)))*NaN;
projdat_CC_RN2O=zeros(3,length(years(1:12:end)))*NaN;
projdat_CC_RF=zeros(no_iterations,3,length(years(1:12:end)))*NaN;
projdat_CC_To=zeros(no_iterations,3,length(years(1:12:end)))*NaN;
projdat_CC_Ts=zeros(no_iterations,3,length(years(1:12:end)))*NaN;

% load model parameters and initial conditions
for n=1:CparameterSet
%     [A,AF,AM,Alk, Cat,Clo,Co,Cs,Cup,Ksp,L,N,OM,P,So,VF,a2,alb0,alphaRF_ch4,alpha_ch4,d,dr0,e,eps,g,k1,k2,ka,kd,kh,m,tau_ch4_pi,vtau]=load_constants(n,fnC);
    [A,AF,AM,Alk, Cat,Clo,Co,Cs,Cup,Ksp,L,N,OM,P,So,VF,a2,alb0,alphaRF_ch4,alpha_ch4,d,dr0,e,eps,g,k1,k2,ka,kd,kh,m,tau_ch4_pi,vtau,tau_N20,alphaRF_N2O,tau_CFC12,alphaRF_CFC12]=load_constants(n,fnC);%<<UPDATE 2>>
    A_cm(n)=A;
    a2_cm(n)=a2;
    d_cm(n)=d;
    dr0_cm(n)=dr0;
    e_cm(n)=e;
    ka_cm(n)=ka;
    kd_cm(n)=kd;
    m_cm(n)=m;
    %Initial conditions - OPTIMISED
    Clo_cm(n)=Clo;
    Cup_cm(n)=Cup;
    Cat_cm(n)=Cat;
    N_cm(n)=N;
    P_cm(n)=P;
    So_cm(n)=So;
end

clear conc_*  emission*_*
for myscen=1:3
    switch myscen
        case 1
            scenario='RCP85';
        case 2
            scenario='RCP85';
        case 3
            scenario='RCP85';
    end
    [years,meinhausenRF_OZO_(myscen,:),meinhausenRF_AER_(myscen,:),TSI,alb,emissionCH4_(myscen,:),emissionCO2_(myscen,:),emissionN2O_(myscen,:),emissionSO2_(myscen,:),emission_volc,radf_other,...
        radf_halo,radf_N2O,radf_CO2,radf_CH4,conc_CO2,conc_CH4,conc_N2O(myscen,:),mTSI,conc_Hal(myscen,:)]=load_RCP(scenario,startyear,endyear);%<<UPDATE 2>>
end

for NI=1:no_iterations
    disp(['ITERATION....',num2str(NI),'/',num2str(no_iterations)])
    
    %switches for forcing
    F1=1; %CO2
    F2=1; %CH4
    F3=1; %SO2
    F4=1; %Volcanics
    F5=1; %Solar
    F6=1; %ozone
    
    %     default_constants_optimised
    
    g=datJD_g(NI);
    Cs=datJD_Cs(NI);
    Co=datJD_Co(NI);
    L=datJD_L(NI);
    sf=datJD_sf(NI);
    Cparam=datJD_Cparam(NI);
    
    A=A_cm(Cparam);
    a2=a2_cm(Cparam);
    d=d_cm(Cparam);
    dr0=dr0_cm(Cparam);
    e=e_cm(Cparam);
    ka=ka_cm(Cparam);
    kd=kd_cm(Cparam);
    m=m_cm(Cparam);
    %Initial conditions - OPTIMISED
    Clo=Clo_cm(Cparam);
    Cup=Cup_cm(Cparam);
    Cat=Cat_cm(Cparam);
    N=N_cm(Cparam);
    P=P_cm(Cparam);
    So=So_cm(Cparam);
    
    
    for myscen=1:3
        switch myscen
             case 1
                scenario='COM20cc';const_year=2020;
            case 2
                scenario='COM30cc';const_year=2030;
            case 3
                scenario='COM40cc';const_year=2040;
        end
        emissionSO2=emissionSO2_(myscen,:);
        emissionCH4=emissionCH4_(myscen,:);
        emissionCO2=emissionCO2_(myscen,:);
        
        R_OZO_=meinhausenRF_OZO_(myscen,:);
        R_AER=meinhausenRF_AER_(myscen,:); %<<UPDATE>>
        R_AER=(R_AER-mean(R_AER(ind1870)))*sf/mRF_AER_2006_2018 + mean(R_AER(ind1870)); %<<UPDATE>>
        R_AER=R_AER-R_AER(1);  %<<UPDATE 2>>
        R_OZO_=R_OZO_-R_OZO_(1);
        
        
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
        
        %<<UPDATE 2>>
        C_N2Opi=275.4251; %UPDATE5 % in 1750 275.4 in 1850; % PI concentration
        C_N2O(1)=0;
        R_N2O(1)=0;
        
        R_Halo(1)=0;
        C_Halo(1)=0;
        

        
        %% Loop through timesteps
        ind_const_year=find(years==const_year);
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
            if years(y)<=const_year
                C_CO2(y)=Cat(y)/2.13;
            else
                C_CO2(y)=C_CO2(ind_const_year);
            end
            R_CO2(y)=5.35*log((C_CO2(y))/C_CO2pi);
            
            % convert CH4 emissions to RF
            T=tau_ch4_pi*( C_CH4pi/(C_CH4(y-1)+C_CH4pi) )^alpha_ch4;
            if years(y)<=const_year
                C_CH4(y)=C_CH4(y-1)  + DT*( F2*emissionCH4(y-1)/2.78 - (1/T)*C_CH4(y-1) );
            else
                C_CH4(y)=C_CH4(ind_const_year);
            end
            R_CH4(y)=alphaRF_ch4*(sqrt(C_CH4pi+C_CH4(y))-sqrt(C_CH4pi)); % based on RF proportional sqrt(CH4)-sqrt(CH4pi) see IPCC TAR chp6
            
            % anthropogenic aerosols
            if years(y)<=const_year
                R_SO2(y)=R_AER(y);%<<UPDATED>>
            else
                R_SO2(y)=R_SO2(ind_const_year);
            end
            
            % N2O %<<UPDATE 2>>
            if years(y)<=const_year
                C_N2O(y)=conc_N2O(myscen,y)-C_N2Opi;
            else
                C_N2O(y)=C_N2O(ind_const_year);
            end
            R_N2O(y)=alphaRF_N2O*(sqrt(C_N2Opi+C_N2O(y))-sqrt(C_N2Opi));
            
            if years(y)<=9999999999; %const_year
                C_Halo(y)=conc_Hal(myscen,y);
            else
                C_Halo(y)=conc_Hal(myscen,ind_const_year);
            end
            R_Halo(y)=alphaRF_CFC12/1000*C_Halo(y); %(convert to ppb)
            
            
            % volcanic aerosols
            OpTkhkV(y)=OpTkhkV(y-1) + DT*( emission_volc(y) - OpTkhkV(y-1)/vtau );
            R_volc(y)= VF*OpTkhkV(y);
            
            % Ozone
            if years(y)<=const_year
            R_OZO(y)=R_OZO_(y);
            else
                R_OZO(y)=R_OZO(ind_const_year);
            end
            
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
        X(NI,myscen)=deltaT;
        
        if mod(NI,1000)==0 %only need to run this a few times to make sure I sample Cparam 
            projdat_CC_CCO2(Cparam,myscen,:)=C_CO2(1:12:end);
            projdat_CC_RCO2(Cparam,myscen,:)=R_CO2(1:12:end);
        end
        if NI==1
            %sanity check
            figure(1);
            subplot(2,2,1)
            plot(years,C_CO2);hold on; title(['CO2 delta=',num2str(mean(C_CO2(ind2012))-mean(C_CO2(ind1870)))])
            text(1870,myscen*70+300,num2str(mean(C_CO2(ind2012))-mean(C_CO2(ind1870))))
            subplot(2,2,2)
            plot(years,C_CH4);hold on; title(['CH4 delta=',num2str(mean(C_CH4(ind2012))-mean(C_CH4(ind1870)))])
            subplot(2,2,3)
            plot(years,C_N2O);hold on; title(['N2O delta=',num2str(mean(C_N2O(ind2012))-mean(C_N2O(ind1870)))])
            subplot(2,2,4)
            plot(years,C_Halo);hold on; title(['Halo delta=',num2str(mean(C_Halo(ind2012))-mean(C_Halo(ind1870)))])
            
            figure(2);
            subplot(3,3,1)
            plot(years,R_CO2);hold on; title(['CO2 delta=',num2str(mean(R_CO2(ind2012))-mean(R_CO2(ind1870)))])
            text(1870,myscen*.7,num2str(mean(R_CO2(ind2012))-mean(R_CO2(ind1870))))
            subplot(3,3,2)
            plot(years,R_CH4);hold on; title(['CH4 delta=',num2str(mean(R_CH4(ind2012))-mean(R_CH4(ind1870)))])
            subplot(3,3,3)
            plot(years,R_N2O);hold on; title(['N2O delta=',num2str(mean(R_N2O(ind2012))-mean(R_N2O(ind1870)))])
            subplot(3,3,4)
            plot(years,R_Halo);hold on; title(['Halo delta=',num2str(mean(R_Halo(ind2012))-mean(R_Halo(ind1870)))])
            subplot(3,3,5)
            plot(years,R_OZO);hold on; title(['Ozone delta=',num2str(mean(R_OZO(ind2012))-mean(R_OZO(ind1870)))])
            subplot(3,3,6)
            plot(years,R_sol);hold on; title(['Solar delta=',num2str(mean(R_sol(ind2012))-mean(R_sol(ind1870)))])
            subplot(3,3,7)
            plot(years,R_volc);hold on; title(['Volcanic delta=',num2str(mean(R_volc(ind2012))-mean(R_volc(ind1870)))])
            subplot(3,3,8)
            plot(years,R_SO2);hold on; title('Aerosol')
            
            figure(3)
            plot(years,Ts);hold on; title('Temperature')
            plot(years,To)
            
            projdat_CC_CCH4(myscen,:)=C_CH4(1:12:end);
            projdat_CC_CN2O(myscen,:)=C_N2O(1:12:end);
            projdat_CC_CHalo(myscen,:)=C_Halo(1:12:end);
            
            projdat_CC_RCH4(myscen,:)=R_CH4(1:12:end);
            projdat_CC_RN2O(myscen,:)=R_N2O(1:12:end);
            projdat_CC_RHalo(myscen,:)=R_Halo(1:12:end);
            projdat_CC_RCOZO(myscen,:)=R_OZO(1:12:end);
            projdat_CC_RVOLC(myscen,:)=R_volc(1:12:end);
            projdat_CC_RSOL(myscen,:)=R_sol(1:12:end);
        end
       
        projdat_CC_RCSO2(NI,myscen,:)=R_SO2(1:12:end);
        projdat_CC_RF(NI,myscen,:)=RF_(1:12:end);
        projdat_CC_To(NI,myscen,:)=To(1:12:end);
        projdat_CC_Ts(NI,myscen,:)=Ts(1:12:end);
        
        %%CCO2_Cparam(Cparam,myscen,:)=C_CO2(1:12:end);
        
        Cup(2:end)=[];
        Cat(2:end)=[]; 
        Clo(2:end)=[]; 
        P(2:end)=[]; 
        So(2:end)=[]; 
        N(2:end)=[];
        clear R_* C_* OpTkhkV To Ts Ts deltaT RF_ aragonite_saturation 
    end

end

yearall=years(1:12:end);



