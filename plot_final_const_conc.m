% plot_final_const_conc.m

figure(1006);clf
cols={'g','b','k','r','b--','k--','r--'}
for myscen=1:3
    figure(1006);
    plot(yearall,squeeze(projdat_CC_RCO2(1,myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('CO2');pbaspect([1.3 1 1])
    plot(yearall,squeeze(projdat_CC_RCO2(2,myscen,:)),cols{myscen})
    plot(yearall,squeeze(projdat_CC_RCO2(3,myscen,:)),cols{myscen})
    plot(yearall,squeeze(projdat_CC_RCO2(4,myscen,:)),cols{myscen})
    figure(1007);
    plot(yearall,squeeze(projdat_CC_RCH4(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('CH4');pbaspect([1.3 1 1]);ylim([0 1])
    figure(1008);
    plot(yearall,squeeze(projdat_CC_RN2O(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('N2O');pbaspect([1.3 1 1]);ylim([0 1])
    figure(1009);
    plot(yearall,squeeze(projdat_CC_RCOZO(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('O3');pbaspect([1.3 1 1]);ylim([0 1])
    figure(1010);
    plot(yearall,squeeze(projdat_CC_RHalo(myscen,:)),cols{myscen});hold on;set(gca,'fontsize',15);title('Halogens');pbaspect([1.3 1 1]);ylim([0 1])
    figure(1011);
    plot(yearall,squeeze(prctile(projdat_CC_RCSO2(:,myscen,:),[ 50 ])),cols{myscen});hold on;set(gca,'fontsize',15);title('Aerosol (median)');pbaspect([1.3 1 1]);ylim([-1 0])
end

print -f1006 -depsc 'figures/all_radiative_CC_forcingCO2'
print -f1007 -depsc 'figures/all_radiative_CC_forcingCH4'
print -f1008 -depsc 'figures/all_radiative_CC_forcingN2O'
print -f1009 -depsc 'figures/all_radiative_CC_forcingO3'
print -f1010 -depsc 'figures/all_radiative_CC_forcingHal'
print -f1011 -depsc 'figures/all_radiative_CC_forcingAer'

ind=find(datJD_Cparam==1);
figure(1);clf
for s=1:3
    plot(yearall,squeeze(prctile(projdat_CC_Ts(ind,s,:),[5 50 95])),cols{s});hold on
end

ind2020=find(yearall==2020);
projdat_CC_Ts_2020=zeros(length(ind),3,length(yearall));
for s=1:3
    for n=1:length(ind)
        projdat_CC_Ts_2020(n,s,:)=squeeze(projdat_CC_Ts(ind(n),s,:))-squeeze(projdat_CC_Ts(ind(n),s,ind2020));
    end
end

figure(2);clf
for s=1:3
    plot(yearall,squeeze(prctile(projdat_CC_Ts_2020(ind,s,:),[5 50 95])),cols{s});hold on
end
xlim([2010 2100])

SATcommit=projdat_Ts(:,1:3,:)*NaN;
% SATmax=projdat_Ts(:,1:3)*NaN;
% SATmax_year=projdat_Ts(:,1:3)*NaN;
for n=1:length(projdat_Ts)
    n
    for s=5:7
        SATcommit(n,s-4,:)= (projdat_Ts(n,s,:)-projdat_Ts(n,s,ind2020));
%         SATmax(n,s-4)= max(projdat_Ts(n,s,:));
%         ind=find(projdat_Ts(n,s,:)==SATmax(n,s-4));SATmax_year(n,s-4)=yearall(ind);
%         SATmax2020(n,s-4)= max(projdat_Ts(n,s,:))-projdat_Ts(n,s,ind2020);
    end
end

plot(yearall,squeeze(prctile(SATcommit(ind,1,:),[50])))
plot(yearall,squeeze(prctile(SATcommit(ind,2,:),[50])))
plot(yearall,squeeze(prctile(SATcommit(ind,3,:),[50])))
set (gca,'fontsize',16)
pbaspect([1.7 1 1])
print -depsc -f2 'figures/committment_constantConc'

