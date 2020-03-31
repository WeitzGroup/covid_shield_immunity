%% Impact of asymptomatic fraction (p) -- no impact of age on p -- R0 as funciton of p (beta-s constant)
%% David Demory -- March 24 2020
% In this scenario p is constant with age and R0 is dynamic as function of
% p with beta-s constant.

%% Set-up figures
clf
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'fig_asymptomatic_pcst_lowR0cst';
tmpfilebwname = sprintf('%s_v0',tmpfilename);
tmpfilenoname = sprintf('%s_v0',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);
set(gcf,'position', [514 145 500 1000]);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto','color','white');

%% Set-up environment
% Reset
clear stats
clear statsb
clear pars
clear agepars
clear population
clear outbreak

% Population
agepars.meanage=5:10:95;
agepars.highage=[9:10:99];  % Age groups
agepars.lowage=[0:10:90];  % Age groups
population.N=10*10^6;
population.agefrac = [0.12 0.14 0.14 0.13 0.13 0.13 0.10 0.06 0.04 0.01]; 
population.meanage = sum(agepars.meanage.*population.agefrac);
pars.Itrigger = 10000/population.N; % Trigger at 10000 total cases, irrespective of type

% Parameters
pars.gamma_e=1/4;   % Transition to infectiousness
pars.gamma_a=1/4;   % Resolution rate for asymptomatic 
pars.gamma_s=1/4;  % Resolution rate for symptomatic
pars.gamma_h=1/10;  % Resolution rate in hospitals
pars.beta_a=3/10;   % Transmission for asymptomatic
pars.beta_s=6/10;      % Transmission for symptomatic
%pars.p=0.6955;         % Fraction asymptomatic
% pars.p=[0.99 0.99 0.95 0.9 0.8 0.7 0.6 0.5 0.5 0.5];         % Fraction asymptomatic
pars.p=[0.95 0.95 0.9 0.8 0.7 0.5 0.4 0.2 0.2 0.2];         % Fraction asymptomatic
pars.overall_p=sum(pars.p.*population.agefrac);
pars.Itrigger = 10000/population.N; % Trigger at 5000 total cases, irrespective of type

% Epi parameters
pars.Ra=pars.beta_a/pars.gamma_a;
pars.Rs=pars.beta_s/pars.gamma_s;
pars.R0=pars.p*pars.Ra+(1-pars.p)*pars.Rs;
pars.R0=sum(pars.p.*population.agefrac*pars.Ra+(1-pars.p).*population.agefrac*pars.Rs);
% R0 = 2.1

% Constante for R0 from Ia and Is:
cstR0a = pars.overall_p*pars.beta_a/pars.gamma_a;
cstR0s = (1-pars.overall_p)*pars.beta_s/pars.gamma_s;
R0verif_baseline = cstR0a+cstR0s

% Age-stratification
agepars.meanage=5:10:95;
agepars.highage=[9:10:99];  % Age groups
agepars.lowage=[0:10:90];  % Age groups
agepars.hosp_frac=[0.1 0.3 1.2 3.2 4.9 10.2 16.6 24.3 27.3 27.3]/100;
agepars.hosp_crit=[5 5 5 5 6.3 12.2 27.4 43.2 70.9 70.9]/100;
agepars.crit_die= 0.5*ones(size(agepars.meanage));
agepars.num_ages = length(agepars.meanage);
N=agepars.num_ages;
agepars.S_ids=1:N;
agepars.E_ids=(N+1):2*N;
agepars.Ia_ids=(2*N+1):3*N;
agepars.Is_ids=(3*N+1):4*N;
agepars.Ihsub_ids=(4*N+1):5*N;
agepars.Ihcri_ids=(5*N+1):6*N;
agepars.R_ids=(6*N+1):7*N;
agepars.D_ids=(7*N+1):8*N;
agepars.Slock_ids=(8*N+1):9*N;
agepars.ageleave = ones(1,10);
agepars.IFR = sum((1-pars.p).*population.agefrac.*agepars.hosp_frac.*agepars.hosp_crit.*agepars.crit_die);



P = 0.5:0.05:0.95;
Ih    = [];Ib    = [];I    = [];
Ddayh = [];Ddayb = [];Dday = [];
Hdayh = [];Hdayb = [];Hday = [];
DDh   = [];DDb   = [];DD   = [];
RR00 = [];

for i = 1:length(P)
    
    pars.p = P(i);  % Fraction asymptomatic
    pars.beta_a = (cstR0a*pars.gamma_a)/pars.p;
    pars.beta_s = (cstR0s*pars.gamma_s)/(1-pars.p);
    
    % Init the population - baseline
    % Open plus hospitals
    % SEIaIS (open) and then I_ha I_hs and then R (open) and D (cumulative) and
    % then S (lockdown)- 9 categories in total, all age-stratified
    % Here, we ignore the lockdown category
    tmpzeros = zeros(size(agepars.meanage));
    outbreak.y0=[population.agefrac tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros];
    % Initiate an outbreak
    outbreak.y0=population.N*outbreak.y0;
    outbreak.y0(3)=outbreak.y0(3)-1;
    outbreak.y0(13)=1;
    outbreak.y0=outbreak.y0/population.N;
    outbreak.pTime=365;
    
    % Options
    opts=odeset('reltol',1e-8,'maxstep',0.1,'events',@intervene_trigger);
    
    % Sims - Baseline
    % no shielding
    pars.alpha=0;  % Shielding
    [t,y,te,ye,ie]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);
    [tb,yb]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime],ye,opts,pars,agepars);
    statsb.I = yb(:,agepars.Ia_ids)+yb(:,agepars.Is_ids);
    statsb.R=yb(:,agepars.R_ids);
    statsb.D=yb(:,agepars.D_ids);
    statsb.Htot=yb(:,agepars.Ihsub_ids)+yb(:,agepars.Ihcri_ids);
    statsb.Hacu=yb(:,agepars.Ihcri_ids);
    statsb.Hacu_day=sum(statsb.Hacu,2);
    tempDb = sum(statsb.D,2);
    tempItotb = sum(statsb.D,2) + sum(statsb.R,2);
    statsb.finalsize_D(i) = tempDb(end);
    statsb.finalsize_I(i) = tempItotb(end);
    Ib = [Ib,sum(statsb.I,2)];
    DDb = [DDb;statsb.D(end,:)];
    Ddayb = [Ddayb,diff(tempDb)];
    Hdayb = [Hdayb,statsb.Hacu_day];
    
    % Sims - Intervene
    pars.alpha=2;  % Shielding
    [t,y,te,ye,ie]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);
    [t,y]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime],ye,opts,pars,agepars);
    stats.I  = y(:,agepars.Ia_ids)+y(:,agepars.Is_ids);
    stats.R=y(:,agepars.R_ids);
    stats.D=y(:,agepars.D_ids);
    stats.Htot=y(:,agepars.Ihsub_ids)+y(:,agepars.Ihcri_ids);
    stats.Hacu=y(:,agepars.Ihcri_ids);
    stats.Hacu_day=sum(stats.Hacu,2);
    tempD = sum(stats.D,2);
    tempItot = sum(stats.D,2) + sum(stats.R,2);
    stats.finalsize_D(i) = tempD(end);
    stats.finalsize_I(i) = tempItot(end);
    I  = [I,sum(stats.I,2)];
    DD = [DD;stats.D(end,:)];
    Dday  = [Dday,diff(tempD)];
    Hday = [Hday,stats.Hacu_day];
    
    % Sims - Intervene High
    pars.alpha=20;  % Shieldingshg
    [t,y,te,ye,ie]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);
    [th,yh]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime],ye,opts,pars,agepars);
    statsh.I = yh(:,agepars.Ia_ids)+yh(:,agepars.Is_ids);
    statsh.R=yh(:,agepars.R_ids);
    statsh.D=yh(:,agepars.D_ids);
    statsh.Htot=yh(:,agepars.Ihsub_ids)+yh(:,agepars.Ihcri_ids);
    statsh.Hacu=yh(:,agepars.Ihcri_ids);
    statsh.Hacu_day=sum(statsh.Hacu,2);
    tempDh = sum(statsh.D,2);
    tempItoth = sum(statsh.D,2) + sum(statsh.R,2);
    statsh.finalsize_D(i) = tempDh(end);
    statsh.finalsize_I(i) = tempItoth(end);
    Ih = [Ih,sum(statsh.I,2)];
    Ddayh = [Ddayh,diff(tempDh)];
    DDh = [DDh;statsh.D(end,:)];
    Hdayh = [Hdayh,statsh.Hacu_day];
    
    % Epi parameters
    pars.Ra=pars.beta_a/pars.gamma_a;
    pars.Rs=pars.beta_s/pars.gamma_s;
    pars.R0=pars.p*pars.Ra+(1-pars.p)*pars.Rs;
    pars.R0=sum(pars.p.*population.agefrac*pars.Ra+(1-pars.p).*population.agefrac*pars.Rs);
    disp(['p = ',num2str(pars.p),' -- \beta_a = ',num2str(pars.beta_a),' -- \beta_s = ',num2str(pars.beta_s),' -- R0 = ',num2str(pars.R0)])
    disp(['Ia0 = ',num2str(sum(pars.p*10000.*population.agefrac)),' Is0 = ',num2str(sum((1-pars.p)*10000.*population.agefrac))])
    disp(' ')
    
    RR00 = [RR00, pars.R0];
    
end

%% Calcul of reduction
% Death
baseline = statsb.finalsize_D*100000;
(baseline(1)-baseline(end))/baseline(1)
a2 = stats.finalsize_D*100000;
(a2(1)-a2(end))/a2(1)
a20 = statsh.finalsize_D*100000;
(a20(1)-a20(end))/a20(1)
% ICU
baseline = max(Hdayb)*100000;
(baseline(1)-baseline(end))/baseline(1)
a2 = max(Hday)*100000;
(a2(1)-a2(end))/a2(1)
a20 = max(Hdayh)*100000;
(a20(1)-a20(end))/a20(1)

%% Plot

subplot(2,1,1) %Total deaths

yyaxis left
hold on
plot(P,statsb.finalsize_D*100000,'k-','LineWidth',3)%,'MarkerFaceColor','k')
plot(P,stats.finalsize_D*100000,'k--','LineWidth',3)%,'MarkerFaceColor','k')
plot(P,statsh.finalsize_D*100000,'k:','LineWidth',3)%,'MarkerFaceColor','k')

set(gca,'fontsize',16);
xlabel('Fraction of asymptomatic cases ($p$)','fontsize',16,'verticalalignment','top','interpreter','latex','color','k');
ylabel('Total deaths per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex','color','k');
title({'Asymptomatic fraction $p$ constant with age';'${\cal{R}}_0$ = 1.56'},'fontsize',18,'interpreter','latex','color','k')
xlim([0.5 0.95])
ylim([0 600])
set(gca,'ycolor',[0.25 0.25 0.25]);

yyaxis right
plot(P,RR00,'o-','linewidth',3,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
set(gca,'ycolor',[0.25 0.25 0.25]);
ylabel('${\cal{R}}_0$','fontsize',16,'verticalalignment','top','interpreter','latex','Interpreter','Latex');
legend({'Baseline','2:1 Shielding','20:1 Shielding','${\cal{R}}_0$'},'location','Best','Interpreter','latex','NumColumns',1);
legend('boxoff');
ylim([1 2.5])
%set(gca,'ycolor',[0.25 0.25 0.25]);

subplot(2,1,2) %peak ICU
yyaxis left
hold on
plot(P,max(Hdayb)*100000,'k-','LineWidth',3)%,'MarkerFaceColor','k')
plot(P,max(Hday)*100000,'k--','LineWidth',3)%,'MarkerFaceColor','k')
plot(P,max(Hdayh)*100000,'k:','LineWidth',3)%,'MarkerFaceColor','k')

set(gca,'fontsize',16);
xlabel('Fraction of asymptomatic cases ($p$)','fontsize',16,'verticalalignment','top','interpreter','latex','color','k');
ylabel('ICU beds per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex','color','k');
xlim([0.5 0.95])
ylim([0 250])
set(gca,'ycolor',[0.25 0.25 0.25]);

yyaxis right
plot(P,RR00,'o-','linewidth',3,'color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
legend({'Baseline','2:1 Shielding','20:1 Shielding','${\cal{R}}_0$'},'location','Best','Interpreter','latex','NumColumns',1);
legend('boxoff');
%legend({'Baseline','2:1 Shielding','20:1 Shielding','${\cal{R}}_0$'},'Interpreter','latex');
%legend('boxoff');
ylim([1 2.5])
set(gca,'ycolor',[0.25 0.25 0.25]);

%% Print
% for writing over the top
% coordinates are normalized again to (0,1.0)
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
% first two points are normalized x, y positions
% text(,,'','Fontsize',14);

% automatic creation of postscript
% without name/date
psprintc(tmpfilenoname);
psprint(tmpfilebwname);

%tmpt = pwd;
%tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
%text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
%datenamer(1.1,.05,90);
% datename(.5,.05);
% datename2(.5,.05); % 2 rows

% automatic creation of postscript
psprintc(tmpfilename);
print([num2str(tmpfilename),'.eps'],'-depsc')
% set following on if zooming of
% plots is required
% may need to get legend up as well
%axes(tmpa1)
%axes(tmplh)
%clear tmp*
