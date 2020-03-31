% Here we extended the core model to account for the waning of immunity over time by adding 
% two explicit classes of shields: the newly recruited early shield ($RH_1$) and a late stage shield ($RH_2$).
% The recovered population would be recruited into the early shield with a time constant of about 1 week. 
% The early shields would then progress into the late shields before returning to the susceptible class. 
% By using two different temporal classes of shields, we approximate a gamma distribution of immunity duration
% which is more realistic than a simple exponential process 

% The simulations below correspond to the high R0 scenario
% FigS7-left

clear
clc 
close all

% main data goes here
% loglog(,, '');
% RE-STIR Model
% Structure has 3 layers
% Layer 1 - Free to Move
% Layer 2 - Hospitals
% Layer 3 - Shelter in Place 

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

% OLD Parameters
% pars.gamma_e=1/4;   % Transition to infectiousness
% pars.gamma_a=1/6;   % Resolution rate for asymptomatic 
% pars.gamma_s=1/10;  % Resolution rate for symptomatic
% pars.gamma_h=1/10;  % Resolution rate in hospitals
% pars.beta_a=3.5/10;   % Transmission for asymptomatic
% pars.beta_s=7/10;      % Transmission for symptomatic
% pars.p=0.9;         % Fraction asymptomatic
% % Could be age structured
% %pars.p=[0.99 0.99 0.95 0.9 0.8 0.7 0.6 0.5 0.5 0.5];         % Fraction asymptomatic
% pars.overall_p=sum(pars.p.*population.agefrac);
% pars.r=1/7;


% Parameters (03/30/2020)
pars.gamma_e=1/4;   % Transition to infectiousness
pars.gamma_a=1/6;   % Resolution rate for asymptomatic 
pars.gamma_s=1/6;  % Resolution rate for symptomatic
pars.gamma_h=1/10;  % Resolution rate in hospitals
pars.beta_a=3/10;   % Transmission for asymptomatic
pars.beta_s=6/10;      % Transmission for symptomatic
% pars.p=0.9;         % Fraction asymptomatic
% Could be age structured
pars.p=[0.95 0.95 0.9 0.8 0.7 0.6 0.4 0.2 0.2 0.2];         % Fraction asymptomatic
pars.overall_p=sum(pars.p.*population.agefrac);
pars.Itrigger = 10000/population.N; % Trigger at 10000 total cases, irrespective of type

%pars.r=1/7;

% Shield Parameters
avg_immunityduration = 2; % Mean (shield) immunity duration, months
avg_immunityduration = avg_immunityduration*30; 
rate_immunity = 2/avg_immunityduration; % Transition rates from H1 (early shield) -> H2 (late shield) AND H2 -> S (susceptible) are the same
                                         % hence, 1/r_1 + 1/r_2 = avg_immunityduration
indx_shield =  3:6; % Indexes of the age groups (20-60) where immune shields are drawn from
pars.r_shield = zeros(1, length(agepars.meanage)); 
pars.r_shield(indx_shield) = 1/7; % Testing time, and transition from Recovered to Shield, ~1 week
pars.r_1 = zeros(1, length(agepars.meanage)); 
pars.r_1(indx_shield) = rate_immunity; % Transition rate from early (RH1) to late shield (RH2)
pars.r_2 = pars.r_1; % Transition rate from late shield (RH2) to Susceptible (S)


% Epi parameters
pars.Ra=pars.beta_a/pars.gamma_a;
pars.Rs=pars.beta_s/pars.gamma_s;
pars.R0=pars.p*pars.Ra+(1-pars.p)*pars.Rs;
pars.R0=sum(pars.p.*population.agefrac*pars.Ra+(1-pars.p).*population.agefrac*pars.Rs);
% R0 = 2.33

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
agepars.RH1_ids=(9*N+1):10*N;
agepars.RH2_ids=(10*N+1):11*N;
agepars.ageleave = ones(1,10);
% IFR of 0.9% (in line with global)
agepars.IFR = sum((1-pars.p).*population.agefrac.*agepars.hosp_frac.*agepars.hosp_crit.*agepars.crit_die);



% Init the population - baseline
% Open plus hospitals
% SEIaIS (open) and then I_ha I_hs and then R (open) and D (cumulative) and 
% then S (lockdown)- 11 categories in total, all age-stratified
% Here, we ignore the lockdown category
% Joey & Rogelio incorporate two more categories, Early shield RH1 and Late shield (RH2)
tmpzeros = zeros(size(agepars.meanage));
%outbreak.y0=[population.agefrac tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros];
outbreak.y0=[population.agefrac tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros];
% Initiate an outbreak
outbreak.y0=population.N*outbreak.y0;
outbreak.y0(3)=outbreak.y0(3)-1; 
outbreak.y0(13)=1;
outbreak.y0=outbreak.y0/population.N;
outbreak.pTime=365;

% First run without functional shields (alpha = 0) until cases reach 10,000
pars.alpha=0; 
opts=odeset('reltol',1e-8,'maxstep',0.1,'events',@intervene_trigger);
[t,y,te,ye,ie]=ode45(@stir_model_youngshields_explicit,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);


% Then, run Sims - Intervene (alpha = 2)
opts=odeset('reltol',1e-8,'maxstep',0.1);
pars.alpha=2;  % Shielding
[t,y]=ode45(@stir_model_youngshields_explicit,[0:1:outbreak.pTime],ye,opts,pars,agepars);


% Stats
stats.R=y(:,agepars.R_ids);
stats.D=y(:,agepars.D_ids);
stats.Htot=y(:,agepars.Ihsub_ids)+y(:,agepars.Ihcri_ids);
stats.Hacu=y(:,agepars.Ihcri_ids);
stats.Dday_age= stats.D(2:end,:)-stats.D(1:end-1,:);
stats.Dday=sum(stats.Dday_age');
stats.Hacu_day=sum(stats.Hacu');
stats.lock = y(:,agepars.Slock_ids);
stats.shield_early = y(:,agepars.RH1_ids);
stats.shield_late = y(:,agepars.RH2_ids);
stats.S = y(:,agepars.S_ids);
stats.I = y(:,agepars.Ia_ids) + y(:,agepars.Is_ids);


% Sims - Intervene High (alpha = 20)
pars.alpha=20;  % Shielding
[th,yh]=ode45(@stir_model_youngshields_explicit,[0:1:outbreak.pTime],ye,opts,pars,agepars);


% Stats
statsh.R=yh(:,agepars.R_ids);
statsh.D=yh(:,agepars.D_ids);
statsh.Htot=yh(:,agepars.Ihsub_ids)+yh(:,agepars.Ihcri_ids);
statsh.Hacu=yh(:,agepars.Ihcri_ids);
statsh.Dday_age=statsh.D(2:end,:)-statsh.D(1:end-1,:);
statsh.Dday=sum(statsh.Dday_age');
statsh.Hacu_day=sum(statsh.Hacu');
statsh.lock=yh(:,agepars.Slock_ids);
statsh.shield_early = yh(:,agepars.RH1_ids);
statsh.shield_late = yh(:,agepars.RH2_ids);

% Sims - Baseline no shielding (alpha = 0)
pars.alpha=0;
[tb,yb]=ode45(@stir_model_youngshields_explicit,[0:1:outbreak.pTime],ye,opts,pars,agepars);

% Stats
statsb.R=yb(:,agepars.R_ids);
statsb.D=yb(:,agepars.D_ids);
statsb.Htot=yb(:,agepars.Ihsub_ids)+yb(:,agepars.Ihcri_ids);
statsb.Hacu=yb(:,agepars.Ihcri_ids);
statsb.Dday_age=statsb.D(2:end,:)-statsb.D(1:end-1,:);
statsb.Dday=sum(statsb.Dday_age');
statsb.Hacu_day=sum(statsb.Hacu');
statsb.lock=yb(:,agepars.Slock_ids);
statsb.shield_early = yb(:,agepars.RH1_ids);
statsb.shield_late = yb(:,agepars.RH2_ids);


%semilogy(t,y);
%legend('S','E','I1','I2','R','D');

subplot(3,1,1);
tmph=plot(t(2:end),statsb.Dday*100000,'k-');
set(tmph,'linewidth',3);
hold on
tmph=plot(t(2:end),stats.Dday*100000,'k--');
set(tmph,'linewidth',3);
tmph=plot(t(2:end),statsh.Dday*100000,'k:');
set(tmph,'linewidth',3);
set(gca,'fontsize',16);
gc= get(gca);
gc.XAxis.LineWidth = 2;
gc.YAxis.LineWidth = 2;
xlabel('Time, days','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Deaths per day per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
title({'COVID-19 Epidemic - High Scenario - Shields Ages 20-60';'${\cal{R}}_0=2.33$, Mean immunity duration: 2 months'},'fontsize',18,'interpreter','latex');
%title({'COVID-19 Epidemic - High Scenario - Shields Ages 20-60';'Asymptomatic incidence $p=0.9$, ${\cal{R}}_0=2.1$, Mean immunity duration: 2 months'},'fontsize',18,'interpreter','latex')
ylim([0 max(statsb.Dday*100000)+1]);
xlim([0 365]);
%ylim([0 3]);
legend('Baseline','2:1 Shielding','20:1 Shielding');
legend('boxoff');
subplot(3,1,2);
tmph=plot(t,statsb.Hacu_day*100000,'b-');
set(tmph,'linewidth',3);
hold on
tmph=plot(t,stats.Hacu_day*100000,'b--');
set(tmph,'linewidth',3);
tmph=plot(t,statsh.Hacu_day*100000,'b:');
set(tmph,'linewidth',3);
set(gca,'fontsize',16);
xlabel('Time, days','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('ICU beds per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
% title('','fontsize',24)
%ylim([0 60]);
xlim([0 365]);
ylim([0 max(statsb.Hacu_day*100000)+10]);
gc= get(gca);
gc.XAxis.LineWidth = 2;
gc.YAxis.LineWidth = 2;
legend('Baseline','2:1 Shielding','20:1 Shielding');
legend('boxoff');
subplot(3,1,3);
[tmpax tmph1 tmph2]=plotyy(agepars.meanage,statsb.D(end,:)*100000,agepars.meanage,population.agefrac);
set(tmph1,'linewidth',3,'color','k','marker','o','linestyle','-','markerfacecolor','k');
set(tmph2,'linewidth',3,'color',[0.5 0.5 0.5],'marker','d');
axes(tmpax(1));
set(gca,'fontsize',16);
hold on
tmpp=plot(agepars.meanage,stats.D(end,:)*100000,'ko--');
set(tmpp,'linewidth',3);
tmpp=plot(agepars.meanage,statsh.D(end,:)*100000,'ko:');
set(tmpp,'linewidth',3);
xlabel('Age range','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Cumulative deaths per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
set(gca,'ycolor','k');
%ylim([0 55]);
ylim([0 max(statsb.D(end,:)*100000)+10]);
get(gca, 'ytick')
set(gca,'ytick',[0:50:300]);
gc= get(gca);
gc.XAxis.LineWidth = 2;
gc.YAxis.LineWidth = 2;
tmpl=legend('Baseline','2:1 Shielding','20:1 Shielding','Age Structure');
set(tmpl,'location','SouthWest');
legend('boxoff');

axes(tmpax(2));
ylabel('Population age structure','fontsize',16,'verticalalignment','top','interpreter','latex','color','k');
set(gca,'ycolor',[0.25 0.25 0.25]);
ylim([0 0.2]);

set(gcf, 'position', [349    76   726   7295]);
saveas(gcf, './Fig_immune_2m_high', 'fig')
saveas(gcf, './Fig_immune_2m_high','epsc')

clear tmp*

%% Plot Shields dynamics through the year, for baseline, alpha = 2 and alpha = 20 interventions

% figure(2)
% plot(sum(statsb.shield_early,2)*100000, 'linewidth', 2.5, 'color', [0, 0.4470, 0.7410] , 'linestyle', '-')
% hold on
% plot(sum(statsb.shield_late,2)*100000, 'linewidth', 2.5, 'color',	[0, 51/255, 153/255] , 'linestyle', '-')
% plot(sum(stats.shield_early,2)*100000, 'linewidth', 2.5, 'color', [0.8500, 0.3250, 0.0980],'linestyle', '--')
% plot(sum(stats.shield_late,2)*100000, 'linewidth', 2.5, 'color',[179/255, 89/255, 0] ,'linestyle', '--')
% plot(sum(statsh.shield_early,2)*100000, 'linewidth', 2.5, 'color', [0.9290, 0.6940, 0.1250]	,'linestyle', ':')
% plot(sum(statsh.shield_late,2)*100000, 'linewidth', 2.5, 'color',[153/255, 128/255, 0.1250]	 ,'linestyle', ':')
% hold off
% legend('Early Shield baseline', 'Late Shield baseline', 'Early Shield, $\alpha$ = 2', 'Late Shield $\alpha$ = 2', 'Early Shield $\alpha$ = 20', 'Late Shield $\alpha$ = 20', 'interpreter','latex')
% legend box off
% xlabel('Time, days','fontsize',16,'verticalalignment','top','interpreter','latex');
% ylabel('Shields per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
% title({'COVID-19 Epidemic - High Scenario - Shields Ages 20-60';'Asymptomatic incidence $p=0.9$, ${\cal{R}}_0=1.85$, Immunity for 2 months'},'fontsize',18,'interpreter','latex')
% xlim([0 365]);
% set(gca, 'fontweight', 'bold', 'fontsize', 12)





