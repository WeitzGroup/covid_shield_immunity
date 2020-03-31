% Here, we show how the Cumulative Deaths, Max ICU beds (per 100,000), and
% Total Cases change as a function of the waning immunity term (i.e., mean immunity
% duration).

% The simulations below correspond to the high R0 scenario
% FigS8-left

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
pars.Itrigger = 10000/population.N; % Trigger at 5000 total cases, irrespective of type

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
agepars.cases_ids=(11*N+1):12*N; % We add a mock equation to track the total exposed cases at the end of the simulation
agepars.ageleave = ones(1,10);
% IFR of 0.9% (in line with global)
agepars.IFR = sum((1-pars.p).*population.agefrac.*agepars.hosp_frac.*agepars.hosp_crit.*agepars.crit_die);


% Init the population - baseline
% Open plus hospitals
% SEIaIS (open) and then I_ha I_hs and then R (open) and D (cumulative) and 
% then S (lockdown)- 11 categories in total, all age-stratified
% Here, we ignore the lockdown category
% Joey & Rogelio also incorporate two more categories, Early shield (RH1) and Late shield (RH2)
tmpzeros = zeros(size(agepars.meanage));
outbreak.y0=[population.agefrac tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros tmpzeros];
% Initiate an outbreak
outbreak.y0=population.N*outbreak.y0;
outbreak.y0(3)=outbreak.y0(3)-1; 
outbreak.y0(13)=1;
outbreak.y0=outbreak.y0/population.N;
outbreak.pTime=365;

% Shield parameters
indx_shield =  3:6; % Indexes of the ages (20-60) where immune shields are drawn from
pars.r_shield = zeros(1, length(agepars.meanage)); 
pars.r_shield(indx_shield) = 1/7; % Testing time, from Recovered to Shield, ~1 week

mean_immuneduration = 60:7:365; % immunity duration vector, days
ye_mat = zeros(length(mean_immuneduration), 120);
opts=odeset('reltol',1e-8,'maxstep',0.1,'events',@intervene_trigger);
% Find the starting point for intervention, i.e., when cases reach 10,000
for i = 1:length(mean_immuneduration)
    % Shield parameters
    avg_immunityduration = mean_immuneduration(i); 
    rate_immunity = 2/avg_immunityduration; % Transition rates from H1 (early shield) -> H2 (late shield) & H2 -> S (susceptible) are the same
                                             %hence, 1/r_1 + 1/r_2 = avg_immunityduration
    pars.r_1 = zeros(1, length(agepars.meanage)); 
    pars.r_1(indx_shield) = rate_immunity; % Transition rate from early (RH1) to late shield (RH2)
    pars.r_2 = pars.r_1; % Transition rate from late shield (RH2) to Susceptible (S)

    % save 'ye' for later use
    pars.alpha=0;  % Shielding
    [tb,yb,te,ye,ie]=ode45(@stir_model_youngshields_exposedcases,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);
    ye_mat(i,:)= ye;
end

save ye_mat.mat ye_mat
load ye_mat.mat

% We test different mean immunity durations ranging from 2-12 months
mean_immuneduration = 60:7:365; % duration vector, days
alpha_vals = [0 2 20];
opts=odeset('reltol',1e-8,'maxstep',0.1);
tic
for i = 1:length(alpha_vals)
    
    for n =1:length(mean_immuneduration)

        % Shield parameters
        avg_immunityduration = mean_immuneduration(n); 
        rate_immunity = 2/avg_immunityduration; % Transition rates from H1 (early shield) -> H2 (late shield) & H2 -> S (susceptible) are the same
                                                 %hence, 1/r_1 + 1/r_2 = avg_immunityduration
        pars.r_1 = zeros(1, length(agepars.meanage)); 
        pars.r_1(indx_shield) = rate_immunity; % Transition rate from early (RH1) to late shield (RH2)
        pars.r_2 = pars.r_1; % Transition rate from late shield (RH2) to Susceptible (S)


        % Sims - Intervene
        pars.alpha= alpha_vals(i);  % Shielding interventions, 1) for alpha = 0, 2) for alpha = 2, 3) for alpha = 20
        [t,y]=ode45(@stir_model_youngshields_exposedcases, [0:1:outbreak.pTime],ye_mat(n,:),opts,pars,agepars);

        % Stats
        stats(n).R=y(:,agepars.R_ids);
        stats(n).D=y(:,agepars.D_ids);
        stats(n).Htot=y(:,agepars.Ihsub_ids)+y(:,agepars.Ihcri_ids);
        stats(n).Hacu=y(:,agepars.Ihcri_ids);
        stats(n).Dday_age= stats(n).D(2:end,:)-stats(n).D(1:end-1,:);
        stats(n).Dday=sum(stats(n).Dday_age');
        stats(n).Hacu_day=sum(stats(n).Hacu');
        stats(n).lock = y(:,agepars.Slock_ids);
        stats(n).shield_early = y(:,agepars.RH1_ids);
        stats(n).shield_late = y(:,agepars.RH2_ids);

        stats(n).S = y(:,agepars.S_ids);
        stats(n).I = y(:,agepars.Ia_ids) + y(:,agepars.Is_ids);
        stats(n).E = y(:,agepars.cases_ids);
    end
    meta_stats(i).stats = stats;
end
toc

% Save results for different alpha values
stats_alpha0 = meta_stats(1).stats;
stats_alpha2 = meta_stats(2).stats;
stats_alpha20 = meta_stats(3).stats;

save stats_alpha0.mat stats_alpha0
save stats_alpha2.mat stats_alpha2
save stats_alpha20.mat stats_alpha20


%% Plot the results (high R0 scenario)
% cumulative deaths, Max ICU beds and total cases
load stats_alpha2.mat
load stats_alpha20.mat
load stats_alpha0.mat
mean_immuneduration = 60:7:365; %days

deaths_alpha2 = zeros(size(stats_alpha2));
deaths_alpha20 = zeros(size(stats_alpha20));
deaths_alpha0 = zeros(size(stats_alpha0));

maxicu_alpha2 = zeros(size(stats_alpha2));
maxicu_alpha20 = zeros(size(stats_alpha20));
maxicu_alpha0 = zeros(size(stats_alpha0));

cases_alpha0 = zeros(size(stats_alpha0));
cases_alpha2 = zeros(size(stats_alpha2));
cases_alpha20 = zeros(size(stats_alpha20));


for n = 1:size(stats_alpha0,2)
    % cumulative deaths
    deaths_alpha2(n) = sum(stats_alpha2(n).Dday)*1e5;
    deaths_alpha20(n) = sum(stats_alpha20(n).Dday)*1e5;
    deaths_alpha0(n) = sum(stats_alpha0(n).Dday)*1e5;
    % max ICU
    maxicu_alpha2(n) = max(stats_alpha2(n).Hacu_day)*1e5;
    maxicu_alpha20(n) = max(stats_alpha20(n).Hacu_day)*1e5;
    maxicu_alpha0(n) = max(stats_alpha0(n).Hacu_day)*1e5;
    % total cases, alpha 2
    E = stats_alpha2(n).E;
    Eday = E';
    Eday = sum(Eday);
    Efinal = Eday(end);
 
    cases_alpha2(n) = Efinal;
    
    % total cases, alpha 20
    E = stats_alpha20(n).E;
    Eday = E';
    Eday = sum(Eday);
    Efinal = Eday(end);
 
    cases_alpha20(n) = Efinal;
    
     % total cases, alpha 0
    E = stats_alpha0(n).E;
    Eday = E';
    Eday = sum(Eday);
    Efinal = Eday(end);
 
    cases_alpha0(n) = Efinal;
end

% cumulative deaths per 100,000 as a function of mean immunity duration
subplot(3,1,1)
plot(mean_immuneduration/30, deaths_alpha0, 'k-', 'linewidth',3)
hold on
plot(mean_immuneduration/30, deaths_alpha2, 'k--', 'linewidth',3)
plot(mean_immuneduration/30, deaths_alpha20, 'k:', 'linewidth',3)
hold off
set(gca,'fontsize',16);
gc = get(gca);
gc.XAxis.LineWidth = 2;
gc.YAxis.LineWidth = 2;
xlabel('Average immunity duration, months','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Cumulative deaths per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
title({'COVID-19 Epidemic - High Scenario - Shields Ages 20-60';'${\cal{R}}_0=2.33$'},'fontsize',18,'interpreter','latex')
xlim([2 mean_immuneduration(end)/30]);
ylim([0 max(deaths_alpha0)+10])
legend('$\alpha = 0$', '$\alpha = 2$','$\alpha = 20$', 'interpreter','latex', 'Location', 'east');
legend('boxoff');

% max ICU beds per 100,000 as a function of average immunity duration
subplot(3,1,2)
plot(mean_immuneduration/30, maxicu_alpha0, 'k-', 'linewidth',3)
hold on
plot(mean_immuneduration/30, maxicu_alpha2, 'k--', 'linewidth',3)
plot(mean_immuneduration/30, maxicu_alpha20, 'k:', 'linewidth',3)
hold off
set(gca,'fontsize',16);
gc = get(gca);
gc.XAxis.LineWidth = 2;
gc.YAxis.LineWidth = 2;
xlabel('Average immunity duration, months','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Max ICU beds per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
xlim([2 mean_immuneduration(end)/30]);
ylim([0 max(maxicu_alpha0)+30])
legend('$\alpha = 0$', '$\alpha = 2$','$\alpha = 20$', 'interpreter','latex', 'Location', 'east');
legend('boxoff');

% Total number of cases from the integration of the exposed class mock eqn.
% 'stats(n).E'
% Formally, the size of outbreak or total number of cases will be 1-S(t_infinity)

subplot(3,1,3)
% max ICU beds per 100,000 as a function of average immunity duration
plot(mean_immuneduration/30, cases_alpha0*1e5, 'k-', 'linewidth',3)
hold on
plot(mean_immuneduration/30, cases_alpha2*1e5, 'k--', 'linewidth',3)
plot(mean_immuneduration/30, cases_alpha20*1e5, 'k:', 'linewidth',3)
hold off
set(gca,'fontsize',16);
gc = get(gca);
gc.XAxis.LineWidth = 2;
gc.YAxis.LineWidth = 2;
xlabel('Average immunity duration, months','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Number of cases per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
xlim([2 mean_immuneduration(end)/30]);
ylim([15000 max(cases_alpha0*1e5)+1000]);
lg = legend('$\alpha = 0$', '$\alpha = 2$','$\alpha = 20$', 'interpreter','latex', 'position', [0.7498    0.1518    0.1533    0.0728]);
legend('boxoff');
set(gcf, 'position', [336    41   547   781]);
saveas(gcf, './Fig_waningimmunity_high', 'fig')
saveas(gcf, './Fig_waningimmunity_high','epsc')