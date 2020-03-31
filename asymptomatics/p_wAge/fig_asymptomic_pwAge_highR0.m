%% Impact of asymptomatic fraction (p) -- p is a funciton of age.
%% David Demory -- March 25 2020

%% set-up figure
clf
figA = figure;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'fig_asymptomatic_pWage_highR0cst';
tmpfilebwname = sprintf('%s_v2',tmpfilename);
tmpfilenoname = sprintf('%s_v2',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);
set(gcf,'position', [514 145 800 699]);

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
pars.gamma_a=1/6;   % Resolution rate for asymptomatic 
pars.gamma_s=1/6;  % Resolution rate for symptomatic
pars.gamma_h=1/10;  % Resolution rate in hospitals
pars.beta_a=3/10;   % Transmission for asymptomatic
pars.beta_s=6/10;      % Transmission for symptomatic
%pars.p=0.6955;         % Fraction asymptomatic
% pars.p=[0.99 0.99 0.95 0.9 0.8 0.7 0.6 0.5 0.5 0.5];         % Fraction asymptomatic
pars.p=[0.95 0.95 0.9 0.8 0.7 0.6 0.4 0.2 0.2 0.2];    % Fraction asymptomatic
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

% p distribution as function of age
p90=[0.99 0.99 0.95 0.93 0.88 0.85 0.825 0.8 0.75 0.72];
p75=[0.88 0.88 0.85 0.80 0.75 0.7 0.65 0.55 0.412 0.30];
p50=[0.8 0.79 0.65 0.55 0.45 0.35 0.24 0.15 0.10 0.1];


%% Dynamics
% Sims - Baseline
% Irrelevant, no shielding
A = [0,2,20]; %alpha (shielding)
P = [p50;p75;p90]; % p distribution (average: 0.5, 0.75, 0.9 -- cf. p_distribution_wAge.m)
opts=odeset('reltol',1e-8,'maxstep',0.1);

for j = 1:3; % loop on p
    for i = 1:length(A); % loop on alpha
        
        % p and beta-s initiation
        pars.p = P(j,:);
        barp = sum(pars.p.*population.agefrac);
        %pars.beta_a = (1.35*pars.gamma_a)./pars.p;
        %pars.beta_s = (0.5*pars.gamma_s)./(1-pars.p);
        pars.beta_a = (cstR0a*pars.gamma_a)./sum(pars.p.*population.agefrac);
        pars.beta_s = (cstR0s*pars.gamma_s)./(1-sum(pars.p.*population.agefrac));
        pars.Ra=pars.beta_a./pars.gamma_a;
        pars.Rs=pars.beta_s./pars.gamma_s;
        pars.R0_age=pars.p.*pars.Ra+(1-pars.p).*pars.Rs;
        %pars.R0_age
        pars.R0=sum(pars.p.*population.agefrac.*pars.Ra+(1-pars.p).*population.agefrac.*pars.Rs);
        
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
        
        % Sims - Get to Crossing
        opts=odeset('reltol',1e-8,'maxstep',0.1,'events',@intervene_trigger);
        pars.alpha=A(i);
        [t,y,te,ye,ie]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);
        % Sims
        [t,y]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime],ye,opts,pars,agepars);
        
        % Stats
        stats.p{j}.alpha{i}.R=y(:,agepars.R_ids);
        stats.p{j}.alpha{i}.D=y(:,agepars.D_ids);
        stats.p{j}.alpha{i}.I = y(:,agepars.Ia_ids)+y(:,agepars.Is_ids);
        stats.p{j}.alpha{i}.Htot = y(:,agepars.Ihsub_ids)+y(:,agepars.Ihcri_ids);
        
        disp(['\bar{p} = ',num2str(sum(pars.p.*population.agefrac)),' -- \alpha = ',num2str(pars.alpha),' -- R0 = ',num2str(pars.R0)])
        %Verif fraction = 1 ode
        sum(sum(y,2))/length(t)
    end
end




%% Plot

subplot(3,1,1);hold on
tmph50=plot(t(2:end),diff(sum(stats.p{1}.alpha{1}.D,2))*100000,'-','Color',[0.1 0.1 0.1]);
tmph75=plot(t(2:end),diff(sum(stats.p{2}.alpha{1}.D,2))*100000,'-','Color',[0.5 0.5 0.5]);
tmph90=plot(t(2:end),diff(sum(stats.p{3}.alpha{1}.D,2))*100000,'-','Color',[0.75 0.75 0.75]);
set(tmph50,'linewidth',3);set(tmph75,'linewidth',3);set(tmph90,'linewidth',3);

tmph50=plot(t(2:end),diff(sum(stats.p{1}.alpha{2}.D,2))*100000,'--','Color',[0.1 0.1 0.1]);
tmph75=plot(t(2:end),diff(sum(stats.p{2}.alpha{2}.D,2))*100000,'--','Color',[0.5 0.5 0.5]);
tmph90=plot(t(2:end),diff(sum(stats.p{3}.alpha{2}.D,2))*100000,'--','Color',[0.75 0.75 0.75]);
set(tmph50,'linewidth',3);set(tmph75,'linewidth',3);set(tmph90,'linewidth',3);

tmph50=plot(t(2:end),diff(sum(stats.p{1}.alpha{3}.D,2))*100000,':','Color',[0.1 0.1 0.1]);
tmph75=plot(t(2:end),diff(sum(stats.p{2}.alpha{3}.D,2))*100000,':','Color',[0.5 0.5 0.5]);
tmph90=plot(t(2:end),diff(sum(stats.p{3}.alpha{3}.D,2))*100000,':','Color',[0.75 0.75 0.75]);
set(tmph50,'linewidth',3);set(tmph75,'linewidth',3);set(tmph90,'linewidth',3);

legB =plot(0,0,'k-','LineWidth',3)
leg =plot(0,0,'k--','LineWidth',3)
legH =plot(0,0,'k:','LineWidth',3)

set(gca,'fontsize',16);
xlabel('Time, days','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Deaths per day per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
title({'Asymptomatic fraction $p$ as a function of age';' High scenario ${\cal{R}}_0=2.33$'},'fontsize',18,'interpreter','latex')
ylim([0 20]);
xlim([0 365])
legend([legB,leg,legH],{'Baseline','2:1 Shielding','20:1 Shielding'});
legend('boxoff');
text(100,18,'\it{p} = 0.5','FontSize',16,'Color',[0.1 0.1 0.1])
text(120,10,'\it{p} = 0.75','FontSize',16,'Color',[0.5 0.5 0.5])
text(140,4,'\it{p} = 0.9','FontSize',16,'Color',[0.75 0.75 0.75])


subplot(3,1,2);hold on
tmph50=plot(t,sum(stats.p{1}.alpha{1}.Htot,2)*100000,'-','Color',[7, 47, 95]/255);
tmph75=plot(t,sum(stats.p{2}.alpha{1}.Htot,2)*100000,'-','Color',[56, 149, 211]/255);
tmph90=plot(t,sum(stats.p{3}.alpha{1}.Htot,2)*100000,'-','Color',[88, 204, 237]/255);
set(tmph50,'linewidth',3);set(tmph75,'linewidth',3);set(tmph90,'linewidth',3);

tmph50=plot(t,sum(stats.p{1}.alpha{2}.Htot,2)*100000,'--','Color',[7, 47, 95]/255);
tmph75=plot(t,sum(stats.p{2}.alpha{2}.Htot,2)*100000,'--','Color',[56, 149, 211]/255);
tmph90=plot(t,sum(stats.p{3}.alpha{2}.Htot,2)*100000,'--','Color',[88, 204, 237]/255);
set(tmph50,'linewidth',3);set(tmph75,'linewidth',3);set(tmph90,'linewidth',3);

tmph50=plot(t,sum(stats.p{1}.alpha{3}.Htot,2)*100000,':','Color',[7, 47, 95]/255);
tmph75=plot(t,sum(stats.p{2}.alpha{3}.Htot,2)*100000,':','Color',[56, 149, 211]/255);
tmph90=plot(t,sum(stats.p{3}.alpha{3}.Htot,2)*100000,':','Color',[88, 204, 237]/255);
set(tmph50,'linewidth',3);set(tmph75,'linewidth',3);set(tmph90,'linewidth',3);

legB =plot(0,0,'-','LineWidth',3,'Color','k')
leg =plot(0,0,'--','LineWidth',3,'Color','k')
legH =plot(0,0,':','LineWidth',3,'Color','k')

set(gca,'fontsize',16);
xlabel('Time, days','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('ICU beds per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
% title('','fontsize',24)
ylim([0 1000]);
xlim([0 365]);
%yticks([0 250 500])
legend([legB,leg,legH],{'Baseline','2:1 Shielding','20:1 Shielding'});
legend('boxoff');
text(100,900,'\it{p} = 0.5','FontSize',16,'Color',[7, 47, 95]/255)
text(120,500,'\it{p} = 0.75','FontSize',16,'Color',[56, 149, 211]/255)
text(140,200,'\it{p} = 0.9','FontSize',16,'Color',[88, 204, 237]/255)


subplot(3,1,3);

yyaxis right
tmpage = plot(agepars.meanage,population.agefrac);
set(tmpage,'linewidth',3,'color',[0.5 0.5 0.5],'marker','d');
ylabel('Population age structure','fontsize',16,'verticalalignment','top','interpreter','latex','color','k');
set(gca,'ycolor',[0.25 0.25 0.25]);
ylim([0 0.2]);

yyaxis left
hold on
tmpb50=plot(agepars.meanage,stats.p{1}.alpha{1}.D(end,:)*100000);
tmpb75=plot(agepars.meanage,stats.p{2}.alpha{1}.D(end,:)*100000);
tmpb90=plot(agepars.meanage,stats.p{3}.alpha{1}.D(end,:)*100000);
set(tmpb50,'linewidth',3,'color',[7, 47, 95]/255,'marker','o','linestyle','-','markerfacecolor',[7, 47, 95]/255);
set(tmpb75,'linewidth',3,'color',[56, 149, 211]/255,'marker','o','linestyle','-','markerfacecolor',[56, 149, 211]/255);
set(tmpb90,'linewidth',3,'color',[88, 204, 237]/255,'marker','o','linestyle','-','markerfacecolor',[88, 204, 237]/255);

tmp50=plot(agepars.meanage,stats.p{1}.alpha{2}.D(end,:)*100000);
tmp75=plot(agepars.meanage,stats.p{2}.alpha{2}.D(end,:)*100000);
tmp90=plot(agepars.meanage,stats.p{3}.alpha{2}.D(end,:)*100000);
set(tmp50,'linewidth',3,'color',[7, 47, 95]/255,'marker','o','linestyle','--','markerfacecolor',[7, 47, 95]/255);
set(tmp75,'linewidth',3,'color',[56, 149, 211]/255,'marker','o','linestyle','--','markerfacecolor',[56, 149, 211]/255);
set(tmp90,'linewidth',3,'color',[88, 204, 237]/255,'marker','o','linestyle','--','markerfacecolor',[88, 204, 237]/255);


tmph50=plot(agepars.meanage,stats.p{1}.alpha{3}.D(end,:)*100000);
tmph75=plot(agepars.meanage,stats.p{2}.alpha{3}.D(end,:)*100000);
tmph90=plot(agepars.meanage,stats.p{3}.alpha{3}.D(end,:)*100000);
set(tmph50,'linewidth',3,'color',[7, 47, 95]/255,'marker','o','linestyle',':','markerfacecolor',[7, 47, 95]/255);
set(tmph75,'linewidth',3,'color',[56, 149, 211]/255,'marker','o','linestyle',':','markerfacecolor',[56, 149, 211]/255);
set(tmph90,'linewidth',3,'color',[88, 204, 237]/255,'marker','o','linestyle',':','markerfacecolor',[88, 204, 237]/255);

xlabel('Age range','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Cumulative deaths per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
set(gca,'ycolor','k');
ylim([0 350]);

ll0=plot(0,0,'k-','LineWidth',3);
ll2=plot(0,0,'k--','LineWidth',3);
ll20=plot(0,0,'k:','LineWidth',3);

tmpl = legend([tmpage,ll0,ll2,ll20],{'Age Structure','Baseline','2:1','Shielding','20:1 Shielding'});
set(tmpl,'location','SouthWest');
legend('boxoff');
text(65,300,'\it{p} = 0.5','FontSize',16,'Color',[7, 47, 95]/255)
text(45,160,'\it{p} = 0.75','FontSize',16,'Color',[56, 149, 211]/255)
text(35,75,'\it{p} = 0.9','FontSize',16,'Color',[88, 204, 237]/255)



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

print([num2str(tmpfilename),'.eps'],'-depsc')
% automatic creation of postscript
psprintc(tmpfilename);

% set following on if zooming of 
% plots is required
% may need to get legend up as well
%axes(tmpa1)
%axes(tmplh)
%clear tmp*