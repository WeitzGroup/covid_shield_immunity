clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figbaseline_shields_high';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);
set(gcf,'position', [514 145 687 799]);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');

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
pars.Itrigger = 10000/population.N; % Trigger at 10000 total cases, irrespective of type

% Parameters
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

% Epi parameters
pars.Ra=pars.beta_a/pars.gamma_a;
pars.Rs=pars.beta_s/pars.gamma_s;
pars.R0=pars.p*pars.Ra+(1-pars.p)*pars.Rs;
pars.R0=sum(pars.p.*population.agefrac*pars.Ra+(1-pars.p).*population.agefrac*pars.Rs);

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
pars.alpha=0;  % Shielding
[t,y,te,ye,ie]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);

% Sims - Intervene w/Shields
opts=odeset('reltol',1e-8,'maxstep',0.1);
pars.alpha=2;  % Shielding
[t,y]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime], ye,opts,pars,agepars);

% Stats
stats.R=y(:,agepars.R_ids);
stats.D=y(:,agepars.D_ids);
stats.Htot=y(:,agepars.Ihsub_ids)+y(:,agepars.Ihcri_ids);
stats.Hacu=y(:,agepars.Ihcri_ids);
stats.Dday_age=stats.D(2:end,:)-stats.D(1:end-1,:);
stats.Dday=sum(stats.Dday_age');
stats.Hacu_day=sum(stats.Hacu');
stats.lock=y(:,agepars.Slock_ids);

% Sims - Intervene High
pars.alpha=20;  % Shielding
[th,yh]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime], ye,opts,pars,agepars);

% Stats
statsh.R=yh(:,agepars.R_ids);
statsh.D=yh(:,agepars.D_ids);
statsh.Htot=yh(:,agepars.Ihsub_ids)+yh(:,agepars.Ihcri_ids);
statsh.Hacu=yh(:,agepars.Ihcri_ids);
statsh.Dday_age=statsh.D(2:end,:)-statsh.D(1:end-1,:);
statsh.Dday=sum(statsh.Dday_age');
statsh.Hacu_day=sum(statsh.Hacu');
statsh.lock=yh(:,agepars.Slock_ids);

% Sims - Baseline
% Irrelevant, no shielding
pars.alpha=0;
[tb,yb]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime], ye,opts,pars,agepars);
% Stats
statsb.R=yb(:,agepars.R_ids);
statsb.D=yb(:,agepars.D_ids);
statsb.Htot=yb(:,agepars.Ihsub_ids)+yb(:,agepars.Ihcri_ids);
statsb.Hacu=yb(:,agepars.Ihcri_ids);
statsb.Dday_age=statsb.D(2:end,:)-statsb.D(1:end-1,:);
statsb.Dday=sum(statsb.Dday_age');
statsb.Hacu_day=sum(statsb.Hacu');
statsb.lock=yb(:,agepars.Slock_ids);

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
xlabel('Time, days','fontsize',16,'verticalalignment','top','interpreter','latex');
ylabel('Deaths per day per 100,000','fontsize',16,'verticalalignment','bottom','interpreter','latex');
title({'Shield Immunity, High Scenario, ${\cal{R}}_0=2.33$'},'fontsize',18,'interpreter','latex')
ylim([0 15]);
legend('Baseline','2:1 Shielding','20:1 Shielding');
legend('boxoff');
xlim([0 365]);
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
ylim([0 300]);
set(gca,'ytick',[0:50:300]);
xlim([0 365]);
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
ylim([0 300]);
set(gca,'ytick',[0:50:300]);
tmpl=legend('Baseline','2:1 Shielding','20:1 Shielding','Age Structure');
set(tmpl,'location','SouthWest');
legend('boxoff');

axes(tmpax(2));
ylabel('Population age structure','fontsize',16,'verticalalignment','top','interpreter','latex','color','k');
set(gca,'ycolor',[0.25 0.25 0.25]);
set(gca,'fontsize',16);
ylim([0 0.2]);
% title('','fontsize',24)
%
%
% Some helpful plot commands
% tmph=plot(x,y,'ko');
% set(tmph,'markersize',10,'markerfacecolor,'k');
% tmph=plot(x,y,'k-');
% set(tmph,'linewidth',2);

% for use with layered plots
% set(gca,'box','off')

% adjust limits
% tmpv = axis;
% axis([]);
% ylim([]);

% change axis line width (default is 0.5)
% set(tmpa1,'linewidth',2)

% fix up tickmarks
% set(gca,'xtick',[1 100 10^4])
% set(gca,'ytick',[1 100 10^4])

% creation of postscript for papers
% psprint(tmpxfigfilename);

% the following will usually not be printed 
% in good copy for papers
% (except for legend without labels)

% legend
% tmplh = legend('stuff',...);
% tmplh = legend('','','');
% remove box
% set(tmplh,'visible','off')
% legend('boxoff');

% title('','fontsize',24)
% 'horizontalalignment','left');

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

tmpt = pwd;
tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
datenamer(1.1,.05,90);
% datename(.5,.05);
% datename2(.5,.05); % 2 rows

% automatic creation of postscript
psprintc(tmpfilename);

% set following on if zooming of 
% plots is required
% may need to get legend up as well
%axes(tmpa1)
%axes(tmplh)
clear tmp*
