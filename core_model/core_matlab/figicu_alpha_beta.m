clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figicu_alpha_beta';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

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

% Parameters
pars.gamma_e=1/2;   % Transition to infectiousness
pars.gamma_a=1/4;   % Resolution rate for asymptomatic 
pars.gamma_s=1/4;  % Resolution rate for symptomatic
pars.gamma_h=1/10;  % Resolution rate in hospitals
pars.beta_a=3/10;   % Transmission for asymptomatic
pars.beta_s=6/10;      % Transmission for symptomatic
% Could be age structured
% pars.p=0.5;         % Fraction asymptomatic
% pars.p=[0.99 0.99 0.95 0.9 0.8 0.7 0.6 0.5 0.5 0.5];         % Fraction asymptomatic
pars.p=[0.95 0.95 0.90 0.8 0.7 0.6 0.4 0.2 0.2 0.2];         % Fraction asymptomatic
pars.overall_p=sum(pars.p.*population.agefrac);
pars.Itrigger = 10000/population.N; % Trigger at 5000 total cases, irrespective of type

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
% IFR of 0.9% (in line with global)
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

% Baseline case
tpre = t;
ypre = y;
[tb,yb]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime],ye,opts,pars,agepars);
stats.deadbase = sum(yb(end,agepars.D_ids));

% Run the intervention iff we haven't done it yet
fid =fopen('alpha_beta_intervene_low.mat');
if (fid <0)
  % Initiate the intervention
  pars.sdistance=0.1:0.02:0.5;
  pars.alpha_range=0:1:20;
  opts=odeset('reltol',1e-8,'maxstep',0.1);
  more off
  for i=1:length(pars.sdistance),
    for j=1:length(pars.alpha_range),
      i
      pars.beta_a=4/10*(1-pars.sdistance(i));   % Transmission for asymptomatic
      pars.beta_s=8/10*(1-pars.sdistance(i));      % Transmission for symptomatic
      pars.alpha=pars.alpha_range(j);
      [t,y]=ode45(@stir_model_youngshields,[0:1:outbreak.pTime],ye,opts,pars,agepars);
      stats.totdead(i,j)=sum(y(end,agepars.D_ids));
      Hacu=y(:,agepars.Ihcri_ids);
      Hacu_day=sum(Hacu');
      stats.Hacu_max(i,j)=max(Hacu_day);
    end
  end
  stats.dead_reduce=(stats.deadbase-stats.totdead)/stats.deadbase;
  save alpha_beta_intervene_low pars agepars population outbreak stats tb yb tpre ypre
else
  load alpha_beta_intervene_low
end

imagesc(pars.alpha_range,pars.sdistance,stats.Hacu_max*10^5);
axis square
set(gca,'YDir','normal');       
hold on
[c,h]=contour(pars.alpha_range,pars.sdistance,stats.Hacu_max*10^5,[25 25]);
set(h(1),'linewidth',3,'color','red')
clabel(c,h)
[c,h]=contour(pars.alpha_range,pars.sdistance,stats.Hacu_max*10^5,[1 2 5 10 15 20 30 40 60]);
set(h(1),'linewidth',2,'color','white')
clabel(c,h)
ylim([0.1 0.4]);
tmpc=colorbar;
%set(tmpc,'ylim',[0 1]);
%set(tmpc,'ticks',[0:0.2:1]);
%set(tmpc,'ticklabels',{'0%','20%','40%','60%','80%','100%'});


xlabel('Shielding, $\alpha$','fontsize',18,'verticalalignment','top','interpreter','latex');
ylabel('Social distancing, reduction in $\beta$','fontsize',18,'verticalalignment','bottom','interpreter','latex');
title({'Peak ICU beds needed, per 100,000';'Low Scenario, ${\cal{R}}_0=1.55$'},'fontsize',18,'interpreter','latex')

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
% xlim([]);

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
