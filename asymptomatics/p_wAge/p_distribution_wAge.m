%% Impact of asymptomatic fraction (p) -- p as function of age -- plot the distributions
%% David Demory -- March 24 2020

%% Set-up figures
clf
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'fig_p_distribution_wAge';
tmpfilebwname = sprintf('%s_v0',tmpfilename);
tmpfilenoname = sprintf('%s_v0',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);
set(gcf,'position', [514 145 450 400]);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto','color','white');

%% Distributions
% p distribution as function of age with average p = 50, 75 and 90.
p90=[0.99 0.99 0.95 0.93 0.88 0.85 0.825 0.8 0.75 0.72];
p75=[0.88 0.88 0.85 0.80 0.75 0.7 0.65 0.55 0.412 0.30];
p50=[0.8 0.79 0.65 0.55 0.45 0.35 0.24 0.15 0.10 0.1];

% age groups
agepars.meanage=5:10:95;
agepars.highage=[9:10:99];  % Age groups
agepars.lowage=[0:10:90];  % Age groups
population.N=10*10^6;
population.agefrac = [0.12 0.13 0.13 0.13 0.13 0.13 0.11 0.06 0.04 0.02];
population.meanage = sum(agepars.meanage.*population.agefrac);

% average p (\bar{p})
barp50 = sum(p50.*population.agefrac)
barp75 = sum(p75.*population.agefrac)
barp90 = sum(p90.*population.agefrac)

% beta_a (R0 = 1.85)
%beta_a_50 = (1.35*pars.gamma_a)./p50;
%beta_s_50 = (0.5*pars.gamma_s)./(1-p50);
%beta_a_75 = (1.35*pars.gamma_a)./p75;
%beta_s_75 = (0.5*pars.gamma_s)./(1-p75);
%beta_a_90 = (1.35*pars.gamma_a)./p90;
%beta_s_90 = (0.5*pars.gamma_s)./(1-p90);

%% plot


hold on
pplot = bar(agepars.meanage,[p50;p75;p90],'EdgeColor','none','BarWidth',1,'FaceColor','flat')
pplot(1).CData = repmat([0.75 0.75 0.75],10,1); 
pplot(2).CData = repmat([0.5 0.5 0.5],10,1);
pplot(3).CData = repmat([0.1 0.1 0.1],10,1);
xlabel('Age','Interpreter','Latex')
ylabel('Asymptomatic fraction ($p$)','Interpreter','Latex')
legend(pplot,{'$\bar{p}= 0.5$','$\bar{p}= 0.75$','$\bar{p}= 0.9$'},'Orientation','horizontal','location','NorthOutside','Interpreter','latex')
legend('boxoff')
set(gca,'Fontsize',16)
%subplot(1,3,2)
%hold on
%bar(agepars.meanage,1./[beta_a_50;beta_a_75;beta_a_90],'EdgeColor','none','BarWidth',1)
%xlabel('Age')
%ylabel('$\beta_a$ (d)','Interpreter','latex')
%subplot(1,3,3)
%hold on
%bar(agepars.meanage,1./[beta_s_50;beta_s_75;beta_s_90],'EdgeColor','none','BarWidth',1)
%xlabel('Age')
%ylabel('$\beta_s$ (d)','Interpreter','latex')

%legend({'$\bar{p}= 0.5$','$\bar{p}= 0.75$','$\bar{p}= 0.9$'},'Orientation','horizontal','location','NorthOutside','Interpreter','latex')
%legend('boxoff')

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

print('Fig_p_distribution_wAge.eps','-depsc')

% automatic creation of postscript
psprintc(tmpfilename);

% set following on if zooming of
% plots is required
% may need to get legend up as well
%axes(tmpa1)
%axes(tmplh)
clear tmp*
