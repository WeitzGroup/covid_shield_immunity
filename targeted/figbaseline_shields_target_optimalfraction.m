agepars.meanage=5:10:95;
population.agefrac = [0.12 0.14 0.14 0.13 0.13 0.13 0.10 0.06 0.04 0.01]; 

figure(1);
subplot(1,2,1);
bar(agepars.meanage,(load('ga_shields_delpoyment_alpha2_highCase.txt'))'./population.agefrac)
ylabel('Shielding Concentration, $\theta_{a}/f_{a}$','fontsize',16,'verticalalignment','bottom','interpreter','latex','color','k');
xlabel('Mean Age','fontsize',16,'verticalalignment','top','interpreter','latex');
set(gca,'fontsize',16);
axis square;
title({'Shield Immunity, High Scenario, ${\cal{R}}_0=2.33$'},'fontsize',18,'interpreter','latex');

subplot(1,2,2);
agepars.meanage=5:10:95;
bar(agepars.meanage,(load('ga_shields_delpoyment_alpha2_lowCase.txt'))'./population.agefrac)
ylabel('Shielding Concentration, $\theta_{a}/f_{a}$','fontsize',16,'verticalalignment','bottom','interpreter','latex','color','k');
xlabel('Mean Age','fontsize',16,'verticalalignment','top','interpreter','latex');
set(gca,'fontsize',16);
axis square;
title({'Shield Immunity, Low Scenario, ${\cal{R}}_0=1.57$'},'fontsize',18,'interpreter','latex')

% output figures in eps 
ff = figure(1);
ff.Units = 'inches';
Width = 30; Height = 15;

ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];

Name = 'fig_baseline_shields_target_optimalfraction';
print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');










