agepars.meanage=5:10:95;
population.agefrac = [0.12 0.14 0.14 0.13 0.13 0.13 0.10 0.06 0.04 0.01]; 

high_alpha2 = load('ga_shields_delpoyment_alpha2_highCase.txt')./population.agefrac';
high_alpha20 = load('ga_shields_delpoyment_alpha20_highCase.txt')./population.agefrac';

low_alpha2 = load('ga_shields_delpoyment_alpha2_lowCase.txt')./population.agefrac';
low_alpha20 = load('ga_shields_delpoyment_alpha20_lowCase.txt')./population.agefrac';

high_bpcombined = [high_alpha2, high_alpha20];
low_bpcombined = [low_alpha2, low_alpha20];

figure(1);
l = cell(1,2);
l{1}='\alpha = 2'; l{2}='\alpha = 20';    
subplot(1,2,1);
h1 = bar(agepars.meanage,high_bpcombined,'grouped');
ylabel('Shielding Concentration, $\theta_{a}/f_{a}$','fontsize',16,'verticalalignment','bottom','interpreter','latex','color','k');
xlabel('Mean Age','fontsize',16,'verticalalignment','top','interpreter','latex');
ylim([0 10])
set(gca,'fontsize',16);
axis square;
title({'Shield Immunity, High Scenario, ${\cal{R}}_0=2.33$'},'fontsize',18,'interpreter','latex');
legend(h1,l, 'Location','northwest');
legend('boxoff');

subplot(1,2,2);
h2 = bar(agepars.meanage,low_bpcombined,'grouped');
ylabel('Shielding Concentration, $\theta_{a}/f_{a}$','fontsize',16,'verticalalignment','bottom','interpreter','latex','color','k');
xlabel('Mean Age','fontsize',16,'verticalalignment','top','interpreter','latex');
ylim([0 10])
set(gca,'fontsize',16);
axis square;
title({'Shield Immunity, Low Scenario, ${\cal{R}}_0=1.57$'},'fontsize',18,'interpreter','latex');
legend(h2,l, 'Location','northwest');
legend('boxoff');


% output figures in eps 
ff = figure(1);
ff.Units = 'inches';
Width = 30; Height = 15;

ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];

Name = 'fig_baseline_shields_target_optimalfraction_v1';
print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');