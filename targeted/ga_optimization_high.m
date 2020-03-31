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

% Sheilds tagert delolyment (optimized given alpha = 2)
agepars.agefrac = [0.12 0.14 0.14 0.13 0.13 0.13 0.10 0.06 0.04 0.01]; % need this pars pass to model
% pars.Rshields_deploy = agepars.agefrac; % for sanity check of code

%% Parameters that may influence the optimization result
% Shielding preference
pars.alpha1 = 2;  
% pars.alpha1 = 20;  

% optimization parameters
pars.B = 200; % max capacity of ICU beds

% control parameters, weights
pars.Wi = 1e-7; pars.Wd = 1; 

% high/low case - whole system parameters


%% --- GA-algorithm ---
% probability simplex feasible domain
% equality constrain

Aeq = ones(1, length(agepars.agefrac));
beq = 1; % Aeq*x = beq => sum of x is 1
% inequality constrain - positive
A = eye(length(agepars.agefrac));
b = zeros(length(agepars.agefrac), 1); 

% define cost/fitness in GA as function of x
fun = @(x) cost(outbreak, pars, agepars, x);


% set max generation
max_iter = 30;
options = optimoptions('ga','MaxGenerations',max_iter,'PlotFcn', @gaplotbestf);
rng default % For reproducibility
x_f = ga(fun,length(agepars.agefrac),-A,b,Aeq,beq,[],[],[],options);

% write file
fid= fopen('ga_shields_delpoyment_alpha2_highCase.txt','w');
%fid= fopen('ga_shields_delpoyment_alpha20_highCase.txt','w');
fprintf(fid,'%.5f\n',x_f);
fclose(fid);

