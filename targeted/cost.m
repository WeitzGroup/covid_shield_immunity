function cst_f = cost(outbreak, pars, agepars, Rshields_deploy)
pars.Rshields_deploy = Rshields_deploy; % set deployment

% Sims - Get to Crossing
opts=odeset('reltol',1e-8,'maxstep',0.1,'events',@intervene_trigger);
pars.alpha = 0;
[t,y,te,ye,ie]=ode45(@stir_model_youngshields_target,[0:1:outbreak.pTime], outbreak.y0,opts,pars,agepars);

% Sims - Intervene w/Shields
pars.alpha = pars.alpha1; % set the level of interest
opts=odeset('reltol',1e-8,'maxstep',0.1);
[t,y]=ode45(@stir_model_youngshields_target,[0:1:outbreak.pTime], ye,opts,pars,agepars);

% ICU over time
ICU_sol = sum(y(:, agepars.Ihcri_ids), 2);

if sum(ICU_sol.*100000 >= pars.B)
    cst_f = 1e7; % arbitrarily  large
else
    % total Deaths
    D_tf = sum(y(end, agepars.D_ids))*100000;
    % return cost
    cst_f = pars.Wi*sum(log10(1./(pars.B - 100000.*ICU_sol))) + pars.Wd*D_tf;
end

end