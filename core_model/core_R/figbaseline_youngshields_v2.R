rm(list=ls())

# (0) About ---------------------------------------------------------------

# 03.30.2020
# Conan Zhao
# czhao98 AT gatech DOT edu
#
# This is a translation of figbaseline_youngshields_v2.m and generates:
# figbaseline_youngshields_v2.pdf

# Model Parms: use ODE45
#            : reltol=1e-8
#            : maxstep=0.1



# (1) Libraries -----------------------------------------------------------

require(deSolve)
require(ggplot2)
require(reshape2)
require(cowplot)

# (2) Model & Util Functions ----------------------------------------------

stir_model_youngshields_v2 = function(t, y, combpars) {
  # usage: dydt = stir_model_youngshields(t,y,pars,agepars)
  # S-E-Ia-Is-R model with I_s further broken down to I_ha I_hs D
  # There are 8 categories in total, all age-stratified
  with(as.list(c(y, combpars)), {
    
    # Assign things
    dydt     = rep(0, length(y));             # concatenated vector of all variables
    Ia       = sum(y[agepars$Ia_ids]);        # Infected, Asymptomatic
    Is       = sum(y[agepars$Is_ids]);        # Infected, Symptomatic
    R        = sum(y[agepars$R_ids]);         # Recovered
    S        = sum(y[agepars$S_ids]);         # Symptomatic
    E        = sum(y[agepars$E_ids]);         # Exposed
    Rshields = sum(y[agepars$R_ids[3:6]])     # Recovered Shields
    Ntot     = S+E+Ia+Is+R;                   # Total # of individuals
    
    dydt[agepars$S_ids]=-pars$beta_asym*y[agepars$S_ids]*Ia/(Ntot+pars$alpha*Rshields)-pars$beta_sym*y[agepars$S_ids]*Is/(Ntot+pars$alpha*Rshields)+agepars$ageleave*y[agepars$Slock_ids];
    dydt[agepars$E_ids]=pars$beta_asym*y[agepars$S_ids]*Ia/(Ntot+pars$alpha*Rshields)+pars$beta_sym*y[agepars$S_ids]*Is/(Ntot+pars$alpha*Rshields)-pars$gamma_e*y[agepars$E_ids];
    dydt[agepars$Ia_ids]=pars$p*pars$gamma_e*y[agepars$E_ids]-pars$gamma_asym*y[agepars$Ia_ids];
    dydt[agepars$Is_ids]=(rep(1,length(pars$p))-pars$p)*pars$gamma_e*y[agepars$E_ids]-pars$gamma_sym*y[agepars$Is_ids];
    dydt[agepars$Ihsub_ids]=agepars$hosp_frac*(1-agepars$hosp_crit)*pars$gamma_sym*y[agepars$Is_ids]-pars$gamma_h*y[agepars$Ihsub_ids];
    dydt[agepars$Ihcri_ids]=agepars$hosp_frac*agepars$hosp_crit*pars$gamma_sym*y[agepars$Is_ids]-pars$gamma_h*y[agepars$Ihcri_ids];
    dydt[agepars$R_ids]=pars$gamma_asym*y[agepars$Ia_ids]+pars$gamma_sym*y[agepars$Is_ids]*(1-agepars$hosp_frac)+pars$gamma_h*y[agepars$Ihsub_ids]+pars$gamma_h*y[agepars$Ihcri_ids]*(1-agepars$crit_die);
    dydt[agepars$D_ids]=pars$gamma_h*y[agepars$Ihcri_ids]*agepars$crit_die;
    dydt[agepars$Slock_ids]=-agepars$ageleave*y[agepars$Slock_ids];
    
    return(list(dydt))
  })
}

run_core_model = function(outbreak_in = outbreak
                          , pars_in = pars
                          , agepars_in = agepars
                          , population_in = population){
  # Runs figbaseline_youngshields_v2 core model
  # Returns (list): 
  #     y     model out
  #     t     time vector
  #     stats summary model statistics
  
  # Set up Baseline
  pars_baseline = pars_in
  pars_baseline$alpha = 0
  
  # Run Baseline Model
  times_baseline=0:outbreak_in$pTime
  model_out = ode(outbreak_in$y0, times_baseline, stir_model_youngshields_v2, c(pars_baseline, agepars_in), rtol=1e-8, hmax=0.1, method='ode45')
  y_baseline = model_out[,-1]
  t_baseline = model_out[,1]
  
  # Back-solve for when <shield_threshold> individuals exposed and begin shielding
  t_shield = min(which((1-rowSums(y_baseline[,agepars_in$S_ids]))*population_in$N > outbreak_in$shield_threshold))-1 # -1 due to 0-indexing

  # Use new starting point for intervention models
  outbreak_shields = outbreak_in
  outbreak_shields$y0=y_baseline[t_shield+1,] # new initial condition
  
  times = 0:outbreak$pTime
  
  # Run Shielding Model
  model_out = ode(outbreak_shields$y0, times, stir_model_youngshields_v2, c(pars_in, agepars_in), rtol=1e-8, hmax=0.1, method='ode45')
  y = model_out[,-1]
  t = model_out[,1]
  
  
  # Stats
  stats=list()
  stats$R=y[,agepars_in$R_ids];
  stats$D=y[,agepars_in$D_ids];
  stats$Htot=y[,agepars_in$Ihsub_ids]+y[,agepars_in$Ihcri_ids];
  stats$Hacu=y[,agepars_in$Ihcri_ids];
  stats$Dday_age=stats$D[2:(nrow(stats$D)),]-stats$D[1:(nrow(stats$D)-1),];
  stats$Dday=rowSums(stats$Dday_age);
  stats$Hacu_day=rowSums(stats$Hacu);
  stats$lock=y[,agepars_in$Slock_ids];
  
  return(list('y' = y, 't' = t, 'stats' = stats))
}

# (3) Inputs --------------------------------------------------------------

# Load Defaults:
#     pars          Epi Parameters
#     population    Population Age Fraction
#     agepars       Age-based Epi Parameters
#     outbreak      Initial

  #source('2020-03-30_input-parameters_high.R')
  source('2020-03-30_input-parameters_low.R')

# (4) Core Model ----------------------------------------------------------

# Sims - baseline
pars$alpha=0
model_res = run_core_model()
yb = model_res$y
tb = model_res$t
statsb = model_res$stats

# Sims - intervene low
pars$alpha=2
model_res = run_core_model()
y = model_res$y
t = model_res$t
stats = model_res$stats

# Sims - intervene high
pars$alpha=20
model_res = run_core_model()
yh = model_res$y
th = model_res$t
statsh = model_res$stats


# (5) Plotting ------------------------------------------------------------

theme_set(theme_gray(base_size = 12))

# DF for deaths per day
df_Dday = data.frame('baseShield' = statsb$Dday * 100000
                     , 'lowShield' = stats$Dday * 100000
                     , 'highShield' = statsh$Dday * 100000
                     , 't' = t[-1])
colnames(df_Dday) = c('Baseline', '2:1 Shielding', '20:1 Shielding', 't')
melt_Dday = melt(df_Dday, 't', variable.name = 'alpha', value.name = 'deaths')
p_Dday = ggplot(melt_Dday, aes(x=t, y=deaths)) + 
  geom_line(size = 1.5, color='black', aes(linetype=alpha)) +
  xlab('Time, days') + 
  ylab('Deaths per 100,000') + 
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

# DF for icu beds per day
df_Hacu_day = data.frame('baseShield' = statsb$Hacu_day * 100000
                         , 'lowShield' = stats$Hacu_day * 100000
                         , 'highShield' = statsh$Hacu_day * 100000
                         , 't' = t)
colnames(df_Hacu_day) = c('Baseline', '2:1 Shielding', '20:1 Shielding', 't')
melt_Hacu_day = melt(df_Hacu_day, 't', variable.name = 'alpha', value.name = 'icu_beds')
p_Hacu_day = ggplot(melt_Hacu_day, aes(x=t, y=icu_beds)) + 
  geom_line(size = 1.5, color='blue', aes(linetype=alpha)) +
  xlab('Time, days') + 
  ylab('ICU beds per 100,000') + 
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted'))

# DF for deaths by age
df_D_byAge = data.frame('baseShield' = statsb$D[nrow(stats$D),]*100000
                        , 'lowShield' = stats$D[nrow(stats$D),]*100000
                        , 'highShield' = statsh$D[nrow(stats$D),]*100000
                        , 'age' = agepars$meanage
                        , 'agefrac' = population$agefrac)
colnames(df_D_byAge) = c('Baseline', '2:1 Shielding', '20:1 Shielding', 'age', 'Population Structure')
melt_D_byAge = melt(df_D_byAge, c('age', 'Population Structure'), variable.name = 'alpha', value.name = 'deaths')

coeff = 6*max(melt_D_byAge$deaths) # For plotting age structure on top of deaths by age
p_D_byage = ggplot(melt_D_byAge, aes(x=age)) + 
  
  geom_line(size = 1.5, aes(y=coeff*`Population Structure`, linetype=alpha, colour =  'Age Fraction')) +
  geom_point(size = 4, shape=5, aes(y=coeff*`Population Structure`, color = 'Age Fraction')) +
  
  geom_line(size = 1.5, aes(y=deaths, linetype=alpha, color = 'Model')) +
  geom_point(size = 4, shape=1, aes(y=deaths, color = 'Model')) +
  scale_y_continuous(name = 'Cumulative Deaths per 100,000'
                     , sec.axis = sec_axis(~.*(1/coeff), name="Population Age Structure")) +
  theme(axis.title.y.right = element_text(color='gray50')) + 
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted')) + 
  scale_color_manual(values = c('grey50', 'black')) +
  labs(color = NULL)


# Combine to single multiplot figure
title = ggdraw() + draw_label(paste("COVID-19 Epidemic - ", Ro_lowhigh, " Scenario - Shields Ages 20-60", sep='', collapse='')
                              , fontface='bold', size=14)
subtitle = ggdraw() + draw_label(paste('Asymptomatic incidence p = ', pars$overall_p, ', R_o = ', pars$R0), fontface='bold', size=14)

p_res = plot_grid(p_Dday, p_Hacu_day, p_D_byage, ncol=1, align='hv')
p_res_titled = plot_grid(title,subtitle,p_res, rel_heights = c(0.1,0.1,3), ncol=1)

# Save Figure
FULL_FIGNAME = paste(substr(Sys.time(), 1, 10), "_", FIGNAME, sep='', collapse='')
ggsave(FULL_FIGNAME, p_res_titled, 'pdf', width=10, height=8)







