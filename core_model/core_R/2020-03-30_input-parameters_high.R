if(exists("agepars")) rm(agepars)
if(exists("population")) rm(population)
if(exists("pars")) rm(pars)
if(exists("outbreak")) rm(outbreak)
if(exists("N")) rm(N)
if(exists("tmpvals")) rm(tmpvals)
if(exists("tmpzeros")) rm(tmpzeros)

# 03.25.2020
# Conan Zhao
# czhao98 AT gatech DOT edu
#
# Here we set up default age & population parameters, as well as epi parameters.

# (0) Setup ---------------------------------------------------------------

# Declares  
agepars=list()    # Population Structure Parameters
population=list() # Population Parameters
pars=list()       # Modeling Parameters
outbreak=list()   # Initial Outbreak

# OUTPUT PDF NAME
FIGNAME = 'figbaseline_youngshields_v2_high.pdf'
Ro_lowhigh = 'High'

# (1) Population Parameters -----------------------------------------------

# Population
# We model a population with an age distribution simlar to that of the UK or US
# Ages binned as (0-9, 10-19, 20-29, .... 80-89, 90-99)
agepars$meanage=c(10*(1:10)-5)
agepars$highage=c(10*(1:10)-1)
agepars$lowage=c(10*(1:10)-10)

population$N=10*10^6  # Total Population Size
population$agefrac=c(0.12,0.14,0.14,0.13,0.13,0.13,0.10,0.06,0.04,0.01)
population$meanage=sum(population$agefrac*agepars$meanage)


# (2) Model Epi Parameters ------------------------------------------------

# Pars (Input)

pars$gamma_e=1/4      # Transition to infectiousness
pars$gamma_asym=1/6      # Resolution rate for asymptomatic
pars$gamma_sym=1/6      # Resolution rate for symptomatic
pars$gamma_h=1/10     # Resolution rate in hospitals

pars$beta_asym=3/10      # Transmission for asymptomatic
pars$beta_sym=6/10      # Transmission for symptomatic

pars$p=c(0.95, 0.95, 0.9, 0.8, 0.7, 0.6, 0.4, 0.2, 0.2, 0.2) # Asymptomatic fraction by age group
pars$r=1/7

pars$overall_p=sum(pars$p*population$agefrac)   # Overall asymptomatic fraction

# Epi Pars (Calculated)
pars$R_asym=pars$beta_asym/pars$gamma_asym
pars$R_sym=pars$beta_sym/pars$gamma_sym
pars$R0=sum(pars$p*population$agefrac*pars$R_asym+(1-pars$p)*population$agefrac*pars$R_sym)

tmpval=pars$beta_asym/pars$beta_sym*pars$p/(1-pars$p)*(pars$r+pars$gamma_sym)/(pars$r+pars$gamma_asym)
pars$q=tmpval/(1+tmpval)


# (3) Population Age Structure Parameters ---------------------------------

# Age Stratification
agepars$num_ages=length(agepars$meanage)  # number of age categories
N=agepars$num_ages                        # temporary variable

agepars$hosp_frac=c(0.1, 0.3, 1.2, 3.2, 4.9, 10.2, 16.6, 24.3, 27.3, 27.3)/100  # Fraction of cases requiring hospitalization
agepars$hosp_crit=c(5, 5, 5, 5, 6.3, 12.2, 27.4, 43.2, 70.9, 70.9)/100          # Fraction of critical hospitalized cases
agepars$crit_die=0.5*rep(1,N)                                                   # Fraction of fatal critical cases

# Indices (For Modeling)
# Our output will be a (t x m) matrix, for t time points and m variables.
# The following are the column indices for each state variable
agepars$S_ids=1:N
agepars$E_ids=(N+1):(2*N)
agepars$Ia_ids=(2*N+1):(3*N)
agepars$Is_ids=(3*N+1):(4*N)
agepars$Ihsub_ids=(4*N+1):(5*N)
agepars$Ihcri_ids=(5*N+1):(6*N)
agepars$R_ids=(6*N+1):(7*N)
agepars$D_ids=(7*N+1):(8*N)
agepars$Slock_ids=(8*N+1):(9*N)
agepars$ageleave=rep(1,1,10)


# (4) Initial Conditions --------------------------------------------------

# Init the population-baseline:
#   -Open plus hospitals
#   -SEIaIS (open) and then I_ha I_hs and then R (open) and D (cumulative) and 
#    then S (lockdown)- 9 categories in total, all age-stratified
#   -Here, we ignore the lockdown category
#   -Initiate single exposed individual age 20-29
#   -Integrate until 10,000 people have been exposed.
tmpzeros=rep(0,length(agepars$meanage))
outbreak$y0=c(population$agefrac, tmpzeros, tmpzeros, tmpzeros, tmpzeros, tmpzeros, tmpzeros, tmpzeros, tmpzeros)

# Initiate an outbreak: a single exposed individual age 20-29
outbreak$y0=population$N*outbreak$y0
outbreak$y0[3]=outbreak$y0[3]-1
outbreak$y0[13]=1
outbreak$y0=outbreak$y0/population$N
outbreak$pTime=365

outbreak$shield_threshold = 10000 # Number of people exposed before shielding measures

print('Loaded: pars, agepars, population, N')
rm(tmpval, tmpzeros, N)