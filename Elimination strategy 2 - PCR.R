
#CLEAR WORKSPACE
rm(list=ls())
ptm <- proc.time()

#SET WORKING DIRECTORY

setwd(dir = "/Users/simonmendelsohn/Dropbox/TT Placement/Plantation Model")
getwd()

#LOAD LIBRARIES
library(deSolve)
library(ggplot2)

elimination = NULL

for (PCR_cov in seq(0.24,1,0.02)) {

  #MODEL DURATION (in years)
  duration <- 50 # years
  time.seq <- seq(0, duration*12, by = 1/30.415) # in months, by 1 day intervals
  
  #MODEL PARAMETERS
  
  #basic model parameters                
  nu=(1/(68.2*12))      # birth rate (1/life expectancy)
  mu=(1/(68.2*12))      # death rate (1/life expectancy)
  gamma=1/(11/30.415)   # rate of movement from liver to asexual blood stage = 1/(average duration of liver stage)
  delta=1/(10/30.415)   # rate of gametocyte maturation  = 1/(average duration of gametocyte maturation)
  p=(0.93)              # proportion of patients who are clinical/symptomatic (no immunity)
  q=(0.50)              # proportion of patients who are clinical/symptomatic (partial immunity)
  tau=0.4               #** Proportion of clinical (symptomatic) patients treated. 50% seek treatment * 90% sensitivity of RDT * 90% complete treatment
  zeta=1/(120/30.415)   # rate of natural recovery from malaria infection = (1/average time to natural clearance of malaria from blood)
  omega=1/(365/30.415)  # rate of loss of immunity = 1/(average duration of immunity)
  
  # variable parameters (force of infection)
  lagm=14/30.415        # time from mosquito biting infectious person  until the mosquito is infectious
  theta_cov=0.45        # coverage of vector control
  theta_eff=0           # efficacy of vector control assuming 100% coverage
  alpha=0.797           #** average number of bites per day
  amp=0.9               # relative amplitude of seasonal forcing
  phi=8                 # month of peak in seasonal forcing
  
  # Migration parameters
  magi=1                 # Amplitude of migration
  eta=0.39               # proportion of workers on plantations who are migrant workers
  phi2=5                 # month of peak in migrant worker immigration
  phi3=12                # month of peak in migrant worker emigration
  ksi=0.00               # prevalence of malaria in area from which migrants come
  
  #PCR Paramaters
  PCR_interv_start= 12*3                # duration in months until MDA intervention begins
  PCR_interv_end= (12*3) + 1/30.415      # duration in months PCR campaign
  PCR_cov= 0                             # coverage of PCR campaign
  PCR_eff= 1.0*0.9                       # reduction in transmission given full PCR coverage
  
  # population parameters
  worker_pop <- 10000
  initN <- worker_pop*(1-eta)          # inital population of plantation(s)
  M <- worker_pop*(eta)                # migrant workers per year 
  prev  <- 0.0284                      # prevalence of asymptomatic chronic infections on plantation(s)
  
  # MODEL INITIAL CONDITIONS: S -> L1 -> A1 -> G1 -> R ( -> S) -> L2 -> A2 -> G2 (-> R)
  initR  <- 0                          # inital treated compartment (partial immunity) 
  initL1 <- 0                          # inital liver stage compartment (no immunity) = incidence malaria on plantation (per day) x average duration of liver stage
  initA1 <- 0                          # inital asexual blood-stage compartment (no immunity) = incidence malaria on plantation (per day) x duration of asexual blood stage 
  initG1 <- 0                          # inital sexual gametocyte blood-stage compartment (no immunity) = initA1 x proportion of initA1 not treated x duration of Gametocyte blood stage 
  initL2 <- 0                          # inital liver stage compartment (with partial immunity)
  initA2 <- 0                          # inital asexual blood-stage compartmen (with partial immunity)
  initG2 <- prev*initN                 # inital sexual gametocyte blood-stage compartment (with partial immunity)
  
  initS  <- initN-initL1-initA1-initG1-initR-initL2-initA2-initG2      # initial susceptible compartment
  
  state  <- c(S = initS, L1 = initL1, A1 = initA1, G1 = initG1, R = initR, L2=initL2, A2 = initA2, G2 = initG2)
  
  #MODEL FUNCTION
  malaria<-function(t, state, parameters) 
  {
    
    # Parameters
    
    parameters <- c(
      # basic model parameters
      nu=nu,          # birth rate (1/life expectancy)
      mu=mu,          # death rate (1/life expectancy)
      gamma=gamma,    # rate of movement from liver to asexual blood stage = 1/(average duration of liver stage)
      delta=delta,    # rate of gametocyte maturation  = 1/(average duration of gametocyte maturation)
      p=p,            # proportion of patients who are clinical/symptomatic (no immunity)
      q=q,            # proportion of patients who are clinical/symptomatic (partial immunity)
      tau=tau,        # Proportion of clinical (symptomatic) patients treated
      zeta=zeta,      # rate of natural recovery from malaria infection = (1/average time to natural clearance of malaria from blood)
      omega=omega,    # rate of loss of immunity = 1/(average duration of immunity)
      
      # variable parameters (force of infection)
      lagm=lagm,           # time from mosquito biting infectious person, until the mosquito is infectious
      alpha=alpha,         # average nunber of bites per day
      theta_cov=theta_cov, # coverage of vector control
      theta_eff=theta_eff, # efficacy of vector control assuming 100% coverage
      amp=amp,             # relative amplitude of seasonal forcing
      phi=phi,             # week of peak in seasonal forcing
      
      #PCR Paramaters
      PCR_interv_start=PCR_interv_start,     # number of weeks until MDA intervention begins
      PCR_interv_end=PCR_interv_end,         # duration in weeks of MDA drug effect
      PCR_cov=PCR_cov,                       # coverage of MDA
      PCR_eff=PCR_eff,                       # reduction in transmission given full MDA coverage
      
      # Migration parameters
      magi=magi,           # speed of migration 
      eta=eta,             # proportion of workers on plantations who are migrant workers
      phi2=phi2,           # month of peak in migrant worker immigration
      phi3=phi3,           # month of peak in migrant worker emigration
      ksi=ksi
    )
    
    with(as.list(c(state, parameters)),
         {
           
           # Define lag in Gt (G total) -- the number of people in compartments G1 and G2 at the time mosquitoes were infected "lagm" weeks prior (ie. time it takes from mosquito bite until the mosquito is infectious)
           if (t<=lagm) 
             Gt <- initG1 + initG2
           
           else
             Gt <- lagvalue(t=t-lagm,nr=4) + lagvalue(t=t-lagm,nr=8)
           
           # variables
           N    <- (S+L1+A1+G1+R+L2+A2+G2)
           seas <- 1+amp*cos(2*pi*(t-phi)/12)
           beta <- (1-theta_cov*theta_eff)*alpha*(Gt/N)
           lam  <- beta*seas
           
           # Migration variables
           iota <- 0 #magi*(1+cos(2*pi*(t-phi2)/12))/12 # rate of immigration
           eps <-  0 #magi*(1+cos(2*pi*(t-phi3)/12))/12 # rate of emigration
           ion<-ceiling(t %% 12)== phi2
           eon<-ceiling(t %% 12)== phi3
           
           #PCR
           PCR  <- (PCR_cov*PCR_eff)*(t>(PCR_interv_start))*(t<=(PCR_interv_start+PCR_interv_end))  # PCR campaign
           
           
           # rate of change
           dS  <-  nu*N-mu*S-lam*S+omega*R-eps*M*S/N+iota*M*(1-ksi)-eon*(S/N)*M+(1-ksi)*ion*M
           dL1 <- -mu*L1+lam*S-gamma*L1-eps*M*L1/N-eon*(L1/N)*M
           dA1 <- -mu*A1+gamma*L1-(delta/(1-p*tau))*A1-eps*M*A1/N-eon*(A1/N)*M-PCR*A1
           dG1 <- -mu*G1+delta*A1-(zeta/(1-p*tau))*G1-eps*M*G1/N-eon*(G1/N)*M-PCR*G1
           dR  <-  ((delta*p*tau)/(1-p*tau))*A1+(zeta/(1-p*tau))*G1+((delta*q*tau)/(1-q*tau))*A2+(zeta/(1-q*tau))*G2-(mu+lam+omega)*R-eps*M*R/N-eon*(R/N)*M+PCR*(A1+G1)+PCR*(A2+G2)
           dL2 <- -mu*L2+lam*R-gamma*L2-eps*M*L2/N-eon*(L2/N)*M
           dA2 <- -mu*A2+gamma*L2-(delta/(1-q*tau))*A2-eps*M*A2/N-eon*(A2/N)*M-PCR*A2
           dG2 <- -mu*G2+delta*A2-(zeta/(1-q*tau))*G2-eps*M*G2/N+iota*M*(ksi)-eon*(G2/N)*M+ksi*ion*M-PCR*G2

           # return the rate of change
           list(c(dS, dL1, dA1, dG1, dR, dL2, dA2, dG2))
         }
    ) 
    
  }
  
  #MODEL OUTPUT
  out  <- dede(y = state, times = time.seq, func = malaria, parms = parameters)
  pop  <- out[,"S"]+out[,"L1"]+out[,"A1"]+out[,"G1"]+out[,"R"]+out[,"L2"]+out[,"A2"]+out[,"G2"]
  time <- out[,"time"]
  OUT  <- data.frame(out)
  OUT$Prevalence.pop <- (OUT$A1 + OUT$G1 + OUT$A2 + OUT$G2) #  Population prevalence (asymptomatic and symptomatic people in blood stage)
  
  OUT$year <- time/12


  
  elimination = rbind(elimination, (data.frame(PCR_cov,((OUT$year[which(((OUT$Prevalence.pop < 0.01) & (OUT$Prevalence.pop > 0.009))== T)])[1])-3)))
}

colnames(elimination) <- c("PCRcoverage","EliminationYears")
elimination$PCRcoverage <- elimination$PCRcoverage*100
plot(elimination$PCRcoverage,elimination$EliminationYears, col = "dark blue",lwd=2, type='l',xlab="PCR coverage (%)",ylab="Years to Elimination", xlim= c(min=20, max=100), main="Years to Elimination after a single PCR campaign")

proc.time() - ptm

#plot(OUT$year, OUT$Prevalence.pop/100, col = "dark blue",lwd=2, type='l')
