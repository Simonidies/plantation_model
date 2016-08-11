### CAMBODIAN PLANTATION MALARIA TRANSMISSION MODEL: R CODE (10 AUGUST 2016)  ###
### AUTHOR: SIMON MENDELSOHN, NUFFIELD DEPARTMENT OF MEDICINE, OXFORD UNIVERSITY ###

#CLEAR WORKSPACE
rm(list=ls())

#LOAD LIBRARIES
library(deSolve)
library(ggplot2)
  
  #MODEL DURATION (in years)
  duration <- 10 # years
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
  eta=0.37               # proportion of workers on plantations who are migrant workers
  phi2=5                 # month of peak in migrant worker immigration
  phi3=12                # month of peak in migrant worker emigration
  ksi=0.00               # prevalence of malaria in area from which migrants come
  
  #PCR Paramaters
  PCR_interv_start= 12*3                 # duration in months until MDA intervention begins
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
  OUT$Gt <- OUT$G1 + OUT$G2
  OUT$N <- OUT$S + OUT$L1 + OUT$A1 + OUT$G1 + OUT$R + OUT$L2 + OUT$A2 + OUT$G2
  OUT$Prevalence.perc <- ((OUT$A1 + OUT$G1 + OUT$A2 + OUT$G2)/OUT$N)*100 # Percentage population prevalence (asymptomatic and symptomatic people in blood stage)
  OUT$Prevalence.pop <- (OUT$A1 + OUT$G1 + OUT$A2 + OUT$G2) #  Population prevalence (asymptomatic and symptomatic people in blood stage)
  OUT$Incidence <- ((delta*p*tau)/(1-p*tau))*OUT$A1 + zeta*OUT$G1 + ((delta*q*tau)/(1-q*tau))*OUT$A2 # Monthly reported incidence = everyone in blood stages who is symptomatic and received treatment)
  OUT$Clinical <- (OUT$A1 + OUT$G1)/(OUT$A2 + OUT$G2) #Proportion of cases with no immunity / partial immunity
  OUT$year <- time/12

  # function to round up to neaerest hundred
  roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
  }
  



#############

#MODEL PLOTS
  
  #Change format of dataframe so that all compartments/states are in one columnn, count in 2nd column, year in 3rd column
  cbPalette <- c("light blue", "dark blue", "red", "maroon","green","dark green","black","yellow")
  Stab<-cbind(cbind(rep('S',length(time.seq)),as.numeric(OUT$S)),OUT$year)
  L1tab<-cbind(cbind(rep('L1',length(time.seq)),as.numeric(OUT$L1)),OUT$year)
  A1tab<-cbind(cbind(rep('A1',length(time.seq)),as.numeric(OUT$A1)),OUT$year)
  G1tab<-cbind(cbind(rep('G1',length(time.seq)),as.numeric(OUT$G1)),OUT$year)
  Rtab<-cbind(cbind(rep('R',length(time.seq)),as.numeric(OUT$R)),OUT$year)
  L2tab<-cbind(cbind(rep('L2',length(time.seq)),as.numeric(OUT$L2)),OUT$year)
  A2tab<-cbind(cbind(rep('A2',length(time.seq)),as.numeric(OUT$A2)),OUT$year)
  G2tab<-cbind(cbind(rep('G2',length(time.seq)),as.numeric(OUT$G2)),OUT$year)
  Prev.perc.tab<-cbind(cbind(rep('Prevalence',length(time.seq)),as.numeric(OUT$Prevalence.perc)),OUT$year)
  Prev.pop.tab<-cbind(cbind(rep('Prevalence',length(time.seq)),as.numeric(OUT$Prevalence.pop)),OUT$year)
  Inctab<-cbind(cbind(rep('Incidence',length(time.seq)),as.numeric(OUT$Incidence)),OUT$year)
  MMWInctab<-cbind(cbind(rep('Incidence',length(time.seq)),as.numeric(OUT$Incidence)),OUT$year)
  MMWInctab[,2] <- as.integer(MMWInctab[,2]) # proportion of people who are treated by the MMW
  MMWInctab[,2] <- (as.integer(MMWInctab[,2])/initN)*1000 # incidence per 1000 population

  # Plot 1: All compartments
  dd<-as.data.frame(rbind(Stab,L1tab,A1tab,G1tab,Rtab,L2tab,A2tab,G2tab))
  colnames(dd)<-c('State','Count','Year')
  dd$Count<-as.numeric(as.character(dd$Count))
  dd$Year<-as.numeric(as.character(dd$Year))
  
  p1 <- ggplot(data=dd,aes(x=Year,y=Count,group=State,color=State)) +
    geom_line(size=1) +
    theme_bw()+
    scale_colour_manual(values=cbPalette,breaks=c("S","L1","A1","G1","R","L2","A2","G2")) +
    labs(title = "Compartmental model of malaria transmission", x="Year", y="Plantation inhabitants") +
    theme(legend.title=element_blank())  +
    theme(plot.title = element_text(lineheight=.8, face="bold")) +
    theme(legend.key = element_blank())
     p1
  
  # Plot 2: All compartments, except susceptible population
  dd2<-as.data.frame(rbind(L1tab,A1tab,G1tab,Rtab,L2tab,A2tab,G2tab))
  colnames(dd2)<-c('State','Count','Year')
  dd2$Count<-as.numeric(as.character(dd2$Count))
  dd2$Year<-as.numeric(as.character(dd2$Year))  
  
  p2 <- ggplot(data=dd2,aes(x=Year,y=Count,group=State,color=State)) +
    geom_line(size=1) +
    theme_bw()+
    scale_colour_manual(values=cbPalette,breaks=c("L1","A1","G1","R","L2","A2","G2")) +
    labs(title = "Compartmental model of malaria transmission", x="Year", y="Plantation inhabitants") +
    theme(legend.title=element_blank())  +
    theme(plot.title = element_text(lineheight=.8, face="bold")) +
    theme(legend.key = element_blank())
   p2
  
  # Plot 3: Infected compartments
  dd3<-as.data.frame(rbind(L1tab,A1tab,G1tab,L2tab,A2tab,G2tab))
  colnames(dd3)<-c('State','Count','Year')
  dd3$Count<-as.numeric(as.character(dd3$Count))
  dd3$Year<-as.numeric(as.character(dd3$Year))  
  
  p3 <- ggplot(data=dd3,aes(x=Year,y=Count,group=State,color=State)) +
    geom_line(size=1) +
    theme_bw()+
    scale_colour_manual(values=cbPalette,breaks=c("L1","A1","G1","R","L2","A2","G2")) +
    labs(title = "Compartmental model of malaria transmission", x="Year", y="Plantation inhabitants") +
    theme(legend.title=element_blank())  +
    theme(plot.title = element_text(lineheight=.8, face="bold")) +
    theme(legend.key = element_blank())
    p3
  
  # Plot 4: Prevalence
  dd4<-as.data.frame(rbind(Prev.perc.tab))
  colnames(dd4)<-c('State','Count','Year')
  dd4$Count<-as.numeric(as.character(dd4$Count))
  dd4$Year<-as.numeric(as.character(dd4$Year))  
  
  p4 <- ggplot(data=dd4,aes(x=Year,y=Count, color="red")) +
    geom_line(size=1) +
    theme_bw()+
    theme(plot.title = element_text(lineheight=.8, face="bold")) +
    labs(title = "Malaria prevalence on plantations", x="Year", y="Prevalence (%)") +
    theme(legend.position="none")
    p4
  
  dd4$Date <- seq.Date(as.Date("01/01/07","%d/%m/%y"),by = "day", length.out = 3650)
  plot(dd4[2373:3102,4],dd4[2373:3102,2],type = 'l',main="Model vs data for plantation prevalence\nin Mondulkiri Province" ,col="black",lwd=2, ylim=c(0,5), xlab="Month", ylab ="Prevalence %")
  points(dd4[2709,4],0,col= "red", pch = 19) #MK 0-2.84, ST 0.53-1.68
  arrows(dd4[2709,4],0,dd4[2709,4],0.67,length=0.05, angle=90, code=3,col="red")
  points(dd4[2831,4],2.84,col= "red", pch = 19)
  arrows(dd4[2831,4],1.59,dd4[2831,4],5,length=0.05, angle=90, code=3,col="red")
  legend("bottomleft", inset=.05,
         c("Model   ","Data   "), fill=c("black","red"), horiz=TRUE)

  dd4$Date <- seq.Date(as.Date("01/01/07","%d/%m/%y"),by = "day", length.out = 3650)
  plot(dd4[2373:3102,4],dd4[2373:3102,2],type = 'l',main="Model vs data for plantation prevalence\nin Stung Treng Province" ,col="black",lwd=3, ylim=c(0,3.7), xlab="Month", ylab ="Prevalence %")
  points(dd4[2709,4],0.53,col= "red", pch = 19) #MK 0-2.84, ST 0.53-1.68
  arrows(dd4[2709,4],0.18,dd4[2709,4],1.56,length=0.05, angle=90, code=3,col="red")
  points(dd4[2831,4],1.68,col= "red", pch = 19)
  arrows(dd4[2831,4],0.77,dd4[2831,4],3.61,length=0.05, angle=90, code=3,col="red")
  legend("bottomleft", inset=.05,
         c("Model   ","Data   "), fill=c("black","red"), horiz=TRUE)
  
  # Plot 5: Incidence
  MMWInctab <- as.data.frame(rbind(MMWInctab))
  colnames(MMWInctab)<-c('State','Count','Year')
  MMWInctab$Count<-as.numeric(as.character(MMWInctab$Count))
  MMWInctab$Year<-as.numeric(as.character(MMWInctab$Year))  
  
  p5 <- ggplot(data=MMWInctab,aes(x=Year,y=Count,color="red")) +
    geom_line(size=1) +
    theme_bw()+
    labs(title = "Malaria incidence on plantations", x="Year", y="Treated symptomatic cases \nper 1000 plantation inhabitants") +
    theme(legend.title=element_blank())  +
    theme(plot.title = element_text(lineheight=.8, face="bold")) +
    theme(legend.position="none")
    p5
  
    MMWInctab$Date <- seq.Date(as.Date("01/01/07","%d/%m/%y"),by = "day", length.out = 3650)
    plot(MMWInctab[1095:3650,4],MMWInctab[1095:3650,2]*1/4,type = 'l', col="black",lwd=2, ylab = "Incidence per 1000 population", xlab="")
    
    
  # Plot 6: Phase plot
  plot(MMWInctab[,2],Prev.perc.tab[,2], type = 'l', xlab = "Incidence of symptomatic cases per 1000 plantation inhabitants", ylab = "Prevalence (%)", main = "Phase plot of incidence vs prevalance", col="red", lwd =2)

  # Plot 7: Total population N
  plot(time,OUT$N,type="l",ylim=c(0,11000), ylab = "N")
