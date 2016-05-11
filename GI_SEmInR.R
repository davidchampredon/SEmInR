# ===================================================================
# Calculate generation interval distribution
# for a SEmInR model. 
#
# Method explained in:
# Champredon D, Dushoff J. 
# Intrinsic and realized generation intervals in infectious-disease transmission. 
# Proceedings of the Royal Society B: Biological Sciences 2015; 282: 20152026.
# ===================================================================

# source("SEmInR_deterministic.R")
library("deSolve")
library("ggplot2")
library("gridExtra")


GI.intrinsic <- function(tau, df){
  res = 0
  N = length(df$time)
  dt = df$time[2]-df$time[1]
  JJ = df$Jall[1:N]
  integ  = sum(JJ)*dt
  res = df$Jall[tau]/integ
  return(res)
}

GI.fwd <- function(tau, s, df){
  res = 0
  N = length(df$time)
  dt = df$time[2]-df$time[1]
  
  JJ = df$Jall[1:(N-s)]
  SS = df$S[(s+1):N]
  integ  = sum(JJ*SS)*dt
  res = df$Jall[tau]*df$S[tau+s]/integ
  return(res)
}


GI.bck <- function(tau, t, df){
  res = 0
  
  if(tau<t){
    N = length(df$time)
    dt = df$time[2]-df$time[1]
    
    JJ = df$Jall[1:(t-1)]
    II = df$inc[(t-1):1]
    integ = sum(JJ*II)*dt
    res = df$Jall[tau]*df$inc[t-tau]/integ
  }
  return(res)
}


calc.GI.SEmInR <- function(latent_mean,
						   infectious_mean,
						   nE,nI,
						   n.points.GI.crv,
						   horizon){
	### Parameters
	tvec <- seq(0,horizon,0.2)
	sigma <- 1/latent_mean
	gamma <- 1/infectious_mean
	R0 = 1.234567 # <-- doesn't matter here - to do : recode properly eq (3.2) in Champredon RSPB 2015 
	beta <- R0*gamma
	f <- 0.0
	mu <- 0.00
	
	# Initial conditions
	I.init <- 1e-5 
	S.init <- 1 - I.init
	E.init <- 0
	
	params.SEmInR <- c(mu=mu, 
					   beta=beta,
					   sigma=sigma,
					   gamma=gamma,
					   f=f,
					   nE=nE, nI=nI)
	
	### Inital conditions
	# W A R N I N G : 
	# initial infectious in I[1] compartment (not E[1])
	EIvec <- c(E= rep(0,ifelse(nE==Inf,1,nE)),
			   I= c(I.init,rep(0,ifelse(nI==Inf,0,nI-1))) )
	YJvec <- c(Y=c(1,rep(0,ifelse(nE==Inf,1,nE-1))),
			   J=rep(0,ifelse(nI==Inf,0,nI) ))
	inits.SEmInR <- c(S=1-I.init, EIvec, R=0, YJvec, Z=I.init)
	
	#### Solutions ####
	
	SEmInR <- as.data.frame(lsoda(inits.SEmInR, tvec, SEmInR.FX, 
								  parms=params.SEmInR))
	SEmInR <- calc.Iall(dat = SEmInR)
	SEmInR <- calc.Jall(dat = SEmInR)
	SEmInR <- calc.Yall(dat = SEmInR)
	# incidence (from solved cumulative incidence)
	SEmInR$inc <- c(SEmInR$Z[1],diff(SEmInR$Z))
	
	##########################
	#### Generation times ####
	##########################
	
	NT = nrow(SEmInR)
	
	# -- mean forward & backward generation interval
	tt = vector()
	f.bar = vector()
	g.bar = vector()
	
	loop.idx = round(seq(1,NT-1,length.out=n.points.GI.crv))
	dt = SEmInR$time[2]-SEmInR$time[1]
	
	# calculate gi.intrinsic
	gi.intrinsic = vector()
	for(tau in 1:NT) gi.intrinsic[tau] <- GI.intrinsic(tau,SEmInR)
	# theoretical mean: 
	theo.mean.gii <- latent_mean + infectious_mean*(nI+1)/2/nI
	# numerical mean:
	mean.gii <- sum(SEmInR$time*gi.intrinsic*dt)
	# variance:
	var.gii <- sum(SEmInR$time^2*gi.intrinsic*dt) - mean.gii^2
	
	cnt = 1
	
	for(s in loop.idx){
		# calculate gi.fwd:
		gi.fwd = vector()
		for(tau in 1:(NT-s))  gi.fwd[tau] <- GI.fwd(tau, s, SEmInR)
		# calculate expectation of gi.fwd:
		tt[cnt] = SEmInR$time[s]
		tvec = SEmInR$time[c(1:(NT-s))]
		f.bar[cnt] = sum(gi.fwd*tvec*dt)
		
		# calculate gi.bck:
		gi.bck = vector()
		for(tau in 1:(s-1)){
			gi.bck[tau] <- GI.bck(tau, s, SEmInR)
		}  
		# calculate expectation of gi.bck:
		tvec.bck = SEmInR$time[c(1:(s-1))]
		g.bar[cnt] = sum(gi.bck*tvec.bck*dt)
		
		cnt = cnt+1
	}
	return(list(intrinsic = gi.intrinsic,
				fwd = gi.fwd,
				bck = gi.bck,
				times = SEmInR$time))
}


calc.theoretical.GI.shiny <- function(latent_mean,
                                      infectious_mean,
                                      R0,
                                      nE,nI,
                                      n.points.GI.crv,
                                      horizon,
                                      do.plot = TRUE){
  ### Parameters
  dt <- seq(0,horizon,0.2)
  sigma <- 1/latent_mean
  gamma <- 1/infectious_mean
  beta <- R0*gamma
  f <- 0.0
  mu <- 0.00
  
  # Initial conditions
  I.init <- 1e-5 #prm[prm$V1=="init_I1",2]/N
  S.init <- 1 - I.init
  E.init <- 0
  
  params.SEmInR <- c(mu=mu, 
                     beta=beta,
                     sigma=sigma,
                     gamma=gamma,
                     f=f,
                     nE=nE, nI=nI)
  
  ### Inital conditions
  # W A R N I N G : 
  # initial infectious in I[1] compartment (not E[1])
  EIvec <- c(E= rep(0,ifelse(nE==Inf,1,nE)),
             I= c(I.init,rep(0,ifelse(nI==Inf,0,nI-1))) )
  YJvec <- c(Y=c(1,rep(0,ifelse(nE==Inf,1,nE-1))),
             J=rep(0,ifelse(nI==Inf,0,nI) ))
  inits.SEmInR <- c(S=1-I.init, EIvec, R=0, YJvec, Z=I.init)
  
  #### Solutions ####
  
  SEmInR <- as.data.frame(lsoda(inits.SEmInR, dt, SEmInR.FX, 
                                parms=params.SEmInR))
  SEmInR <- calc.Iall(dat = SEmInR)
  SEmInR <- calc.Jall(dat = SEmInR)
  SEmInR <- calc.Yall(dat = SEmInR)
  
  # incidence (from solved cumulative incidence)
  SEmInR$inc <- c(SEmInR$Z[1],diff(SEmInR$Z))
  
  
  ##########################
  #### Generation times ####
  ##########################
  
  NT = nrow(SEmInR)
  
  # -- mean forward & backward generation interval
  tt = vector()
  f.bar = vector()
  g.bar = vector()
  
  loop.idx = round(seq(1,NT-1,length.out=n.points.GI.crv))
  dt = SEmInR$time[2]-SEmInR$time[1]
  
  # calculate gi.intrinsic
  gi.intrinsic = vector()
  for(tau in 1:NT) gi.intrinsic[tau] <- GI.intrinsic(tau,SEmInR)
  # theoretical mean: 
  theo.mean.gii <- latent_mean + infectious_mean*(nI+1)/2/nI
  # numerical mean:
  mean.gii <- sum(SEmInR$time*gi.intrinsic*dt)
  # variance:
  var.gii <- sum(SEmInR$time^2*gi.intrinsic*dt) - mean.gii^2
  
  cnt = 1
  
  for(s in loop.idx){
    # calculate gi.fwd:
    gi.fwd = vector()
    for(tau in 1:(NT-s))  gi.fwd[tau] <- GI.fwd(tau, s, SEmInR)
    
    # calculate expectation of gi.fwd:
    tt[cnt] = SEmInR$time[s]
    
    tvec = SEmInR$time[c(1:(NT-s))]
    f.bar[cnt] = sum(gi.fwd*tvec*dt)
    
    # calculate gi.bck:
    gi.bck = vector()
    for(tau in 1:(s-1)){
      gi.bck[tau] <- GI.bck(tau, s, SEmInR)
    }  
    
    # calculate expectation of gi.bck:
    tvec.bck = SEmInR$time[c(1:(s-1))]
    g.bar[cnt] = sum(gi.bck*tvec.bck*dt)
    
    cnt = cnt+1
  }
  
  GI.ODE <- data.frame(time=tt,
                       GI.fwd.mean = f.bar,
                       GI.bck.mean = g.bar)
  
  
  #############################
  ###### PLOTS TO CHECK #######
  #############################
  
  if(do.plot){
    par(mfrow=c(2,3),cex.main=2,cex.axis=2,cex.lab=2)
    
    TIME = SEmInR$time
    
    ### PREVALENCE & INCIDENCE
    
    #     plot(x=TIME, y=SEmInR$Iall, typ="l", col="red",lwd=6, main="Prevalence")
    plot(x=TIME, y=SEmInR$inc, typ="l", 
         col="blue",lwd=6, main="Incidence")
    
    ### PROBABILITIES
    
    # time since infection
    t.inf <- TIME[which(TIME<3*(latent_mean+infectious_mean))]
    n.t.inf <- length(t.inf)
    
    plot(x=t.inf, y=SEmInR$Yall[1:n.t.inf], 
         typ="l",lwd=6, 
         xlab="Time since infection", ylab="",las=1,
         main="Probability of being in E[k]")
    col.Y <- which(grepl("Y",names(SEmInR)))
    
    for(i in 1:length(col.Y)){
      col.i = 1-i/length(col.Y)
      lines(x=t.inf, y=SEmInR[1:n.t.inf,col.Y[i]], 
            col=rgb(col.i,col.i,col.i))
    }
    abline(v=latent_mean,lty=3)
    
    plot(x=t.inf, y=SEmInR$Jall[1:n.t.inf], 
         typ="l",lwd=6, 
         xlab="Time since infection", ylab="",las=1,
         main="Probability of being in I[k]")
    col.J <- which(grepl("J",names(SEmInR)))
    
    for(i in 1:length(col.J)){
      col.i = 1-i/length(col.J)
      lines(x=t.inf, y=SEmInR[1:n.t.inf,col.J[i]], 
            col=rgb(col.i,col.i,col.i))  
    }
    abline(v=latent_mean+infectious_mean*(nI+1)/2/nI,lty=2)
    
    ### GENERATION INTERVALS
    
    # Forward:
    
    plot(x=tt,y=f.bar,pch=16,typ="o",
         main="Mean Forward GI",
         xlab = "Calendar time",
         ylim=c(0,max(f.bar,latent_mean+infectious_mean,na.rm = T)),
         las=1, lwd=6)
    abline(h=mean.gii,lty=2)
    
    # Backward:
    
    plot(x=tt,y=g.bar,pch=16,typ="o",
         main="Mean Backward GI",
         xlab = "Calendar time",
         ylim=c(0,max(g.bar,latent_mean+infectious_mean,na.rm = T)),
         las=1, lwd=6)
    abline(h=mean.gii,lty=2)
  }
  
  # Intrinsic:
  
  # equivalent gamma distribution
  shape <- mean.gii^2/var.gii
  rate <- mean.gii/var.gii
  equiv.gamma <- dgamma(x=t.inf, shape=shape, rate=rate)
  
  title=paste0("Intrinsic GI \n",
               "lat.mean=",latent_mean,
               " ; infec.mean=",infectious_mean,
               " ; nE=",nE, " ; nI=",nI)
  
  plot(x=t.inf, y=gi.intrinsic[1:n.t.inf], 
       typ="l", 
       main = title,
       xlab="Time Since Infection",
       ylab="",las=1,
       ylim=range(gi.intrinsic[1:n.t.inf],equiv.gamma),
       lwd=6)
  lines(x=t.inf,y=equiv.gamma,
        col=rgb(1,0,0,0.6),lwd=6,lty=3)
  abline(v=mean.gii,lty=2,lwd=3)
  abline(v=theo.mean.gii,lty=1,col="green",lwd=3)
  
  # quantiles
  qq=qgamma(p=c(0.025,0.975),shape =shape,rate=rate)
  xxq = t.inf[t.inf>qq[1] & t.inf<qq[2]]
  yyq=dgamma(xxq,shape =shape,rate=rate)
  polygon(x=c(xxq,rev(xxq)),border = NA,
          y=c(yyq,rep(0,length(yyq))),
          col=rgb(0,0,0,0.2))
  
  qq=qgamma(p=c(0.25,0.75),shape =shape,rate=rate)
  xxq = t.inf[t.inf>qq[1] & t.inf<qq[2]]
  yyq=dgamma(xxq,shape =shape,rate=rate)
  polygon(x=c(xxq,rev(xxq)),border = NA,
          y=c(yyq,rep(0,length(yyq))),
          col=rgb(0,0,0,0.2))
  
  grid()
  legend(x = "topright", legend = c("GI integrated","Gamma(m,v)"),cex = 2,
         lwd=8, col=c("black","red"),lty=c(1,3))
  
  return(list(GI.ODE=GI.ODE, 
              SEmInR=SEmInR, 
              GI.intrinsic=gi.intrinsic)
  )
}