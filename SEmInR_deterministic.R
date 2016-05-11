library("deSolve")
library(snowfall)

SEmInR.det <- function(t, x, parms)
{
	### (Code based on Ben Bolker's)
	
	with(c(as.list(parms),as.list(x)), {
		
		S = x[1]
		E = x[2:(nE+1)]
		I = x[(nE+2):(nE+nI+1)]
		R = x[nE+nI+2]
		C = x[(nE+2)]
		D = x[nE+nI+2]
		K = x[(nE+2)] # cumul symptoms onset
		
		sigma2 <- sigma*nE ## transition rate through E boxes
		gamma2 <- gamma*nI ## transition rate through I boxes
		
		infrate = beta*S*sum(I)
		
		dS = mu * (1  - S)
		dS = dS - infrate
		
		## Exposed
		if (nE>1) {
			dE = sigma2*(c(0,E[1:(nE-1)])-E)
		} else dE=-E*sigma
		dE[1] = dE[1] + infrate
		dE = dE - mu*E
		
		## Infected
		if (nI>1) {
			dI = gamma2*(c(0,I[1:(nI-1)])-I)
		} else dI=-I*gamma
		onset = sigma2*E[nE]
		
		dI[1] = dI[1]+onset
		dI = dI - mu*I
		recovery = gamma2*I[nI]
		
		## Recovered 
		dR = (1-f)*recovery - mu*R
		
		## Cumulative Incidence
		dC = infrate
		
		## Cumulative symptoms onset
		dK = onset
		
		## Cumulative RECOVERY
		dD = recovery
		
		res=c(dS, dE, dI, dR, dC, dD, dK)
		list(res)
	})
}

calc.Iall <- function(dat)
{
	col.I <- which(grepl("I",names(dat)))
	if(length(col.I)>1) dat$Iall <- rowSums(dat[,col.I])
	if(length(col.I)==1) dat$Iall <- dat[,col.I]
	return(dat)
}

calc.Eall <- function(dat)
{
	col.I <- which(grepl("E",names(dat)))
	if(length(col.I)>1) dat$Eall <- rowSums(dat[,col.I])
	if(length(col.I)==1) dat$Eall <- dat[,col.I]
	return(dat)
}

simul.SEmInR.det <- function(prm.to.fit, prm.fxd){
	### Solve the SEmInR
	### ODEs equations
	
	# Parameters
	prm <- c(prm.to.fit, prm.fxd)
	horizon         <- prm[["horizon"]]
	nE              <- prm[["nE"]]
	nI              <- prm[["nI"]]
	latent_mean     <- prm[["latent_mean"]]
	infectious_mean <- prm[["infectious_mean"]]
	N               <- prm[["popSize"]] 
	R0              <- prm[["R0"]]
	I.init          <- prm[["init_I1"]]/N
	n.time.steps    <- prm[["n.time.steps"]]
	per.capita      <- prm[["per.capita"]]
	
	dt    <- seq(0, horizon, horizon/n.time.steps)
	sigma <- 1/latent_mean
	gamma <- 1/infectious_mean
	beta  <- R0*gamma
	f     <- 0.0
	mu    <- 0.0
	
	# Initial conditions
	S.init <- 1 - I.init
	E.init <- 0
	
	params.SEmInR <- c(mu=mu, 
					   beta=beta,
					   sigma=sigma,
					   gamma=gamma,
					   f=f,
					   nE=nE, nI=nI)
	
	### Inital conditions
	EIvec <- c(E = rep(0,ifelse(nE==Inf,1,nE)),
			   I = c(I.init,rep(0,ifelse(nI==Inf,0,nI-1))))
	inits.SEmInR <- c(S=1-I.init,EIvec,R=0,C=I.init,D=0,K=I.init)
	
	#### Solutions
	SEmInR <- as.data.frame(lsoda(inits.SEmInR, dt, 
								  SEmInR.det, 
								  parms = params.SEmInR))
	SEmInR <- calc.Iall(SEmInR)  # <- global prevalence
	SEmInR$inc <- c(I.init,diff(SEmInR$C))  # <- global incidence
	
	if(!per.capita){
		SEmInR <- SEmInR * N
		SEmInR$time <- SEmInR$time/N
	}
	return(list(ts=SEmInR, R0=R0))
}


calc.Jall <- function(dat)
{
	col.I <- which(grepl("J",names(dat)))
	if(length(col.I)>1) dat$Jall <- rowSums(dat[,col.I])
	if(length(col.I)==1) dat$Jall <- dat[,col.I]
	return(dat)
}

calc.Yall <- function(dat)
{
	col.I <- which(grepl("Y",names(dat)))
	if(length(col.I)>1) dat$Yall <- rowSums(dat[,col.I])
	if(length(col.I)==1) dat$Yall <- dat[,col.I]
	return(dat)
}

SEmInR.FX <- function(t, x, parms)
{
	### (Code based on Ben Bolker's)
	
	with(c(as.list(parms),as.list(x)), {
		
		S = x[1]
		E = x[2:(nE+1)]
		I = x[(nE+2):(nE+nI+1)]  # prevalence
		R = x[nE+nI+2]
		Y = x[(nE+nI+3):(2*nE+nI+2)]  # proba to be in Ek
		J = x[(2*nE+nI+3):(2*nE+2*nI+2)] # proba to be in Ik
		Z = x[2*nE+2*nI+3]  # cumulative incidence
		
		
		sigma2=sigma*nE ## transition rate through E boxes
		gamma2=gamma*nI ## transition rate through I boxes
		
		infrate = beta*S*sum(I)
		
		dS = mu * (1  - S)
		dS = dS - infrate
		
		## Exposed
		if (nE>1) {
			dE = sigma2*(c(0,E[1:(nE-1)])-E)
		} else dE=-E*sigma
		dE[1] = dE[1] + infrate
		dE = dE - mu*E
		
		## Infected
		if (nI>1) {
			dI = gamma2*(c(0,I[1:(nI-1)])-I)
		} else dI=-I*gamma
		onset = sigma2*E[nE]
		
		dI[1] = dI[1]+onset
		dI = dI - mu*I
		recovery = gamma2*I[nI]
		
		## Recovered 
		dR = (1-f)*recovery - mu*R
		
		
		## Proba being in Ek
		if (nE>1) {
			dY = sigma2*(c(0,Y[1:(nE-1)])-Y)
		} else dY=-Y*sigma
		dY[1] = dY[1] + 0*infrate
		
		## Proba being in Ik
		if (nI>1) {
			dJ = gamma2*(c(0,J[1:(nI-1)])-J)
		} else dJ=-J*gamma
		onset.proba = sigma2*Y[nE]
		
		dJ[1] = dJ[1]+onset.proba
		
		dZ = infrate
		
		res=c(dS, dE, dI, dR, dY, dJ, dZ)
		list(res)
	})
}

calc.intrinsic<- function(seminr.obj){
	R0 <- seminr.obj[["R0"]]
	seminr <- seminr.obj[["ts"]]
	
	delta.t <- seminr$time[2]-seminr$time[1]
	message(paste("R0 =",R0))
	message(paste("delta.t =",delta.t))
	
	g = vector()
	g[1] <- 0  #seminr$inc[1]/delta.t/R0/seminr$S[1]/seminr$inc[1]
	
	for(k in 2:(length(seminr$inc)-1)){
		tmp1 <- seminr$inc[k]/delta.t/R0/seminr$S[k]/seminr$inc[1]
		tmp2 = 0 
		for(q in 1:(k-1)) tmp2 = tmp2 + g[q]*seminr$inc[k-q+1]/seminr$inc[1]
		g[k] <- min(1,max(0,tmp1-tmp2))
	}
	par(mfrow=c(1,2))
	plot(x=seminr$time[-1], g, typ="l",main="GI")
	zoom<-40
	zoom.idx = which(seminr$time<zoom)
	plot(x=seminr$time[zoom.idx], g[zoom.idx],typ="o",main="GI")
	lines(dexp(x = seminr$time[zoom.idx],rate = 1/5),col="red")
	return(g)
}

