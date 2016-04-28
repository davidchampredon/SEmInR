library("deSolve")
library(snowfall)

SEmInR.beta.gamma <- function(t, x, parms)
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

simul.SEmInR <- function(prm.to.fit, prm.fxd){
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
								  SEmInR.beta.gamma, 
								  parms = params.SEmInR))
	SEmInR <- calc.Iall(SEmInR)  # <- global prevalence
	SEmInR$inc <- c(I.init,diff(SEmInR$C))  # <- global incidence
	
	if(!per.capita){
		SEmInR <- SEmInR * N
		SEmInR$time <- SEmInR$time/N
	}
	return(list(ts=SEmInR, R0=R0))
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

llk.pois <- function(prm.to.fit, prm.fxd, t.obs, inc.obs) {
	### RETURNS POISSON LIKELIHOOD GIVEN INCIDENCE DATA
	### (assume Poisson distributed observation errors)
	
	sim <- simul.SEmInR(prm.to.fit, prm.fxd)
	df  <- sim$ts
	sim.inc <- df$inc[df$time %in% t.obs]
	llk <- -sum(dpois(x = inc.obs, lambda=sim.inc, log=TRUE))
	return(llk)
}


fit.mle.SEmInR <- function(prm.to.fit, prm.fxd, 
						   t.obs, inc.obs,
						   method,
						   upper=Inf, lower=-Inf, maxit = 500){
	### FIT A SEmInR MODEL TO INCIDENCE OBSERVATIONS
	### (observation error ~ Poisson)
	param.fit <- optim(par     = prm.to.fit,
					   fn      = llk.pois,
					   prm.fxd = prm.fxd,
					   t.obs   = t.obs,
					   inc.obs = inc.obs,
					   method  = method, #'SANN',"L-BFGS-B",
					   lower   = lower, 
					   upper   = upper,
					   control = list(trace = 2, maxit = maxit))
	prm.fitted <- param.fit$par
	llkmin <- param.fit$value
	return(list(prm.fitted=prm.fitted, llkmin=llkmin))
}

sf.outer <- function(x,y,fun, ...) {
	### PARALLEL VERSION OF 'outer' FUNCTION
	### USING PACKAGE 'snowfall'
	
	sfInit(parallel = TRUE, cpu = parallel::detectCores())
	sfLibrary(deSolve)
	
	fun2 <- function(i,m,xx,yy, ...){
		for(j in 1:length(yy)){
			m[i,j] <- fun(xx[i],yy[j], ...)
		}
		return(m)
	}
	m <- matrix(nrow=length(x),ncol=length(y))
	
	sfExportAll()
	zz <- sfSapply(x=1:length(x), fun = fun2, 
				   m=m, xx=x, yy=y, ...,
				   simplify = FALSE)
	sfStop()
	
	res <- list()
	for(k in 1:length(zz)){
		res[[k]] <- zz[[k]][k,]
	}
	res <- do.call("rbind",res)
	return(res)
}


llk.surf <- function(grid.size, prm.to.fit,
					 boundaries,
					 prm.fxd,
					 t.obs,
					 cases.obs,
					 prm.fitted = NULL,
					 paral = TRUE){
	### CALCULATE THE LIKELIHOOD SURFACE (GIVEN A GRID)
	###
	wrapfct <- function(x,y,xynames,
						prm.fxd,
						t.obs,cases.obs){
		z <- c(x, y)
		names(z) <- xynames
		llk.pois(z,
				 prm.fxd, 
				 t.obs,
				 inc.obs = cases.obs)
	}
	wrapfct.vec <- Vectorize(wrapfct, list("x","y"))
	
	lsurf <- list()
	lsurf.idx <- list()
	pval <- list()
	nf <- length(prm.to.fit)
	cnt <- 1
	for(i in 1:nf)
		for(j in 1:nf){
			if(j>i){
				print(paste("likelihood surface:", i,j))
				pi <- seq(boundaries[i,1], boundaries[i,2],length.out = grid.size)
				pj <- seq(boundaries[j,1], boundaries[j,2],length.out = grid.size)
				pnames <- names(prm.to.fit)[c(i,j)]
				prm.fxd2 <- c(prm.fxd, prm.to.fit[-c(i,j)])
				
				if(!paral){
					lsurf[[cnt]] <- outer(X = pi,Y = pj, FUN = wrapfct.vec,
										  xynames = pnames,
										  prm.fxd = prm.fxd2,
										  t.obs, cases.obs)
				}
				if(paral){
					lsurf[[cnt]] <- sf.outer(x = pi,y = pj, fun = wrapfct.vec, 
											 xynames = pnames,
											 prm.fxd = prm.fxd2,
											 t.obs, cases.obs)
				}
				# - - - - 
				
				lsurf.idx[[cnt]] <- c(i,j)
				pval[[i]] <- pi
				if(j==nf) pval[[nf]] <- pj
				
				contour(pi,pj,lsurf[[cnt]],
						nlevels = 20,
						main = "Negative Log-Likelihood Contours",
						xlab=pnames[1], ylab=pnames[2])
				if(!is.null(prm.fitted)){
					points(x=prm.fitted[i],y=prm.fitted[j],
						   col="red",pch=16,cex=2)	
				}
				cnt <- cnt + 1
			}
		}
	return(list(lsurf = lsurf,
				lsurf.idx = lsurf.idx,
				pval = pval))
}

CI.llk.sample <- function(CIlevel, nsample,  
						  prm.fitted, llkmin,
						  prm.fxd, 
						  t.obs, inc.obs,
						  maxit = 50,
						  prop.search = 0.5){
	
	### Return parameter values that _should_ be in the CI region
	
	nf    <- length(prm.fitted)
	thres <- llkmin + qchisq(CIlevel,nf)/2
	
	radius <- 0.1 * prm.fitted
	
	cnt <- 1
	lk <- thres-1
	vp <- list()
	incr.radius <- 1.1
	
	# Try to find approximately the CI region
	while(mean(lk>thres)<prop.search & cnt<maxit){
		print(paste('Searching CI region',cnt,'<',maxit))
		lk <- vector()
		for(k in 1:2^nf){
			# going in random directions:
			vp[[k]] <- prm.fitted +radius*rbinom(n=nf,1,0.5)
			lk[k]   <- llk.pois(prm.to.fit = vp[[k]], prm.fxd, t.obs, inc.obs )
		}
		radius <- incr.radius*radius
		cnt <- cnt +1
	}
	if(cnt>=maxit) warning('CI.llk.sample: Search CI reached maximum iterations!')
	# Within the approximate region found above,
	# sample parameter values.
	vp.inside <- list()
	for(i in 1:nsample){
		print(paste('sampling',i,'/',nsample))
		radius2 <- radius/incr.radius*runif(1,0,1)
		vp.inside[[i]] <- prm.fitted +radius2*rbinom(n=nf,1,0.5)
	}
	
	# vp.inside is a list of vectors where more than half 
	# should be (?!) within the CIlevel, 
	# as defined by the llk ratio test (chi squared)
	return(vp.inside)
}

#- - - - - - - - - -
# - - - - - - - - - -
# - - - - - - - - - -
	
runtest <- FALSE

if(runtest){
	
	### Example code using the functions 
	### defined above: 
	### - Make synthetic data
	### - fit SEmInR
	### - forecast future incidence with SEmInR
	
	set.seed(1234)
	t1 <- as.numeric(Sys.time())
	
	# Simple simulation:
	prm.to.fit <- c(
		latent_mean     = 4,
		infectious_mean = 2,
		popSize         = 10000,
		R0              = 1.9)
	
	prm.fxd <-  c(
		horizon = 100,
		nE = 2,
		nI = 2,
		init_I1 = 1,
		n.time.steps = 500,
		per.capita = FALSE
	)
	# Synthetic data
	sim <- simul.SEmInR(prm.to.fit, prm.fxd)
	df  <- sim$ts
	true.inc <- df$inc
	tt <- df$time
	t.obs <- seq(5,40,by=5)
	round.inc <- ceiling(true.inc[tt %in% t.obs])
	cases.obs <- rpois(n=length(t.obs), lambda = round.inc)
	
	# Test likelihood function:
	llk <- llk.pois(prm.to.fit,
					prm.fxd, 
					t.obs,
					inc.obs = cases.obs)
	print(llk)
	
	# Fit SEmInR model to incidence data:
	boundaries <- matrix(data=c(0.5,5,
								0.6,4,
								0.88,7),
						 ncol = 2, byrow = TRUE)
	
	init.prm.fit <- prm.to.fit*runif(n=length(prm.to.fit),0.5,2)
	
	FIT <- fit.mle.SEmInR(prm.to.fit = init.prm.fit, 
						  prm.fxd, 
						  t.obs, 
						  inc.obs = cases.obs,
						  # lower=boundaries[,1], 
						  # upper=boundaries[,2], 
						  method = "SANN",#"SANN",# "CG",# "Nelder-Mead",#'SANN',#"L-BFGS-B",
						  maxit = 80)
	prm.fitted <- FIT[['prm.fitted']]
	llkmin     <- FIT[['llkmin']]    
	print(prm.fitted)
	
	# Try to approximate CI value:
	cival <- CI.llk.sample(CIlevel = 0.95, 
						   nsample = 100, 
						   prm.fitted, 
						   llkmin, prm.fxd, t.obs,
						   inc.obs = cases.obs,
						   prop.search = 0.5)
	
	
	# Simulate forward with fitted data (point estimate):
	sim.fit <- simul.SEmInR(prm.fitted, prm.fxd)
	df.fit  <- sim.fit$ts
	
	# Simulate forward with fitted data (CI envelop):
	inc.CI <- list()
	for(s in 1:length(cival)){
		print(paste('simulating CI incidence',s,'/',length(cival)))
		sim.CI <- simul.SEmInR(cival[[s]],prm.fxd)
		inc.CI[[s]] <- sim.CI$ts$inc
	}
	
	m.CI <- matrix(unlist(inc.CI),ncol=length(inc.CI[[1]]),byrow = TRUE)
	inc.ci.lo <- apply(m.CI,MARGIN = 2,FUN=min)
	inc.ci.hi <- apply(m.CI,MARGIN = 2,FUN=max)
	
	# Likelihood surfaces (for each parameter to fit pairs):
	if(FALSE){  # <-- take some time and not used, just for example.
	lsurf <- llk.surf(grid.size = 20, prm.to.fit ,boundaries,
					  prm.fxd ,t.obs ,cases.obs, prm.fitted = prm.fitted )
	}
	
	
	#  - - - Plots - - - #
	
	inc.best <- df.fit$inc
	plot(tt,true.inc,typ='l',lty=2,
		 ylim=1+range(inc.best,inc.ci.lo,inc.ci.hi),
		 log='')
	lines(df.fit$time,inc.best, 
		  lwd = 2,
		  typ='l')
	for(s in 1:length(cival)){
		lines(df.fit$time,inc.CI[[s]], col=rgb(0,0,0,0.1))
	}
	lines(df.fit$time, inc.ci.lo,col='blue')
	lines(df.fit$time, inc.ci.hi,col='blue')
	points(t.obs,cases.obs,pch=16,col='orange',cex=2)
	abline(v=max(t.obs), lwd=2, lty=2, col='orange')
	grid()
	t2 <- as.numeric(Sys.time())
	print(prm.fitted)
	message(paste("completed in",round((t2-t1)/60,1),'minutes'))
}
