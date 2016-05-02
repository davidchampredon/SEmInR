###
###  MAXIMUM LIKELIHOOD FIT FOR A DETERMINISTIC SEmInR
###


llk.pois <- function(prm.to.fit, prm.fxd, t.obs, inc.obs) {
	### RETURNS POISSON LIKELIHOOD GIVEN INCIDENCE DATA
	### (assume Poisson distributed observation errors)
	
	sim <- simul.SEmInR.det(prm.to.fit, prm.fxd)
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
		if(i%%10==0) print(paste('sampling',i,'/',nsample))
		radius2 <- radius/incr.radius*runif(1,0,1)
		vp.inside[[i]] <- prm.fitted +radius2*rbinom(n=nf,1,0.5)
	}
	
	# vp.inside is a list of vectors where more than half 
	# should be (?!) within the CIlevel, 
	# as defined by the llk ratio test (chi squared)
	return(vp.inside)
}

