library(snowfall)
# source("SEmInR_deterministic.R")

### Summary statistics
sum.stat <- function(target, x){
	n <- length(target)
	
	# Summary stat #1
	s1 <- sqrt(sum(  (target-x[1:n])^2) )
	
	# Summary stat #2
	t.obs <- 1:n
	d.target <- data.frame(t.obs,logtarget=log(1+target))
	lm.target <- lm(data = d.target, 
					formula = logtarget ~ poly(t.obs,degree = 3, raw=TRUE))
	lm.target.coef <-  lm.target$coefficients
	ftarget <- lm.target.coef[1] + lm.target.coef[2]*t.obs + lm.target.coef[3]*t.obs^2 + lm.target.coef[4]*t.obs^3
	
	d.x <- data.frame(t.obs,logx=log(1+x))
	lm.x <- lm(data = d.x, 
			   formula = logx ~ poly(t.obs,degree = 3, raw=TRUE))
	lm.x.coef <-  lm.x$coefficients
	
	fx <- lm.x.coef[1] + lm.x.coef[2]*t.obs + lm.x.coef[3]*t.obs^2 + lm.x.coef[4]*t.obs^3
	s2 <- sqrt(sum( (fx-ftarget)^2 ))
	
	if(F){ # DEBUG
		plot(t.obs,log(1+target), main=s2,typ='o',ylim=range(log(1+target),fx))
		lines(t.obs,log(1+x),col='red',typ='o',pch=15)
		lines(t.obs,ftarget,col='blue',typ='l',pch=15)
		lines(t.obs,fx,col='red',typ='l',pch=16)
	}
	return(s2)
}

### Sample parameter values from prior distributions
sample.priors <- function(priors, n){
	
	# number of parameters that will be sampled:
	n.priors <- length(priors)
	
	# matrix that will hold all prior samples
	M <- matrix(nrow = n, ncol = n.priors)
	for(j in 1:n.priors){
		z <- priors[[j]]
		M[,j] <- do.call(what = paste0("r",z[[1]]) , args = c(n, z[[2]]))
	}
	colnames(M) <- names(priors)
	return(M)
}

### Wrap function for snowfall (parallel execution)
fct.wrap <- function(i, t.obs, inc.obs,  M, prm.fxd, nABC){
	sim <- simul.SEmInR.det(prm.to.fit = M[i,], prm.fxd)
	df  <- sim$ts
	inc <- df$inc
	tt  <- df$time
	
	dt  <- tt[2]-tt[1]
	idx <- vector()
	for(k in 1:length(t.obs)) { idx[k] <- which(abs(tt-t.obs[k]) <= dt/1.99)[1] }
	t <- tt[idx]
	
	return(sum.stat(target = inc.obs, x = inc[idx]))
}

## Fit SEmInR with ABC
fit.ABC.SEmInR <- function(t.obs, inc.obs, 
						   prm.fxd, priors.prm.to.fit, 
						   nABC, post.prop){
	
	M <- sample.priors(priors.prm.to.fit, nABC)
	
	# Parallel execution:
	sfInit(parallel = TRUE, cpu = parallel::detectCores())
	sfLibrary(deSolve)
	sfExportAll()
	ss <- sfSapply(x = 1:nABC, fun = fct.wrap, 
				   t.obs=t.obs, inc.obs=inc.obs,
				   M=M, prm.fxd=prm.fxd, nABC=nABC, 
				   simplify = FALSE)
	sfStop()
	
	M <- cbind(M, ss=unlist(ss))
	M <- M[order(M[,'ss']),]
	thres.row <- round(nABC*post.prop)
	Mpost <- M[1:thres.row,]
	return(list(posteriors=Mpost, priors=priors))
}


### Sample incidence from posterior distribution
### and return median and credible interval
post.incidence <- function(Mpost,t,CI){
	Mpost2 <- Mpost[,-ncol(Mpost)]
	n <- nrow(Mpost2)
	inc.post <- matrix(ncol=nrow(Mpost2),nrow=length(t))
	for(i in 1:n){
		if(i%%10==0) cat('post incidence for SEmInR ',i,'/',n,'\n')
		sim <- simul.SEmInR.det(prm.to.fit=Mpost2[i,], prm.fxd)
		df  <- sim$ts
		inc.post[,i] <- df$inc
	}
	inc.md <- apply(X = inc.post,MARGIN = 1, FUN = quantile, probs= 0.5)
	inc.lo <- apply(X = inc.post,MARGIN = 1, FUN = quantile, probs= 0.5-CI/2)
	inc.hi <- apply(X = inc.post,MARGIN = 1, FUN = quantile, probs= 0.5+CI/2)
	
	return(list(inc.md=inc.md, inc.lo=inc.lo, inc.hi=inc.hi))
}







