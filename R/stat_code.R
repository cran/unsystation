unsys.station.test <- function(x, M=2000, sig.lev=.05, max.scale=NULL, m=NULL, B=200, eps=5, use.all=FALSE, do.parallel=FALSE){

	T <- length(x)
	if(is.null(max.scale)) max.scale <- round(log(log(T, 2), 2))
	if(is.null(m)) m <- round(sqrt(T))
	if(do.parallel){ cl <- parallel::makeCluster(4); doParallel::registerDoParallel(cl) }

	top.cand0 <- bottom.cand0 <- NULL
	y.mat <- matrix(0, ncol=max.scale, nrow=T)
	for(k in 1:max.scale){
	  y.mat[, k] <- y <- func_coef(x, -k)^2
	  ref <- sort(y, decreasing=TRUE, index.return=TRUE)
	  top.cand0 <- c(top.cand0, setdiff(ref$ix[1:(eps + 2*2^k)], c(1:(2^k), (T-2^k+1):T))[1:eps])
	  bottom.cand0 <- c(bottom.cand0, setdiff(ref$ix[T:(T-eps-2*2^k+1)], c(1:(2^k), (T-2^k+1):T))[1:eps])
	}
	top.cand <- sample(top.cand0, M, replace=TRUE); bottom.cand <- sample(bottom.cand0, M, replace=TRUE)

	fr <- funcRes(y.mat, M, m, rep(1/T, T), top.cand, bottom.cand, apply(y.mat, 2, mean))

	ind <- which((!duplicated(fr$res[, 5])) & fr$res[, 5] > 0)
	if(use.all) ref <- ind else{
	  ref <- ind[sort(fr$res[ind, 4+max.scale+1], decreasing=TRUE, index.return=TRUE)$ix[1:(10*M)]]
	}
	se.mat <- fr$res[ref, 1:4]
	I <- length(ind); R <- length(ref)

	arx <- ar(x, order.max=log(T), method='yw')
	coef <- arx$ar
	ep <- arx$resid[!is.na(arx$resid)];
	sig <- 1.4826*median(abs(ep-median(ep)))
	ep <- ep[abs(ep-median(ep)) < sig*qt(1-.005, 10)]
	ep <- ep-mean(ep)
	if(length(coef)==0) boot.x <- matrix(sample(ep, B*T, replace=TRUE), ncol=B) else boot.x <- funcSimX(coef, matrix(sample(ep, (T+length(coef))*B, replace=TRUE), ncol=B))

	b <- 0
	if(do.parallel){
		null.stat <- foreach::foreach(b=iterators::iter(1:B), .combine=rbind, .packages=c('Rcpp', 'RcppArmadillo', 'unsystation')) %dopar% {
		  bx <- boot.x[, b]
		  by.mat <- matrix(0, ncol=max.scale, nrow=T)
		  for(k in 1:max.scale){
				by.mat[, k] <- func_coef(bx, -k)^2
			}
			tmp <- funcResVar(by.mat, se.mat, apply(by.mat, 2, mean))
			c(tmp)
		}
	} else{
		null.stat <- foreach::foreach(b=iterators::iter(1:B), .combine=rbind, .packages=c('Rcpp', 'RcppArmadillo', 'unsystation')) %do% {
		  bx <- boot.x[, b]
		  by.mat <- matrix(0, ncol=max.scale, nrow=T)
		  for(k in 1:max.scale){
				by.mat[, k] <- func_coef(bx, -k)^2
			}
			tmp <- funcResVar(by.mat, se.mat, apply(by.mat, 2, mean))
			c(tmp)
		}
	}
	stat <- abs(fr$res[ref, 4+1:max.scale])/funcApplyVar(null.stat, max.scale, R)

	k <- which.max(apply(stat, 1, max))
	intervals <- fr$res[ref[k], 1:4]
	test.stat <- max(stat[k,])
	test.criterion <- qnorm(1-sig.lev/2/I/max.scale)
	test.res <- test.stat > test.criterion

	if(do.parallel) parallel::stopCluster(cl)

	return(list(intervals=intervals, test.stat=test.stat, test.criterion=test.criterion, test.res=test.res))

}
