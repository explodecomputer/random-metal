DL = function(y,sW){
k          = length(y)
w = sW^2
sum.w      = sum(w)
mu.hat     = sum(y*w)/sum.w
W          = sum.w - (sum(w^2)/sum.w)
Q          = sum(w*(y-mu.hat)^2)
typS       = sum(w*(k-1))/(sum.w^2 - sum(w^2))
Qp         = 1-pchisq(Q,k-1)

tausq.hat  = max(0,(Q-(k-1))/W)
w.tauhat   = 1/((1/w)+tausq.hat)

return(tausq.hat)

}



eff <- c(0.1, 0.2, 0.21, 0.5, 0.3)
se <- c(0.02, 0.3, 0.03, 0.07, 0.1)

a <- rma(yi=eff, sei=se, method="FE")
b <- rma(yi=eff, sei=se, method="DL")




DL = function(y=logHR, s = se.logHR)
{
	k          = length(y)
	w          = 1/s^2
	sum.w      = sum(w)
	mu.hat     = sum(y*w)/sum.w
	W          = sum.w - (sum(w^2)/sum.w)
	Q          = sum(w*(y-mu.hat)^2)
	typS       = sum(w*(k-1))/(sum.w^2 - sum(w^2))
	Qp         = 1-pchisq(Q,k-1)

	# fixed effects estimate

	Var.muhat = 1/sum(w)
	muFE.CI   = mu.hat + c(-1.96,1.96)*sqrt(Var.muhat)
	Zfe       = (mu.hat)/sqrt(Var.muhat)
	FEp       = 2*(1 - pnorm(abs(Zfe)))

	# DL random effects estimate

	tausq.hat  = max(0,(Q-(k-1))/W)
	w.tauhat   = 1/(s^2+tausq.hat)
	mu.tau.hat = sum(w.tauhat*y)/sum(w.tauhat)
	Var.mutauhat = 1/sum(w.tauhat)
	muRE.CI    = mu.tau.hat + c(-1.96,1.96)*sqrt(Var.mutauhat)
	Isq        = tausq.hat/(tausq.hat+typS)
	Zre        = (mu.tau.hat)/sqrt(Var.mutauhat)
	REp        = 2*(1 - pnorm(abs(Zre)))

	sum.wre = sum(w.tauhat)

	summary = c(k,Q,Qp,(100*Isq),mu.hat,muFE.CI,
	            FEp,mu.tau.hat,muRE.CI,REp)

	return(summary)
}

DL(eff, se)




