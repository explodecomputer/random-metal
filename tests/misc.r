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
	k            = length(y)
	w            = 1/s^2
	sum.w        = sum(w)
	mu.hat       = sum(y*w)/sum.w
	W            = sum.w - (sum(w^2)/sum.w)
	Q            = sum(w*(y-mu.hat)^2)
	typS         = sum(w*(k-1))/(sum.w^2 - sum(w^2))
	Qp           = 1-pchisq(Q,k-1)

	# fixed effects estimate

	Var.muhat    = 1/sum(w)
	muFE.CI      = mu.hat + c(-1.96,1.96)*sqrt(Var.muhat)
	Zfe          = (mu.hat)/sqrt(Var.muhat)
	FEp          = 2*(1 - pnorm(abs(Zfe)))

	# DL random effects estimate

	tausq.hat    = max(0,(Q-(k-1))/W)
	w.tauhat     = 1/(s^2+tausq.hat)
	mu.tau.hat   = sum(w.tauhat*y)/sum(w.tauhat)
	Var.mutauhat = 1/sum(w.tauhat)
	muRE.CI      = mu.tau.hat + c(-1.96,1.96)*sqrt(Var.mutauhat)
	Isq          = tausq.hat/(tausq.hat+typS)
	Zre          = (mu.tau.hat)/sqrt(Var.mutauhat)
	REp          = 2*(1 - pnorm(abs(Zre)))

	sum.wre      = sum(w.tauhat)

	summary      = c(k,Q,Qp,(100*Isq),mu.hat,muFE.CI,
	                 FEp,mu.tau.hat,muRE.CI,REp)

	out = list(
		fe=mu.hat,
		fe_se = sqrt(Var.muhat),
		re=mu.tau.hat,
		re_se = sqrt(Var.mutauhat)
	)

	return(out)
}

DL(eff, se)




rs1100405       t       c       -0.0178 0.0132  0.177   +--     0.0     0.388   2       0.8237  0.01514304      0.02311397      0.9879



CHR	POS	SNP	N	EFFECT_ALLELE	NON_EFFECT_ALLELE	EFFECT_ALLELE_FREQ	BETA	SE	r2	r2hat	P_VAL
7	43517495	rs1100405	1467	4	2	0.342195	0.00101	0.03945	4.47e-07	0.996	0.9798

CHR POS SNP STRAND EFFECT_ALLELE NON_EFFECT_ALLELE FREQ_EFFECT N BETA SE r2 PVALUE r2hat
7 43517495 rs1100405 + C T 0.665 1233 0.023 0.016 0.00152 0.1674 0.9958

SNP	CHR	POS	Rsq	AL1	AL2	FREQ1	TRAIT	EFFECT	SE	H2	LOD	PVALUE
rs1100405	7	43517495	0.9986	C	T	0.666	glucose_ND	0.011	0.029	0.008	0.032	0.7029



library(metafor)

rma(y=c(-0.00101, 0.023, 0.011), sei=c(0.03945, 0.016, 0.029), method="FE")
rma(y=c(-0.00101, 0.023, 0.011), sei=c(0.03945, 0.016, 0.029), method="DL")

DL(c(-0.00101, 0.023, 0.011), c(0.03945, 0.016, 0.029))







rs2433681       a       g       -0.0149 0.0225  0.5072  -+-     68.2    6.297   2       0.04292 -0.03706755     0.04442108      0.9704



CHR	POS	SNP	N	EFFECT_ALLELE	NON_EFFECT_ALLELE	EFFECT_ALLELE_FREQ	BETA	SE	r2	r2hat	P_VAL
2	169498916	rs2433681	1467	1	3	0.0978187	-0.1559	0.06133	0.004391	0.946	0.0118


CHR POS SNP STRAND EFFECT_ALLELE NON_EFFECT_ALLELE FREQ_EFFECT N BETA SE r2 PVALUE r2hat
2 169498916 rs2433681 + G A 0.926 1233 -0.014 0.029 0.00018 0.6296 0.9918


SNP	CHR	POS	Rsq	AL1	AL2	FREQ1	TRAIT	EFFECT	SE	H2	LOD	PVALUE
rs2433681	2	169498916	0.967	G	A	0.894	glucose_ND	0.009	0.044	0.002	0.008	0.8466


rma(y=c(-0.1559, 0.014, -0.009), sei=c(0.06133, 0.029, 0.044), method="FE")
rma(y=c(-0.1559, 0.014, -0.009), sei=c(0.06133, 0.029, 0.044), method="DL")

