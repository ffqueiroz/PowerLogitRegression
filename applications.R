##### Power logit regression for modeling bounded data
##### Queiroz, F.F and Ferrari, S.L.P.

require("PLreg")
require("betareg")
require("EnvStats")
require("FlexReg")
require("gamlss")

### Application 1: Employment in non-agricultural sectors
y = read.table("agro.txt")$x

# fit beta
fit.beta <- betareg(y~1)
summary(fit.beta)

a_hat = fit.beta$coefficients$precision*( exp(fit.beta$coefficients$mean)/(1+exp(fit.beta$coefficients$mean))  )
b_hat = fit.beta$coefficients$precision - a_hat

fda_hat = pbeta(y[order(y)],a_hat,b_hat )
residuos = qnorm(fda_hat)
Upsilom_beta = (length(y)^(-1))*sum(abs(qnorm(fda_hat)-evNormOrdStats(n = length(y))))

# Beta envelope
fda_hat = pbeta(y[order(y)],a_hat,b_hat )
residuos = qnorm(fda_hat)
rep=100
residuos_env=matrix(rep(0,length(y)*rep),ncol=rep)
for(i in 1:rep){
  y_env=rbeta(length(y), shape1 = a_hat, shape2 = b_hat)
  fit_env = betareg(y_env~1)

  a_hat.env = fit_env$coefficients$precision*( exp(fit_env$coefficients$mean)/(1+exp(fit_env$coefficients$mean))  )
  b_hat.env = fit_env$coefficients$precision - a_hat.env

  fda_hat.env = pbeta(y_env,a_hat.env,b_hat.env )

  residuos_env[,i]=sort(qnorm(fda_hat.env))
  print(i)
}

resid_env=residuos_env
e1=numeric(length(resid_env[,1]))
e2=numeric(length(resid_env[,1]))
med=numeric(length(resid_env[,1]))

for(i in 1:length(y)){
  eo=resid_env[i,]
  e1[i]=quantile(eo,0.025, na.rm = TRUE)
  e2[i]=quantile(eo,0.975, na.rm = TRUE)
  med[i]=median(eo, na.rm = TRUE)
}

fda_hat = pbeta(y, a_hat, b_hat )
residuos = qnorm(fda_hat)

resRid <- residuos
e1Rid <- sort(e1)
med1Rid <- sort(med)
e2Rid <- sort(e2)

faixaRid <- c(-3.5,3.5)
qqnorm(resRid, pch="+", las=1, ylim=faixaRid, main="")
qqnormInt <- function(y,IDENTIFY = TRUE){
  qqnorm(y,pch="+", las=1, ylim=faixaRid, xlab="Quantile N(0,1)", ylab="quantile residual", main="", cex = 0.8) -> X
  if(IDENTIFY) return(identify(X,cex=0.8))
  invisible(X)
}
qqnormInt(resRid)

par(new=T)
qqnorm(e1Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(e2Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(med1Rid, axes=F, lty=2, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)


# fit GJS
fit.GJSN <- PLreg(y ~ 1, family = "NO", control = PLreg.control(lambda = 1))
summary(fit.GJSN)
fitted.values(fit.GJSN)
envelope(fit.GJSN, type = "quantile", ylim=c(-3.5,3.5), ylab = "quantile residual")

fit.GJSPE <- PLreg(y ~ 1, family = "PE", zeta = 1.44,  control = PLreg.control(lambda = 1))
#extra.parameter(fit.GJSPE, 0.2, 3)
summary(fit.GJSPE)
fitted.values(fit.GJSPE)
envelope(fit.GJSPE, type = "quantile", ylim=c(-3.5,3.5), ylab = "quantile residual")

fit.GJST <- PLreg(y ~ 1, family = "TF", zeta = 4.56,  control = PLreg.control(lambda = 1))
#extra.parameter(fit.GJST, 0.2, 10)
summary(fit.GJST)
fitted.values(fit.GJST)

fit.GJSSN <- PLreg(y ~ 1, family = "SN", zeta = 0.53,  control = PLreg.control(lambda = 1))
#extra.parameter(fit.GJSSN, 0.1, 4)
summary(fit.GJSSN)
fitted.values(fit.GJSSN)


# fit PL

fit.PLN <- PLreg(y ~ 1, family = "NO")
summary(fit.PLN)
fitted.values(fit.PLN)
envelope(fit.PLN, type = "quantile", ylim=c(-3.5,3.5), ylab = "quantile residual")

fit.PLPE <- PLreg(y ~ 1, family = "PE", zeta = 2.69)
#extra.parameter(fit.PLPE, 0.2, 3)
summary(fit.PLPE)
fitted.values(fit.PLPE)
envelope(fit.PLPE, type = "quantile", ylim=c(-3.5,3.5), ylab = "quantile residual")

fit.PLT <- PLreg(y ~ 1, family = "T", zeta = 100)
#extra.parameter(fit.PLT, 0.2, 15)
summary(fit.PLT)
fitted.values(fit.PLT)

fit.PLSN <- PLreg(y ~ 1, family = "SN", zeta = 0.97)
#extra.parameter(fit.PLSN, 0.1, 4)
summary(fit.PLSN)
fitted.values(fit.PLSN)


### Tabela

# estimativas
est = rbind(c(fit.beta$fitted.values[1], fit.beta$coefficients$precision, 0),
            c(fit.GJSN$fitted.values[1], exp(fit.GJSN$coefficients$dispersion), 0),
            c(fit.GJST$fitted.values[1], exp(fit.GJST$coefficients$dispersion), 0),
            c(fit.GJSPE$fitted.values[1], exp(fit.GJSPE$coefficients$dispersion), 0),
            c(fit.GJSSN$fitted.values[1], exp(fit.GJSSN$coefficients$dispersion), 0),
            c(fit.PLN$fitted.values[1], exp(fit.PLN$coefficients$dispersion), fit.PLN$coefficients$skewness),
            c(fit.PLT$fitted.values[1], exp(fit.PLT$coefficients$dispersion), fit.PLT$coefficients$skewness),
            c(fit.PLPE$fitted.values[1], exp(fit.PLPE$coefficients$dispersion), fit.PLPE$coefficients$skewness),
            c(fit.PLSN$fitted.values[1], exp(fit.PLSN$coefficients$dispersion), fit.PLSN$coefficients$skewness))

# Plots - Figure 5

hist(y, nclass=15, main="", las=1, ylab="Density", prob=TRUE, col = "white", ylim = c(0, 6), xlim = c(0.3, 1))
curve(dbeta(x, a_hat, b_hat), add = TRUE, lty = 2)
curve(dPL(x, est[6, 1], est[6, 2], est[6, 3], family = "NO"), 0 , 1, add = TRUE, lwd = 2, col= "blue")
curve(dPL(x, est[2, 1], est[2, 2], 1, family = "NO"), 0 , 1, add = TRUE, lty =3)

parametros=c(expression(PL-N),
             expression("beta"),
             expression(GJS-N))
legend(x=0.3,y=6,
       legend=parametros,
       col=c("blue", "black","black"),lty=c(1,2,3), lwd = c(2,1,1),cex=0.9,box.col = "white", bty="n")


plot(ecdf(y),verticals = TRUE,  lwd = 2, col = "gray", main = "", ylab="Cumulative density function",
     xlab = "y", las=1, xlim = c(0.3, 1.04))
stripchart(y, add = TRUE, at = 0, col = "black", pch=3)

curve(pbeta(x, a_hat, b_hat), add = TRUE, lty = 2)
curve(pPL(x, est[6, 1], est[6, 2], est[6, 3], family = "NO"), 0 , 1, add = TRUE, lwd = 2, col= "blue")
curve(pPL(x, est[2, 1], est[2, 2], 1, family = "NO"), 0 , 1, add = TRUE, lty =3)

parametros=c(expression(PL-N),
             expression("beta"),
             expression(GJS-N))
legend(x=0.3,y=0.95,
       legend=parametros,
       col=c("blue", "black","black"),lty=c(1,2,3), lwd = c(2,1,1),cex=0.9,box.col = "white", bty="n")

# Quantile relative discrepancies

p = seq(0.1, 0.99, 0.005)

quantis_empiricos = quantile(y, p)

q_PLNO  = qPL(p, est[6, 1], est[6, 2], est[6, 3], family = "NO")
q_PLPE  = qPL(p, est[8, 1], est[8, 2], est[8, 3], family = "PE", zeta = 2.69 )
q_beta  = qbeta(p, a_hat, b_hat)
q_GJSPE = qPL(p, est[4, 1], est[4, 2], 1, family = "PE", zeta = 1.44)
q_GJSSN = qPL(p, est[5, 1], est[5, 2], 1, family = "SN", zeta = 0.53 )
q_GJSN = qPL(p, est[2, 1], est[2, 2], 1, family = "NO")

D_PLNO  = (q_PLNO   -  quantis_empiricos)/quantis_empiricos
D_PLPE  = (q_PLPE   - quantis_empiricos)/quantis_empiricos
D_beta  = (q_beta   - quantis_empiricos)/quantis_empiricos
D_GJSPE = (q_GJSPE  - quantis_empiricos)/quantis_empiricos
D_GJSSN = (q_GJSSN  - quantis_empiricos)/quantis_empiricos
D_GJSN = (q_GJSN  - quantis_empiricos)/quantis_empiricos

plot(quantis_empiricos,D_PLNO, type = "l", col ="blue", lwd=2, ylim = c(-0.2,0.2), ylab = "Relative quantile discrepancy",
     xlab = "Empirical quantiles", las=1)
lines(quantis_empiricos,D_beta, type = "l", lty = 2)
lines(quantis_empiricos,D_GJSN, type = "l", lty =3)
abline(h = 0, lty = 2, col="grey")

parametros=c(expression(PL-N),
             expression("beta"),
             expression(GJS-N))
legend(x=0.61,y=0.2,
       legend=parametros,
       col=c("blue", "black","black"),lty=c(1,2,3), lwd = c(2,1,1),cex=0.9,box.col = "white", bty="n")

# Extra distributions
# fit logitSHASH 
gen.Family("SHASH", type = "logit")

fit.SHASH <- gamlss(y~1,family="logitSHASH", 
                    control = gamlss.control(n.cyc = 200, c.crit = 0.005))
AIC(fit.SHASH)

mu_hat.SHASH    <- fit.SHASH$mu.coefficients
sigma_hat.SHASH <- exp(fit.SHASH$sigma.coefficients)
nu_hat.SHASH    <- exp(fit.SHASH$nu.coefficients)
tau_hat.SHASH   <- exp(fit.SHASH$tau.coefficients)

hist(y, nclass=15, main="", las=1, ylab="Density", prob=TRUE, col = "white", ylim = c(0, 6), xlim = c(0.3, 1))
curve(dlogitSHASH(x, mu = mu_hat.SHASH, sigma = sigma_hat.SHASH, nu = nu_hat.SHASH,
                  tau = tau_hat.SHASH), add = TRUE, lty = 2)

fda_hatSHASH = plogitSHASH(y[order(y)], mu = mu_hat.SHASH, sigma = sigma_hat.SHASH, nu = nu_hat.SHASH,
                           tau = tau_hat.SHASH)
residuos = qnorm(fda_hatSHASH)
Upsilom.SHASH = (length(y)^(-1))*sum(abs(qnorm(fda_hatSHASH)-evNormOrdStats(n = length(y))))

# logitSHASH envelope
residuos = qnorm(fda_hatSHASH)
rep=100
residuos_env=matrix(rep(0,length(y)*rep),ncol=rep)
for(i in 1:rep){
  y_env=rlogitSHASH(length(y), mu = mu_hat.SHASH, sigma = sigma_hat.SHASH, nu = nu_hat.SHASH,
                    tau = tau_hat.SHASH)
  fit_env = gamlss(y_env ~ 1,family="logitSHASH", 
                   control = gamlss.control(n.cyc = 200, c.crit = 0.005))
  
  mu_hat_env    <- fit_env$mu.coefficients
  sigma_hat_env <- exp(fit_env$sigma.coefficients)
  nu_hat_env    <- exp(fit_env$nu.coefficients)
  tau_hat_env   <- exp(fit_env$tau.coefficients)
  
  fda_hat.env = plogitSHASH(y_env, mu = mu_hat_env, sigma = sigma_hat_env, 
                            nu = nu_hat_env, tau = tau_hat_env)
  residuos_env[,i]=sort(qnorm(fda_hat.env))
  print(i)
}

resid_env=residuos_env
e1=numeric(length(resid_env[,1]))
e2=numeric(length(resid_env[,1]))
med=numeric(length(resid_env[,1]))

for(i in 1:length(y)){
  eo=resid_env[i,]
  e1[i]=quantile(eo,0.025, na.rm = TRUE)
  e2[i]=quantile(eo,0.975, na.rm = TRUE)
  med[i]=median(eo, na.rm = TRUE)
}



resRid <- residuos
e1Rid <- sort(e1)
med1Rid <- sort(med)
e2Rid <- sort(e2)

faixaRid <- c(-3.5,3.5)
qqnorm(resRid, pch="+", las=1, ylim=faixaRid, main="")
qqnormInt <- function(y,IDENTIFY = TRUE){
  qqnorm(y,pch="+", las=1, ylim=faixaRid, xlab="Quantile N(0,1)", ylab="quantile residual", main="", cex = 0.8) -> X
  if(IDENTIFY) return(identify(X,cex=0.8))
  invisible(X)
}
qqnormInt(resRid)

par(new=T)
qqnorm(e1Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(e2Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(med1Rid, axes=F, lty=2, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)


# fit GB1 

fit.GB1 <- gamlss(y~1,family="GB1", 
                  control = gamlss.control(n.cyc = 200, c.crit = 0.005))
AIC(fit.GB1)

mu_hat.GB1    <- exp(fit.GB1$mu.coefficients)/(1 + exp(fit.GB1$mu.coefficients))
sigma_hat.GB1 <- exp(fit.GB1$sigma.coefficients)/(1 + exp(fit.GB1$sigma.coefficients))
nu_hat.GB1    <- exp(fit.GB1$nu.coefficients)
tau_hat.GB1   <- exp(fit.GB1$tau.coefficients)

hist(y, nclass=15, main="", las=1, ylab="Density", prob=TRUE, col = "white", ylim = c(0, 6), xlim = c(0.3, 1))
curve(dGB1(x, mu = mu_hat.GB1, sigma = sigma_hat.GB1, nu = nu_hat.GB1,
           tau = tau_hat.GB1), 0.01, 0.99, add = TRUE, lty = 2)


fda_hatGB1 = pGB1(y[order(y)], mu = mu_hat.GB1, sigma = sigma_hat.GB1, nu = nu_hat.GB1,
                  tau = tau_hat.GB1)
residuos = qnorm(fda_hatGB1)
Upsilom_GB1 = (length(y)^(-1))*sum(abs(qnorm(fda_hatGB1)-evNormOrdStats(n = length(y))))

# GB1 envelope
rep=100
residuos_env=matrix(rep(0,length(y)*rep),ncol=rep)
for(i in 1:rep){
  y_env=rGB1(length(y), mu = mu_hat.GB1, sigma = sigma_hat.GB1, nu = nu_hat.GB1,
             tau = tau_hat.GB1)
  fit_env = gamlss(y_env ~ 1,family="GB1", 
                   control = gamlss.control(n.cyc = 200, c.crit = 0.005))
  
  mu_hat_env    <- exp(fit_env$mu.coefficients)/(1 + exp(fit_env$mu.coefficients))
  sigma_hat_env <- exp(fit_env$sigma.coefficients)/(1 + exp(fit_env$sigma.coefficients))
  nu_hat_env    <- exp(fit_env$nu.coefficients)
  tau_hat_env   <- exp(fit_env$tau.coefficients)
  
  fda_hat.env = pGB1(y_env, mu = mu_hat_env, sigma = sigma_hat_env, 
                     nu = nu_hat_env, tau = tau_hat_env)
  residuos_env[,i]=sort(qnorm(fda_hat.env))
  print(i)
}

resid_env=residuos_env
e1=numeric(length(resid_env[,1]))
e2=numeric(length(resid_env[,1]))
med=numeric(length(resid_env[,1]))

for(i in 1:length(y)){
  eo=resid_env[i,]
  e1[i]=quantile(eo,0.025, na.rm = TRUE)
  e2[i]=quantile(eo,0.975, na.rm = TRUE)
  med[i]=median(eo, na.rm = TRUE)
}

resRid <- residuos
e1Rid <- sort(e1)
med1Rid <- sort(med)
e2Rid <- sort(e2)

faixaRid <- c(-3.5,3.5)
qqnorm(resRid, pch="+", las=1, ylim=faixaRid, main="")
qqnormInt <- function(y,IDENTIFY = TRUE){
  qqnorm(y,pch="+", las=1, ylim=faixaRid, xlab="Quantile N(0,1)", ylab="quantile residual", main="", cex = 0.8) -> X
  if(IDENTIFY) return(identify(X,cex=0.8))
  invisible(X)
}
qqnormInt(resRid)

par(new=T)
qqnorm(e1Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(e2Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(med1Rid, axes=F, lty=2, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)


# fit FB distribution 

loglik <- function(par, y){
  mu1 <- exp(par[1])/(1 + exp(par[1]))
  mu2 <- exp(par[2])/(1 + exp(par[2]))
  phi <- exp(par[3])
  p <- exp(par[4])/(1 + exp(par[4]))
  
  f <- p*dbeta(y, mu1*phi, (1-mu1)*phi) + (1 - p)*dbeta(y, mu2*phi, (1-mu2)*phi)
  sum(log(f))
}

# starting values
FB <- summary(flexreg(y ~ 1, type="FB", n.iter=1000))

mu_hat <- exp(FB$Summary.mu[1])/(1+exp(FB$Summary.mu[1]))
phi_hat <- FB$Summary.phi[1]
w_hat <- FB$Summary.add[2,1]
p_hat <- FB$Summary.add[1,1]

mu1.start <- mu_hat + (1 - p_hat)*(w_hat*min(mu_hat/p_hat, (1-mu_hat)/(1-p_hat)))
mu2.start <- mu_hat - p_hat*(w_hat*min(mu_hat/p_hat, (1-mu_hat)/(1-p_hat)))

start <- c(log(mu1.start/(1 - mu1.start)), log(mu2.start/(1 - mu2.start)), 
           log(phi_hat), log(p_hat/(1 - p_hat)) )
fit.FB <- optim(start, loglik, y = y, control = list(fnscale = -1), hessian = TRUE)

mu_hat1.FB <- exp(fit.FB$par[1])/(1+exp(fit.FB$par[1]))
mu_hat2.FB <- exp(fit.FB$par[2])/(1+exp(fit.FB$par[2]))
phi_hat.FB <- exp(fit.FB$par[3])
p_hat.FB <- exp(fit.FB$par[4])/(1+exp(fit.FB$par[4]))

mu_hat <- mu_hat1.FB*p_hat.FB + (1-p_hat.FB)*mu_hat2.FB
w <- (mu_hat1.FB - mu_hat2.FB)/min(mu_hat/p_hat.FB, (1-mu_hat)/(1-p_hat.FB))

hist(y, nclass=15, main="", las=1, ylab="Density", prob=TRUE, col = "white", ylim = c(0, 6), xlim = c(0.3, 1))
curve(p_hat.FB*dbeta(x, mu_hat1.FB*phi_hat.FB, (1-mu_hat1.FB)*phi_hat.FB) + 
        (1 - p_hat.FB)*dbeta(x, mu_hat2.FB*phi_hat.FB, (1-mu_hat2.FB)*phi_hat.FB),
      0.01, 1, add = TRUE, lty = 2)

fda_hat.FB = p_hat.FB*pbeta(y[order(y)], mu_hat1.FB*phi_hat.FB, (1-mu_hat1.FB)*phi_hat.FB) + 
  (1 - p_hat.FB)*pbeta(y[order(y)], mu_hat2.FB*phi_hat.FB, (1-mu_hat2.FB)*phi_hat.FB)
residuos = qnorm(fda_hat.FB)

Upsilom.FB = (length(y)^(-1))*sum(abs(qnorm(fda_hat.FB)-evNormOrdStats(n = length(y))))

# FB envelope
rep=100
residuos_env=matrix(rep(0,length(y)*rep),ncol=rep)
for(i in 1:rep){
  u <- rbinom(length(y), 1, p_hat.FB)
  beta1 <- rbeta(length(y), mu_hat1.FB*phi_hat.FB, (1-mu_hat1.FB)*phi_hat.FB)
  beta2 <- rbeta(length(y), mu_hat2.FB*phi_hat.FB, (1-mu_hat2.FB)*phi_hat.FB)
  y_env <- ifelse(u == 1, beta1, beta2)
  
  fit.FB_env <- optim(start, loglik, y = y_env, control = list(fnscale = -1))
  
  mu_hat1.FB_env <- exp(fit.FB_env$par[1])/(1+exp(fit.FB_env$par[1]))
  mu_hat2.FB_env <- exp(fit.FB_env$par[2])/(1+exp(fit.FB_env$par[2]))
  phi_hat.FB_env <- exp(fit.FB_env$par[3])
  p_hat.FB_env <- exp(fit.FB_env$par[4])/(1+exp(fit.FB_env$par[4]))
  
  fda_hat.env = p_hat.FB_env*pbeta(y_env, mu_hat1.FB_env*phi_hat.FB_env, (1-mu_hat1.FB_env)*phi_hat.FB_env) + 
    (1 - p_hat.FB_env)*pbeta(y_env, mu_hat2.FB_env*phi_hat.FB_env, (1-mu_hat2.FB_env)*phi_hat.FB_env)
  residuos_env[,i]=sort(qnorm(fda_hat.env))
  print(i)
}

resid_env=residuos_env
e1=numeric(length(resid_env[,1]))
e2=numeric(length(resid_env[,1]))
med=numeric(length(resid_env[,1]))

for(i in 1:length(y)){
  eo=resid_env[i,]
  e1[i]=quantile(eo,0.025, na.rm = TRUE)
  e2[i]=quantile(eo,0.975, na.rm = TRUE)
  med[i]=median(eo, na.rm = TRUE)
}

resRid <- residuos
e1Rid <- sort(e1)
med1Rid <- sort(med)
e2Rid <- sort(e2)

faixaRid <- c(-3.5,3.5)
qqnorm(resRid, pch="+", las=1, ylim=faixaRid, main="")
qqnormInt <- function(y,IDENTIFY = TRUE){
  qqnorm(y,pch="+", las=1, ylim=faixaRid, xlab="Quantile N(0,1)", ylab="quantile residual", main="", cex = 0.8) -> X
  if(IDENTIFY) return(identify(X,cex=0.8))
  invisible(X)
}
qqnormInt(resRid)

par(new=T)
qqnorm(e1Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(e2Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(med1Rid, axes=F, lty=2, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)



# fit logitskewnormal type 2


gen.Family("SN2", type = "logit")

fit.SN2 <- gamlss(y~1,family="logitSN2", 
                  control = gamlss.control(n.cyc = 200, c.crit = 0.005))
AIC(fit.SN2)

mu_hat.SN2    <- fit.SN2$mu.coefficients
sigma_hat.SN2 <- exp(fit.SN2$sigma.coefficients)
nu_hat.SN2    <- exp(fit.SN2$nu.coefficients)

hist(y, nclass=15, main="", las=1, ylab="Density", prob=TRUE, col = "white", ylim = c(0, 6), xlim = c(0.3, 1))
curve(dlogitSN2(x, mu = mu_hat.SN2, sigma = sigma_hat.SN2, nu = nu_hat.SN2), add = TRUE, lty = 2)


fda_hatSN2 = plogitSN2(y[order(y)], mu = mu_hat.SN2, sigma = sigma_hat.SN2, nu = nu_hat.SN2)
residuos = qnorm(fda_hatSN2)
Upsilom.SN2 = (length(y)^(-1))*sum(abs(qnorm(fda_hatSN2) - evNormOrdStats(n = length(y))))


# logitSHASH envelope
residuos = qnorm(fda_hatSN2)
rep=100
residuos_env=matrix(rep(0,length(y)*rep),ncol=rep)
for(i in 1:rep){
  y_env=rlogitSN2(length(y), mu = mu_hat.SN2, sigma = sigma_hat.SN2, nu = nu_hat.SN2)
  fit_env = gamlss(y_env ~ 1,family="logitSN2", 
                   control = gamlss.control(n.cyc = 200, c.crit = 0.005))
  
  mu_hat_env    <- fit_env$mu.coefficients
  sigma_hat_env <- exp(fit_env$sigma.coefficients)
  nu_hat_env    <- exp(fit_env$nu.coefficients)
  
  fda_hat.env = plogitSN2(y_env, mu = mu_hat_env, sigma = sigma_hat_env, 
                          nu = nu_hat_env)
  residuos_env[,i]=sort(qnorm(fda_hat.env))
  print(i)
}

resid_env=residuos_env
e1=numeric(length(resid_env[,1]))
e2=numeric(length(resid_env[,1]))
med=numeric(length(resid_env[,1]))

for(i in 1:length(y)){
  eo=resid_env[i,]
  e1[i]=quantile(eo,0.025, na.rm = TRUE)
  e2[i]=quantile(eo,0.975, na.rm = TRUE)
  med[i]=median(eo, na.rm = TRUE)
}

resRid <- residuos
e1Rid <- sort(e1)
med1Rid <- sort(med)
e2Rid <- sort(e2)

faixaRid <- c(-3.5,3.5)
qqnorm(resRid, pch="+", las=1, ylim=faixaRid, main="")
qqnormInt <- function(y,IDENTIFY = TRUE){
  qqnorm(y,pch="+", las=1, ylim=faixaRid, xlab="Quantile N(0,1)", ylab="quantile residual", main="", cex = 0.8) -> X
  if(IDENTIFY) return(identify(X,cex=0.8))
  invisible(X)
}
qqnormInt(resRid)

par(new=T)
qqnorm(e1Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(e2Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(med1Rid, axes=F, lty=2, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)




# fit logitskewt3 
gen.Family("ST3", type = "logit")

fit.ST3 <- gamlss(y~1,family="logitST3", 
                  control = gamlss.control(n.cyc = 200, c.crit = 0.005))
AIC(fit.ST3)

mu_hat.ST3    <- fit.ST3$mu.coefficients
sigma_hat.ST3 <- exp(fit.ST3$sigma.coefficients)
nu_hat.ST3    <- exp(fit.ST3$nu.coefficients)
tau_hat.ST3    <- exp(fit.ST3$tau.coefficients)

hist(y, nclass=15, main="", las=1, ylab="Density", prob=TRUE, col = "white", ylim = c(0, 6), xlim = c(0.3, 1))
curve(dlogitST3(x, mu = mu_hat.ST3, sigma = sigma_hat.ST3, nu = nu_hat.ST3, tau = tau_hat.ST3), add = TRUE, lty = 2)

fda_hatST3 = plogitST3(y[order(y)], mu = mu_hat.ST3, sigma = sigma_hat.ST3, nu = nu_hat.ST3, tau = tau_hat.ST3)
residuos = qnorm(fda_hatST3)
Upsilom.ST3 = (length(y)^(-1))*sum(abs(qnorm(fda_hatST3) - evNormOrdStats(n = length(y))))


# logitST3 envelope
residuos = qnorm(fda_hatST3)
rep=100
residuos_env=matrix(rep(0,length(y)*rep),ncol=rep)
for(i in 1:rep){
  y_env=rlogitST3(length(y), mu = mu_hat.ST3, sigma = sigma_hat.ST3, nu = nu_hat.ST3, tau = tau_hat.ST3)
  fit_env = gamlss(y_env ~ 1,family="logitST3", 
                   control = gamlss.control(n.cyc = 200, c.crit = 0.005))
  
  mu_hat_env    <- fit_env$mu.coefficients
  sigma_hat_env <- exp(fit_env$sigma.coefficients)
  nu_hat_env    <- exp(fit_env$nu.coefficients)
  tau_hat_env   <- exp(fit_env$tau.coefficients)
  
  fda_hat.env = plogitST3(y_env, mu = mu_hat_env, sigma = sigma_hat_env, 
                          nu = nu_hat_env, tau = tau_hat_env)
  residuos_env[,i]=sort(qnorm(fda_hat.env))
  print(i)
}

resid_env=residuos_env
e1=numeric(length(resid_env[,1]))
e2=numeric(length(resid_env[,1]))
med=numeric(length(resid_env[,1]))

for(i in 1:length(y)){
  eo=resid_env[i,]
  e1[i]=quantile(eo,0.025, na.rm = TRUE)
  e2[i]=quantile(eo,0.975, na.rm = TRUE)
  med[i]=median(eo, na.rm = TRUE)
}



resRid <- residuos
e1Rid <- sort(e1)
med1Rid <- sort(med)
e2Rid <- sort(e2)

faixaRid <- c(-3.5,3.5)
qqnorm(resRid, pch="+", las=1, ylim=faixaRid, main="")
qqnormInt <- function(y,IDENTIFY = TRUE){
  qqnorm(y,pch="+", las=1, ylim=faixaRid, xlab="Quantile N(0,1)", ylab="quantile residual", main="", cex = 0.8) -> X
  if(IDENTIFY) return(identify(X,cex=0.8))
  invisible(X)
}
qqnormInt(resRid)

par(new=T)
qqnorm(e1Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(e2Rid, axes=F, lty=1, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)
par(new=T)
qqnorm(med1Rid, axes=F, lty=2, ylim=faixaRid, type="l", xlab="", ylab="", main="", lwd=1)

# Plots new distributions

hist(y, nclass=15, main="", las=1, ylab="Density", prob=TRUE, col = "white", ylim = c(0, 6), xlim = c(0.3, 1))
curve(dlogitSHASH(x, mu = mu_hat.SHASH, sigma = sigma_hat.SHASH, nu = nu_hat.SHASH,
                  tau = tau_hat.SHASH), add = TRUE, lty = 2)
curve(dGB1(x, mu = mu_hat.GB1, sigma = sigma_hat.GB1, nu = nu_hat.GB1,
           tau = tau_hat.GB1), 0.01, 0.99, add = TRUE, lty = 3)
curve(p_hat.FB*dbeta(x, mu_hat1.FB*phi_hat.FB, (1-mu_hat1.FB)*phi_hat.FB) + 
        (1 - p_hat.FB)*dbeta(x, mu_hat2.FB*phi_hat.FB, (1-mu_hat2.FB)*phi_hat.FB),
      0.01, 1, add = TRUE)
curve(dlogitSN2(x, mu = mu_hat.SN2, sigma = sigma_hat.SN2, nu = nu_hat.SN2), add = TRUE, lty =4)


parametros=c(expression(FB),
             expression(logitSHASH),
             expression(GB1),
             expression(logitSN2))
legend(x=0.3,y=6,
       legend=parametros,
       col=c("black", "black","black", "black"),lty=c(1,2,3, 4), lwd = c(1,1,1, 1),cex=0.9,box.col = "white", bty="n")


plot(ecdf(y),verticals = TRUE,  lwd = 2, col = "gray", main = "", ylab="Cumulative density function",
     xlab = "y", las=1, xlim = c(0.3, 1.04))
stripchart(y, add = TRUE, at = 0, col = "black", pch=3)

curve(plogitSHASH(x, mu = mu_hat.SHASH, sigma = sigma_hat.SHASH, nu = nu_hat.SHASH,
                  tau = tau_hat.SHASH), 0.0001,0.9999, add = TRUE, lty = 2)
curve(p_hat.FB*pbeta(x, mu_hat1.FB*phi_hat.FB, (1-mu_hat1.FB)*phi_hat.FB) + 
        (1 - p_hat.FB)*pbeta(x, mu_hat2.FB*phi_hat.FB, (1-mu_hat2.FB)*phi_hat.FB), 0 , 1, add = TRUE, lwd = 1)
curve(pGB1(x, mu = mu_hat.GB1, sigma = sigma_hat.GB1, nu = nu_hat.GB1,
           tau = tau_hat.GB1), 0.00001 , 0.9999, add = TRUE, lty =3)
curve(plogitSN2(x, mu = mu_hat.SN2, sigma = sigma_hat.SN2, nu = nu_hat.SN2), 0.000001, 0.999999, add = TRUE, lty = 4)

legend(x=0.3,y=0.95,
       legend=parametros,
       col=c("black", "black","black", "black"),lty=c(1,2,3, 4), lwd = c(1,1,1, 1),cex=0.9,box.col = "white", bty="n")



### Application 2: Firm cost data

data("Firm")
help("Firm", package = "PLreg")

# Fitting PL slash model with zeta = 2

fit <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm,
             family = "SLASH", zeta = 2)
extra.parameter(fit, 1,2.5, grid = 30) # zeta = 1.88

# Model with zeta = 1.88
fit <- PLreg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm,
             family = "SLASH", zeta = 1.88)
summary(fit)

# Fitting PL slash with constant dispersion
fit <- PLreg(firmcost ~ sizelog + indcost, data = Firm,
             family = "SLASH", zeta = 2)
extra.parameter(fit, 1,2.5, grid = 30) # zeta = 2.29

fit <- PLreg(firmcost ~ sizelog + indcost, data = Firm,
             family = "SLASH", zeta = 2.29)
summary(fit)


# GJS slash with constant dispersion
fit <- PLreg(firmcost ~ sizelog + indcost, data = Firm,
             family = "SLASH", zeta = 1.52, control = PLreg.control(lambda=1))
summary(fit)


# Plots - Figure 6
resid = residuals(fit, type = "standardized")
plot(resid, xlab = "index", ylab = "standardized residual", pch ="+", ylim = c(-5,7), las = 1)
identify(1:73, resid, cex = 0.9)
abline(h = -2.5, lty = 2)
abline(h = 2.5, lty = 2)
abline(h = 0, col = "gray")

envelope(fit, type = "standardized",ylab = "standardized residual", xlim = c(-2.5, 2.5), ylim = c(-7, 7))

influence_measures = influence(fit, graph = FALSE)
plot(influence_measures$case.weights, xlab = "index", las = 1, type = "h", ylim = c(0,0.8), ylab = expression( group("|", h[max], "|") ))
identify(1:73, influence_measures$case.weights, cex = 0.9)

plot(Firm$indcost, influence_measures$GL, xlab = "indcost", pch = "+",
     las = 1, ylab = expression(GL[ii]), ylim = c(0, 0.5))
identify(Firm$indcost, influence_measures$GL, cex = 0.9)

plot(resid, fit$v, xlab = "standardized residual", pch = "+",
     las = 1, ylab = expression(v(z)))
identify(resid, fit$v, cex = 0.9)


# Fitted lines

# beta model
fitbeta      <- betareg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm)
fitbeta_wo16 <- betareg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm[-16,])
fitbeta_wo15 <- betareg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm[-15,])
fitbeta_wo72 <- betareg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm[-72,])
fitbeta_wo10 <- betareg(firmcost ~ sizelog + indcost | sizelog + indcost, data = Firm[-10,])

coef.mu      <- coef(fitbeta)[1:3]
coef.mu_wo16 <- coef(fitbeta_wo16)[1:3]
coef.mu_wo15 <- coef(fitbeta_wo15)[1:3]
coef.mu_wo72 <- coef(fitbeta_wo72)[1:3]
coef.mu_wo10 <- coef(fitbeta_wo10)[1:3]

plot(Firm$indcost, Firm$firmcost,
     pch = "+",
     xlab = "indcost",
     ylab = "firm cost",
     las = 1)
#identify(Firm$indcost, Firm$firmcost, cex = 0.8)
covariate = matrix(c(rep.int(1, 1000),
                     rep(median(Firm$sizelog), 1000),
                     seq(0, 1.22, length.out = 1000)),
                   ncol = 3)
lines(covariate[,3],
      as.vector(fit$link$median$linkinv(covariate%*%coef.mu)),
      type = "l")
lines(covariate[,3],
      as.vector(fit$link$median$linkinv(covariate%*%coef.mu_wo15)),
      type = "l", lty = 2, col = "blue")
lines(covariate[,3],
      as.vector(fit$link$median$linkinv(covariate%*%coef.mu_wo16)),
      type = "l", lty = 3, col = "red")
lines(covariate[,3],
      as.vector(fit$link$median$linkinv(covariate%*%coef.mu_wo72)),
      type = "l", lty = 4, col = "green")
lines(covariate[,3],
      as.vector(fit$link$median$linkinv(covariate%*%coef.mu_wo10)),
      type = "l", lty = 5, col = "violet")

parameters = c("mle",
               "mle w/o #15",
               "mle w/o #16",
               "mle w/o #72",
               "mle w/o #10")
legend(x = 0.5,
       y = 0.8,
       legend = parameters,
       col = c("black", "blue", "red", "green", "violet"),
       lty = c(1, 2, 3, 4, 5),
       cex = 0.75)

# PL-slash model

fitPL <- PLreg(firmcost ~  sizelog + indcost, data = Firm,
               family = "SLASH", zeta = 2.29)

fitPL_wo16 <- PLreg(firmcost ~ sizelog + indcost,
                    data = Firm[-16,],
                    family = "SLASH",
                    zeta = 2.29)

fitPL_wo15 <- PLreg(firmcost ~ sizelog + indcost,
                    data = Firm[-15,],
                    family = "SLASH",
                    zeta = 2.29)

fitPL_wo72 <- PLreg(firmcost ~ sizelog + indcost,
                    data = Firm[-72,],
                    family = "SLASH",
                    zeta = 2.29)

fitPL_wo10 <- PLreg(firmcost ~ sizelog + indcost,
                    data = Firm[-10,],
                    family = "SLASH",
                    zeta = 2.29)

coef.mu      <- coef(fitPL)[1:3]
coef.mu_wo16 <- coef(fitPL_wo16)[1:3]
coef.mu_wo15 <- coef(fitPL_wo15)[1:3]
coef.mu_wo72 <- coef(fitPL_wo72)[1:3]
coef.mu_wo10 <- coef(fitPL_wo10)[1:3]

plot(Firm$indcost, Firm$firmcost,
     pch = "+",
     xlab = "indcost",
     ylab = "firm cost",
     las = 1)
#identify(Firm$indcost, Firm$firmcost, cex = 0.8)
covariate = matrix(c(rep.int(1, 1000),
                     rep(median(Firm$sizelog), 1000),
                     seq(0, 1.22, length.out = 1000)),
                   ncol = 3)
lines(covariate[,3],
      as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu)),
      type = "l")
lines(covariate[,3],
      as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo15)),
      type = "l", lty = 2, col = "blue")
lines(covariate[,3],
      as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo16)),
      type = "l", lty = 3, col = "red")
lines(covariate[,3],
      as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo72)),
      type = "l", lty = 4, col = "green")
lines(covariate[,3],
      as.vector(fitPL$link$median$linkinv(covariate%*%coef.mu_wo10)),
      type = "l", lty = 5, col = "violet")

parameters = c("pmle",
               "pmle w/o #15",
               "pmle w/o #16",
               "pmle w/o #72",
               "pmle w/o #10")
legend(x = 0.5,
       y = 0.8,
       legend = parameters,
       col = c("black", "blue", "red", "green", "violet"),
       lty = c(1, 2, 3, 4, 5),
       cex = 0.75)

