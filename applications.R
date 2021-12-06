##### Power logit regression for modeling bounded data
##### Queiroz, F.F and Ferrari, S.L.P.

require("PLreg")
require("betareg")
require("EnvStats")

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

fda_hat = pbeta(y,a_hat,b_hat )
residuos = qnorm(fda_hat)

resRid <- residuos
e1Rid <- sort(e1)
med1Rid <- sort(med)
e2Rid <- sort(e2)

faixaRid <- c(-3.5,3.5)
qqnorm(resRid, pch="+", las=1, ylim=faixaRid, xlab="quantis te?ricos", ylab="res?duo quant?lico", main="")
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



### Application 3: Body fat of little brown bats

data("bodyfat_Aeolus")
help("bodyfat_Aeolus", package = "PLreg")

# Plots - Figure 8
boxplot(bodyfat_Aeolus$percentfat ~ bodyfat_Aeolus$sex, ylab = "y", xlab = "sex", las = 1, pch = 16,
        names=c("0","1"))
boxplot(bodyfat_Aeolus$percentfat ~ bodyfat_Aeolus$year, ylab = "y", xlab = "year", las = 1, pch = 16,
        names=c("0","1"))
plot(bodyfat_Aeolus$days, bodyfat_Aeolus$percentfat, pch = 16, ylab = "y", xlab = "days", las = 1)
identify(bodyfat_Aeolus$days, bodyfat_Aeolus$percentfat, cex = 0.8)
summary(bodyfat_Aeolus$percentfat[bodyfat_Aeolus$year == 2016])
summary(bodyfat_Aeolus$percentfat[bodyfat_Aeolus$year == 2009])

## Fitting PL-N, PL-PE, PL-Hyp, and PL-SN models

PL_NO <- PLreg(percentfat ~ days + sex + year | days + sex + year, data = bodyfat_Aeolus,
               family = "NO")
summary(PL_NO)

# PL-PE

#Initial model with zeta = 2
PL_PE <- PLreg(percentfat ~ days + sex + year | days + sex + year, data = bodyfat_Aeolus,
               family = "PE", zeta = 2)
# Choosing the best value for zeta
extra.parameter(PL_PE, lower = 1, upper = 3, grid = 10) # zeta = 1.6

PL_PE <- PLreg(percentfat ~ days + sex + year |  days + sex + year , data = bodyfat_Aeolus,
               family = "PE", zeta = 1.6)
summary(PL_PE)


# PL-Hyp

#Initial model with zeta = 2
PL_Hyp <- PLreg(percentfat ~ days + sex + year | days + sex + year, data = bodyfat_Aeolus,
                family = "Hyp", zeta = 2)
# Choosing the best value for zeta
extra.parameter(PL_Hyp, lower = 1, upper = 10, grid = 15) # zeta = 5.5

PL_Hyp <- PLreg(percentfat ~ days + sex + year | days + sex + year , data = bodyfat_Aeolus,
                family = "Hyp", zeta = 5.5)
summary(PL_Hyp)


# PL-SN

#Initial model with zeta = 2
PL_SN <- PLreg(percentfat ~ days + sex + year | days + sex + year, data = bodyfat_Aeolus,
               family = "SN", zeta = 2)
extra.parameter(PL_SN, lower = 0.1, upper = 2, grid = 10) # zeta = 0.52

PL_SN <- PLreg(percentfat ~ days + sex + year |  days + sex + year, data = bodyfat_Aeolus,
               family = "SN", zeta = 0.52)
summary(PL_SN)


## Fitting loglog-N, loglog-PE, loglog-Hyp, and loglog-SN models

loglog_NO <- PLreg(percentfat ~ days + sex + year | days + sex + year, data = bodyfat_Aeolus,
               family = "NO", control = PLreg.control(lambda = 0))
summary(loglog_NO)

# loglog-PE

#Initial model with zeta = 2
loglog_PE <- PLreg(percentfat ~ days + sex + year | days + sex + year, data = bodyfat_Aeolus,
               family = "PE", zeta = 2, control = PLreg.control(lambda = 0))
# Choosing the best value for zeta
extra.parameter(loglog_PE, lower = 1, upper = 3, grid = 10) # zeta = 1.6

loglog_PE <- PLreg(percentfat ~ days + sex + year |  days + sex + year , data = bodyfat_Aeolus,
               family = "PE", zeta = 1.6, control = PLreg.control(lambda = 0))
summary(loglog_PE)


# loglog-Hyp

#Initial model with zeta = 2
loglog_Hyp <- PLreg(percentfat ~ days + sex + year | days + sex + year, data = bodyfat_Aeolus,
                family = "Hyp", zeta = 2, control = PLreg.control(lambda = 0))
# Choosing the best value for zeta
extra.parameter(loglog_Hyp, lower = 1, upper = 10, grid = 15) # zeta = 5.5

loglog_Hyp <- PLreg(percentfat ~ days + sex + year | days + sex + year , data = bodyfat_Aeolus,
                family = "Hyp", zeta = 5.5, control = PLreg.control(lambda = 0))
summary(loglog_Hyp)


# loglog-SN

#Initial model with zeta = 2
loglog_SN <- PLreg(percentfat ~ days + sex + year | days + sex + year, data = bodyfat_Aeolus,
               family = "SN", zeta = 2, control = PLreg.control(lambda = 0))
extra.parameter(loglog_SN, lower = 0.1, upper = 2, grid = 10) # zeta = 0.52

loglog_SN <- PLreg(percentfat ~ days + sex + year |  days + sex + year, data = bodyfat_Aeolus,
               family = "SN", zeta = 0.52, control = PLreg.control(lambda = 0))
summary(loglog_SN)


# Fitting the loglog-NO model


#Initial model
loglog_NO <- PLreg(percentfat ~ days  + year | sex , data = bodyfat_Aeolus,
               family = "NO", control = PLreg.control(lambda= 0 ))
summary(loglog_NO)


plot(loglog_NO)
envelope(loglog_NO, type = "standardized")
influence(loglog_NO)


# Plots - Figure 9
resid = residuals(loglog_NO, type = "standardized")
plot(resid, xlab = "index", ylab = "standardized residual", pch ="+", ylim = c(-3,3), las = 1)
abline(h = -2.5, lty = 2)
abline(h = 2.5, lty = 2)
abline(h = 0, col = "gray")

envelope(loglog_NO, type = "standardized",ylab = "standardized residual")

influence_measures = influence(loglog_NO, graph = FALSE)
plot(influence_measures$case.weights, xlab = "index", las = 1, type = "h", ylim = c(0,0.8), ylab = expression( group("|", h[max], "|") ))
identify(1:159, influence_measures$case.weights, cex = 0.9)

plot(influence_measures$GL, xlab = "index", las = 1, type = "h", ylim = c(0,0.1), ylab = expression(GL[ii]))

