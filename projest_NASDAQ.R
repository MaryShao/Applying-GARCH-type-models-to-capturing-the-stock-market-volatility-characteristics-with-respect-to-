## install.packages('rugarch')
# @Manual{Ghalanos_2014,
#   author       = {Alexios Ghalanos},
#   title        = {{rugarch}: Univariate GARCH models.},
#   year         = {2014},
#   note         = {R package version 1.3-5.},}

library(rugarch)
## install.packages('quantmod')
library(quantmod)
# install.packages('fBasics')
library(fBasics)
#install.packages('forecast')
library(forecast)
#install.packages('aTSA')
library(aTSA)
#install.packages('PerformanceAnalytics')
library(PerformanceAnalytics) 
library(tseries)
library(stats)

# improt data
# period1 2004.01.01-2006.12.31: ^IXIC.csv
getSymbols("^IXIC", from='2004-01-01', to= '2006-12-31')
rt1 <- dailyReturn(Cl(IXIC), type='log')
rt1 <- 100*rt1
# period2 2007.01.01-2009.12.31: ^IXIC-2.csv
getSymbols("^IXIC", from='2007-01-01', to= '2009-12-31')
rt2 <- dailyReturn(Cl(IXIC), type='log')
rt2 <- 100*rt2
# period3 2010.01.01-2012.12.31: ^IXIC-3.csv
getSymbols("^IXIC", from='2010-01-01', to= '2012-12-31')
rt3 <- dailyReturn(Cl(IXIC), type='log')
rt3 <- 100*rt3

# plot daily nominal return series
plot(rt2,type='l',main='NASDAQ') 
par(new = TRUE)
plot(rt1,type='l',main='NASDAQ',ylim=range(c(rt2,rt1)), axes = FALSE)
par(new = TRUE)
plot(rt3,type='l',main='NASDAQ',ylim=range(c(rt2,rt3)), axes = FALSE)
# plot daily actual volatility--assessed by daily squared returns
plot(rt2^2,type='l',main='NASDAQ') 
par(new = TRUE)
plot(rt1^2,type='l',main='NASDAQ',ylim=range(c(rt2^2,rt1^2)), axes = FALSE)
par(new = TRUE)
plot(rt3^2,type='l',main='NASDAQ',ylim=range(c(rt2^2,rt3^2)), axes = FALSE)

# data analysis
skew1<-skewness(rt1)
skew2<-skewness(rt2, na.rm = FALSE, type = 3)
skew3<-skewness(rt3, na.rm = FALSE, type = 3)
kur1<-mean(((rt1-mean(rt1))^4))/sd(rt1)^4
kur2<-mean(((rt2-mean(rt2))^4))/sd(rt2)^4
kur3<-mean(((rt3-mean(rt3))^4))/sd(rt3)^4
jb1<-jarque.bera.test(rt1)
jb2<-jarque.bera.test(rt2)
jb3<-jarque.bera.test(rt3)
n1<-length(rt1)
n2<-length(rt2)
n3<-length(rt3)
#w1<-Box.test(rt1, lag = 1, type = c("Box-Pierce", "Ljung-Box"), fitdf = 0)


Observations <- c(n1,n2,n3)
Mean <- c(mean(rt1),mean(rt2),mean(rt3)) 
Median <- c(median(rt1),median(rt2),median(rt3)) 
Maximum <- c(max(rt1),max(rt2),max(rt3)) 
Minimum <- c(min(rt1),min(rt2),min(rt3)) 
Std.Dev. <- c(sd(rt1),sd(rt2),sd(rt3)) 
Skewness <- c(skew1,skew2,skew3)
Kurtosis <- c(kur1,kur2,kur3)
JBtest <- c(jb1[1],jb2[1],jb3[1])
names(Observations) <- c('2004-2006','2007-2009','2010-2012')
rbind(Observations,Mean,Median,Maximum,Minimum,Std.Dev.,Skewness,Kurtosis,JBtest)
# H0:non-stationary time series (pvalue<0.01 strong evidence reject)
adf.test(rt1)
adf.test(rt2)
adf.test(rt3)

# ARMA
m1<-auto.arima(rt1,stationary=TRUE)
o1<-arimaorder(m1)
c(o1[1],o1[2],o1[3])
arch.test(arima(rt1,order=c(o1[1],o1[2],o1[3])))
m2<-auto.arima(rt2,stationary=TRUE)
o2<-arimaorder(m2)
c(o2[1],o2[2],o2[3])
arch.test(arima(rt2,order=c(o2[1],o2[2],o2[3])))
m3<-auto.arima(rt3,stationary=TRUE)
o3<-arimaorder(m3)
c(o3[1],o3[2],o3[3])
arch.test(arima(rt1,order=c(o3[1],o3[2],o3[3])))

#GARCH(1,1)
garch11.spec = ugarchspec(variance.model = list(model="sGARCH",  
                                                garchOrder=c(1,1)), 
                          mean.model = list(armaOrder=c(0,0)))
garch11_1.fit = ugarchfit(spec=garch11.spec, data=rt1)
garch11_1.fit 
coef(garch11_1.fit) 
garch11_2.fit = ugarchfit(spec=garch11.spec, data=rt2)
garch11_2.fit 
coef(garch11_2.fit)
garch11_3.fit = ugarchfit(spec=garch11.spec, data=rt3)
garch11_3.fit 
coef(garch11_3.fit)

# EGARCH(1,1)
egarchsnp.spec = ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)), 
                            mean.model=list(armaOrder=c(0,0))) 
egarchsnp_1.fit = ugarchfit(spec=egarchsnp.spec, data=rt1) 
egarchsnp_1.fit
coef(egarchsnp_1.fit) 
egarchsnp_2.fit = ugarchfit(spec=egarchsnp.spec, data=rt2) 
egarchsnp_2.fit
coef(egarchsnp_2.fit) 
egarchsnp_3.fit = ugarchfit(spec=egarchsnp.spec, data=rt3) 
egarchsnp_3.fit
coef(egarchsnp_3.fit) 
#residuals(egarchsnp.fit) 


#GJR-GARCH(1,1)
gjrgarchsnp.spec = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                              mean.model=list(armaOrder=c(2,1), include.mean=TRUE), 
                              distribution.model="std")
gjrgarchsnp_1.fit = ugarchfit(spec=gjrgarchsnp.spec, data=rt1) 
gjrgarchsnp_1.fit
coef(gjrgarchsnp_1.fit) 
gjrgarchsnp_2.fit = ugarchfit(spec=gjrgarchsnp.spec, data=rt2) 
gjrgarchsnp_2.fit
coef(gjrgarchsnp_2.fit) 
gjrgarchsnp_3.fit = ugarchfit(spec=gjrgarchsnp.spec, data=rt3) 
gjrgarchsnp_3.fit
coef(gjrgarchsnp_3.fit) 
#residuals(gjrgarchsnp.fit) 　 

# t.test
options(warn=-1)
t.test(coef(garch11_1.fit),coef(garch11_2.fit))$p.value
t.test(coef(egarchsnp_1.fit),rt1)$p.value
t.test(coef(gjrgarchsnp_1.fit),rt1)$p.value
t.test(coef(garch11_2.fit),rt2)$p.value
t.test(coef(egarchsnp_2.fit),rt2)$p.value
t.test(coef(gjrgarchsnp_2.fit),rt2)$p.value
t.test(coef(garch11_3.fit),rt3)$p.value
t.test(coef(egarchsnp_3.fit),rt3)$p.value
t.test(coef(gjrgarchsnp_3.fit),rt3)$p.value

# forecast
# get new value
getSymbols("^IXIC", from='2007-01-01', to= '2007-01-31')
d1 <- dailyReturn(Cl(IXIC), type='log')
getSymbols("^IXIC", from='2010-01-01', to= '2010-01-31')
d2 <- dailyReturn(Cl(IXIC), type='log')
getSymbols("^IXIC", from='2018-01-01', to= '2018-01-31')
d3 <- dailyReturn(Cl(IXIC), type='log')

#GARCH(1,1)
garch_f1<-ugarchforecast(garch11_1.fit, data=NULL,n.ahead=length(d1))
accuracy(ts(fitted(garch_f1)),d1)
garch_f2<-ugarchforecast(garch11_2.fit, data=NULL,n.ahead=length(d2))
accuracy(ts(fitted(garch_f2)),d2)
garch_f3<-ugarchforecast(garch11_3.fit, data=NULL,n.ahead=length(d3))
accuracy(ts(fitted(garch_f3)),d3)
#EGARCH(1,1)
egarch_f1<-ugarchforecast(egarchsnp_1.fit, data=NULL,n.ahead=length(d1))
accuracy(ts(fitted(egarch_f1)),d1)
egarch_f2<-ugarchforecast(egarchsnp_2.fit, data=NULL,n.ahead=length(d2))
accuracy(ts(fitted(egarch_f2)),d2)
egarch_f3<-ugarchforecast(egarchsnp_3.fit, data=NULL,n.ahead=length(d3))
accuracy(ts(fitted(egarch_f3)),d3)
#GJR_GARCH(1,1)
gjrgarch_f1<-ugarchforecast(gjrgarchsnp_1.fit, data=NULL,n.ahead=length(d1))
accuracy(ts(fitted(gjrgarch_f1)),d1)
gjrgarch_f2<-ugarchforecast(gjrgarchsnp_2.fit, data=NULL,n.ahead=length(d2))
accuracy(ts(fitted(gjrgarch_f2)),d2)
gjrgarch_f3<-ugarchforecast(gjrgarchsnp_3.fit, data=NULL,n.ahead=length(d3))
accuracy(ts(fitted(gjrgarch_f3)),d3)

