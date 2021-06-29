library(tseries)
library(zoo)
library(fGarch)

setwd("C:/Users/LENOVO/Desktop/ST Series Temporais/Projeto")

data_read2<-read.csv("ts2.csv",header = TRUE,sep = ",")

ts2<-as.ts(data_read2$Adj.Close)
head(data_read2)
tail(data_read2)
plot.zoo(ts2)

#log-diff the data
log_diff<-diff(log(ts2))
plot(log_diff)

#clear evidence of "clustered" volatility

mean(log_diff)
wilcox.test(log_diff,mu = 0)
acf(log_diff)
pacf(log_diff)

#model selection->garch parameters

modelaux<-garchFit(formula =~garch(1,0),data = log_diff,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula =~garch(2,0),data = log_diff,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula =~garch(1,1),data = log_diff,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula =~garch(2,1),data = log_diff,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula =~garch(1,2),data = log_diff,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula =~garch(2,2),data = log_diff,include.mean = T,trace = F)
summary(modelaux)

#model selection->aparch parameters
modelaux<-garchFit(formula = ~aparch(1,0), data=log_diff, include.delta=TRUE,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula = ~aparch(1,0), data=log_diff, include.delta=TRUE,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula = ~aparch(2,0), data=log_diff, include.delta=TRUE,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula = ~aparch(1,1), data=log_diff, include.delta=TRUE,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula = ~aparch(2,1), data=log_diff, include.delta=TRUE,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula = ~aparch(1,2), data=log_diff, include.delta=TRUE,include.mean = T,trace = F)
summary(modelaux)

modelaux<-garchFit(formula = ~aparch(2,2), data=log_diff, include.delta=TRUE,include.mean = T,trace = F)
summary(modelaux)


#min bic models are garch(1,1) and aparch(1,1)

#distribution estimation

model_0<-garchFit(formula =~garch(1,1),data = log_diff,include.mean = T,trace = F)
summary(model_0)
plot(model_0,which=13)

model_1<-garchFit(formula = ~garch(1,1),data = log_diff,cond.dist = "std",include.mean = T,trace = F)
summary(model_1)
plot(model_1,which=13)



model_a_0<-garchFit(formula = ~aparch(1,1), data=log_diff, include.delta=TRUE,include.mean = T,trace = F)
summary(model_a_0)
plot(model_a_0,which=13)

model_a_1<-garchFit(formula = ~aparch(1,1), data=log_diff, include.delta=TRUE, cond.dist ="std",include.mean = T,trace = F)
summary(model_a_1)
plot(model_a_1,which=13)

#model_1 and model_a_1 with the better fit

par(mfrow=c(1,3))

hist(residuals(object=model_1,standardize=T))
acf(residuals(object=model_1,standardize=T))
pacf(residuals(object=model_1,standardize=T))

hist(residuals(object=model_a_1,standardize=T))
acf(residuals(object=model_a_1,standardize=T))
pacf(residuals(object=model_a_1,standardize=T))

#VaR backstest
par(mfrow=c(1,1))
plot(log_diff)

#GARCH(1,1) with t-student residuals
summary(model_1)
coef<-coef(model_1)

#var95 garch
lines(coef["mu"]+volatility(model_1)*qstd(1-0.95,nu = coef["shape"]),col="blue")
# number of observations under the predicted VaR
sum(coef["mu"]+volatility(model_1)*qstd(1-0.95,nu=coef["shape"])>log_diff)/length(log_diff)
# unexpected losses
sum((log_diff-coef["mu"]+volatility(model_1)*qstd(1-0.95,nu=coef["shape"]))[coef["mu"]+volatility(model_1)*qstd(1-0.95,nu=coef["shape"])>log_diff])


#var99 garch
lines(coef["mu"]+volatility(model_1)*qstd(1-0.99,nu=coef["shape"]),col="red")
# number of observations under the predicted VaR
sum(coef["mu"]+volatility(model_1)*qstd(1-0.99,nu=coef["shape"])>log_diff)/length(log_diff)
# unexpected losses
sum((log_diff-coef["mu"]+volatility(model_1)*qstd(1-0.99,nu=coef["shape"]))[coef["mu"]+volatility(model_1)*qstd(1-0.99,nu=coef["shape"])>log_diff])



#APARCH(1,1) with t-student residuals
par(mfrow=c(1,1))
plot(log_diff)
coef1<-coef(model_a_1)

#var95 
lines(coef1["mu"]+volatility(model_a_1)*qstd(1-0.95,nu=coef1["shape"]),col="blue")
# number of observations over the var
sum(coef1["mu"]+volatility(model_a_1)*qstd(1-0.95,nu=coef1["shape"])>log_diff)/length(log_diff)
# unexpected losses
sum((log_diff-coef1["mu"]+volatility(model_a_1)*qstd(1-0.95,nu=coef1["shape"]))[coef1["mu"]+volatility(model_a_1)*qstd(1-0.95,nu=coef1["shape"])>log_diff])

#var99 

lines(coef1["mu"]+volatility(model_a_1)*qstd(1-0.99,nu=coef1["shape"]),col="red")
# number of observations over the var
sum(coef1["mu"]+volatility(model_a_1)*qstd(1-0.99,nu=coef1["shape"])>log_diff)/length(log_diff)
# unexpected losses
sum((log_diff-coef1["mu"]+volatility(model_a_1)*qstd(1-0.99,nu=coef1["shape"]))[coef1["mu"]+volatility(model_a_1)*qstd(1-0.99,nu=coef1["shape"])>log_diff])

#model_a_1 with the better results beeing the one choosen for building the forecasts

#forecast

#var 1 day
pred<-predict(model_a_1,n.ahead=30,mse="cond")

z_95<-as.numeric(qstd(1-0.95,nu=coef1["shape"]))
z_99<-as.numeric(qstd(1-0.99,nu=coef1["shape"]))

print("Var 95-1 day")
print(pred[1,1]+pred[1,2]*z_95)
print(exp(pred[1,1]+pred[1,2]*z_95)-1)


print("Var 99-1 day")
print(pred[1,1]+pred[1,2]*z_99)
print(exp(pred[1,1]+pred[1,2]*z_99)-1)


print("Var 95-7 day")
print(pred[7,1]+pred[7,2]*z_95)
print(exp(pred[7,1]+pred[7,2]*z_95)-1)

print("Var 99-7 day")
print(pred[7,1]+pred[7,2]*z_99)
print(exp(pred[7,1]+pred[7,2]*z_99)-1)

print("Var 95-30 day")
print(pred[30,1]+pred[30,2]*z_95)
print(exp(pred[30,1]+pred[30,2]*z_95)-1)

print("Var 99-30 day")
print(pred[30,1]+pred[30,2]*z_99)
print(exp(pred[30,1]+pred[30,2]*z_99)-1)

#note-> pred function uses the long term distribution of X to build forecast with more than 1 step


