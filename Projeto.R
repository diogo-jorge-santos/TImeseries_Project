library(tseries)
library(zoo)
library(forecast)


setwd("C:/Users/LENOVO/Desktop/ST Series Temporais/Projeto")

data_read1<-read.csv("ts1.csv",header = TRUE,sep = ";")

head(data_read1)
tail(data_read1)
dim(data_read1)
sum(is.na(data_read1))


#compute daily mean, with existing hourly data, i.e, if it has only 22 hours of observed data
#compute the mean of those 22 days, instead of 24
ts<-aggregate(ts(data_read1$Av.da.Lib...µg.m3., freq = 24), 1, mean,na.rm = T)
#replace null values with the last non-null value
ts<-na.locf(ts)


sum(is.na(ts))

#replace null values with previous non-NA observation (in order to not create bias)

ts1<-ts(ts)


par(mfrow=c(2,1))
plot(ts1)
hist(ts1)
par(mfrow=c(2,1))
acf(ts1)
pacf(ts1)

#non stable variance -> box transformation

#estimation for the value of lambda on a boxcax transformation
BoxCox.lambda(ts1)
ts1<-BoxCox(ts1,lambda = 0.05655603)

par(mfrow=c(2,1))

plot(ts1)
hist(ts1)

adf.test(ts1)

#rej H0, do not diff
#but mean does not equal to 0-> introduce mean term to sarima model



#non-normal data
shapiro.test(ts1)
jarque.bera.test(ts1)


#box-transformed data can be modeled bay a arima(p,0,p) model


#model selection->arma
out <- matrix(NA, nrow=25, ncol=3)
aux<-1

for(p in 0:4){
  for(q in 0:4){
    model<-arima(ts1,order=c(p,0,q),include.mean = T)
    out[aux,]<-c(p,q,BIC(model))
    aux<-aux+1
    }
}
aux_df<-as.data.frame(out[order(out[,3]),])
colnames(aux_df)<-c("p","q","BIC")
print(aux_df)


#model diagnostics -> arma
model1<-arima(ts1,order=c(2,0,2),include.mean = T,method = "ML")
summary(model1)

par(mfrow=c(1,3))
hist(model1$residuals)
plot(model1$residuals)
qqnormPlot(model1$residuals)

par(mfrow=c(1,2))
acf(model1$residuals)
pacf(model1$residuals)

Box.test(model1$residuals,type="Ljung-Box",lag = 10,fitdf = 2+2)
Box.test(model1$residuals,type="Ljung-Box",lag = 20,fitdf = 2+2)
Box.test(model1$residuals,type="Ljung-Box",lag = 30,fitdf = 2+2)
Box.test(model1$residuals,type="Ljung-Box",lag = 50,fitdf = 2+2)

shapiro.test(model1$residuals)
jarque.bera.test(model1$residuals)

model2<-arima(ts1,order=c(3,0,1),include.mean = T,method = "ML")
summary(model2)

par(mfrow=c(1,3))
hist(model2$residuals)
plot(model2$residuals)
qqnormPlot(model2$residuals)
par(mfrow=c(1,2))
acf(model2$residuals)
pacf(model2$residuals)

Box.test(model2$residuals,type="Ljung-Box",lag = 10,fitdf = 3+1)
Box.test(model2$residuals,type="Ljung-Box",lag = 20,fitdf = 3+1)
Box.test(model2$residuals,type="Ljung-Box",lag = 30,fitdf = 3+1)
Box.test(model2$residuals,type="Ljung-Box",lag = 50,fitdf = 3+1)


shapiro.test(model2$residuals)
jarque.bera.test(model2$residuals)

#both do not pass whiteness test 
#fit sarma model with m=28


#model selection-> sarma
out2 <- matrix(NA, nrow=9, ncol=3)
aux<-1
#note-> these for loops take a long time to finish (1-2 minutes)
for(p in 0:2){
  for(q in 0:2){
    model<-arima(ts(ts1,freq=28),order=c(2,0,2),seasonal = c(p,0,q),include.mean = T,method = "ML")
    out2[aux,]<-c(p,q,BIC(model))
    aux<-aux+1
  }
}
aux_df1<-as.data.frame(out2[order(out2[,3]),])
colnames(aux_df1)<-c("P","Q","BIC")
print(aux_df1)

#model diagnostics->sarma
model3<-arima(ts(ts1,freq=28),order=c(2,0,2),seasonal = c(1,0,1),include.mean = T,method = "ML")
summary(model3)


par(mfrow=c(1,3))
hist(model3$residuals)
plot(model3$residuals)
qqnormPlot(model3$residuals)

shapiro.test(model3$residuals)
jarque.bera.test(model3$residuals)

par(mfrow=c(1,2))
acf(model3$residuals)
pacf(model3$residuals)

Box.test(model3$residuals,type="Ljung-Box",lag = 10,fitdf = 2+2+1+1)
Box.test(model3$residuals,type="Ljung-Box",lag = 15,fitdf = 2+2+1+1)
Box.test(model3$residuals,type="Ljung-Box",lag = 20,fitdf = 2+2+1+1)
Box.test(model3$residuals,type="Ljung-Box",lag = 30,fitdf = 2+2+1+1)
Box.test(model3$residuals,type="Ljung-Box",lag = 100,fitdf = 2+2+1+1)
Box.test(model3$residuals,type="Ljung-Box",lag = 365,fitdf = 2+2+1+1)

#out of sample testing

#5-lags out of sample
ts_test<-ts(ts1,freq=28)
train<-subset(ts_test,end=2190-5)
Real_value<-subset(ts_test,start=2190-5+1)

model_test<-arima(train,order=c(2,0,2),seasonal = c(1,0,1),include.mean = T,method = "ML")
summary(model_test)
summary(model3)

forecast_test<-forecast(model_test,h=5,level=0.95)
forecast_test_bt<-forecast(model_test,h=5,level=0.95,bootstrap = T)
print(forecast_test)
autoplot(forecast_test,include=20) + autolayer(Real_value)

print(forecast_test_bt)
autoplot(forecast_test_bt,include=20) + autolayer(Real_value)


#100-lags out of sample
train<-subset(ts_test,end=2190-100)
Real_value<-subset(ts_test,start=2190-100+1)

model_test<-arima(train,order=c(2,0,2),seasonal = c(1,0,1),include.mean = T)

forecast_test_l<-forecast(model_test,h=100,level=0.95)
forecast_test1_l<-forecast(model_test,h=100,level=0.95,bootstrap = T)
autoplot(forecast_test_l,include=20) + autolayer(Real_value)
autoplot(forecast_test1_l,include=20) + autolayer(Real_value)

sum(forecast_test_l$lower>Real_value)+sum(Real_value>forecast_test_l$upper)
sum(forecast_test1_l$lower>Real_value)+sum(Real_value>forecast_test_l$upper)



#forecast

model_forecast<-forecast(model3,h = 5,level = 0.95)
#undo the box-cox transformations previosly applied-> using the invariant property of maximum likelihood estimates
model_forecast$mean<-InvBoxCox(model_forecast$mean,lambda = 0.05655603)
model_forecast$lower<-InvBoxCox(model_forecast$lower,lambda = 0.05655603)
model_forecast$upper<-InvBoxCox(model_forecast$upper,lambda = 0.05655603)
model_forecast$x<-InvBoxCox(model_forecast$x,lambda = 0.05655603)
model_forecast$fitted<-InvBoxCox(model_forecast$fitted,lambda = 0.05655603)
model_forecast$residuals<-InvBoxCox(model_forecast$residuals,lambda = 0.05655603)

print(model_forecast)

autoplot(model_forecast,include = 100)

#forecast bootstrap

model_forecast_b<-forecast(model3,h = 5,bootstrap = T,level = 0.95)

model_forecast_b$mean<-InvBoxCox(model_forecast_b$mean,lambda = 0.05655603)
model_forecast_b$lower<-InvBoxCox(model_forecast_b$lower,lambda = 0.05655603)
model_forecast_b$upper<-InvBoxCox(model_forecast_b$upper,lambda = 0.05655603)
model_forecast_b$x<-InvBoxCox(model_forecast_b$x,lambda = 0.05655603)
model_forecast_b$fitted<-InvBoxCox(model_forecast_b$fitted,lambda = 0.05655603)
model_forecast_b$residuals<-InvBoxCox(model_forecast_b$residuals,lambda = 0.05655603)
print(model_forecast_b)

autoplot(model_forecast_b,include = 100)



