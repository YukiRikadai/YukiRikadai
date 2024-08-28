library("forecast")
library("TTR")
library("statmod")
library("tseries")
library("cpm")
library("changepoint")
library("Kendall")
library("trend")
library("wbs")
library("ecp")
library("vars")
library("BVAR")
library("tidyverse")
#library("pcalg")
library("KFAS")
####################ここから###########################
data<-read.csv("/Users/ohashi/Desktop/統数研/ohashi_業務/Site2_(Lower_site).csv",header = T)

#変数の設定
water_level <- data$Site2.Lower.site.Water.level.m. ## Ground Water Level（地下水位）
rain <- data$Site2.Rain.mm. ## Rain　（降雨量）
Time <- data$Sampling.Time ## Time　（時間）

## Inclinometer （３つの計測器の値）
inc18 <- data$Site2.Lower.site.Inclinometer.GL.1.8m.X.
inc38 <- data$Site2.Lower.site.Inclinometer.GL.3.8m.X.
inc58 <- data$Site2.Lower.site.Inclinometer.GL.5.8m.X.

#displacement (mm) <- data$Displacement*1000　
dis1 <- 1.8*sin(inc18*pi/180)/cos(32.5*pi/180)*1000　#地表面変位
vel <- c(0,diff(dis1))

#データの階差
n1=1
n2=9000
Time1<-Time[n1:n2]
rain1<-diff(rain[n1:n2])
water_level1<-diff(water_level[n1:n2])
dis11<-diff(dis1[n1:n2])
####################ここまで###########################

#1期先の予測AR
prelv1<-c()
for(i in 1:228){
  nn1=1    　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(6h分)
  ar1<-ar(water_level[n3:n4],method = "yule-walker") #AR
  prelv0<-predict(ar1,n.ahead = nn1)
  prelv1[i]<-prelv0$pred[nn1]
}

#1期先の予測ARIMA
prelv1_arima <- c()
for (i in 1:228) {
  nn1 <- 1     # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり(1日分)
  # データのサブセットを取得
  water_subset <- water_level[n3:n4]
  # ARIMAモデルの自動選択
  arima_model <- auto.arima(water_subset)
  # stepステップ先を予測
  prelv0_arima<- forecast(arima_model, h = nn1)
  # 予測値を保存
  prelv1_arima[i] <- prelv0_arima$mean[nn1]
}


#1期先の予測VAR
prelv1_var <- c()
for (i in 1:228) {
  nn1 <- 1      # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり(1日分)
  
  # VARモデルのためのデータフレームを作成
  data1 <- data.frame(dis11[n3:n4], water_level1[n3:n4])
  # VARモデルを適用
  data1_var <- VAR(data1, type = "both", ic = "AIC")
  # 予測を実施
  predata <- predict(data1_var, n.ahead = nn1)
  # 予測結果を保存 (water_levelの予測値)
  prelv1_var[i] <- predata$fcst$dis11[1, 1]
}

#比較
x1<-1073:1300
z<-1:228
plot(x1, water_level1[x1], type = "l", xlim = c(), ylim = c(), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
plot(x1, dis11[x1], type = "l", xlim = c(), ylim = c(0,2.5), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, prelv1_var, type = "l", xlim = c(), ylim = c(0,2.5), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
par(new = TRUE)
plot(z, prelv1_var, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "yellow", axes = FALSE)
legend("topleft", legend = c("true value", "AR value", "ARIMA value","VAR value"), col = c("green", "red", "blue","yellow"), lty = 1)

##########地面変位を予測##############
predata1_var<-c()
for(i in 1:827){
  nn1=1     　　#予測数
  n3=2500+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(1日分)
  
  data1=data.frame(rain[n3:n4],dis1[n3:n4])
  data1_var<-VAR(data1,type="both",ic="AIC")  #VAR
  predata<-predict(data1_var,n.ahead = nn1)
  predata1_var[i]<-predata$fcst$dis1[,1]
}

x<-2574:3400
z<-1:827
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,410),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
plot(z,predata1_var,type="l",xlim=c(),ylim=c(180,410),lty=2,xlab="",ylab="",col="purple",axes="F")
legend("topleft", legend = c("ture value","predict value"), col = c("red","purple"), lty =c(1,2))

#例えば、3000地点までを読みこんで以降を予測
predata1_var<-c()
for(i in 1:827){
  nn1=1     　　#予測数
  n3=2500+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(1日分)
  if(i < 500) {
    data1 <- data.frame(rain[n3:n4], dis1[n3:n4])
    data1_var <- VAR(data1, type = "both", ic = "AIC")  # VAR
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$dis1[,1]
  } else {
    previous_predictions <- predata1_var[(i - 73):(i - 1)]  # 直前の 72 個の予測値を取得
    data1 <- data.frame(rain[n3:n4], previous_predictions)
    data1_var <- VAR(data1, type = "both", ic = "AIC")
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$previous_predictions[,1]
  }
}

x<-2574:3400
z<-1:827
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,410),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
plot(z,predata1_var,type="l",xlim=c(),ylim=c(180,410),lty=2,xlab="",ylab="",col="purple",axes="F")
legend("topleft", legend = c("ture value","predict value"), col = c("red","purple"), lty =c(1,2))

x<-2574:3400
z<-1:827
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,410),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
plot(z,predata1_var,type="l",xlim=c(),ylim=c(180,410),lty=2,xlab="",ylab="",col="purple",axes="F")
legend("topleft", legend = c("ture value","predict value"), col = c("red","purple"), lty =c(1,2))

################displacementデータ########
x<-1:length(dis1)
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,600),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
dev.off()
x<-2200:3400
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,500),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)

######2200〜2400を使って2400〜2450を予測し、#########
######適合しているか調べる#########
predata1_var<-c()
for(i in 1:250){
  nn1=1     　　#予測数
  n3=2200+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(1日分)
  if(i < 200) {
    data1 <- data.frame(rain[n3:n4], dis1[n3:n4])
    data1_var <- VAR(data1, type = "both", ic = "AIC")  # VAR
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$dis1[,1]
  } else {
    previous_predictions <- predata1_var[(i - 73):(i - 1)]  # 直前の 72 個の予測値を取得
    data1 <- data.frame(rain[n3:n4], previous_predictions)
    data1_var <- VAR(data1, type = "both", ic = "AIC")
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$previous_predictions[,1]
  }
}

dev.off()
x<-2274:2523
z<-1:250
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,250),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
plot(z,predata1_var,type="l",xlim=c(),ylim=c(180,250),lty=2,xlab="",ylab="",col="purple",axes="F")
legend("topleft", legend = c("ture value","predict value"), col = c("red","purple"), lty =c(1,2))

##仮説検定##
hani_x<- 2473:2523
hani_y<-200:250
x1 <-dis1[hani_x]
y1<-predata1_var[hani_y]
residuals <- x1 - y1
mse <- mean(residuals^2)
mae <- mean(abs(residuals))
cat("MSE:", mse, "\n")
cat("MAE:", mae, "\n")
#MSE=0.01455129, MAE=0.09510347 





###########次は2474から2673で2723を予測######
predata1_var<-c()
for(i in 1:250){
  nn1=1     　　#予測数
  n3=2400+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(1日分)
  if(i < 200) {
    data1 <- data.frame(rain[n3:n4], dis1[n3:n4])
    data1_var <- VAR(data1, type = "both", ic = "AIC")  # VAR
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$dis1[,1]
  } else {
    previous_predictions <- predata1_var[(i - 73):(i - 1)]  # 直前の 72 個の予測値を取得
    data1 <- data.frame(rain[n3:n4], previous_predictions)
    data1_var <- VAR(data1, type = "both", ic = "AIC")
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$previous_predictions[,1]
  }
}

dev.off()
x<-2474:2723
z<-1:250
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,250),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
plot(z,predata1_var,type="l",xlim=c(),ylim=c(180,250),lty=2,xlab="",ylab="",col="purple",axes="F")
legend("topleft", legend = c("ture value","predict value"), col = c("red","purple"), lty =c(1,2))

##仮説検定##
hani_x<- 2673:2723
hani_y<-200:250
x1 <-dis1[hani_x]
y1<-predata1_var[hani_y]
residuals <- x1 - y1
mse <- mean(residuals^2)
mae <- mean(abs(residuals))
cat("MSE:", mse, "\n")
cat("MAE:", mae, "\n")
#MSE=5787.835 , MAE=49.00124


###########次は2574から2773で2823を予測######
predata1_var<-c()
for(i in 1:250){
  nn1=1     　　#予測数
  n3=2500+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(1日分)
  if(i < 200) {
    data1 <- data.frame(rain[n3:n4], dis1[n3:n4])
    data1_var <- VAR(data1, type = "both", ic = "AIC")  # VAR
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$dis1[,1]
  } else {
    previous_predictions <- predata1_var[(i - 73):(i - 1)]  # 直前の 72 個の予測値を取得
    data1 <- data.frame(rain[n3:n4], previous_predictions)
    data1_var <- VAR(data1, type = "both", ic = "AIC")
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$previous_predictions[,1]
  }
}

dev.off()
x<-2574:2823
z<-1:250
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,300),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
plot(z,predata1_var,type="l",xlim=c(),ylim=c(180,300),lty=2,xlab="",ylab="",col="purple",axes="F")
legend("topleft", legend = c("ture value","predict value"), col = c("red","purple"), lty =c(1,2))

##仮説検定##
hani_x<- 2773:2823
hani_y<-200:250
x1 <-dis1[hani_x]
y1<-predata1_var[hani_y]
residuals <- x1 - y1
mse <- mean(residuals^2)
mae <- mean(abs(residuals))
cat("MSE:", mse, "\n")
cat("MAE:", mae, "\n")
#MSE=1.364934  , MAE=1.046811 

#>>帰無仮説採択


########ある地点から帰無仮説棄却→採択になれば、#######
########簡易的ではあるが、予測可能地点がわかるのでは？####

for(k in 1:1000){
predata1_var<-c()
shoki = 2400+k-1
for(i in 1:250){
  nn1=1     　　#予測数
  n3=shoki+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(1日分)
  if(i < 200) {
    data1 <- data.frame(rain[n3:n4], dis1[n3:n4])
    data1_var <- VAR(data1, type = "both", ic = "AIC")  # VAR
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$dis1[,1]
  } else {
    previous_predictions <- predata1_var[(i - 73):(i - 1)]  # 直前の 72 個の予測値を取得
    data1 <- data.frame(rain[n3:n4], previous_predictions)
    data1_var <- VAR(data1, type = "both", ic = "AIC")
    predata <- predict(data1_var, n.ahead = nn1)
    predata1_var[i] <- predata$fcst$previous_predictions[,1]
  }
}
##仮説検定##
hani_x<- (2673+k):(2723+k)
hani_y<-200:250
x1 <-dis1[hani_x]
y1<-predata1_var[hani_y]
test_statistic <- mean(y1 - x1)
n <- length(x1) # 標本サイズ
alpha <- 0.05 #有意水準
degrees_of_freedom <- n - 1 #自由度
critical_value <- qt(1 - alpha, df = degrees_of_freedom) #棄却点
if (test_statistic < critical_value) {
  cat("帰無仮説受容ポイント\n")
  print(shoki+200)
  break
}
}
#よくよく考えてみたらこれARIMAでやった方が理屈にあうじゃん
