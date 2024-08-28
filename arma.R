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
library("pcalg")
library("KFAS")

#csvファイルの読み込み
#file.choose("Site2_(Lower_site).csv")
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


# water_level,rain,displacementのプロット
N<-length(Time)
x<-1:N

par(mar = c(5, 8, 3, 4))


#データの階差
n1=1
n2=9000
Time1<-Time[n1:n2]
rain1<-diff(rain[n1:n2])
water_level1<-diff(water_level[n1:n2])
dis11<-diff(dis1[n1:n2])

#AR(1)モデルで地下水位の一期先の予測
predata1<-c()
for(i in 1:228){
  nn1=1     　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(1日分)
  
  data1=data.frame(rain1[n3:n4],water_level1[n3:n4])
  ar1<-ar(water_level1[n3:n4],method = "yule-walker",order.max = 1) #AR
  predata0<-predict(ar1,n.ahead = 1)
  predata1[i]<-predata0$pred[1]
}

#ARIMAモデル
predata1_arima<- c()
for (i in 1:228) {
  nn1 <- 1      # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり(1日分)
  # データのサブセットを取得
  water_subset <- water_level1[n3:n4]
  # ARIMAモデルの自動選択。今回は(0,2,1)が最適。結局ARMAモデルではある
  arima_model <- auto.arima(water_subset)
  # 1ステップ先を予測
  predata0_arima <- forecast(arima_model, h = 1)
  # 予測値を保存
  predata1_arima[i] <- predata0_arima$mean
}

#ARとARIMAを比較
x<-1073:1300
z<-1:228
plot(x, water_level1[x], type = "l", xlim = c(), ylim = c(-0.2, 0.5), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
#plot(z, predata1, type = "l", xlim = c(), ylim = c(-0.2, 0.5), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
#par(new = TRUE)
plot(z, predata1_arima, type = "l", xlim = c(), ylim = c(-0.2, 0.5), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
legend("topleft", legend = c("true value", "AR value", "ARIMA value"), col = c("green", "red", "blue"), lty = 1)






#ARIMAの特性を考えて、water_levelで予測
#water_levelをAR
prelv1<-c()
for(i in 1:228){
  nn1=1     　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(1日分)
  ar1<-ar(water_level[n3:n4],method = "yule-walker",order.max = 1) #AR
  prelv0<-predict(ar1,n.ahead = 1)
  prelv1[i]<-prelv0$pred[1]
}

#water_levelをARIMA
prelv1_arima <- c()
step <- 1  # 1ステップ先を予測
for (i in 1:228) {
  nn1 <- 1      # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり(1日分)
  # データのサブセットを取得
  water_subset <- water_level[n3:n4]
  # ARIMAモデルの自動選択
  arima_model <- auto.arima(water_subset)
  # stepステップ先を予測
  prelv0_arima<- forecast(arima_model, h = step)
  # 予測値を保存
  prelv1_arima[(step*(i-1) + 1):(step*i)] <- prelv0_arima$mean
}

#ARとARIMAを比較
plot(x, water_level[x], type = "l", xlim = c(), ylim = c(0,8), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
plot(z, prelv1, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, prelv1_arima, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
legend("topleft", legend = c("true value", "AR value", "ARIMA value"), col = c("green", "red", "blue"), lty = 1)



#n期先の予測データをAR(1)とARIMAで比較
#一期先の予測AR
prelv1<-c()
for(i in 1:228){
  nn1=1     　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(6h分)
  ar1<-ar(water_level[n3:n4],method = "yule-walker") #AR
  prelv0<-predict(ar1,n.ahead = nn1)
  prelv1[i]<-prelv0$pred[1]
}
#一期先の予測ARIMA
prelv1_arima <- c()
for (i in 1:228) {
  nn1 <- 1      # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり(1日分)
  # データのサブセットを取得
  water_subset <- water_level[n3:n4]
  # ARIMAモデルの自動選択
  arima_model <- auto.arima(water_subset)
  # stepステップ先を予測
  prelv0_arima<- forecast(arima_model, h = nn1)
  # 予測値を保存
  prelv1[i] <- prelv0_arima$mean
}
#比較
x1<-1073:1300
z<-1:228
plot(x, water_level[x], type = "l", xlim = c(), ylim = c(0,8), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
plot(z, prelv1, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, prelv1_arima, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
legend("topleft", legend = c("true value", "AR value", "ARIMA value"), col = c("green", "red", "blue"), lty = 1)



#二期先の予測AR
prelv2<-c()

for(i in 1:228){
  nn1=2    　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(6h分)
  ar1<-ar(water_level[n3:n4],method = "yule-walker") #AR
  prelv0<-predict(ar1,n.ahead = nn1)
  prelv2[i]<-prelv0$pred[nn1]
}
#二期先の予測ARIMA
prelv2_arima <- c()
for (i in 1:228) {
  nn1 <- 2     # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり(1日分)
  # データのサブセットを取得
  water_subset <- water_level[n3:n4]
  # ARIMAモデルの自動選択
  arima_model <- auto.arima(water_subset)
  # stepステップ先を予測
  prelv0_arima<- forecast(arima_model, h = nn1)
  # 予測値を保存
  prelv2_arima[i] <- prelv0_arima$mean[nn1]
}
#比較
x2<-1073:1300
z<-1:228
plot(x2, water_level[x], type = "l", xlim = c(), ylim = c(0,8), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
plot(z, prelv2, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, prelv2_arima, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
legend("topleft", legend = c("true value", "AR value", "ARIMA value"), col = c("green", "red", "blue"), lty = 1)



#三期先の予測AR
prelv3<-c()

for(i in 1:228){
  nn1=3    　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(6h分)
  ar1<-ar(water_level[n3:n4],method = "yule-walker") #AR
  prelv0<-predict(ar1,n.ahead = nn1)
  prelv3[i]<-prelv0$pred[nn1]
}
#三期先の予測ARIMA
prelv3_arima <- c()
for (i in 1:228) {
  nn1 <- 3     # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり(1日分)
  # データのサブセットを取得
  water_subset <- water_level[n3:n4]
  # ARIMAモデルの自動選択
  arima_model <- auto.arima(water_subset)
  # stepステップ先を予測
  prelv0_arima<- forecast(arima_model, h = nn1)
  # 予測値を保存
  prelv3_arima[i] <- prelv0_arima$mean[nn1]
}
#比較
x3<-1073:1300
z<-1:228
plot(x3, water_level[x], type = "l", xlim = c(), ylim = c(0,8), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
plot(z, prelv3, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, prelv3_arima, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
legend("topleft", legend = c("true value", "AR value", "ARIMA value"), col = c("green", "red", "blue"), lty = 1)



#四期先の予測AR
prelv4<-c()

for(i in 1:228){
  nn1=4    　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり
  ar1<-ar(water_level[n3:n4],method = "yule-walker") #AR
  prelv0<-predict(ar1,n.ahead = nn1)
  prelv4[i]<-prelv0$pred[nn1]
}
#四期先の予測ARIMA
prelv4_arima <- c()
for (i in 1:228) {
  nn1 <- 4     # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり
  # データのサブセットを取得
  water_subset <- water_level[n3:n4]
  # ARIMAモデルの自動選択
  arima_model <- auto.arima(water_subset)
  # stepステップ先を予測
  prelv0_arima<- forecast(arima_model, h = nn1)
  # 予測値を保存
  prelv4_arima[i] <- prelv0_arima$mean[nn1]
}
#比較
x4<-1073:1300
z<-1:228
plot(x4, water_level[x], type = "l", xlim = c(), ylim = c(0,8), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
plot(z, prelv4, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, prelv4_arima, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
legend("topleft", legend = c("true value", "AR value", "ARIMA value"), col = c("green", "red", "blue"), lty = 1)


#五期先の予測AR
prelv5<-c()

for(i in 1:228){
  nn1=5    　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり
  ar1<-ar(water_level[n3:n4],method = "yule-walker") #AR
  prelv0<-predict(ar1,n.ahead = nn1)
  prelv5[i]<-prelv0$pred[nn1]
}
#五期先の予測ARIMA
prelv5_arima <- c()
for (i in 1:228) {
  nn1 <- 5     # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり
  # データのサブセットを取得
  water_subset <- water_level[n3:n4]
  # ARIMAモデルの自動選択
  arima_model <- auto.arima(water_subset)
  # stepステップ先を予測
  prelv0_arima<- forecast(arima_model, h = nn1)
  # 予測値を保存
  prelv5_arima[i] <- prelv0_arima$mean[nn1]
}
#比較
x5<-1073:1300
z<-1:228
plot(x5, water_level[x], type = "l", xlim = c(), ylim = c(0,8), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
plot(z, prelv5, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, prelv5_arima, type = "l", xlim = c(), ylim = c(0,8), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
legend("topleft", legend = c("true value", "AR value", "ARIMA value"), col = c("green", "red", "blue"), lty = 1)





#各期ごとの誤差の絶対値の変化のグラフとデータ
#1期目#-water_level[1073:1300]ここの括弧内の数字正しいのは？
shin_waterlv<-water_level[1073:1300]
gosa1<-c()
gosa1_arima<-c()
for (i in 1:length(shin_waterlv)) {
  shin = shin_waterlv[i]
  gosa1[i] <- abs(shin-prelv1[i])/sd(shin_waterlv)
  gosa1_arima[i] <- abs(shin-prelv1_arima[i])/sd(shin_waterlv)
}
x<-1073:1300
z<-1:228

# プロット
plot(z,gosa1, type = "l", col = "red", xlim = c(), ylim = c(0,0.8), xlab = "", ylab = "", lty = 1)
lines(z, gosa1_arima, col = "blue", lty = 1)
# 凡例
legend("topleft", legend = c("AR_err", "ARIMA_err"), col = c("red", "blue"), lty = 1)
#平均2乗誤差
mean_sq1 <- mean(gosa1^2)
mean_sq1_arima <- mean(gosa1_arima^2)
print(c(mean_sq1,mean_sq1_arima))

#二期目
gosa2<-c()
gosa2_arima<-c()
for (i in 1:length(shin_waterlv)) {
  shin = shin_waterlv[i]
  gosa2[i] <- abs(shin-prelv2[i])/sd(shin_waterlv)
  gosa2_arima[i] <- abs(shin-prelv2_arima[i])/sd(shin_waterlv)
}
# プロット
plot(z, gosa2, type = "l", col = "red", xlim = c(), ylim = c(0,1), xlab = "", ylab = "", lty = 1)
lines(z, gosa2_arima, col = "blue", lty = 1)
#平均2乗誤差
mean_sq2 <- mean(gosa2^2)
mean_sq2_arima <- mean(gosa2_arima^2)
print(c(mean_sq2,mean_sq2_arima))


#三期目
gosa3<-c()
gosa3_arima<-c()
for (i in 1:length(shin_waterlv)) {
  shin = shin_waterlv[i]
  gosa3[i] <- abs(shin-prelv3[i])
  gosa3_arima[i] <- abs(shin-prelv3_arima[i])
}
# プロット
plot(z, gosa3, type = "l", col = "red", xlim = c(), ylim = c(0,1), xlab = "", ylab = "", lty = 1)
lines(z, gosa3_arima, col = "blue", lty = 1)
#平均2乗誤差
mean_sq3 <- mean(gosa3^2)
mean_sq3_arima <- mean(gosa3_arima^2)
print(c(mean_sq3,mean_sq3_arima))


#四期目
gosa4<-c()
gosa4_arima<-c()
for (i in 1:length(shin_waterlv)) {
  shin = shin_waterlv[i]
  gosa4[i] <- abs(shin-prelv4[i])
  gosa4_arima[i] <- abs(shin-prelv4_arima[i])
}
# プロット
plot(z, gosa4, type = "l", col = "red", xlim = c(), ylim = c(0,1), xlab = "", ylab = "", lty = 1)
lines(z, gosa4_arima, col = "blue", lty = 1)
#平均2乗誤差
mean_sq4 <- mean(gosa4^2)
mean_sq4_arima <- mean(gosa4_arima^2)
print(c(mean_sq4,mean_sq4_arima))


#五期目の誤差
gosa5<-c()
gosa5_arima<-c()
for (i in 1:length(shin_waterlv)) {
  shin = shin_waterlv[i]
  gosa5[i] <- abs(shin-prelv5[i])/sd(shin_waterlv)
  gosa5_arima[i] <- abs(shin-prelv5_arima[i])/sd(shin_waterlv)
}
# プロット
plot(z, gosa5, type = "l", col = "red", xlim = c(), ylim = c(0,1.3), xlab = "", ylab = "", lty = 1)
lines(z, gosa5_arima, col = "blue", lty = 1)
#平均2乗誤差
mean_sq5 <- mean(gosa5^2)
mean_sq5_arima <- mean(gosa5_arima^2)
print(c(mean_sq5,mean_sq5_arima))
#ARIMA(今回はパラメータ的にARMA)はいい予測精度を持っている。


plot(x5, water_level1[x], type = "l", xlim = c(), ylim = c(0, 0.5), xlab = "Sample Time (5min)", ylab = "Ground Water level_diff (m)", col = "green")
par(new = TRUE)
plot(z, gosa5, type = "l", xlim = c(), ylim = c(0, 1.5), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, gosa5_arima, type = "l", xlim = c(), ylim = c(0, 1.5), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
# y軸のラベルを追加
axis(4)  # 右側のy軸を表示
# 凡例
legend("topleft", legend = c("true_diff", "AR_err", "ARIMA_err"), col = c("green", "red", "blue"), lty = 1)








#差分取って予測した時との誤差の比較
#一期先の予測AR
predif1<-c()
for(i in 1:228){
  nn1=1     　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(6h分)
  ar1<-ar(water_level1[n3:n4],method = "yule-walker") #AR
  predif0<-predict(ar1,n.ahead = nn1)
  predif1[i]<-predif0$pred[1]
}
#一期先の予測ARIMA
predif1_arima <- c()
for (i in 1:228) {
  nn1 <- 1      # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり(1日分)
  # データのサブセットを取得
  water_subset <- water_level1[n3:n4]
  # ARIMAモデルの自動選択
  arima_model <- auto.arima(water_subset)
  # stepステップ先を予測
  predif0_arima<- forecast(arima_model, h = nn1)
  # 予測値を保存
  predif1_arima[i] <- predif0_arima$mean
}
#比較
x1<-1073:1300
z<-1:228
plot(x, water_level1[x], type = "l", xlim = c(), ylim = c(0,0.6), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
plot(z, predif1, type = "l", xlim = c(), ylim = c(0,0.6), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, predif1_arima, type = "l", xlim = c(), ylim = c(0,0.6), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
legend("topleft", legend = c("true value", "AR value", "ARIMA value"), col = c("green", "red", "blue"), lty = 1)


#1期の誤差
shin_waterlv1<-water_level1[1073:1300]
gosadif1<-c()
gosadif1_arima<-c()
for (i in 1:length(shin_waterlv)) {
  shin = shin_waterlv1[i]
  gosadif1[i] <- abs(shin-predif1[i])/sd(shin_waterlv1)
  gosadif1_arima[i] <- abs(shin-predif1_arima[i])/sd(shin_waterlv1)
}
x<-1073:1300
z<-1:228

plot(z, gosadif1, type = "l", col = "green", xlim = c(), ylim = c(0,5), xlab = "", ylab = "", lty = 1)
#lines(z, gosadif1, col = "green", lty = 1)
lines(z, gosadif1_arima, col = "red", lty = 1)
#lines(z, predif1_arima, col = "blue", lty = 1)
# 凡例
legend("topleft", legend = c("AR_err", "ARIMA_err"), col = c("green", "red", "blue"), lty = 1)
mean_dif1 <- mean(gosadif1^2)
mean_dif1_arima <- mean(gosadif1_arima^2)
print(c(mean_dif1,mean_dif1_arima))
print(c(mean_sq1,mean_sq1_arima))


plot(z, gosadif1, type = "l", col = "red", xlim = c(), ylim = c(0,8), xlab = "", ylab = "", lty = 1)
lines(z, gosadif1_arima, col = "blue", lty = 1)
lines(z, gosa1_arima, col = "green", lty = 1)
lines(z, predif1_arima, col = "purple", lty = 1)
# 凡例
legend("topleft", legend = c("ARdif_err", "ARIMAdif_err","AR_err","ARIMA_err"), col = c( "red", "blue","green","purple"), lty = 1)





#5期先の予測AR
predif5<-c()
for(i in 1:228){
  nn1=5    　　#予測数
  n3=999+i 　　 #観測始め
  n4=n3+72　 　#観測終わり(6h分)
  ar1<-ar(water_level1[n3:n4],method = "yule-walker") #AR
  predif0<-predict(ar1,n.ahead = nn1)
  predif5[i]<-predif0$pred[nn1]
}
#5期先の予測ARIMA
predif5_arima <- c()
for (i in 1:228) {
  nn1 <- 5      # 予測数
  n3 <- 999 + i # 観測始め
  n4 <- n3 + 72 # 観測終わり(1日分)
  # データのサブセットを取得
  water_subset <- water_level1[n3:n4]
  # ARIMAモデルの自動選択
  arima_model <- auto.arima(water_subset)
  # stepステップ先を予測
  predif0_arima<- forecast(arima_model, h = nn1)
  # 予測値を保存
  predif5_arima[i] <- predif0_arima$mean[nn1]
}

x1<-1073:1300
z<-1:228
plot(x, water_level1[x], type = "l", xlim = c(), ylim = c(0,0.6), xlab = "Sample Time (5min)", ylab = "Ground Water level (m)", col = "green")
par(new = TRUE)
plot(z, predif5, type = "l", xlim = c(), ylim = c(0,0.6), lty = 1, xlab = "", ylab = "", col = "red", axes = FALSE)
par(new = TRUE)
plot(z, predif5_arima, type = "l", xlim = c(), ylim = c(0,0.6), lty = 1, xlab = "", ylab = "", col = "blue", axes = FALSE)
#par(new = TRUE)
#plot(z, gosadif5, type = "l", xlim = c(), ylim = c(0,0.6), lty = 1, xlab = "", ylab = "", col = "purple", axes = FALSE)
#par(new = TRUE)
#plot(z, gosadif5_arima, type = "l", xlim = c(), ylim = c(0,0.6), lty = 1, xlab = "", ylab = "", col = "yellow", axes = FALSE)


legend("topleft", legend = c("true value", "AR value", "ARIMA value"), col = c("green", "red", "blue"), lty = 1)



#5期の誤差
shin_waterlv1<-water_level1[1073:1300]
gosadif5<-c()
gosadif5_arima<-c()
for (i in 1:length(shin_waterlv)) {
  shin = shin_waterlv1[i]
  gosadif5[i] <- abs(shin-predif5[i])/sd(shin_waterlv1)
  gosadif5_arima[i] <- abs(shin-predif5_arima[i])/sd(shin_waterlv1)
}
x<-1073:1300
z<-1:228

plot(z, gosadif5, type = "l", col = "green", xlim = c(), ylim = c(0,0.6), xlab = "", ylab = "", lty = 1)
lines(z, gosadif5_arima, col = "red", lty = 1)
# 凡例
legend("topleft", legend = c("AR_err", "ARIMA_err"), col = c("green", "red", "blue"), lty = 1)
mean_dif5 <- mean(gosadif5^2)
mean_dif5_arima <- mean(gosadif5_arima^2)
print(c(mean_dif5,mean_dif5_arima))
print(c(mean_sq5,mean_sq5_arima))


plot(z, gosadif5, type = "l", col = "red", xlim = c(), ylim = c(0,6), xlab = "", ylab = "", lty = 1)
lines(z, gosadif5_arima, col = "blue", lty = 1)
lines(z, gosa5, col = "green", lty = 1)
lines(z, gosa5_arima, col = "purple", lty = 1)
# 凡例
legend("topleft", legend = c("ARdif_err", "ARIMAdif_err","AR_err","ARIMA_err"), col = c( "red", "blue","green","purple"), lty = 1)
#分散標準偏差で割ったら、直がけの方が成績いい？















#################ARIMAモデルで2274から2473で2523までを予測######
predis<- c()
for (i in 1:250) {
  nn1=1     　　#予測数
  n3=2200+i 　　 #観測始め
  n4=n3+72　
  # データのサブセットを取得
  if(i<200){
    dis_sub <- dis1[n3:n4]
    # ARIMAモデルの自動選択。今回は(0,2,1)が最適。結局ARMAモデルではある
    arima_model <- auto.arima(dis_sub)
    # 1ステップ先を予測
    predis0 <- forecast(arima_model, h = 1)
    # 予測値を保存
    predis[i] <- predis0$mean
  }else{
    dis_sub <- predis[(i-73):(i-1)]
    arima_model <- auto.arima(dis_sub)
    predis0 <- forecast(arima_model, h = 1)
    predis[i] <- predis0$mean
  }
}

dev.off()
x<-2274:2523
z<-1:250
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,250),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
plot(z,predis,type="l",xlim=c(),ylim=c(180,250),lty=2,xlab="",ylab="",col="purple",axes="F")
legend("topleft", legend = c("ture value","predict value"), col = c("red","purple"), lty =c(1,2))

#2乗誤差
hani_x<- 2473:2523
hani_y<-200:250
x1 <-dis1[hani_x]
y1<-predis[hani_y]
residuals <- x1 - y1
mse <- mean(residuals^2)
mae <- mean(abs(residuals))
cat("MSE:", mse, "\n")  #MSE:0.01526046 
cat("MAE:", mae, "\n")  #MAE:0.09754634 

####################ARIMAモデルで2474から2673で2723までを予測######
predis<- c()
for (i in 1:250) {
  nn1=1     　　#予測数
  n3=2400+i 　　 #観測始め
  n4=n3+72　
  # データのサブセットを取得
  if(i<200){
    dis_sub <- dis1[n3:n4]
    # ARIMAモデルの自動選択。今回は(0,2,1)が最適。結局ARMAモデルではある
    arima_model <- auto.arima(dis_sub)
    # 1ステップ先を予測
    predis0 <- forecast(arima_model, h = 1)
    # 予測値を保存
    predis[i] <- predis0$mean
  }else{
    dis_sub <- predis[(i-73):(i-1)]
    arima_model <- auto.arima(dis_sub)
    predis0 <- forecast(arima_model, h = 1)
    predis[i] <- predis0$mean
  }
}

dev.off()
x<-2474:2723
z<-1:250
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,250),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
plot(z,predis,type="l",xlim=c(),ylim=c(180,250),lty=2,xlab="",ylab="",col="purple",axes="F")
legend("topleft", legend = c("ture value","predict value"), col = c("red","purple"), lty =c(1,2))

#2乗誤差
hani_x<- 2673:2723
hani_y<-200:250
x1 <-dis1[hani_x]
y1<-predis[hani_y]
residuals <- x1 - y1
mse <- mean(residuals^2)
mae <- mean(abs(residuals))
cat("MSE:", mse, "\n")  #MSE: 39.15869 
cat("MAE:", mae, "\n")  #MAE: 4.886543 


####################ARIMAモデルで2574から2773で2823までを予測######
predis<- c()
for (i in 1:250) {
  nn1=1     　　#予測数
  n3=2500+i 　　 #観測始め
  n4=n3+72　
  # データのサブセットを取得
  if(i<200){
    dis_sub <- dis1[n3:n4]
    # ARIMAモデルの自動選択。今回は(0,2,1)が最適。結局ARMAモデルではある
    arima_model <- auto.arima(dis_sub)
    # 1ステップ先を予測
    predis0 <- forecast(arima_model, h = 1)
    # 予測値を保存
    predis[i] <- predis0$mean
  }else{
    dis_sub <- predis[(i-73):(i-1)]
    arima_model <- auto.arima(dis_sub)
    predis0 <- forecast(arima_model, h = 1)
    predis[i] <- predis0$mean
  }
}

dev.off()
x<-2574:2823
z<-1:250
plot(x,dis1[x],type="l",xlim=c(),ylim=c(180,300),xlab="Sample time (5min)",ylab="Displacement (mm)",col="red")
par(new=T)
plot(z,predis,type="l",xlim=c(),ylim=c(180,300),lty=2,xlab="",ylab="",col="purple",axes="F")
legend("topleft", legend = c("ture value","predict value"), col = c("red","purple"), lty =c(1,2))

#2乗誤差
hani_x<- 2774:2823
hani_y<-200:250
x1 <-dis1[hani_x]
y1<-predis[hani_y]
residuals <- x1 - y1
mse <- mean(residuals^2)
mae <- mean(abs(residuals))
cat("MSE:", mse, "\n")  #MSE: 12.28631  
cat("MAE:", mae, "\n")  #MAE: 1.303257 


