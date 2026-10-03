### EpiPredFunctionsRevised20260819.R ###
##
## EpiPredFunctionsRevised20240913.Rからの更新内容 ## 
##  関数param_estim_by_lmの戻り値にR squaredを追加
##  関数calculate_Rsquaredを追加
##  関数vi.param_estim_by_emの内容変更--最適化関数をL-BFGS-BからCGへ変更
##
## EpiPredFunctionsRevised20241007.Rからの更新内容 ##
##  Shedding Profile Modelを追加
##
## EpiPredFunctionsRevised20250506.Rからの更新内容 ##
##  関数adjusted_Rsquaredを追加（決定係数と自由度調整済決定係数）
##  関数Rsquared_in_lmを追加:
##   calculate_Rsquaredの戻り値に自由度調整済決定係数を追加したもの
##  パラメータ推定関数の戻り値に決定係数と自由度調整済決定係数を追加
##
## EpiPredFunctionsRevised20250616.Rからの更新内容 ##
##  関数shedding_profileを修正
##
## EpiPredFunctionsRevised20250707.Rからの更新内容 ##
##  residual varianceを修正
##
## EpiPredFunctionsRevised20250718.Rからの更新内容 ##
##  最適化とその他軽微なミスを修正
##  関数Rsquared_in_lmの戻り値変更　(SSE,SSR,SST,決定係数1=1-(SSE/SST),決定係数2=SSR/SST)
##  関数adjusted_Rsquaredの戻り値変更　(SSE,SSR,SST,決定係数1=1-(SSE/SST),決定係数2=SSR/SST)
##　パラメータ推定関数の戻り値を追加　(SSE,SSR,SST,決定係数1=1-(SSE/SST),決定係数2=SSR/SST)
##
## EpiPredFunctionsRevised20260426.Rからの更新内容 ##
##  関数infected_prediction_by_SPMを修正
##  感染者数を対数値に変換して最適化と平滑化を実行するように変更


## required libraries 
## library("dplyr") # for mutate 
## library("magrittr") # for %<>%
## library(nleqslv)　# for nleqslv
## library(nnls) # for nnls　# 2026/03/31追加

# find shape and scale parameters used in weibull distribution from mean and standard deviation of a data
find_me_from_sd <-function(xmean,xsd){
  
  g <- function(z){
  
    m=(z[1])^2  # shape parameter
    e=(z[2])^2  # scale parameter
    
    f1=e*gamma(1+(1/m))-xmean  # mean
    f2=sqrt( (e^2)*( gamma(1+(2/m)) - gamma(1+(1/m))^2 ) ) -xsd  
    
    c(f1,f2)
  } 
  
  init <- c(1,1)  # initial value
  
  gsol <- nleqslv(init, g) 
  
 
  return(c((gsol$x[1])^2, (gsol$x[2])^2)) # return(m:shape,e:scale)
  
}


#back projection of disease onset counts from the number of reported cases using Monte Carlo method (Notifiable disease surveillance)
generate_onset <- function(reporteddata,xmu,xsd){  # xmu:mean, xsd:SD
  
  # find shape parameter m and scale parameter e from mean and SD, respectively.
  m <- find_me_from_sd(xmu,xsd)[1] # shape
  e <- find_me_from_sd(xmu,xsd)[2] # scale
  
  colnames(reporteddata) <- c("date","reported")

  rep_delay = list(shape = m, scale = e) 
  generated_onset_day_list <- list()
  
  for (d in 1:length(reporteddata$date)) {
    rep.day <-  reporteddata$date[d]
    n <-  reporteddata$reported[d]
    if (n > 0) {
      r <- rweibull(n, shape = rep_delay$shape, scale = rep_delay$scale)
      generated_onset_day_list <- c(generated_onset_day_list, rep.day-r)
    }
  }
  generated_onset_day_list <- as.Date(unlist(generated_onset_day_list))
  
  df_gen_onset <- as.data.frame(table(generated_onset_day_list))
  colnames(df_gen_onset) <- c("date", "onset")
  df_gen_onset %<>% mutate(df_gen_onset, date=as.Date(date))
  
  df_rep_ons <- merge(reporteddata, df_gen_onset, by="date", all.x=T)
  df_rep_ons [is.na(df_rep_ons)] = 0 # NA in column onset is replaced by 0
  df_rep_ons$onset <- replace(df_rep_ons$onset,nrow(df_rep_ons),df_rep_ons$onset[nrow(df_rep_ons)-1]) 
  
  
  return(df_rep_ons[,-2])
}

# generate back-projected disease onset counts from sentinel data -- 2024/08/21
generate_onset_sentinel <- function(sentineldata,xmu,xsd){  # xmu:mean, xsd:SD
  
  # find shape parameter m and scale parameter e from mean and SD, respectively.
  m <- find_me_from_sd(xmu,xsd)[1] # shape
  e <- find_me_from_sd(xmu,xsd)[2] # scale
  
  # multiply the number of reported by 100 or 1000 for following calculations
  if(100*min(sentineldata[,2])<1.0){ord<-1000}else{ord<-100}
  
  mulrep <- ord*sentineldata[,2]
  reporteddata <- data.frame(date=sentineldata[,1],reported=mulrep)
  
  rep_delay = list(shape = m, scale = e) 
  generated_onset_day_list <- list()
  
  for (d in 1:length(reporteddata$date)) {
    rep.day <-  reporteddata$date[d]
    n <-  reporteddata$reported[d]
    if (n > 0) {
      r <- rweibull(n, shape = rep_delay$shape, scale = rep_delay$scale)
      generated_onset_day_list <- c(generated_onset_day_list, rep.day-r)
    }
  }
  generated_onset_day_list <- as.Date(unlist(generated_onset_day_list))
  
  df_gen_onset <- as.data.frame(table(generated_onset_day_list))
  colnames(df_gen_onset) <- c("date", "onset")
  df_gen_onset$onset <- df_gen_onset$onset/ord # restore the original order 
  df_gen_onset %<>% mutate(df_gen_onset, date=as.Date(date))
  
  df_rep_ons <- merge(reporteddata, df_gen_onset, by="date", all.x=T)
  df_rep_ons [is.na(df_rep_ons)] = 0 # NA in column onset is replaced by 0
  df_rep_ons$onset <- replace(df_rep_ons$onset,nrow(df_rep_ons),df_rep_ons$onset[nrow(df_rep_ons)-1]) 
  
  return(df_rep_ons[,-2])
}



# 週の定点データを7で割って線形補間する関数（日次データ変換）# 2026/03/31追加
lerp_data <- function(timeseries){　
  
  xest <- seq(min(timeseries$date),max(timeseries$date),by=1)
  
  ipdata <- data.frame(
    date = xest,
    value = approx(timeseries$date, timeseries[,2]/7, xout = xest)$y
  )
  return(ipdata)
  
} # end of function



# 時系列データの信頼区間　- 2026/03/31追加
calc_ts_CI <- function(obs, pred, n_param = 3, conf_level = 0.95) {
  # 1. 残差の計算（対数スケールを想定）
  res <- obs - pred
  n <- length(res)
  
  # 2. 不偏分散 Ve の計算 (自由度: n - n_param)
  ve <- sum(res^2) / (n - n_param)
  
  # 3. 有効標本サイズ n_eff の計算 (ラグ1の自己相関rho1を使用)
  rho1 <- acf(res, plot = FALSE)$acf[2]
  n_eff <- n * (1 - rho1) / (1 + rho1)
  
  # 4. 有効自由度 (df_eff) の設定
  # 有効なデータ数からパラメータ数を引く
  df_eff <- n_eff - n_param
  
  # 安全策：自由度が1を下回らないように補正（計算エラー防止）
  if (df_eff < 1) df_eff <- 1
  
  # 5. t分布の分位点 (t-score) を取得
  t_score <- qt(1 - (1 - conf_level) / 2, df = df_eff)
  
  # # 6. 信頼区間 (CI) と 予測区間 (PI) の幅を計算（有効標本サイズを使用する場合）
  # ci_width <- t_score * sqrt(ve / n_eff)
  # pi_width <- t_score * sqrt( ve*(1+1/n_eff) ) 
  # # pi_width <- 1.96 * sqrt(ve) # サンプルサイズが大きい時
  
  # 6. 信頼区間 (CI) と 予測区間 (PI) の幅を計算（元の標本サイズを使用する場合）
  ci_width <- t_score * sqrt(ve / n) # 1.96 * sqrt(ve / n)
  pi_width <- t_score * sqrt( ve*(1+1/n) ) # 1.96 * sqrt(ve)
  # pi_width <- 1.96 * sqrt(ve) # サンプルサイズが大きい時
  
  # 結果をリストで返す
  return(list(
    ve = ve,
    n_eff = n_eff,
    rho1 = rho1,
    ci_width = ci_width,
    pi_width = pi_width,
    summary = paste0("元のn: ", n, " / 有効標本サイズ n_eff: ", round(n_eff, 2))
  ))
} # end of function 



# CI PI付き予測グラフ用データフレーム作成 added at 2026/03/24
plot_result_cipi <- function(dataframe, ci,pi){ # dataframe("date","predicted",ci,pi)
  
  #　注意：予測値の対数・非対数と合わせること
  
  date<- dataframe[,1]
  predicted <- dataframe[,2]
  
  ## 対数スケール用 ##
  log_df <-  data.frame(date=date, predicted.log=predicted, CI_lower=predicted-ci, CI_upper=predicted+ci, PI_lower=predicted-pi, PI_upper=predicted+pi)
  
  return( log_df )
  
}# end of function

# CI PI付き予測グラフ用データフレーム作成(非対数) added at 2026/03/24　　*入力は対数データ
plot_result_cipi_nonlog <- function(dataframe, ci,pi){ # dataframe("date","predicted",ci,pi)
  
  #　注意：予測値の対数・非対数と合わせること
  
  date<- dataframe[,1]
  predicted <- dataframe[,2]
  
    ## 非対数スケール用 ##
  df <-  data.frame(date=date, predicted=10^predicted, CI_lower=10^predicted * 10^(-ci), CI_upper=10^predicted * 10^(+ci), PI_lower=10^predicted * 10^(-pi), PI_upper=10^predicted * 10^(+pi))

  return( df )
  
}# end of function


# 日次予測感染者数を実測データ日付に合わせて集計する -- added at 2026/03/24
weekly_prediction <- function(weeklyobserved,dailypredicted){ 
  
  # predI_df: data.frame(date, predicted)  日次予測 # predI.period1
  predI_df <- data.frame(date=dailypredicted$date, predicted=dailypredicted$estimated_cases)
  # predI_aligned_df : data.frame(date)    週次の観測日 # teitendata
  predI_aligned_df <- data.frame(date=weeklyobserved$date) 
  
  # 日付型の確認
  predI_df$date <- as.Date(predI_df$date)
  predI_aligned_df$date  <- as.Date(predI_aligned_df$date)
  
  # 集計
  predI_aligned_df$pred_sum <- sapply(predI_aligned_df$date, function(d) {
    sum(
      predI_df$pred[
        predI_df$date >= (d - 6) & predI_df$date <= d
      ],
      na.rm = TRUE
    )
  })
  
  return(predI_aligned_df)
  
} # end of function


# R-squared in a simple linear regression model
Rsquared_in_lm<- function(xdata,ydata,pa,pb){
  
  xy.data <- merge(xdata, ydata, by=1)
  xy.data <- na.omit(xy.data)
  
  nn <- nrow(xy.data)
  
  colnames(xy.data) <- c("date","x","y")
  
  edat.log <- data.frame(lx=log10(xy.data$x), ly=log10(xy.data$y))
  
  
  lx <- log10(xdata[,2])
  
  # predicted y ... this prediction needs regression parameters estimated from other data.
  pre.y <- pa + pb*edat.log$lx
  
  
  ### R-squared ###
  
  # observed y
  obs.y <- edat.log$ly
  
  # mean y
  ymean <- mean(obs.y) # mean of y
  
  # sum of squared error = residual sum of squares
  SSE <- sum ( ( obs.y - pre.y )^2 )
  
  #  regression sum of squares
  SSR <- sum ( ( pre.y - ymean)^2 )
  
  # total sum of squared
  SST <- sum( (obs.y - ymean)^2 )
  
  # R-squared.1
  rsq1 <-   1-(SSE/SST)
  
  # R-squared.2
  rsq2 <-   SSR/SST
  
  
  # adjusted R-squared 
  #adj.rsq <- 1-(SSE/SST)*( (nn-1)/(nn-2) )
  
  #return( c(rsq1,rsq2) )
  return( c(SSE,SSR,SST,rsq1,rsq2) )
  
}


Rsquared <- function(observeddata, predicteddata){　
  
  # The number of data (sample size)
  # nn <- length(observeddata)
  
  if(length(observeddata) != length(predicteddata))
    stop("length mismatch")
  
  # mean of observed onset data  
  meano <-  mean( observeddata )
  
  # sse: sum of squared errors = rss
  #sse <- rss
  sse <- sum( ( observeddata - predicteddata )^2 ) 
  
  # ssr: sum of squares of predicted values around mean (valid as SSR only for OLS with intercept)
  ssr <- sum( (predicteddata - meano)^2 )
  
  # sst: sum of squares total
  sst <- sum( (observeddata - meano)^2 )
  
  if(sst == 0)
    stop("SST is zero: R-squared undefined")
  
  # R-squared.1
  rsq1 <- 1-(sse/sst)
  
  # adjusted R-squared.1 
  #adj.rsq1 <- 1-(sse/sst)*( (nn-1)/(nn-2) )
  
  # R-squared.2
  rsq2 <- ssr/sst
  
  #return(c(rsq1,adj.rsq1,rsq2,sse,ssr,sst))
  return(c(sse,ssr,sst,rsq1,rsq2))
  
  
} # end of function


############# パラメータ推定と予測関数 #############

# parameter estimation by linear regression model: updated on 2025/08/19
param_estim_by_LM <- function(xdata,ydata){
  
  
  xy.data <- merge(xdata, ydata, by=1) 
  xy.data <- na.omit(xy.data)
  
  colnames(xy.data) <- c("date","x","y")
  
  edat.log <- data.frame(lx=log10(xy.data$x), ly=log10(xy.data$y)) 
  tr.lm <- lm(edat.log$ly ~ edat.log$lx, data=edat.log)
  
  tr.a <- as.numeric(coef(tr.lm)[1]) # intercept
  tr.b <- as.numeric(coef(tr.lm)[2]) # slope
  
  est_value <- function(x){ 
    a <- tr.a # intercept
    b <- tr.b # slope
    
    y <- a + b*x
    return(y)
  }
  
  pre.y <- unlist( lapply(edat.log$lx,est_value) )
  obs.y <- edat.log$ly
  
  
  num.d <- length(pre.y) # the number of data
  if(num.d>=3){nn<-num.d}else{nn<-3} # degree of freedom > 0
  
  # t-distribution
  tval <- qt(0.025, nn-2 ,lower.tail = FALSE ) # t(N-2,alpha=0.025)
  
  obs.x <- edat.log$lx # observed x
  xmean <- mean(obs.x) # mean of x
  
  sxx <- sum( (obs.x - xmean)^2 )
  
  # residual sum of squares
  rss <- sum( (obs.y - pre.y)^2 )
  
  # unbiased variance
  uv <-  rss/(nn-2)
  
  ### R-squared ### 
  # mean y
  ymean <- mean(obs.y) # mean of y
  
  # sum of squared error = residual sum of squares
  SSE <- rss 
  
  #  regression sum of squares 
  SSR <- sum ( ( pre.y - ymean)^2 )
  
  # total sum of squared
  SST <- sum( (obs.y - ymean)^2 )
  
  # R squared.1
  rsq1 <-   1-(SSE/SST) 
  
  # R squared.2
  rsq2 <-   SSR/SST 
  
  # adjusted R-squared 
  # adj.rsq <- 1-(SSE/SST)*( (nn-1)/(nn-2) )　

  
  #return(c(tr.a,tr.b,tval,num.d,xmean,sxx,uv,rsq1,rsq2))
  return(c(tr.a,tr.b,tval,num.d,xmean,sxx,uv,SSE,SSR,SST,rsq1,rsq2))
  
} # end of function



# prediction of disease onset counts by linear regression model: updated on 2024/10/07
epi.prediction_by_LM <- function(xdata,pa,pb,tval,num.d,xmean,sxx,uv,rsq){
  
  lx <- log10(xdata[,2])
  
  # predicted y
  pred.y <- pa + pb*lx
  
  # confidence interval
  ci.minus <- pa + pb*lx - tval*sqrt( (1/num.d + ((lx-xmean)^2/sxx))*uv )
  ci.plus <- pa + pb*lx + tval*sqrt( (1/num.d + ((lx-xmean)^2/sxx))*uv )
  
  # prediction interval
  pi.minus <- pa + pb*lx - tval*sqrt( (1 + (1/num.d) + ((lx-xmean)^2/sxx))*uv )
  pi.plus <- pa + pb*lx + tval*sqrt( (1 + (1/num.d) + ((lx-xmean)^2/sxx))*uv )
  
  result.x <- data.frame(date=xdata[,1], xdata=lx, prediction.y=pred.y, ci.up=ci.plus, ci.lw=ci.minus, pi.up=pi.plus, pi.lw=pi.minus)
  
  return(result.x)
  
} # end of function



# I->V: parameters estimation by Shedding Profile Model (Forward Prediction) updated at 2026/03/24
param_estim_by_SPM <- function(epidemicdata,sewagedata){ # I:説明変数, V:目的変数
  
  # パラメータ推定の予測方向はI->V(Forward Prediction)で統一
  # 戻り値: v, m, sigma, ci_width, pi_width, R2 (予測Vと観測Vの決定係数)
  # 週次定点データを入力 → 7で割って線形補間 → 週次ウイルスデータと共通日付でマージ = 共通週次データ
  # ウイルス量予測関数畳み込みは日にち単位で行う
  
  ## 日付データ確認
  epidemicdata[,1] <- as.Date(epidemicdata[,1])
  sewagedata[,1]   <- as.Date(sewagedata[,1])
  
  ## NA削除
  epidemicdata <- na.omit(epidemicdata)
  sewagedata <- na.omit(sewagedata)
  
  ## 日次データ判定
  date_diffs <- diff( epidemicdata[,1] )
  median_diff <- as.numeric(median(date_diffs))
  
  ## weekly epidemic data must be changed to daily data ##
  if(median_diff>5)epidemicdata <- lerp_data(epidemicdata) # weekly --> daily
  
  colnames(epidemicdata) <- c("date","infected") # daily
  colnames(sewagedata) <- c("date","virus")
  
  mergeddata <- merge(epidemicdata,sewagedata,by="date") # データ数カウント用
  
  ## Shedding curve ###
  shedding_profile <- function(x, a, m, sigma) {
    return( ((10^a)*exp( -(x-m)^2/(2*sigma^2) )) ) # revised
  }
  
  pI <- 2/3 # fraction of symptomatic infections
  T <- nrow(epidemicdata) # daily data
  
  ## 予測ウイルス量 = 一人あたり排出ウイルス量 x 発症者数　（畳み込みは日にち単位）
  Vt_pred_profile <- function (t, a, m, sigma) {
    vprd <- sum( shedding_profile(t-(1:T), a, m, sigma)*epidemicdata$infected*(1/pI) ) 
    return(vprd)
  }
  
  ## 残差平方和（対数）SSE
  eval_F_log_resid <- function(z) {
    
    ## 予測ウイルス量データ（疫学データに基づく）
    pv <- sapply( 1:T, function(t){Vt_pred_profile(t, z[1],z[2],z[3])} ) # 予測ウイルスデータ
    p.virus <- data.frame(date=epidemicdata$date, pred.virus=pv)
    
    ## ウイルス量実測データとの統合
    virusdata <- merge( mergeddata, p.virus, by="date") #c("date","virus","pred.virus")
    
    ### predicted ###
    p <- virusdata$pred.virus
    p[p <= 0] <-  1e-6
    plog <- log10(p)
    
    ### observed ###
    q <- virusdata$virus
    q[q <= 0] <-  1e-6
    qlog <- log10(q)
    
    return ( sum( (plog - qlog)^2 ) ) 
  }
  
  ## パラメータ初期値
  vv <- max( sewagedata$virus )  # vv is a initial value
  z0 <- c(log10(vv),1,1) # c(a,m,sigma)
  
  ## Optimization using L-BFGS-B ###
  res_log_resid <- optim(
    z0,
    eval_F_log_resid,
    method = "L-BFGS-B",
    lower = c(1e-5, -10, 0.01),  # lower: a, m, sigma
    upper = c(10, 30, 100)  # upper: a, m, sigma
  )
  
  ## optimized parameters and predictor
  reg_log <- c(res_log_resid$par[1],res_log_resid$par[2],res_log_resid$par[3])  #c(a,m,sigma)
  
  ## predicted virus concentration
  pre.y <- sapply( 1:T, function(t){Vt_pred_profile(t, reg_log[1],reg_log[2],reg_log[3])} )
  # pre.logy <-log10(pre.y)
  
  ## predicted virus data
  predicteddata <- data.frame(date=epidemicdata$date, pred.virus=pre.y) #c("date","pred.virus")
  
  ## merge data by common date
  op.virus <- merge(sewagedata, predicteddata, by="date")　#c("date","virus","pred.virus")
  
  observed.v <- op.virus$virus
  predicted.v <- op.virus$pred.virus
  
  ## calculate CI and PI
  res_calcCI <- calc_ts_CI(log10(observed.v), log10(predicted.v)) 
  
  ## number of data (num.d)
  nn <- nrow(mergeddata)
  if(nn>=3){num.d<-nn}else{num.d<-3}
  
  ## truncating small values
  #trncdata <- op.virus[ ! (op.virus$virus<10 | op.virus$pred.virus<10), ]
  
  ## R-squared
  #R2 <- Rsquared(log10(op.virus$pred.virus),log10(op.virus$virus)) # ウイルス量予測精度としての
  R2 <- Rsquared(log10(op.virus$virus), log10(op.virus$pred.virus)) # ウイルス量予測精度としての
  
  result_list <- list(
    v     = 10^res_log_resid$par[1], # v = 10^a
    m     = res_log_resid$par[2],
    sigma = res_log_resid$par[3],
    # num = num.d, #　データ数（標本サイズ）
    # n_eff = res_calcCI$n_eff, #有効標本サイズ
    ci_width = res_calcCI$ci_width,
    pi_width = res_calcCI$pi_width,
    R2 =R2[4] # 1-(sse/sst)
  )
  
  output <- list(
    params = result_list, # 推定パラメータ
    predicted = predicteddata # 推定期間の予測ウイルス量
  )
  
  return(  output  )
  
} # end of function



# ### COVIVISでは使用しない ### 2026/03/24 追加
# I->V: prediction of virus concentration by Shedding Profile Model (Forward Prediction)
virus_prediction_by_SPM <- function(epidemicdata,v,m,sigma){
    
    # パラメータ変換
    a <- log10(v)
    
    ## 日次データかどうかの判定
    date_diffs <- diff( epidemicdata[,1] )
    median_diff <- as.numeric(median(date_diffs))
    
    ## weekly epidemic data must be changed to daily data ##
    if(median_diff>5)epidemicdata <- lerp_data(epidemicdata) # weekly --> daily
    
    colnames(epidemicdata) <- c("date","infected")
    
    ## Shedding curve ###
    shedding_profile <- function(x, a, m, sigma) {
      return( ((10^a)*exp( -(x-m)^2/(2*sigma^2) )) ) # revised
    }
    
    pI <- 2/3 # Fraction of symptomatic infections
    T <- nrow(epidemicdata) # daily data
    
    ## 予測ウイルス量 = 一人あたり排出ウイルス量 x 発症者数
    Vt_pred_profile <- function (t, a, m, sigma) {
      vprd <- sum( shedding_profile(t-(1:T), a, m, sigma)*epidemicdata$infected*(1/pI) ) 
      return(vprd)
    }
    
    pv_val <- sapply( 1:T, function(t){Vt_pred_profile(t, a, m, sigma)} ) # 非logスケールV予測値
    pv_val <- pmax(pv_val, 1e-6) # 0回避
    pv_log_val <- log10(pv_val) # logスケールV予測値
    
    ### 予測値 ###
    # return(data.frame(date=epidemicdata$date, pred.virus=pv_val)) # non-log
    return(data.frame(date=epidemicdata$date, pred.virus=pv_log_val)) # log
    
  } # end of function



# V->I: prediction of the number of infected cases by Shedding Profile Model (Backward Prediction)
### the latest version 2026-08-19 ###
infected_prediction_by_SPM <- function(sewagedata, est_v, est_m, est_sigma, lambda=10^2) {
  
    # Preprocessing data
    sewagedata <- na.omit(sewagedata)
    sewagedata <- sewagedata[order(sewagedata$date), ]
    offset_days <- floor(est_m)
    start_date <- min(sewagedata$date) - offset_days
    end_date <- max(sewagedata$date)
    rep_dates <- seq(start_date, end_date, by = "day")
    tV <- as.numeric(sewagedata$date - start_date)
    tR <- as.numeric(rep_dates - start_date)
    NV <- length(tV) # length of virus data
    NR <- length(tR) # length of reported cases data
    
    
    # Fraction of symptomatic infections
    pI <- 2/3 
    
    # Matrix G
    G <- matrix(0, nrow = NV, ncol = NR)
    for (i in 1:NV) {
      diff_t <- tV[i] - tR  # diff_t = (t^V_i - t^R_j)_{j=1,2,..., NR}
      z <- (diff_t - est_m) / est_sigma
      val <- est_v * exp(-(z^2 / 2)) 
      val[abs(z) > 10] <- 0
      G[i, ] <- val  # (G_{i,j})_{j=1,2,...,NR} = val
    }
    
    
    if (NR < 3) {
      stop("At least three time points are required.")
    }
    # (NR-2)xNR dimensional matrix of zero
    D2 <- matrix(0, NR - 2, NR)
    # 2階差分行列 D2 の作成
    for (i in 1:(NR - 2)) {
      D2[i, i]     <-  1
      D2[i, i + 1] <- -2
      D2[i, i + 2] <-  1
    }
    
    # V = log v: log scaled virus concentrations
    eps <- 1e-8 # <- set 0 to 1e-8 in virus concentration
    V <- log(sewagedata$virus + eps)
    
    # Objection function (X) -- X: control variable
    # X = log(x): log scaled reported cases
    # V = log(v): log scaled virus concentrations
    
    obj_fun <- function(X, G, V, D2, lambda) {
      
      
      x <- exp(X) # Given reported cases x = exp(X)
      # Predicted log scaled virus concentration \hat V_i at day t^V_i:
      # \hat V_i = log( \sum_j G_{ij} exp(X_j) )
      eps_hat <- 1e-12
      V_hat <- as.numeric( log (G %*% x + eps_hat))
      
      # Error sum of squares of the observed log scaled
      # concentrations V from the expected ones V_hat
      loss <- sum((V - V_hat)^2)
      
      # The penalty for the deviation of smoothness
      # of curve X (log scaled reported cases)
      penalty <- 0.5 * lambda * sum((D2 %*% X)^2)
      
      return (loss + penalty)
    }
    
    # Gradient of objective function (X)
    # X = log(x): log scaled reported cases
    # V = log(v): log scaled virus concentrations
    grad_fun <- function(X, G, V, D2, lambda) {
      
      
      x <- exp(X) # Reported cases in normal scale
      virus_hat <- as.numeric(G %*% x) # predicted virus conc
      
      eps_hat <- 1e-12
      V_hat <- as.numeric( log (virus_hat + eps_hat)) # log scaled
      
      diff <- V_hat - V
      grad_loss <- x * colSums(G * (2 * diff / (virus_hat + eps_hat)))
      
      grad_penalty <- lambda * as.numeric(t(D2) %*% (D2 %*% X))
      
      return (grad_loss + grad_penalty)
    }
    
    x_init <- rep(0.1, NR)
    X_init <- log(x_init) # Initial value of X
    
    # Optimization
    res <- nloptr(
      x0 = X_init, # vector of the initial values
      eval_f = obj_fun, # objective function
      eval_grad_f = grad_fun, # its gradient
      # x_i = exp(X_i) > 0 is always satisfied. Hence no lower bound necessary.
      # lb = rep(0, NR), 
      opts = list(
        "algorithm" = "NLOPT_LD_LBFGS",
        "xtol_rel" = 1.0e-6, # Changed (AS)
        "maxeval" = 5000
      ),
      G = G, V = V, D2 = D2, lambda = lambda
    ) 
    
    # nloptr() returns solution, objective, etc
    est_X <- res$solution # Optimized log scaled reported
    estimated_x <- exp( est_X )
    
    return(
      list(
        resI = data.frame(
          date = rep_dates,
          estimated_cases = estimated_x*pI  #予測感染者数*有症率
        ),
        resV = as.numeric(G %*% (estimated_x))
      )
    )
  }





