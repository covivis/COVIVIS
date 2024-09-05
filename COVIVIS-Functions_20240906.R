### EpiPredFunctionsRevised20240906.R ###
## revised at 2024/09/06 ##

## required libraries 
## library("dplyr") # for mutate 
## library("magrittr") # for %<>%
## library(nleqslv)　# for nleqslv


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



# V-I and I-V parameters estimation by linear regression model: updated on 2024/05/22
param_estim_by_lm <- function(xdata,ydata){
  
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
  
  return(c(tr.a,tr.b,tval,num.d,xmean,sxx,uv))
  
}



# VI and IV predictions of disease onset counts by linear regression model: updated on 2024/05/22
epi.prediction_by_lm <- function(xdata,pa,pb,tval,num.d,xmean,sxx,uv){
  
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
  
}


# V-I parameters estimation by Shedding profile model (epidemic model)
vi.param_estim_by_em <- function(sewagedata,onsetdata){
  
  mi <- min(onsetdata$date)
  ma <- max(onsetdata$date)
  
  medata <- merge(onsetdata,sewagedata,by=1, all=T) # daily data with NA
  colnames(medata) <- c("date","onset","sewage") 
  xy.data <- medata[(mi<=medata$date & medata$date<=ma),]
  
  sval <- xy.data$sewage
  oval <- xy.data$onset
  
  T <- nrow(xy.data)
  
  Ot_pred <- function (TT, v, omega, gamma) { 
    
    pI <- 2/3 # Fraction of symptomatic infections
    
    xlist <- c()
    t <- 1
    
    while(t<TT){
      xlist <- c( xlist, exp(-gamma*((t-1):0) ) ) # amount of virus shedding after disease onset
      xlist <- c( xlist, exp(-omega*(1:(TT-t)) )) # amount of virus shedding before disease onset
      
      t <- t+1
    }
    xlist <- c( xlist, exp(-gamma*((TT-1):0) ) )  
    
    shedd.mat <- matrix(xlist,TT,TT,byrow = T)
    
    observed.v <- sval
    
    if( anyNA(observed.v) ){ 
      na.val <- which(is.na(observed.v)) 
      shedd.mat <- shedd.mat[-c(na.val),-c(na.val)]
    }
    
    vna <-observed.v[!is.na(observed.v)]
    
    solve.ons <- solve(shedd.mat,(pI/v)*vna) 
    
    return (solve.ons)
    
  }
  
  
  eval_F_log_resid_Ot <- function(z) {
    
    p <- Ot_pred(T,z[1],z[2],z[3]) #predicted onset 
    if( anyNA(sval) ){ 
      q <- oval[-which(is.na(sval))]
    }else{q <- oval} #observed onset 
    
    p <- replace(p,which(p<=0),1) 
    plog <- log10(p)
    qlog <- log10(q)
    
    return ( sum( (plog - qlog)^2 ) ) 
  }
  
  # an initial value of v in z0
  mendata <- merge(onsetdata,sewagedata,by=1) # daily data without NA
  smax <- max(mendata$sewage)
  ons_of_the_day <- mendata[mendata$sewage==smax,]$onset
  vv <- smax/ons_of_the_day
  
  # minimizing with log transformed variables
  z0 <- c(vv, 1, 1) #(v, omega, gamma)
  res_log_resid_Ot <- optim(z0, eval_F_log_resid_Ot)
  
  # optimized parameters and predictor
  reg_log <- c(res_log_resid_Ot$par[1],res_log_resid_Ot$par[2],res_log_resid_Ot$par[3])
  
  # the number of data used with calculations
  numofdata <- T-sum(is.na(sval)) #the number of data without NA
  
  # residual sum of squares
  rss <- eval_F_log_resid_Ot(reg_log)
  
  # residual variance
  myv <-  rss/(numofdata-2)
  
  # residual standard deviation
  rsd <- (myv)^(1/2)
  
  reg_log2 <- list(
    v     = res_log_resid_Ot$par[1],
    omega = res_log_resid_Ot$par[2],
    gamma = res_log_resid_Ot$par[3],
    rsd = rsd
  )
  
  return(reg_log2)
  
}


# V-I prediction of disease onset counts by Shedding profile model (epidemic model)
vi.prediction_by_em <- function(sewagedata,v,omega,gamma,pr){
  
  dmin <- min(sewagedata$date)
  dmax <- max(sewagedata$date)
  
  prd <- as.integer(dmax-dmin) # period
  daydata <- as.Date( sapply(c(0:prd),function(d){ dmin+d }) )
  vdatawna <- merge(daydata,sewagedata,by=1,all=T) # virus data with NA
  colnames(vdatawna) <- c("date","sewage")
  
  sday <- vdatawna$date
  sval <- vdatawna$sewage # with NA
  
  nd <- length(sval)
  
  Ot_pred_new <- function (TT, v, omega, gamma) { 
    
    pI <- 2/3 # Fraction of symptomatic infections
    
    xlist <- c()
    t <- 1
    
    while(t<TT){
      xlist <- c( xlist, exp(-gamma*((t-1):0) ) ) # amount of virus shedding after disease onset
      xlist <- c( xlist, exp(-omega*(1:(TT-t)) )) # amount of virus shedding before disease onset
      
      t <- t+1
    }
    xlist <- c( xlist, exp(-gamma*((TT-1):0) ) )  
    
    shedd.mat <- matrix(xlist,TT,TT,byrow = T)
    
    observed.v <- sval
    
    
    if( anyNA(observed.v) ){ 
      na.val <- which(is.na(observed.v)) 
      shedd.mat <- shedd.mat[-c(na.val),-c(na.val)]
    }
    
    vna <-observed.v[!is.na(observed.v)]
    
    solve.ons <- solve(shedd.mat,(pI/v)*vna) 
    
    return (solve.ons)
    
  }# end of Ot_pred_new
  
  prediction.y <- log10( Ot_pred_new(nd, v, omega, gamma) ) # prediction onset
  
  # confidence interval
  ci.up <- prediction.y + 1.96*pr*( 1/length(prediction.y) )^(1/2) 
  ci.lw <- prediction.y - 1.96*pr*( 1/length(prediction.y) )^(1/2)
  
  # prediction interval
  pi.up <- prediction.y + 1.96*pr*( 1+( 1/length(prediction.y)) )^(1/2) 
  pi.lw <- prediction.y - 1.96*pr*( 1+( 1/length(prediction.y)) )^(1/2)
  
  pred.data <- data.frame("date"=na.omit(vdatawna)$date,"sewage"=log10(na.omit(vdatawna)$sewage),prediction.y, ci.up, ci.lw, pi.up, pi.lw)
  
  return(pred.data)
  
} # end of function


# I-V parameters estimation by Shedding profile model (epidemic model)
iv.param_estim_by_em <- function(onsetdata,sewagedata){
  
  colnames(onsetdata) <- c("date", "onset") 
  colnames(sewagedata) <- c("date", "sewage") 
  
  pI <- 2/3
  tmax <- nrow(onsetdata)
  T <- tmax
  
  n_onsets <- onsetdata$onset
  
  Vt_pred <- function (t, v, omega, gamma) {
    Pvec <- exp( - omega*(1:(T-t)))
    Qvec <- exp( - gamma*(0:(t-1)))
    ret <-       (v/pI)        * sum( n_onsets[(t+1):T] * Pvec )
    ret <- ret + (v/pI)        * sum( n_onsets[t:1]     * Qvec )
    return (ret)
  }
  
  sdate <- sewagedata$date #weekly
  odate <- onsetdata$date #daily
  comdate <- which( odate %in% sdate )
  
  
  # Objective function
  eval_F_log_resid <- function(z) {
    
    vlist <- list()
    for (t in 1:(tmax-1))
      vlist[[t]] <- Vt_pred(t, z[1], z[2], z[3])
    vlist <- unlist(vlist)
    
    xandy <- na.omit( data.frame(x=vlist[comdate],y=sewagedata$sewage) )
    
    xlog <- log10(xandy$x) # predicted virus concentrations
    ylog <- log10(xandy$y) # observed virus concentrations
    
    return ( sum( (xlog - ylog)^2, na.rm=TRUE ) )
  }
  
  
  # an initial value of v in z0
  mendata <- merge(onsetdata,sewagedata,by=1) # daily data without NA
  smax <- max(mendata$sewage)
  ons_of_the_day <- mendata[mendata$sewage==smax,]$onset
  vv <- smax/ons_of_the_day
  
  # Optimization  -- minimizing with log transformed variables 
  z0 <- c(vv, 1, 1) #(v, omega, gamma)
  res_log_resid <- optim(z0, eval_F_log_resid)
  
  reg_log <- c(res_log_resid$par[1], res_log_resid$par[2], res_log_resid$par[3])
  
  
  # residual sum of squares
  rss <- eval_F_log_resid(reg_log)
  numofdata <- length(comdate)-sum(is.na(sewagedata$sewage)) #the number of data without NA
  
  # residual variance
  myv <-  rss/(numofdata-2)
  
  # residual standard deviation
  rsd <- (myv)^(1/2)
  
  
  reg_log2 <- list(
    v     = res_log_resid$par[1], 
    omega = res_log_resid$par[2],
    gamma = res_log_resid$par[3],
    rsd = rsd
  )
  
  return(reg_log2)
  
}


# I-V prediction of virus concentrations from cases with disease onset by Shedding profile model (epidemic model)
iv.prediction_by_em <- function(onsetdata,v,omega,gamma,pr){
  
  colnames(onsetdata) <- c("date", "onset") 
  
  pI <- 2/3
  tmax <- nrow(onsetdata)
  T <- tmax
  
  n_onsets <- onsetdata$onset # a list of the number of onsets
  odate <- onsetdata$date #daily # a list of date
  
  Vt_pred <- function (t, v, omega, gamma) {
    Pvec <- exp( - omega*(1:(T-t)))
    Qvec <- exp( - gamma*(0:(t-1)))
    ret <-       (v/pI)        * sum( n_onsets[(t+1):T] * Pvec ) 
    ret <- ret + (v/pI)        * sum( n_onsets[t:1]     * Qvec )
    return (ret)
  }
  
  nolog_prediction.v <- sapply(1:tmax,v=v, omega=omega, gamma=gamma, Vt_pred)
  prediction.v <- log10(nolog_prediction.v)
  
  
  # confidence interval
  ci.up <- prediction.v + 1.96*pr*( 1/length(prediction.v) )^(1/2) 
  ci.lw <- prediction.v - 1.96*pr*( 1/length(prediction.v) )^(1/2)
  
  # prediction interval
  pi.up <- prediction.v + 1.96*pr*( 1+( 1/length(prediction.v)) )^(1/2) 
  pi.lw <- prediction.v - 1.96*pr*( 1+( 1/length(prediction.v)) )^(1/2)
  
  logons <- log10(onsetdata$onset)
  logons.na <- replace(logons, which(is.infinite(logons) & (logons < 0)), 0)
  
  pred.data <- na.omit( data.frame("date"=onsetdata$date,"onset"=logons.na,"predicted.v"=prediction.v, ci.up, ci.lw, pi.up, pi.lw) )
  
  return(pred.data)
  
} 
