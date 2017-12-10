## ---------------------------------------------------------------------------
##
## Template for DM assignment Cognitive Modeling Course
##
## Hedderik van Rijn, 091119
##
## ---------------------------------------------------------------------------

##Noise
actr.noise <- function(s,n=1) {
  rand <- runif(n,min=0.0001,max=0.9999)
  s * log((1 - rand ) / rand)
}

#Pulse to time, t == start pulse time,
#a == a, b == b, n == number of pulses
pulse.to.time <-function(t,a,b, n) {
  tn <- t
  tt <- t
  pulse <- 1
  #while number of pulse not met, compute formula
  while(pulse < n){
    sd <-  b * a * tn
    noise <- actr.noise(sd)
    tn <- a * tn + noise
    tt <- tt + tn
    pulse <- pulse +1
  }
  #return total time needed for given pulses(n)
  return(tt)
}

#Time to pulse, t == start pulse time,
#a == a, b == b, n == time in miliseconds
time.to.pulse <-function(t,a,b, n) {
  tn <- t
  tt <- t
  pulse <- 1
  #while total time has not exceeded n, compute formula 
  while(tt <= n){
    sd <-  b * a * tn
    noise <- actr.noise(sd)
    tn <- a * tn + noise
    tt <- tt + tn
    pulse <- pulse +1
  }
  #return total pulses in given time(n)
  return(pulse)
}

## DM functions

create.dm <- function(chunks,encounters) {
  if (chunks > 52) {
    stop("Only up to 52 chunks allowed.")
  }
  DM <- array(NA,c(chunks,encounters))
  row.names(DM) <- c(1:chunks)
  DM
}

add.encounter <- function(DM,chunk,time) {
  tmp <- DM[chunk,]
  DM[chunk,sum(!is.na(tmp))+1] <- time
  DM
}

get.encounters <- function(DM,chunk) {
  tmp <- DM[chunk,]
  tmp[!is.na(tmp)]
}


actr.B <- function(encounters,curtime) {
  if (length(curtime)>1) {
    sapply(curtime,function(X) { actr.B(encounters,X)})
  } else {
    if (curtime < min(encounters)) {
      return(NA)
    } else {
      log(sum((curtime - encounters[encounters<curtime])^-params$d))	
    }
  }
}

actr.B.optimized <- function(n,Time,curtime=NULL) {
  if (is.null(curtime)) {
    log(n/(1-params$d)) - params$d * log(Time)
  } else {
    n <- curtime/Time * n
    log(n/(1-params$d)) - params$d * log(curtime)
  }
}

## ---------------------------------------------------------------------------

load("dataJS.Rdat");
DatJS <- with(datJS,aggregate(list(Tp=Tp),list(Ts=Ts,Cond=Cond),mean))
params <- list()
cond = c(1,2,3)
participants = c(1,2,3,4,5,6)

params$d <- .5
params$num.chunks <- 20
params$max.num.encounters <- 1500
params$curtime <- 0
params$addtime <- 800
params$trials <- 1
sum_y = 0

blend_e <- data.frame(Subj = numeric(0), Cond = numeric(0), Ts = numeric(0), Tp = numeric(0))

## List with parameter values:

create.train.DM <-function (cond,participant) {
  sub_DM <- subset(datJS, datJS$Cond == cond & datJS$Subj == participant)
  return(sub_DM)
} 

create.train.data <- function(cond, participant){
  sub_DM <<- create.train.DM(cond, participant)
  DM <<- create.dm(params$num.chunks,params$max.num.encounters)
  train_data = c()
  for (i in 1:20) {
    train_data = c(train_data, time.to.pulse(100,1.02,0.015,sub_DM$Ts[i]))
  }
  for(i in 1:length(train_data)) {
   DM <<- add.encounter(DM, train_data[i], params$curtime)
   params$curtime <<- params$curtime + params$addtime
  }
  return(DM)
}

blend.data <- function() {
  y <- c()
  z <- c()
  for(i in 1:params$num.chunks) {
    y <- c(y, actr.B(get.encounters(DM,i),params$curtime))
    if(is.na(y[i]) == FALSE) {
      sum_y = sum_y + y[i]
    }
  }
  for(i in 1:length(y)) {
      z = c(z, (y[i]/sum_y))
  }
    
  sum_z = 0
  for(i in 1:length(z)) {
     if(is.na(z[i]) == FALSE) {
       z[i] = (z[i] * i)
        sum_z = round(sum_z + (z[i]))
     }
  }
  return(y)
}


create.test.data <-function() {
  for(i in 1:6) {
    params$curtime <<- 0
    DM <<- create.train.data(1,i)
    for(j in 1:3) {
      blend_c = c()
      cat("Processing condition ", j, "of", "3", "of participant", i,"\n")
      subset_plot = c()
      subset = create.train.DM(j,i)
      for(x in 1:500) {

        DM <<- add.encounter(DM, time.to.pulse(100,1.02,0.015,subset$Ts[x]), params$curtime)
        params$curtime <<- params$curtime + params$addtime
        blend_c = c(blend_c, blend.data())
        newrow = data.frame(Subj = i, Cond =  j, Tp = pulse.to.time(100,1.02,0.015,blend_c[x]), Ts = subset$Ts[x])
        blend_e <<- rbind(blend_e, newrow)
      }
    }
  }
  plot.data(blend_e)
}


plot.data <- function(blend_e) {
  brown <- "#8b4513";
  red <- "#ff1100";
  black <- "#000000";
  brownT <- "#8b451322";
  redT <- "#ff110022";
  blackT <- "#00000022";
  
  
  datJS <- blend_e
  ## ---
  
  par(mfrow=c(1,1))
  
  plotDatJS <- with(datJS,aggregate(list(Tp=Tp),list(Ts=Ts,Cond=Cond),mean))
  
  yrange <- range(plotDatJS$Ts)*c(.95,1.05)
  
  with(plotDatJS[plotDatJS$Cond==3,],plot(Ts,Tp,type="b",col=red,lwd=2,ylim=yrange,xlim=yrange,main="J&S All"))
  with(plotDatJS[plotDatJS$Cond==2,],lines(Ts,Tp,type="b",col=brown,lwd=2,ylim=yrange,xlim=yrange))
  with(plotDatJS[plotDatJS$Cond==1,],lines(Ts,Tp,type="b",col=black,lwd=2,ylim=yrange,xlim=yrange))
  
  lines(c(yrange[1],yrange[2]),c(yrange[1],yrange[2]),col="darkgrey",lty=2)
  
  with(datJS[datJS$Cond==3,],points(jitter(Ts),Tp,col=redT,pch=".",cex=3))
  with(datJS[datJS$Cond==2,],points(jitter(Ts),Tp,col=brownT,pch=".",cex=3))
  with(datJS[datJS$Cond==1,],points(jitter(Ts),Tp,col=blackT,pch=".",cex=3))
  #legend("topleft", legend = c("Short", "Intermediate", "Long"), col=c("#ff110022", "#8b451322", "#00000022"),pch=".")
}

DM <- create.train.DM(1,1)
y = blend.data()
