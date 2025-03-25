library("ggplot2")
library("patchwork")
library("earlywarnings")
library("ggpubr")
library("RColorBrewer")
library("Bhat")
library("doParallel")
library("doSNOW")
library("changepoint")
library('rgl')
library('car')
library('plotly')
library('akima')
library('mgcv')
library('deldir')
library('reticulate')
library('webshot')
library('webshot2')
library('png')
library('grid')
library('gridExtra')
library("trend")



###funs########
generatets <- function(eq = 0,
                       ntotal = 5000, 
                       rr,
                       wnsd = 0.5,
                       disfre = 0.01,
                       disintense = c(-0.5),
                       occurange = c(0.4,0.6),
                       return.tsdis = F){
  
  delta0 = 0.01
  
  ts.out = rep(NA, ntotal)
  
  ts.wn = rnorm(n = ntotal,mean = 0,sd = wnsd)
  disloc =  sort(sample(seq(from = round(occurange*ntotal)[1],
                            to = round(occurange*ntotal)[2],
                            by = 1),round(disfre*ntotal)))
  
  ts.dis = rep(0,ntotal)
  ts.dis[disloc] = sample(disintense, length(disloc), replace = T)
  
  if(length(rr) == 1){
    ts.rr = rep(rr,ntotal)
  }else{
    ts.rr = seq(from = rr[1], to =rr[2], length.out = ntotal)
  }
  
  t.before = eq
  
  for (i in 1:ntotal) {
    
    ti = t.before + ts.rr[i]*(eq - t.before) * delta0 +
      ts.wn[i] * sqrt(delta0) + ts.dis[i] * sqrt(delta0)
    
    ts.out[i] = ti
    t.before = ti
  }
  ts.out = ts(ts.out,start = 1)
  
  if(return.tsdis){
    return(c(list(ts.out),list(ts.dis)))
  }else{
    return(ts.out)
  }
}


generatets.multi <- function(n.ts,
                             eq = 0,
                             ntotal = 5000, 
                             rr,
                             wnsd = 0.5,
                             disfre = 0.01,
                             disintense = c(-0.5),
                             occurange = c(0.4,0.6)){
  
  result.matrix = matrix(NA, nrow = n.ts, ncol = ntotal)
  for (ri in 1:n.ts) {
    result.matrix[ri,] = generatets(eq =eq,ntotal = ntotal,
                                    rr = rr, wnsd = wnsd,
                                    disfre = disfre,
                                    disintense = disintense,
                                    occurange = occurange)
  }
  
  return(result.matrix)
  
}




generatets.ass <- function(
    ntotal = 500,
    k=10,
    r=1,
    vgh=1,
    gi=c(1,3),
    nsd =0.2,
    disfre = 0.05,
    disintense = c(-0.2)*3,
    occurange = c(0.35,0.65)){
  
  
  delta0 = 0.01
  
  ts.out = rep(NA, ntotal)
  
  ts.wn = rnorm(n = ntotal,mean = 0,sd = nsd)
  
  ts.gi = seq(gi[1],gi[2],length.out = ntotal)
  
  
  disloc =  sort(sample(seq(from = round(occurange*ntotal)[1],
                            to = round(occurange*ntotal)[2],
                            by = 1),round(disfre*ntotal)))
  
  ts.dis = rep(0,ntotal)
  
  ts.dis[disloc] = sample(disintense, length(disloc), replace = T)
  
  
  
  v0 = k-gi[1]
  
  for (i in 1:ntotal) {
    v1 = v0 + delta0 *(r*v0*(1-v0/k)-ts.gi[i]*v0^2/(v0^2+vgh^2))+
      sqrt(delta0) * ts.wn[i] + sqrt(delta0) * ts.dis[i]
    
    if(v1<0){
      v1 = 0
    }
    ts.out[i] = v1
    v0 = v1
  }
  ts.out = ts(ts.out,start = 1)
  return(ts.out)
}


generatets.ass.multi <- function(n.ts,
                                 ntotal = 500,
                                 k=10,
                                 r=1,
                                 vgh=1,
                                 gi=c(1,3),
                                 nsd =0.2,
                                 disfre = 0.05,
                                 disintense = c(-0.2)*3,
                                 occurange = c(0.35,0.65)){
  
  result.matrix = matrix(NA, nrow = n.ts, ncol = ntotal)
  for (ri in 1:n.ts) {
    result.matrix[ri,] = generatets.ass(ntotal = ntotal,
                                        k=k,
                                        r=r,
                                        vgh=vgh,
                                        gi=gi,
                                        nsd =nsd,
                                        disfre = disfre,
                                        disintense = disintense,
                                        occurange = occurange)
  }
  return(result.matrix)
}




AR1.bts <- function(ts0, max.bp = 20, min.length = 0.05){
  
  cpt.p = cpt.mean(ts0,penalty= "MBIC",method="PELT",
                   minseglen = round(min.length*length(ts0)))
  
  attr.cpt =  attributes(cpt.p)
  bp.locs = attr.cpt$cpts
  bp.locs = bp.locs[-(length(bp.locs))]
  
  stage.means =  attr.cpt$param.est$mean
  rank.bp = order(abs(diff(stage.means)),decreasing = T)
  if(length(bp.locs)>max.bp){
    bp.locs = bp.locs[which(rank.bp<=max.bp)]
  }
  bp.locs = c(1,bp.locs,length(ts0))
  
  df.cen = c()
  for (i in 1:(length(bp.locs)-1)) {
    ts.seg = ts0[bp.locs[i]:bp.locs[i+1]]
    
    # det.df = data.frame(x = 1:length(ts.seg),y = as.vector(ts.seg))
    # det.mod = lm(y~x,data = det.df)
    # ts.seg = as.vector(det.mod$residuals)
    
    x0 = ts.seg[1:(length(ts.seg)-1)]
    y0 = ts.seg[2:length(ts.seg)]
    segorder = rep(i,length(x0))
    x0 = (x0-mean(ts.seg))/(sd(ts.seg))
    y0 = (y0-mean(ts.seg))/(sd(ts.seg))
    df0 = data.frame(x0,y0,segorder)
  }
  mod.lm = lm(y0 ~ x0,df.cen)
  result0 = mod.lm$coefficients[["x0"]]
  return(result0)
}


gentews_noplot <- function (timeseries) {
  timeseries <- data.matrix(timeseries)
  Y <- timeseries
  timeindex <- length(Y)
  
  nsmY <- Y
  mw <- length(Y)
  omw <- length(nsmY) - mw + 1
  low <- 2
  high <- mw
  
  nMR <- Y
  
  nSD <- sd(nMR, na.rm = T)
  
  nYR <- ar.ols(nMR, aic = FALSE, order.max = 1, dmean = TRUE, 
                intercept = FALSE)
  nARR <- nYR$ar
  nSK <- abs(moments::skewness(nMR, na.rm = TRUE))
  nKURT <- moments::kurtosis(nMR, na.rm = TRUE)
  nCV <- nSD/mean(nMR)
  ACF <- acf(nMR, lag.max = 1, type = c("correlation"), 
             plot = FALSE)
  nACF <- ACF$acf[2]
  spectfft <- spec.ar(nMR, n.freq = length(Y)+1, plot = F, 
                      order = 1)
  nSPECT <- spectfft$spec
  nDENSITYRATIO <- log10(spectfft$spec[low]/spectfft$spec[high])
  
  nRETURNRATE = 1/nARR
  
  AR1_bts = as.numeric(AR1.bts(ts0 = as.vector(timeseries)))
  
  out <- data.frame(timeindex, nARR, nSD, 
                    nSK, nKURT, nCV, nRETURNRATE, 
                    nDENSITYRATIO, nACF, AR1_bts)
  colnames(out) <- c("timeindex", "AR1", "SD", 
                     "SKE", "KURT", "cv", "returnrate", 
                     "SDR", "acf1", "AR1_bts")
  return(out)
}


Bh.co <- function(x1,x2) {
  
  np = round(diff(range(c(x1,x2)))*20/sd(c(x1,x2)))
  
  presta = hist(c(x1,x2),breaks = np, plot = F)
  
  box1 = hist(x1, breaks = presta$breaks, plot = F)
  box2 = hist(x2, breaks = presta$breaks, plot = F)
  
  bhco = sum(sqrt(box1$density*box2$density))*diff(presta$breaks[1:2])
  return(bhco)
}


sen_mk <- function(tt){
  if(is.na(mean(tt))){
    result2 = c(NA,NA,NA)
  }else{
    result1 = sens.slope(tt)
    result2 = c(result1$estimates,result1$statistic,result1$p.value)
  }
  names(result2)<-c("SenSlope","z","p_value")
  return(result2)
}






##################

setwd("E:/ESM_Research/IceCoreTest")
icecoredata00 = read.csv("NGRIPdata.csv",header = T)



###############################

icecoredata = icecoredata00[1:3000,]

ggplot(data = icecoredata,mapping = aes(x = Age, y = NGRIP_d18O))+
  geom_line()

datalength = nrow(icecoredata)
boxsize = 1000
stepsize = 20

irtimes = floor((datalength-boxsize)/stepsize)

resultform = c()

for (i in 1:irtimes) {
  s1 = 1+(i-1)*stepsize
  s2 = (i-1)*stepsize+boxsize
  
  subdata = icecoredata[s1:s2,]
  
  
  
  det.mod = lm(NGRIP_d18O~Age,data = subdata)
  subdata$NGRIP_d18O = as.vector(det.mod$residuals)
  
  
  result0 = gentews_noplot(subdata$NGRIP_d18O)
  
  result0 = cbind(Age = subdata$Age[1],Age_end = max(subdata$Age),result0)
  resultform = rbind(resultform,result0)
}


plotmargin = margin(.5,.2,.2,.3,unit = 'cm')
shifttext = 'Glacial-Holocene Transition'


pic.ts = ggplot(data = icecoredata,mapping = aes(x = Age, y = NGRIP_d18O))+
  geom_line(color = "gray35",linewidth = 0.2)+
  geom_vline(xintercept = 11500, color = "black",
             linetype='dashed',linewidth = 1)+
  annotate(geom = 'text',x = 13000,y = -35,label = shifttext,hjust = 1)+
  scale_x_reverse()+
  labs(tag = letters[1])+
  xlab("Age (years before 2000)")+
  ylab(expression('NGRIP'~' '~ delta^18~'O'))+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        #plot.tag.position = c(0.05,1.05),
        #plot.margin = plotmargin,
        plot.title = element_text(size = 9),
        axis.line = element_blank(),
        #axis.ticks = element_blank(),
        plot.tag = element_text(face = "bold"),
        axis.title = element_text(color="black",size = 10))



regressionform = resultform[resultform$Age>=11500,]

mm1 = lm(AR1~Age,regressionform)
mm2 = lm(AR1_bts~Age,regressionform)

sen_mk(regressionform$AR1)
sen_mk(regressionform$AR1_bts)

locy.ar1 = max(resultform$AR1)-0.85*diff(range(resultform$AR1))

locy.ar1bts = max(resultform$AR1_bts)-0.85*diff(range(resultform$AR1_bts))

pic.ar1 = ggplot(data = resultform)+
  geom_line(aes(x = Age, y = AR1),color = "#D98220",linewidth = 0.6)+
  geom_vline(xintercept = 11500, color = "black",
             linetype='dashed',linewidth = 1)+
  geom_smooth(formula = 'y ~ x',data = regressionform, 
              mapping = aes(x = Age, y = AR1), 
              method = "lm",se = F,color = 'gray40',
              linewidth = 1)+
#  annotate(geom = 'text',x = 40000,y = locy.ar1,label = "p<0.001",hjust = 0)+
  scale_x_reverse()+
  labs(tag = letters[2])+
  xlab("Age (years before 2000)")+
  ylab("AR1")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        #plot.tag.position = c(0.05,1.05),
        #plot.margin = plotmargin,
        plot.title = element_text(size = 9),
        axis.line = element_blank(),
        #axis.ticks = element_blank(),
        plot.tag = element_text(face = "bold"),
        axis.title = element_text(color="black",size = 10))








pic.ar1bts = ggplot(data = resultform)+
  geom_line(aes(x = Age, y = AR1_bts),color = "#187596",linewidth = 0.6)+
  geom_vline(xintercept = 11500, color = "black",
             linetype='dashed',linewidth = 1)+
  geom_smooth(formula = 'y ~ x',data = regressionform, 
              mapping = aes(x = Age, y = AR1_bts), 
              method = "lm",se = F,color = 'gray40',
              linewidth = 1)+
#  annotate(geom = 'text',x = 40000,y = locy.ar1bts,label = "p<0.001",hjust = 0)+
  scale_x_reverse()+
  labs(tag = letters[3])+
  xlab("Age (years before 2000)")+
  ylab("AR1_bts")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        #plot.tag.position = c(0.05,1.05),
        #plot.margin = plotmargin,
        plot.title = element_text(size = 9),
        axis.line = element_blank(),
        #axis.ticks = element_blank(),
        plot.tag = element_text(face = "bold"),
        axis.title = element_text(color="black",size = 10))


pic.group.clipts = pic.ts/(pic.ar1|pic.ar1bts)

pic.group.clipts


outfile = paste0("IceCoreClipts.jpg")
jpeg(outfile, units="cm",
     width= 20,
     height= 11,
     res = 600)

plot.new()

plot(pic.group.clipts)

dev.off()

