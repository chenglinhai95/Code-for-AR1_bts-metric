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




AR1.bts <- function(ts0, max.bp = 5, min.length = 0.1){
  
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





###generate simulate dataset############

setwd("E:/ESM_Research/result")

rrs = seq(from = 0.2, to = 0.8, by = 0.2)
gis = seq(from = 0.5, to = 2, by = 0.5)

fre.range = seq(from = 0, to = 0.015, by = 0.005)
dis.intenf = c(1,2,3)

#T1 = Sys.time()

df.para = c()
for (rr0 in rrs) {
  for (fre0 in fre.range) {
    if (fre0 == 0){
      df.para = rbind(df.para, data.frame(rr0,fre0,intf0 = 0))
      next
    }
    for (intf0 in dis.intenf) {
      df.para = rbind(df.para, data.frame(rr0,fre0,intf0))
    }
  }
}



df.para.ass = c()
for (gi0 in gis) {
  for (fre0 in fre.range) {
    if (fre0 == 0){
      df.para.ass = rbind(df.para.ass, data.frame(gi0,fre0,intf0 = 0))
      next
    }
    for (intf0 in dis.intenf) {
      df.para.ass = rbind(df.para.ass, data.frame(gi0,fre0,intf0))
    }
  }
}



# T2 = Sys.time()
# print(T2-T1)
# print("Generating Dataset complete!")

nseries = nrow(df.para)


cl <- makeCluster(15, type = "SOCK")
registerDoSNOW(cl)

pb = txtProgressBar(max = nseries,style = 3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress = progress)

T1 = Sys.time()

df.simu <- foreach(x = 1:nseries,
                   .combine='rbind',
                   .inorder = TRUE,
                   .packages= c('earlywarnings','changepoint'),
                   .options.snow = opts,
                   .errorhandling = "pass") %dopar%
  {
    rr0 = df.para$rr0[x]
    fre0 = df.para$fre0[x]
    intf0 = df.para$intf0[x]
    mat.ts = generatets.multi(n.ts = 500,
                              eq = 0,
                              ntotal = 50000,
                              rr = rr0,
                              wnsd = 0.2,
                              disfre = fre0,
                              disintense = c(-0.2)*intf0,
                              occurange = c(0.4,0.6))
    simuname = paste0("E:/ESM_Research/data/simusub_",
                      sprintf("%02d", round(rr0*100)),
                      sprintf("%02d", round(fre0*100)),
                      sprintf("%02d", round(intf0*10)),".rda")
    
    saveRDS(mat.ts, file = simuname)

    mat.ts = readRDS(file = simuname)

    df.res = data.table::rbindlist(apply(mat.ts,1,gentews_noplot))

    dis.fre = rep(fre0,nrow(df.res))
    dis.inten = rep(intf0,nrow(df.res))
    rec.rate = rep(rr0,nrow(df.res))

    df.res = cbind(df.res,rec.rate,dis.fre,dis.inten)
    rownames(df.res) = NULL

    calcname = paste0("E:/ESM_Research/data/ESM_",
                      sprintf("%02d", round(rr0*100)),
                      sprintf("%02d", round(fre0*100)),
                      sprintf("%02d", round(intf0*10)),".rda")
    saveRDS(df.res, file = calcname)

    return(calcname)
  }

stopCluster(cl)

T2 = Sys.time()
print(T2-T1)
print("ESM calculation complete!")



######


nseries = nrow(df.para.ass)


cl <- makeCluster(15, type = "SOCK")
registerDoSNOW(cl)

pb = txtProgressBar(max = nseries,style = 3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress = progress)

T1 = Sys.time()

df.simu <- foreach(x = 1:nseries,
                   .combine='rbind',
                   .inorder = TRUE,
                   .packages= c('earlywarnings','changepoint'),
                   .options.snow = opts,
                   .errorhandling = "pass") %dopar%
  {
    gi0 = df.para.ass$gi[x]
    fre0 = df.para.ass$fre0[x]
    intf0 = df.para.ass$intf0[x]
    mat.ts = generatets.ass.multi(n.ts = 500,
                                  ntotal = 50000,
                                  k=10,
                                  r=1,
                                  vgh=1,
                                  gi=c(gi0,gi0),
                                  nsd =0.2,
                                  disfre = fre0,
                                  disintense = c(-0.2)*intf0,
                                  occurange = c(0.4,0.6))
    
    plot(ts(mat.ts[10,]))
    
    simuname = paste0("E:/ESM_Research/data/simusub_ass_",
                      sprintf("%02d", round(gi0*10)),
                      sprintf("%02d", round(fre0*100)),
                      sprintf("%02d", round(intf0*10)),".rda")
    saveRDS(mat.ts, file = simuname)
    
    mat.ts = readRDS(file = simuname)
    
    df.res = data.table::rbindlist(apply(mat.ts,1,gentews_noplot))
    
    dis.fre = rep(fre0,nrow(df.res))
    dis.inten = rep(intf0,nrow(df.res))
    graz.int = rep(gi0,nrow(df.res))
    
    df.res = cbind(df.res,graz.int,dis.fre,dis.inten)
    rownames(df.res) = NULL
    
    calcname = paste0("E:/ESM_Research/data/ESM_ass_",
                      sprintf("%02d", round(gi0*10)),
                      sprintf("%02d", round(fre0*100)),
                      sprintf("%02d", round(intf0*10)),".rda")
    saveRDS(df.res, file = calcname)
    
    return(calcname)
  }

stopCluster(cl)

T2 = Sys.time()
print(T2-T1)
print("ESM calculation complete!")



#saveRDS(df.simu, file = "E:/ESM_Research/data/df_simulate.rda")
#df.simu = readRDS(file = "E:/ESM_Research/data/df_simulate.rda")

df.esmsimu = c()
for (x in 1:nrow(df.para)) {
  rr0 = df.para$rr0[x]
  fre0 = df.para$fre0[x]
  intf0 = df.para$intf0[x]
  calcname = paste0("E:/ESM_Research/data/ESM_",
                    sprintf("%02d", round(rr0*100)),
                    sprintf("%02d", round(fre0*100)),
                    sprintf("%02d", round(intf0*10)),".rda")
  
  df.res = readRDS(file = calcname)
  df.esmsimu = rbind(df.esmsimu,df.res)
}


df.esmsimu.ass = c()
for (x in 1:nrow(df.para.ass)) {
  gi0 = df.para.ass$gi0[x]
  fre0 = df.para.ass$fre0[x]
  intf0 = df.para.ass$intf0[x]
  calcname = paste0("E:/ESM_Research/data/ESM_ass_",
                    sprintf("%02d", round(gi0*10)),
                    sprintf("%02d", round(fre0*100)),
                    sprintf("%02d", round(intf0*10)),".rda")
  
  df.res = readRDS(file = calcname)
  df.esmsimu.ass = rbind(df.esmsimu.ass, df.res)
}




###ana2_ESMability########################################

setwd("E:/ESM_Research/result")

df.result = df.esmsimu[df.esmsimu$dis.fre == 0,]
varnameshorts = c("AR1","SD","SKE","KURT","AR1_bts")
disnames = c("AR1","Standard Deviation","Skewness","Kurtosis","AR1_bts")

clrs = brewer.pal(length(varnameshorts),name = "Dark2")

df.result$rec.rate = as.factor(df.result$rec.rate)


pic.list = list()

for (iii in 1:length(varnameshorts)) {
  varshort = varnameshorts[iii]
  disvar = disnames[iii]
  clr0 = clrs[iii]
  
  x = df.result[["rec.rate"]]
  y = df.result[[varshort]]
  df0 = data.frame(x,y)
  
  x.breaks = sort(unique(x))
  n.notes = length(x.breaks)-1
  
  
  com.list = list()
  bhco.strs = c()
  
  for (nn in 1:n.notes) {
    com.list = c(com.list,list(paste0(c(x.breaks[nn],x.breaks[nn+1]))))
    
    yy1 = df0$y[df0$x == x.breaks[nn]]
    yy2 = df0$y[df0$x == x.breaks[nn+1]]
    bhco.str = sprintf("%.2f", Bh.co(yy1,yy2)) 
    bhco.strs = c(bhco.strs,bhco.str)
    
  }
  
  range.ybox = range(y)

  inter.y = diff(range.ybox)*0.1
  
  labloc.x = rev(seq(from = 1.5, by = 1, length.out = n.notes))
  labloc.y = rev(seq(from = range.ybox[2], by = inter.y, length.out = n.notes))
  
  df.bhco = data.frame(labloc.x, labloc.y, bhco.strs)
  
  plotrange.y = c(range.ybox[1]-inter.y,max(labloc.y))
  
  
  n.xbreak = length(unique(df0$x))
  clr.g = colorRampPalette(c(clr0,"gray85"))(2*n.xbreak)
  clr.g = (clr.g[1:n.xbreak])
  
  if(iii == 5){
    clr.g = colorRampPalette(c("#0D0970","gray85"))(2*n.xbreak)
    clr.g = (clr.g[1:n.xbreak])
  }
  
  p0 =ggplot(data = df0, mapping = aes(x= x, y =y)) +
    geom_boxplot(mapping = aes(fill = x, alpha = 0.95),
                 outlier.size = 0.5,
                 outlier.color = "gray20",
                 lwd=0.3)+
    scale_x_discrete(limits = rev(levels(df0$x)))+
    scale_fill_manual(values = clr.g)+
    stat_summary(fun=median, geom="line", aes(group=1), linewidth = 0.5,
                 color = 'gray20')  + 
    stat_summary(fun=median, geom="point", size = 1.5)+
    stat_compare_means(comparisons=com.list,
                       method="wilcox.test",
                       label="p.signif",
                       label.y = labloc.y)+
    geom_text(data = df.bhco, mapping = aes(x = labloc.x,
                                            y = labloc.y,
                                            label = bhco.strs))+
    font("title", size = 12)+
    labs(tag = paste0(letters[iii]))+
    xlab("recovery rate")+
    ylab(disvar)+
    expand_limits(y = c(labloc.y+inter.y))+
    theme_classic()+
    theme(plot.tag.position = c(0.005,1),
          axis.text = element_text(color = "gray20"),
          axis.line = element_line(color = "gray20",linewidth = 0.3),
          plot.tag = element_text(face = "bold"),
          legend.position = 'none')
  
  pic.list = c(pic.list,list(p0))
  
}



df.result.ass = df.esmsimu.ass[df.esmsimu.ass$dis.fre == 0,]
varnameshorts = c("AR1","SD","SKE","KURT","AR1_bts")
disnames = c("AR1","Standard Deviation","Skewness","Kurtosis","AR1_bts")

clrs = brewer.pal(length(varnameshorts),name = "Dark2")

df.result.ass$graz.int = as.factor(df.result.ass$graz.int)

pic.list.ass = list()

for (iii in 1:length(varnameshorts)) {
  varshort = varnameshorts[iii]
  disvar = disnames[iii]
  clr0 = clrs[iii]
  
  x = df.result.ass[["graz.int"]]
  y = df.result.ass[[varshort]]
  df0 = data.frame(x,y)
  
  x.breaks = sort(unique(x))
  n.notes = length(x.breaks)-1
  
  
  com.list = list()
  bhco.strs = c()
  
  for (nn in 1:n.notes) {
    com.list = c(com.list,list(paste0(c(x.breaks[nn],x.breaks[nn+1]))))
    
    yy1 = df0$y[df0$x == x.breaks[nn]]
    yy2 = df0$y[df0$x == x.breaks[nn+1]]
    bhco.str = sprintf("%.2f", Bh.co(yy1,yy2)) 
    bhco.strs = c(bhco.strs,bhco.str)
    
  }
  
  range.ybox = range(y)
  
  inter.y = diff(range.ybox)*0.1
  
  labloc.x = seq(from = 1.5, by = 1, length.out = n.notes)
  labloc.y = (seq(from = range.ybox[2], by = inter.y, length.out = n.notes))
  
  df.bhco = data.frame(labloc.x, labloc.y, bhco.strs)
  
  plotrange.y = c(range.ybox[1]-inter.y,max(labloc.y))
  
  
  n.xbreak = length(unique(df0$x))
  clr.g = colorRampPalette(c(clr0,"gray85"))(2*n.xbreak)
  clr.g = rev(clr.g[1:n.xbreak])
  
  if(iii == 5){
    clr.g = colorRampPalette(c("#A85D12","gray85"))(2*n.xbreak)
    clr.g = rev(clr.g[1:n.xbreak])
  }
  
  p0 =ggplot(data = df0, mapping = aes(x= x, y =y)) +
    geom_boxplot(mapping = aes(fill = x, alpha = 0.95),
                 outlier.size = 0.5,
                 outlier.color = "gray20",
                 lwd=0.3)+
    scale_fill_manual(values = clr.g)+
    stat_summary(fun=median, geom="line", aes(group=1), linewidth = 0.5,
                 color = 'gray20')  + 
    stat_summary(fun=median, geom="point", size = 2)+
    stat_compare_means(comparisons=com.list,
                       method="wilcox.test",
                       label="p.signif",
                       label.y = labloc.y)+
    geom_text(data = df.bhco, mapping = aes(x = labloc.x,
                                            y = labloc.y,
                                            label = bhco.strs))+
    font("title", size = 12)+
    labs(tag = paste0(letters[iii+4]))+
    xlab("grazing intensity")+
    ylab(disvar)+
    expand_limits(y = c(labloc.y+inter.y))+
    theme_classic()+
    theme(plot.tag.position = c(0.005,1),
          axis.text = element_text(color = "gray20"),
          axis.line = element_line(color = "gray20",linewidth = 0.3),
          plot.tag = element_text(face = "bold"),
          legend.position = 'none')
  
  
  
  pic.list.ass = c(pic.list.ass,list(p0))
  
}






pic.group = pic.list[[1]] + pic.list[[2]] + pic.list[[3]] + 
  pic.list[[4]] + plot_layout(nrow = 1, guides = 'collect')


pic.group.ass = pic.list.ass[[1]] + pic.list.ass[[2]] + pic.list.ass[[3]] + 
  pic.list.ass[[4]] + plot_layout(nrow = 1, guides = 'collect')

pic.group.both =  pic.group/pic.group.ass




outfile = paste0("Fig2_ESMnoPUDs.jpg")
jpeg(outfile, units="cm",
     width= 25,
     height= 15,
     res = 600)

plot.new()

plot(pic.group.both)

dev.off()


pic.ar1bts.performance.mono = pic.list[[5]] + labs(tag = paste0(letters[1]))
pic.ar1bts.performance.bi = pic.list.ass[[5]] + labs(tag = paste0(letters[3]))


###ana3_combined_distributions##############################################

varnameshorts <- c("AR1","SD","SKE","KURT","AR1_bts")
disnames = c("AR1","Standard Deviation","Skewness","Kurtosis","AR1_bts")

pic.list.main <- list()
pic.list.main.ass <- list()


colormap <- colorRampPalette(c("#006400","#EEC900","#8B008B"))(length(unique(df.esmsimu$dis.fre)))

clsfig2 = c("#1B7A9D", "#49A35E", "#DD7520")


for (iii in 1:length(varnameshorts)) {
  varname = varnameshorts[iii]
  disvar = disnames[iii]
  
  df.dens = df.esmsimu[df.esmsimu$rec.rate == rrs[1],]
  range(df.dens$AR1)
  
  

  int.labs = paste0("Dis. Int. \nLevel ", sort(unique(df.dens$dis.inten))," ")
  names(int.labs) = sort(unique(df.dens$dis.inten))
  
  rr.labs = paste0("Rec. Rate ", sort(unique(df.dens$rec.rate))," ")
  names(rr.labs) = sort(unique(df.dens$rec.rate))
  
  
  df.sta = c()

  for (fre0 in unique(df.dens$dis.fre)){
    for (int0 in unique((df.dens$dis.inten))) {
      if(fre0==0&int0!=0){
        df00 = df.dens[df.dens$dis.fre == 0&df.dens$dis.inten ==0,]
      }else{
        df00 = df.dens[df.dens$dis.fre == fre0&df.dens$dis.inten == int0,]
      }


      std0 = sd(df00[[varname]])
      mean0 = mean(df00[[varname]])
      df.sta = rbind(df.sta,data.frame(fre0,int0,std0,mean0))
    }
  }
  
  df.sta = df.sta[df.sta$int0 != 0,]
  
  df.sta$int0 = as.factor(df.sta$int0)
  
  nue.x = 0
  nue.ymean = df.sta$mean0[df.sta$fre0 == nue.x][1]
  nue.ystd = df.sta$std0[df.sta$fre0 == nue.x][1]
  
  pic.m = ggplot(df.sta, aes(x = fre0 , y = mean0, group = int0, 
                     color = int0,fill = int0)) + 
    geom_errorbar(aes(x = fre0 ,ymin=mean0-std0, ymax=mean0+std0), 
                  width=0.001,
                  linewidth = 0.6,
                  alpha = 0.7) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5)+
    scale_fill_manual(values = clsfig2, name = "PUD Intensity")+
    scale_color_manual(values = clsfig2, name = "PUD Intensity")+
    xlab("PUD frequency")+
    ylab(disvar)+
    labs(tag = paste0(letters[iii]))+
    theme_classic()+
    theme(axis.text = element_text(color = "gray20"),
          axis.line = element_line(color = "gray20",linewidth = 0.3),
          plot.tag = element_text(face = "bold"))
  
  #####################################20241108
  
  if(iii <= length(varnameshorts)){
    pic.m = pic.m +
      theme(legend.position = 'none')
  }
  
  
  pic.list.main <- c(pic.list.main, list(pic.m))

  
}


#####


for (iii in 1:length(varnameshorts)) {
  varname = varnameshorts[iii]
  disvar = disnames[iii]
  
  df.dens = df.esmsimu.ass[df.esmsimu.ass$graz.int == rev(gis)[1],]

  df.sta = c()
  
  for (fre0 in unique(df.dens$dis.fre)){
    for (int0 in unique((df.dens$dis.inten))) {
      if(fre0==0&int0!=0){
        df00 = df.dens[df.dens$dis.fre == 0&df.dens$dis.inten ==0,]
      }else{
        df00 = df.dens[df.dens$dis.fre == fre0&df.dens$dis.inten == int0,]
      }
      
      
      std0 = sd(df00[[varname]])
      mean0 = mean(df00[[varname]])
      df.sta = rbind(df.sta,data.frame(fre0,int0,std0,mean0))
    }
  }
  
  df.sta = df.sta[df.sta$int0 != 0,]
  
  df.sta$int0 = as.factor(df.sta$int0)
  
  nue.x = 0
  nue.ymean = df.sta$mean0[df.sta$fre0 == nue.x][1]
  nue.ystd = df.sta$std0[df.sta$fre0 == nue.x][1]
  
  pic.m = ggplot(df.sta, aes(x = fre0 , y = mean0, group = int0, 
                             color = int0,fill = int0)) + 
    geom_errorbar(aes(x = fre0 ,ymin=mean0-std0, ymax=mean0+std0), 
                  width=0.001,
                  linewidth = 0.6,
                  alpha = 0.7) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5)+
    scale_fill_manual(values = clsfig2, name = "PUD Intensity")+
    scale_color_manual(values = clsfig2, name = "PUD Intensity")+
    xlab("PUD frequency")+
    ylab(disvar)+
    labs(tag = paste0(letters[iii+4]))+
    theme_classic()+
    theme(axis.text = element_text(color = "gray20"),
          axis.line = element_line(color = "gray20",linewidth = 0.3),
          plot.tag = element_text(face = "bold"))
  

  
  if(iii < length(varnameshorts)){
    pic.m = pic.m +
      theme(legend.position = 'none')
  }
  
  
  pic.list.main.ass <- c(pic.list.main.ass, list(pic.m))
  
  
}





pg.main1  =  pic.list.main[[1]]+pic.list.main[[2]]+
  pic.list.main[[3]]+pic.list.main[[4]]+
  #plot_annotation(title = "Monostability Model (return rate = 0.1)")+
  plot_layout(nrow = 1)+
  plot_layout( guides = 'collect')&theme(legend.position = 'none')

pg.main2 = pic.list.main.ass[[1]]+pic.list.main.ass[[2]]+
  pic.list.main.ass[[3]]+pic.list.main.ass[[4]]+
  plot_annotation(title = "Bistability Model (grazing intensity = 2)")+plot_layout(nrow = 1)+
  plot_layout( guides = 'collect')&theme(legend.position = 'bottom')

pg.main =  (pg.main1/pg.main2)




outfile = paste0("Fig3_Distri.jpg")
jpeg(outfile, units="cm",
     width= 25,
     height= 15,
     res = 600)

plot.new()

plot(pg.main)

dev.off()




###ana4_AR1bts############

#Run ana 2 first to get piclist 


colormap <- colorRampPalette(c("#006400","#EEC900","#8B008B"))(length(unique(df.esmsimu$dis.fre)))
varname = "AR1_bts"


df.dens = df.esmsimu[df.esmsimu$rec.rate >= rrs[1],]
range(df.dens$AR1)

int.labs = paste0("PUD Int. ", sort(unique(df.dens$dis.inten))," ")
names(int.labs) = sort(unique(df.dens$dis.inten))

rr.labs = paste0("Rec. Rate ", sort(unique(df.dens$rec.rate))," ")
names(rr.labs) = sort(unique(df.dens$rec.rate))

df.dens$dis.inten = as.factor(df.dens$dis.inten)
df.dens$dis.fre = as.factor(df.dens$dis.fre)

pic.bts.pud.mono <- ggplot(data = df.dens)+
  geom_density(mapping = aes_string(x = varname, group = "dis.fre", 
                                    color = "dis.fre"), 
               alpha = 0.6,linewidth = 0.8) +
  facet_grid(dis.inten ~ rec.rate,scales = "free",
             labeller = labeller(dis.inten = int.labs, rec.rate = rr.labs))+
  geom_hline(yintercept = 0)+
  labs(tag = paste0(letters[2]))+
  scale_color_manual(values = colormap)+
  scale_fill_manual(values = colormap)+
  scale_y_continuous(expand=c(0,0))+
  guides(color=guide_legend(title="PUD Fre."))+
  theme(strip.placement = "inside",
        strip.background = element_blank(),
        strip.text.y = element_text(color="black",angle = 0,size = 11),
        strip.text.x = element_text(color="black",angle = 0,size = 11),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        plot.tag = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_line(linewidth = 0.6),
        axis.title = element_text(color="black",size = 11),
        axis.line.x = element_blank(),
        plot.tag.position = c(0.005,1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())


#######



colormap <- colorRampPalette(c("#006400","#EEC900","#8B008B"))(length(unique(df.esmsimu$dis.fre)))
varname = "AR1_bts"



df.dens = df.esmsimu.ass[df.esmsimu.ass$graz.int >= gis[1],]
range(df.dens$AR1)

int.labs = paste0("PUD Int. ", sort(unique(df.dens$dis.inten))," ")
names(int.labs) = sort(unique(df.dens$dis.inten))

rr.labs = paste0("Gra. Int. ", sort(unique(df.dens$graz.int))," ")
names(rr.labs) = sort(unique(df.dens$graz.int))

df.dens$dis.inten = as.factor(df.dens$dis.inten)
df.dens$dis.fre = as.factor(df.dens$dis.fre)

pic.bts.pud.bi <- ggplot(data = df.dens)+
  geom_density(mapping = aes_string(x = varname, group = "dis.fre", 
                                    color = "dis.fre"), 
               alpha = 0.6,linewidth = 0.8) +
  facet_grid(dis.inten ~ graz.int,scales = "free",
             labeller = labeller(dis.inten = int.labs, graz.int = rr.labs))+
  geom_hline(yintercept = 0)+
  labs(tag = paste0(letters[4]))+
  scale_color_manual(values = colormap)+
  scale_fill_manual(values = colormap)+
  scale_y_continuous(expand=c(0,0))+
  guides(color=guide_legend(title="PUD Fre."))+
  theme(strip.placement = "inside",
        strip.background = element_blank(),
        strip.text.y = element_text(color="black",angle = 0,size = 11),
        strip.text.x = element_text(color="black",angle = 0,size = 11),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        plot.tag = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_line(linewidth = 0.6),
        axis.title = element_text(color="black",size = 11),
        axis.line.x = element_blank(),
        plot.tag.position = c(0.005,1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())





pic.ar1bts.performance.mono = pic.ar1bts.performance.mono +
  theme(axis.line = element_line(linewidth = 0.6))
pic.ar1bts.performance.bi = pic.ar1bts.performance.bi +
  theme(axis.line = element_line(linewidth = 0.6))
  
pic.group.btsperform = pic.ar1bts.performance.mono + 
  pic.bts.pud.mono +
  pic.ar1bts.performance.bi + 
  pic.bts.pud.bi + 
  plot_layout(nrow = 2, widths = c(1,3))





outfile = paste0("Fig4_AR1btsperform.jpg")
jpeg(outfile, units="cm",
     width= 25,
     height= 16,
     res = 600)

plot.new()

plot(pic.group.btsperform)

dev.off()


###ana4_supplyment#######################

varnameshorts.tradi <- c("AR1","SD","SKE","KURT")
disnames.tradi = c("AR1","Standard Deviation","Skewness","Kurtosis")
colormap <- colorRampPalette(c("#006400","#EEC900","#8B008B"))(length(unique(df.esmsimu$dis.fre)))

pics.mono = list()
pics.bi = list()

for (i in 1:length(varnameshorts.tradi)) {
  varname = varnameshorts.tradi[i]
  longname = disnames.tradi[i]
  
  
  df.dens = df.esmsimu[df.esmsimu$rec.rate >= rrs[1],]
  range(df.dens$AR1)
  
  int.labs = paste0("PUD Int. ", sort(unique(df.dens$dis.inten))," ")
  names(int.labs) = sort(unique(df.dens$dis.inten))
  
  rr.labs = paste0("Rec. Rate ", sort(unique(df.dens$rec.rate))," ")
  names(rr.labs) = sort(unique(df.dens$rec.rate))
  
  df.dens$dis.inten = as.factor(df.dens$dis.inten)
  df.dens$dis.fre = as.factor(df.dens$dis.fre)
  
  pic.pud.mono <- ggplot(data = df.dens)+
    geom_density(mapping = aes_string(x = varname, group = "dis.fre", 
                                      color = "dis.fre"), 
                 alpha = 0.6,linewidth = 0.8) +
    facet_grid(dis.inten ~ rec.rate,scales = "free",
               labeller = labeller(dis.inten = int.labs, rec.rate = rr.labs))+
    geom_hline(yintercept = 0)+
    labs(tag = paste0(letters[i]))+
    xlab(longname)+
    scale_color_manual(values = colormap)+
    scale_fill_manual(values = colormap)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(n.breaks = 4)+
    guides(color=guide_legend(title="PUD Fre."))+
    theme(strip.placement = "inside",
          strip.background = element_blank(),
          strip.text.y = element_text(color="black",angle = 0,size = 11),
          strip.text.x = element_text(color="black",angle = 0,size = 11),
          axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
          plot.tag = element_text(face = "bold"),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_line(),
          axis.title = element_text(color="black",size = 11),
          axis.line.x = element_blank(),
          plot.tag.position = c(0.005,1),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank())
  
  pics.mono = c(pics.mono,list(pic.pud.mono))
  ###
  
  
  df.dens = df.esmsimu.ass[df.esmsimu.ass$graz.int >= gis[1],]
  range(df.dens$AR1)
  
  int.labs = paste0("PUD Int. ", sort(unique(df.dens$dis.inten))," ")
  names(int.labs) = sort(unique(df.dens$dis.inten))
  
  rr.labs = paste0("Gra. Int. ", sort(unique(df.dens$graz.int))," ")
  names(rr.labs) = sort(unique(df.dens$graz.int))
  
  df.dens$dis.inten = as.factor(df.dens$dis.inten)
  df.dens$dis.fre = as.factor(df.dens$dis.fre)
  
  pic.pud.bi <- ggplot(data = df.dens)+
    geom_density(mapping = aes_string(x = varname, group = "dis.fre", 
                                      color = "dis.fre"), 
                 alpha = 0.6,linewidth = 0.8) +
    facet_grid(dis.inten ~ graz.int,scales = "free",
               labeller = labeller(dis.inten = int.labs, graz.int = rr.labs))+
    geom_hline(yintercept = 0)+
    labs(tag = paste0(letters[i]))+
    xlab(longname)+
    scale_color_manual(values = colormap)+
    scale_fill_manual(values = colormap)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(n.breaks = 4)+
    guides(color=guide_legend(title="PUD Fre."))+
    theme(strip.placement = "inside",
          strip.background = element_blank(),
          strip.text.y = element_text(color="black",angle = 0,size = 11),
          strip.text.x = element_text(color="black",angle = 0,size = 11),
          axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
          plot.tag = element_text(face = "bold"),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_line(),
          axis.title = element_text(color="black",size = 11),
          axis.line.x = element_blank(),
          plot.tag.position = c(0.005,1),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank())
  
  pics.bi = c(pics.bi,list(pic.pud.bi))
  
}



pic.group.mono =   (pics.mono[[1]]+pics.mono[[2]]+pics.mono[[3]]+pics.mono[[4]])+
  plot_layout(ncol = 1, guides = 'collect')


pic.group.bi =   pics.bi[[1]]+pics.bi[[2]]+pics.bi[[3]]+pics.bi[[4]]+
  plot_layout(ncol = 1, guides = 'collect')




outfile = paste0("Figss_ESMperform_mono.jpg")
jpeg(outfile, units="cm",
     width= 22,
     height= 25,
     res = 600)

plot.new()

plot(pic.group.mono)

dev.off()




outfile = paste0("Figss_ESMperform_bi.jpg")
jpeg(outfile, units="cm",
     width= 22,
     height= 25,
     res = 600)

plot.new()

plot(pic.group.bi)

dev.off()


###ana5_contribution########################

setwd("E:/ESM_Research/result")

rrbox = 0.6
disfrebox = c(0,0.015)
disintbox = c(0,3)
n.data = 50


picnamelist = c()
npic = 0



cl <- makeCluster(15, type = "SOCK")
registerDoSNOW(cl)

pb = txtProgressBar(max = nseries,style = 3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress = progress)

T1 = Sys.time()

df.ana <- foreach(x = 1:n.data,
                   .combine='rbind',
                   .inorder = FALSE,
                   .packages= c('earlywarnings','changepoint'),
                   .options.snow = opts,
                   .errorhandling = "pass") %dopar%
  {
    df.result0 =c()
    for (ii in 1:100) {
      rr0 = rrbox
      set.seed(x*1000+ii)
      fre0 = runif(1,min = disfrebox[1], max = disfrebox[2])
      set.seed(x*1000+ii+1)
      intf0 = runif(1,min = disintbox[1], max = disintbox[2])
      
      pars = data.frame(rec.rate = rr0,dis.fre = fre0,dis.inten = intf0)
      
      set.seed(x*1000+ii+2)
      ts.x = generatets( eq = 0,
                         ntotal = 50000,
                         rr = rr0,
                         wnsd = 0.2,
                         disfre = fre0,
                         disintense = c(-0.2)*intf0,
                         occurange = c(0.4,0.6))
      esms = gentews_noplot(ts.x)
      
      df.result0 = rbind(df.result0, cbind(esms,pars))
    }
    
    return(df.result0)
  }
stopCluster(cl)






###

clrs2 = list(c(0,0.5,1),
             c("#187596", "#EEEE8A", "#D98220")) 

znames = c("AR1_bts","AR1")
rangez = c()
zmlist = list()

x = df.ana$dis.fre
y = df.ana$dis.inten

for (i in 1:length(znames)) {
  
  zname = znames[i]
  

  z = df.ana[[zname]]
  
  df3d = data.frame(x,y,z)
  
  mod <- gam(z ~ te(x, y), data = df3d)
  xp = seq(0,0.015,length.out = 15)
  yp = seq(0,3,length.out = 30)
  
  df.plot = c()
  zm = c()
  
  for (x0 in xp) {
    
    zline = c()
    for (y0 in yp) {
      
      z0 = as.numeric(predict(mod,data.frame(x = x0,y = y0)))
      zline = c(zline,z0)
      
      df.plot = rbind(df.plot, data.frame(x0,y0,z0))
    }
    
    zm = cbind(zm,zline)
  }
  
  zmlist = c(zmlist,list(zm))
  rangez =  range(c(zm,rangez))
}

rangez = c(min(rangez)-0.05*diff(rangez),
           max(rangez)+0.05*diff(rangez))




for (i in 1:length(znames)){
  zm = zmlist[[i]]
  zname = znames[i]
  
  # col1 = colorvec100[round((range(zm)[1]-rangez[1])/diff(rangez)*100)]
  # col2 = colorvec100[round((range(zm)[2]-rangez[1])/diff(rangez)*100)]
  
  axx <- list(automargin=T,
              zerolinecolor = 'black',
              zerolinewidth = 2,
              gridcolor = 'gray20',
              title=list(text='PUD frequency', standoff = 40,
                       font=list(color='black',size = 20)))
  axy <- list(automargin=T,
              zerolinecolor = 'black',
              zerolinewidth = 2,
              gridcolor = 'gray20',
              title=list(text='PUD intensity', standoff = 40,
                         font=list(color='black',size = 20)))
  axz <- list(automargin=T,range = rangez, standoff = 25,
              zerolinecolor = 'black',
              zerolinewidth = 2,
              gridcolor = 'gray20',
              title=list(text= paste0("  ",zname,"  "), 
                         font=list(color='black',size = 20)))
  
  pic3d = plot_ly(x = xp, y = yp, z = zmlist[[i]]) %>% 
    add_surface(colorscale = clrs2,
                opacity  = 0.85,
                colorbar=list(xpad = 1, ypad =0.9,
                              thickness = 15,
                              ticklen=6,
                              tickwidth = 2,
                              orientation = 'h',
                  tickfont = list(size = 20))) %>% 
    layout(font = list(size = 15),
           margin = list(l = 0, r = 0, t = 0, b = 0),
           scene = list(
             zaxis = axz,
             xaxis = axx,
             yaxis = axy,
             camera = list(eye = list(x = -1.6, y = -1.6, z = 0.7))))
  
  htmlname = "temprorary.html"
  
  pngname = paste0("GAM3d_rr01_",zname,".png")
  picnamelist = c(picnamelist,pngname)
  npic = npic+1
  
  htmlwidgets::saveWidget(pic3d, file = htmlname)
  webshot2::webshot(htmlname, 
                    pngname, 
                    zoom =1, vwidth = 650, vheight = 650)
  
}


###########ass

gibox = 1
disfrebox = c(0,0.015)
disintbox = c(0,3)




cl <- makeCluster(15, type = "SOCK")
registerDoSNOW(cl)

pb = txtProgressBar(max = nseries,style = 3)
progress = function(n) setTxtProgressBar(pb,n)
opts = list(progress = progress)

T1 = Sys.time()

df.ana.ass <- foreach(x = 1:n.data,
                  .combine='rbind',
                  .inorder = FALSE,
                  .packages= c('earlywarnings','changepoint'),
                  .options.snow = opts,
                  .errorhandling = "pass") %dopar%
  {
    df.result0 =c()
    
    for (ii in 1:100) {
      gi0 = gibox
      set.seed(x*1000+ii)
      fre0 = runif(1,min = disfrebox[1], max = disfrebox[2])
      set.seed(x*1000+ii+1)
      intf0 = runif(1,min = disintbox[1], max = disintbox[2])
      
      pars = data.frame(graz.int = gi0,dis.fre = fre0,dis.inten = intf0)
      set.seed(x*1000+ii+2)
      ts.x = generatets.ass(ntotal = 50000,
                            k=10,
                            r=1,
                            vgh=1,
                            gi=c(gi0,gi0),
                            nsd =0.2,
                            disfre = fre0,
                            disintense = c(-0.2)*intf0,
                            occurange = c(0.4,0.6))
      esms = gentews_noplot(ts.x)
      
      df.result0 = rbind(df.result0, cbind(esms,pars))
    }
    
  
    return(df.result0)
  }
stopCluster(cl)


###


znames = c("AR1_bts","AR1")
rangez = c()
zmlist = list()

x = df.ana.ass$dis.fre
y = df.ana.ass$dis.inten

for (i in 1:length(znames)) {
  
  zname = znames[i]
  
  
  z = df.ana.ass[[zname]]
  
  df3d = data.frame(x,y,z)
  
  mod <- gam(z ~ te(x, y), data = df3d)
  xp = seq(0,0.015,length.out = 15)
  yp = seq(0,3,length.out = 30)
  
  df.plot = c()
  zm = c()
  
  for (x0 in xp) {
    
    zline = c()
    for (y0 in yp) {
      
      z0 = as.numeric(predict(mod,data.frame(x = x0,y = y0)))
      zline = c(zline,z0)
      
      df.plot = rbind(df.plot, data.frame(x0,y0,z0))
    }
    
    zm = cbind(zm,zline)
  }
  
  zmlist = c(zmlist,list(zm))
  rangez =  range(c(zm,rangez))
}

rangez = c(min(rangez)-0.05*diff(rangez),
           max(rangez)+0.05*diff(rangez))


for (i in 1:length(znames)){
  zm = zmlist[[i]]
  zname = znames[i]
  
  # col1 = colorvec100[round((range(zm)[1]-rangez[1])/diff(rangez)*100)]
  # col2 = colorvec100[round((range(zm)[2]-rangez[1])/diff(rangez)*100)]
  
  axx <- list(automargin=T,
              zerolinecolor = 'black',
              zerolinewidth = 2,
              gridcolor = 'gray20',
              title=list(text='PUD frequency', standoff = 40,
                         font=list(color='black',size = 20)))
  axy <- list(automargin=T,
              zerolinecolor = 'black',
              zerolinewidth = 2,
              gridcolor = 'gray20',
              title=list(text='PUD intensity', standoff = 40,
                         font=list(color='black',size = 20)))
  axz <- list(automargin=T,range = rangez, standoff = 25,
              zerolinecolor = 'black',
              zerolinewidth = 2,
              gridcolor = 'gray20',
              title=list(text= paste0("  ",zname,"  "), 
                         font=list(color='black',size = 20)))
  
  pic3d = plot_ly(x = xp, y = yp, z = zmlist[[i]]) %>% 
    add_surface(colorscale = clrs2,
                opacity  = 0.85,
                colorbar=list(xpad = 1, ypad =0.9,
                              thickness = 15,
                              ticklen=6,
                              tickwidth = 2,
                              orientation = 'h',
                              tickfont = list(size = 20))) %>% 
    layout(font = list(size = 15),
           margin = list(l = 0, r = 0, t = 0, b = 0),
           scene = list(
             zaxis = axz,
             xaxis = axx,
             yaxis = axy,
             camera = list(eye = list(x = -1.6, y = -1.6, z = 0.7))))
  
  htmlname = "temprorary.html"
  
  pngname = paste0("GAM3d_gi20_",zname,".png")
  picnamelist = c(picnamelist,pngname)
  npic = npic+1
  
  
  htmlwidgets::saveWidget(pic3d, file = htmlname)
  webshot2::webshot(htmlname, 
                    pngname, 
                    zoom =1, vwidth = 650, vheight = 650)
  
}



grouppic = ggarrange(rasterGrob(readPNG(picnamelist[1])),
                     rasterGrob(readPNG(picnamelist[2])),
                     rasterGrob(readPNG(picnamelist[3])),
                     rasterGrob(readPNG(picnamelist[4])),
                     ncol=2,nrow = 2,
                     labels = letters[1:npic],
                     font.label = list(size = 10, face = "bold"))

ggsave(filename = "Fig5grouped.jpg", grouppic,width = 12,height = 12,units = 'cm',
       dpi = 500)









