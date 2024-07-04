
library("changepoint")

##this function is for calculating AR1_bts of a time series

AR1.bts <- function(ts0,
                    max.bp.ratio = 0.01,
                    min.length = 0.1){
  #ts0 is the original time series for calculation
  #max.bp.ratio restrains the maximum number of breaking points by its ratio to the total length of ts0
  #min.length defines the minimum length of each broken time seires segment by its retio to the total length of ts0
  
  max.bp = round(length(ts0)*max.bp.ratio)
  cpt.p = cpt.meanvar(ts0,penalty= "MBIC",method="PELT",
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
