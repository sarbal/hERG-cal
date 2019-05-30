# helper functions 
if( !require("gplots")){
  install.packages("gplots")
}
if( !require("plyr")){
  install.packages("plyr")
}
if( !require("RColorBrewer")){
  install.packages("RColorBrewer")
}
if( !require("viridis")){
  install.packages("viridis")
}
if( !require("beanplot")){
  install.packages("beanplot")
}
if( !require("beeswarm")){
  install.packages("beeswarm")
}
if( !require("corrgram")){
  install.packages("corrgram")
}
if( !require("vioplot")){
  install.packages("vioplot")
}
if( !require("pheatmap")){
  install.packages("pheatmap")
}
if( !require("vioplot")){
  install.packages("vioplot")
}

if( !require("tidyverse")){
  install.packages("tidyverse")
}

# Colors 
cols = colorpanel(16, "red", "blue")
cols2 = brewer.pal(8, "Spectral")
cols3= rainbow(30)
cols4 = colorpanel(63, "lightgrey", "blue", "darkblue")
cols5 = colorpanel(300, "lightgrey", "red", "darkred")
cols6 = colorpanel(100, "lightgrey", "red", "darkmagenta")
cols7 = c("seagreen", "black", "darkmagenta")
cols8 = viridis(10)
cols9 = colorpanel(100, "white", "red", "darkmagenta")
cols10 = colorpanel(100, "white", "blue", "darkcyan")
cols11 = colorpanel(100, "white", "orange", "deeppink4")
cols12 = magma(100)
cols13 = viridis(100)
cols.rec  =c(magma(5)[2], "deeppink4", "darkcyan" )



# Random functions 
rank_std <- function(x)  { r = rank( x, na.last="keep"); r/max(r, na.rm=T)  }
colSD <- function( data){ return(apply( data, 2, sd, na.rm=T))}
rowSD <- function( data){ return(apply( data, 1, sd, na.rm=T))}
colSE <- function( data){ return( apply( data, 2, sd, na.rm=T)/sqrt(dim(data)[2]))}
rowSE <- function( data){ return( apply( data, 1, sd, na.rm=T)/sqrt(dim(data)[1]))}
se    <- function(x){ sd(x,na.rm=T)/sqrt(length(!is.na(x))) }
rmse  <- function(error){sqrt(mean(error^2, na.rm=T) )}
mae   <- function(error){ mean(abs(error), na.rm=T)}

geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

geo_sd <- function(data) {
  log_data <- log(data)
  gs <- exp(sd(log_data[is.finite(log_data)]))
  return(gs)
}

geo_se <- function(data) {
  gs <- geo_sd(data)
  log_data <- log(data)
  gse <- gs/sqrt(sum(is.finite(log_data)))
  return(gse)
}


lm.studentized  <- function(x,y){
  z = lm(y ~ x )
  z = rstudent(z)
  return( rank(abs(z)) )
}

lm.function  <- function(x,y){
  z = lm(y ~ x )
  return( rank(abs(z$residuals)) )
}


residuals <- function(x,y,A,B,C){ (A*x + B*y + C) }
residuals2 <- function(x,y,A,B,C) { (A*x + B*y + C)/sqrt(A^2+B^2) }
residuals3 <- function(x,y,A,B,C) { abs(A*x + B*y + C)/sqrt(A^2+B^2) }

z_scores <- function(x) {
  mean_x = mean(x, na.rm=T)
  sd_x = sd(x, na.rm=T)
  z =  (x - mean_x) / (sd_x)
  return(z)
}

z_scores_mod <- function(x) {
  med_x = median(x, na.rm=T)
  mad_x = median(abs(x-med_x), na.rm=T)
  z =  0.6745 * (x - med_x) / (mad_x)
  return(z)
}

calc_cpm <-function(X){
  K  = colSums(X)
  X.cpm = sapply(1:length(K), function(k) 10^6*X[,k]/K[k] )
  return(X.cpm)
}


heatmap.3 <- function(mat, ...){
  heatmap.2( mat, ..., density="none", trace="none")
}




# Transparent colors
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


# Given x and two points
get_value <- function( x1, x2, y1,y2, x) {
  m = (y2 - y1) / (x2 - x1 )
  y = y1 + m *( x - x1)
  return(y)
}

# Given y and two points
get_value_x <- function( x1, x2, y1,y2, y) {
  m = (y2 - y1) / (x2 - x1 )
  x = x1 + (y - y1)/m
  return(x)
}


## Formats the density distribution from the histogram function
get_density <- function(hist)
{
  x = sort(rep(hist$breaks,2))
  y = matrix(rbind(hist$density, hist$density))
  y = c(0,y,0)
  
  return(cbind(x,y))
}


## Formats the counts distribution from the histogram function
get_counts <- function(hist)
{
  x = sort(rep(hist$breaks,2))
  y = matrix(rbind(hist$counts, hist$counts))
  y = c(0,y,0)
  
  return(cbind(x,y))
}

# Tic toc functions
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}

convolve_nets <- function(netA,netB,f){
  n <- order(netA)
  temp_netA <- netA[n]
  
  temp_netB = convolve( netB[n], rep(1,f),type="filter")
  convolved = cbind(temp_netA[(f/2):(length(temp_netA)-f/2)],temp_netB/f)
}

gene_set_enrichment <- function(genes, genes.labels, voc){
  
  genes.names = rownames(genes.labels)
  labels.names = colnames(genes.labels)
  genes.counts = rowSums(genes.labels)
  labels.counts = colSums(genes.labels)              			# p
  
  m = match ( genes, genes.names )
  filt.genes  = !is.na(m)
  filt.labels = m[filt.genes]
  
  
  labels.counts.set = rep( sum(filt.genes), length(labels.counts) )	# g
  
  m = match (labels.names, voc[,1])
  v.f = !is.na(m)
  v.g = m[v.f]
  
  universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
  if(  length(filt.labels) == 1 ) { genes.counts.set = genes.labels[filt.labels,] }
  else { genes.counts.set = colSums(genes.labels[filt.labels,]) }             ## does weird things with 0 sets
  
  test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
  pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
  pvals.adj = p.adjust( pvals, method="BH")
  
  results = cbind(voc[v.g,1:2], test[v.f,c(1)]+1, test[v.f,c(2)] , pvals[v.f], pvals.adj[v.f] )
  colnames(results) = c("term", "descrp","p", "q", "pvals", "padj" )
  return (results)
  
}

filter_coexp_median <- function(n.coexp, runid, filtMin=6, medK=0.28, FLAG_PLOT=FALSE){
  temp = n.coexp
  temp[temp >medK] = 1
  temp[temp <= medK] = 0
  
  consTree = hclust(as.dist(temp), method = "average");
  consDend = as.dendrogram(consTree)
  unmergedLabels3 = cutreeDynamic(dendro = consTree, distM = temp, deepSplit = 2, cutHeight = 0.995, minClusterSize = 2, pamRespectsDendro = FALSE );
  nsclust = as.numeric(unmergedLabels3)+1
  unmergedColors3 = magma( max(nsclust))[nsclust]
  if( FLAG_PLOT == TRUE) {
    heatmap.2(temp, density="none", trace="none", col=cols, Rowv=consDend, Colv=consDend, RowSideColors=unmergedColors3, ColSideColors=unmergedColors3,cexRow=0.5, cexCol=0.5, main=runid)
    heatmap.2(n.coexp, density="none", trace="none", col=cols, Rowv=consDend, Colv=consDend, RowSideColors=unmergedColors3, ColSideColors=unmergedColors3,cexRow=0.5, cexCol=0.5, main=runid)
  }
  i.prev = ""
  ki = 1
  ji = 0
  unmergedLabels.mod = as.numeric(unmergedLabels3[consTree$order]) * 0
  for( ii in as.numeric(unmergedLabels3[consTree$order]) ){
    if( ii == i.prev){
      unmergedLabels.mod[ki] = ji
    } else {
      i.prev = ii
      ji = ji + 1
      unmergedLabels.mod[ki] = ji
      
    }
    ki = ki + 1
  }
  
  f.freq = count(unmergedLabels.mod)[,2] < filtMin
  f.keep = count(unmergedLabels.mod)[f.freq,1]
  unmergedLabels.mod2 = unmergedLabels.mod[order(consTree$order)]
  m = match( unmergedLabels.mod2, f.keep)
  f.fm = !is.na(m)
  # deg.filt.list[[i]] = f.fm
  if( FLAG_PLOT == TRUE & sum(f.fm) > 0) { heatmap.2(n.coexp[f.fm,f.fm], density="none", trace="none", cexRow=0.5, cexCol=0.5, main=runid) }
  return(f.fm)
}

calc_fdrs_recur  <- function( data, pp = 0.05) {
  
  temp1 = lapply(1:1000, function(i) shuffle_cols(data) )
  temp2 = sapply(1:1000, function(i) rowSums(temp1[[i]], na.rm=T))
  
  nmax = dim(data)[2] + 1
  
  recur = rowSums(data, na.rm=T)
  ob =  count_recur( as.matrix(recur) , nmax)
  ex = rowSums(sapply(1:1000, function(j) count_recur( as.matrix(temp2[,j] ), nmax ) ) )/1000
  
  test = cbind( ob, ex)
  
  observed = cbind( (rev(cumsum(rev(test[,1])))), 0:(nmax-1))
  expected = cbind( (rev(cumsum(rev(test[,2])))), 0:(nmax-1))
  
  FDR = expected[,1]/observed[,1]
  Pt = expected[ min(which(FDR < pp ) ), 2]
  sig = sum(test[  (which(FDR < pp )) ,1], na.rm=T)
  res = list(FDR, test, Pt, sig, pp)
  names(res) = c("FDR", "test", "Pt", "sig", "pp")
  return( res )
}


count_recur <- function(data,nmax){
  freq = plyr::count(data[,1])
  res = matrix(0, nrow=nmax, ncol=1 )
  rownames(res) = 0:(nmax-1)
  m = match(0:(nmax-1), freq[,1])
  f.r = !is.na(m)
  f.f = m[f.r]
  res[f.r,1] = freq[f.f,2]
  return(res)
}


shuffle_cols <- function( data ) {
  nc = dim(data)[2]
  nr = dim(data)[1]
  return (sapply(1:nc, function(i) data[sample(nr),i] ))
}



require(dynamicTreeCut)
### filter common coExp
run_filtering <- function( subgenesets, flag, net, outputflag=T){
  
  m = match( rownames(network), rownames(subgenesets) )
  f.n = !is.na(m)
  f.ee = m[f.n]
  
  net  = net[f.n,f.n]
  subgenesets = subgenesets[f.ee,]
  subgenesets.filt = subgenesets*0
  
  res.prec = list()
  n = dim(subgenesets)[2]
  
  for(i in 1:n ){
    f.sub =  subgenesets[,i]==1
    subnet = net[f.sub,f.sub]
    genes.sub = subgenesets[f.sub,i]
    #save( subnet, f.sub, genes.sub, file= paste(colnames(subgenesets)[i], nettype, flag,"subnet.Rdata",sep=".") )
    n = dim(subnet)[1]
    if( length(n) == 0 ) { print(n); next; }
    tempfilt =  filter_coexp_median(subnet, colnames(subgenesets)[i])
    tempfilt2 = rep(0,n)
    tempfilt2[tempfilt] = 1:sum(tempfilt)
    names(tempfilt2) = rownames(subnet)
    res.prec[[i]] = tempfilt2
    m = match( names(which(res.prec[[i]]>0)), rownames(subgenesets))
    subgenesets.filt[m,i] = 1
  }
  
  pre.post = cbind(colSums(subgenesets), colSums(subgenesets.filt))
  recur.filt = rowSums(subgenesets.filt)
  recur      = rowSums(subgenesets)
  fdrs = calc_fdrs_recur( subgenesets )
  fdrs.filt = calc_fdrs_recur( subgenesets.filt )
  if(outputflag==T){
    save(res.prec, file=paste(disease,"res.prec", flag, "Rdata", sep="."))
    save(fdrs.filt, subgenesets.filt, file=paste("fdr_calcs", flag, "Rdata", sep="."))
    save(fdrs.filt, subgenesets.filt, file=paste("fdr_calcs.filt", flag, "Rdata", sep="."))
    save(recur, recur.filt, subgenesets.filt, subgenesets, pre.post, file=paste(disease,"recur.filt.genes", flag, "Rdata", sep="."))
  }
  return( list(recur, recur.filt,fdrs, fdrs.filt, subgenesets, subgenesets.filt, pre.post))
}

calc_binom_recur <- function(subgenesets){
  freq = plyr::count( rowSums(subgenesets) ) 
  n = dim(subgenesets)[2]
  N = dim(subgenesets)[1]
  p = max(colSums(subgenesets))/N

  fdrs.bin = p.adjust(pbinom(0:n, n, p, lower.tail=F), n=N )
  cbind(fdrs.bin, 0:n)
  return(fdrs.bin)
}


gene_set_enrichment <- function(genes, genes.labels, voc){
  
  genes.names = rownames(genes.labels)
  labels.names = colnames(genes.labels)
  genes.counts = rowSums(genes.labels)
  labels.counts = colSums(genes.labels)              			# p
  
  m = match (  genes, genes.names)
  filt.genes = !is.na(m) 
  filt.labels  = m[filt.genes]
  
  
  labels.counts.set = rep( sum(filt.genes), length(labels.counts) )	# g
  
  m = match (labels.names, voc[,1])
  v.f = !is.na(m)
  v.g = m[v.f]
  
  universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
  if(  length(filt.labels) == 1 ) { genes.counts.set = genes.labels[filt.labels,] }
  else { genes.counts.set = colSums(genes.labels[filt.labels,]) }             ## does weird things with 0 sets
  
  test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
  pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
  pvals.adj = p.adjust( pvals, method="BH")
  
  results = cbind(voc[v.g,1:2], test[v.f,c(1)]+1, test[v.f,c(2)] , pvals[v.f], pvals.adj[v.f] )
  colnames(results) = c("term", "descrp","p", "q", "pvals", "padj" )
  #results =  results[results[,3] > 0 ,]
  return (results)
  
}

plot_recurrence <- function(xrange, n.sig, p.tests, ...){
  plot( xrange , n.sig , type="l", lwd=2, axes=F, ...)
  axis(1); axis(2);
  par(new=TRUE)
  plot( xrange , -log10(p.tests) , type="l",   axes=F, col="lightgrey", ylab="", xlab="")
  axis(4); mtext("-log10 adjusted P-value of significance threshold",4)
}


run_recur_compare(recurs,pt.all,  "../fig6.panelC.Rdata")

require(plyr)
require(tidyverse)

recur_mat <- function(recurs){
  temp1 = recurs
  temp = plyr::count( temp1)
  temp.mat = spread(temp, key=x.recur.filt, value=freq)
  rownames(temp.mat) = temp.mat[,1]
  temp.mat = temp.mat[,-1]
  temp.mat[is.na(temp.mat)] = 0
  temp.mat = as.matrix(temp.mat )
  
  mat10 =   log10(temp.mat )  + 1
  mat10[!is.finite(mat10)] = 0
  return(mat10)
}

run_recur_compare <- function(recurs, pt.all, filename){
  require(plyr)
  require(tidyverse)
  # Set up upregulated
  temp1 = recurs[,c(4,3)]
  temp = plyr::count( temp1)
  temp.mat = spread(temp, key=x.recur.u.filt, value=freq)
  rownames(temp.mat) = temp.mat[,1]
  temp.mat = temp.mat[,-1]
  temp.mat[is.na(temp.mat)] = 0
  temp.mat = as.matrix(temp.mat )
  
  up.mat10 =   log10(temp.mat )  + 1
  up.mat10[!is.finite(up.mat10)] = 0
  Pt.pre.up = pt.all[3]
  Pt.post.up = pt.all[4]
  
  
  # Set up downregulated
  temp1 = recurs[,c(6,5)]
  temp = plyr::count( temp1)
  temp.mat = spread(temp, key=x.recur.d.filt, value=freq)
  rownames(temp.mat) = temp.mat[,1]
  temp.mat = temp.mat[,-1]
  temp.mat[is.na(temp.mat)] = 0
  temp.mat = as.matrix(temp.mat )
  
  down.mat10 =   log10(temp.mat )  + 1
  down.mat10[!is.finite(down.mat10)] = 0
  Pt.pre.down = pt.all[5]
  Pt.post.down = pt.all[6]
  
  save(up.mat10, Pt.pre.up, Pt.post.up, down.mat10, Pt.pre.down ,Pt.post.down, plot_2D_hist, file=filename)
  
}

plot_2D_hist <- function(mat,pt.x, pt.y, ...){
  image( mat, axes=F, ... )
  n=dim(mat)[2] -1
  ni = diff((0:n/n))[1]
  abline( h= (0:n/n)[pt.y]+ni/2, lwd=3, col="grey")
  axis(2, at=0:n/n, lab=0:n)
  n=dim( mat)[1] -1
  ni = diff((0:n/n))[1]
  axis(1, at=0:n/n, lab=(0:n)  )
  abline( v= (0:n/n)[pt.x]+ni/2, lwd=3, col="grey")
}

