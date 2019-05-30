# Load helper data and libraries 
source("helper.R")

### This file is too large to store here, need other data repository for it 
netfile="agg.rerank.Rdata"
label="agg.rerank"
load(netfile)
network = diag(length(genes.t))
bottom = row(network) > col(network)
colnames(network) = genes.t
rownames(network) = genes.t
network[bottom] = temp
network = network + t(network)
diag(network) = 1


### load in rhythmonome genes and results
load("nanostring_results_dec.Rdata")
m = match( genes, attr$name)
genes = attr[m,]

### plot fig 2 B 


# Set up genesets
m = match( genes$entrezID , colnames(network) )
f.g = !is.na(m)
f.n = m[f.g]
subgeneset = matrix(0, nrow=length(genes.t), ncol=1 )
rownames(subgeneset) = genes.t
subgeneset[f.n,1] = 1

### Sub network for visualization   
sub = network[,f.n]
save(sub, subgeneset, f.n, file=paste("rhyth.sub", label, "Rdata", sep=".") )

### Others
# nvrocs = neighbor_voting_nfold(subgeneset , network ) ;
# save( nvrocs, file=paste("rhyth.nv_rocs", label, "Rdata", sep=".") )

### 
# testn = calculate_functional_effects_networks(subgeneset, network, label,"./" ,"test", nsub)
# save( testn, file=paste("rhyth.enrich", label, "Rdata", sep=".") )

# res.sub = read.table("test.agg.rerank.residuals")
# file = "random.test.agg.rerankresiduals.random"

# load("test.agg.rerank.residuals.Rdata")
# resd.sub = res.sub[,6]
# wilcox.test(resd.sub, resd.rand, alt="g")
# plot(density(resd.rand, bw=0.2), ylim=c(0,0.5))
# lines(density(resd.sub, bw=0.2))
# wilcox.test(resd.sub, resd.rand, alt="g")



### Figure 1 
a = colSums(sub)
b = rowSums(sub)
m = match( names(b), genes$entrezID )
f.m = !is.na(m)
f.g= m[f.m]
a.sub = a[m[f.m]]
b.sub = b[f.m]
plot(a.sub, b.sub, pch=19, xlab="Node degree, global", ylab="Node degree, local")
colnames(sub)  = genes$name
heatmap.3( t(sub[f.m,]) , col=viridis(100) )
n = length(a)
N = length(b)

ids = colnames(sub)
m  = match(rownames(sub),genes$entrez )
f.m = !is.na(m)
f.s = m[f.m]
temp = t(sub[f.m,f.s])
colnames(temp) = genes$name[m2]


clusts.agg = get_cluster_ids(temp)
consDend = clusts.agg[[2]]
clusts.col = clusts.agg[[1]][,4]

pdf("fig1A.pdf"); 
heatmap.3(temp, col=viridis(100), Rowv=consDend, Colv=consDend, RowSideColors=clusts.col, ColSideColors=clusts.col); 
dev.off() 

cac = temp[(rownames(temp)   == "CACNA1C"), (rownames(temp)   != "CACNA1C")]
kcn = temp[(rownames(temp)   == "KCNH2"),(rownames(temp)   != "KCNH2")]

pdf("fig1b.pdf"); 
#heatmap.3( cbind(rank(kcn),rank(cac)), col=inferno(100)); 
heatmap.3( cbind(sort(kcn, decreasing=T),temp3), Rowv=F, Colv=F, col=viridis(100)); 
heatmap.3( cbind(sort(cac, decreasing=T),temp3), Rowv=F, Colv=F, col=viridis(100)); 
dev.off() 


pdf("supfig1.pdf") ; 
scatterhist(a.sub/N, b.sub/n,  col=makeTransparent(1,250), xlab="Node degree, global", ylab="Node degree, local", main="", xlim=c(0,1), ylim=c(0,1)); 
scatterhist(a.sub/N, b.sub/n,  xlab="Node degree, global", ylab="Node degree, local", main="", xlim=c(0,1), ylim=c(0,1)); 
dev.off()

clusts.nano = get_cluster_ids(nano.cor)
consDend = clusts.nano[[2]]
clusts.col = clusts.nano[[1]][,4]
### Figure 2 
pdf("fig2B.pdf"); 
heatmap.3(nano.cor, col=viridis(100), Rowv=consDend, Colv=consDend, RowSideColors=clusts.col, ColSideColors=clusts.col); 
dev.off() 

pdf("fig2A.pdf", width=8)
bp= beeswarm( log10(a) ~ b , ylim=c(-1,6), xlab="", ylab="mRNA count", corral="gutter", axes=F, pch=19, 
              pwcol=makeTransparent(rep(cividis(16),35), 250))  
text(1:35, -0.5, labels = genes$name, xpd=-1, srt=60)
axis(2, lab=10^(-1:5), at=0:6)
dev.off() 

pdf("fig2C.pdf", width=12, height=3.2)

i = 1 
par(mfrow=c(1,4))
for( ft in c("CACNA1C", "CACNA1H", "KCNIP2", "SCN5A")){ 
  f2 = genes$name==ft; 
  x= nanostring[f1,];
  y = nanostring[f2,]
  
  r2 = round(cor(x,y)^2,2)
  plot( x,y, xlab="", ylab="", ylim=c(0,1200), pch=19, col=makeTransparent(cols2[i],250), sub=r2, main=ft)
  z = lm(y~x)
  abline(z,lty=2, col=cols2[i])
  
  i = i + 1 
}
dev.off()




get_cluster_ids <- function( temp, filtMin=2,rowcols=rainbow(35)){
  ids = rownames(temp)
  consTree = hclust(  dist(temp), method="average")
  consDend = as.dendrogram(consTree)
  unmergedLabels3 = cutreeDynamic(dendro = consTree, distM = temp,
                                  deepSplit = 2, cutHeight = max(dist(temp))*0.995 ,
                                  #        deepSplit = 2, cutHeight = max(dist(temp))*0.4 ,
                                  minClusterSize = 2,
                                  pamRespectsDendro = FALSE );
  
  
  
  unmergedColors3 = sample(magma( max(unmergedLabels3) + 1 ))[ as.numeric(unmergedLabels3)+1]
  
  # Label and count modules
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
  
  
  # Filter modules
  f.freq = count(unmergedLabels.mod)[,2] < filtMin
  f.keep = count(unmergedLabels.mod)[f.freq,1]
  unmergedLabels.mod2 = unmergedLabels.mod[order(consTree$order)]
  unmergedColors4 = magma( max(unmergedLabels.mod ) + 1 )[ as.numeric(unmergedLabels.mod )+1]
  m = match( unmergedLabels.mod2, f.keep)
  f.fm = !is.na(m)
  clust =  cbind(rownames(temp)[consTree$order], ids[consTree$order], unmergedLabels.mod, unmergedColors3[consTree$order] )
  heatmap.2(temp, density="none", trace="none", col=viridis(100), Rowv=consDend, Colv=consDend, RowSideColors=rowcols, ColSideColors=unmergedColors3)
  return( list(clust, consDend) ) 
} 


scatterhist <- function(x, y, xlab="", ylab="", main, ...){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y, pch=19, main=main)
  abline(0,1, col=2, lwd=2)
  
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0, col=magma(10)[2])
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, col=magma(10)[6],horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0,
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0,
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}


calculate_functional_effects_networks <- function( gene.sets, network, nettype, dir, studies, nsub ){
  studies.genes = rownames(gene.sets)
  print(studies.genes)
  print(dim(gene.sets) )
  i =1 
  temp = as.matrix( cbind(studies.genes, as.numeric(gene.sets[,i]))[gene.sets[,i]>0,] )
  print(temp) 
  
  genes = lapply(1:length(studies), function(i) as.matrix( cbind(studies.genes, as.numeric(gene.sets[,i]))[gene.sets[,i]>0,] ) )
  labels =  lapply(1:length(studies), function(i) paste(studies[i] ,nettype,sep=".") )
  print(genes)
  r = 1000
  j= 2
  
  file = paste(dir,"random.", nsub, ".",nettype,".residuals.random",sep="" )
  
  an = lapply(1:length(studies), function(i) residual_connectivity_score(network , dir, labels[[i]], genes[[i]] ) )
  an = matrix(unlist(an), ncol =length(studies), nrow=3, byrow=F)
  dn = lapply(1:length(studies), function(i) sub_sample_resd(network , dir,labels[[i]],genes[[i]], nsub, r, file) )
  dn = matrix(unlist(dn), ncol =length(studies), nrow=3, byrow=F)
  
  colnames(an) = studies
  rownames(an) = c("res","rand", "p")
  colnames(dn) = studies
  rownames(dn) = c("res","rand", "p")
  
  return( list( an,dn))
}


residual_connectivity_score <- function(network, dir, label, gene.list){
  print("here")
  if( length(gene.list) < 0 ){ return (list(0,0,1)) }
  m = match ( rownames(network), attr$entrezID)
  f = !is.na(m)
  g = m[f]
  print("there")
  network = network[f,f]
  
  N = dim(network)[1]
  
  m = match( gene.list[,1], attr$entrezID[g] )
  f.r = !is.na(m)
  f.ar = m[f.r]
  print(gene.list[,1]) 
  
  if( sum(f.r) <= 2 ) { return (list(0,0,1)) }
  
  network_sub = network[f.ar,f.ar]
  
  node_degree =  rowSums(network, na.rm=T)
  print( mean(node_degree, na.rm=T) )
  print ( sd(node_degree, na.rm=T) )
  
  node_degree_sub =  rowSums(network_sub, na.rm=T)
  print( mean(node_degree_sub, na.rm=T) )
  print ( sd(node_degree_sub, na.rm=T) )
  
  N = dim(network)[1]
  n = dim(network_sub)[1]
  x = c(0,N)
  y = c( n/N * x[1], n )
  print(N)
  print(n)
  if( n <= 2 ) { return (list(0,0,1)) }
  print("what") 
  resd.sub = residuals( node_degree[f.ar], node_degree_sub, -length( node_degree_sub), N, 0 )
  
  file= paste0(dir, paste("random",label, sep="."),"residuals.random")
  print(file)
  if( file.exists( file )){
    print("oops")
    resd.rand = unlist(read.table(file))
  } else {
    print("phew")
    resd.rand  = unlist(explore_sub_network_list_random3(network, dir, paste("random",label, sep="."),n))
  }
  print("the")
  # resd.rand  = unlist(explore_sub_network_list_random3(network, dir, paste("random",label, sep="."),n))
  resd.rand = as.matrix(resd.rand )
  resd.rand = sort(resd.rand)
  
  X = cbind(node_degree[f.ar], node_degree_sub)
  H = X %*% solve( t(X) %*% X ) %*% t(X)
  X = resd.sub/ sd(resd.sub) * sqrt(1 - diag(H) )
  
  a=round(mean(resd.rand , na.rm=T),3)
  d=round(mean(X , na.rm=T),3)
  test = wilcox.test(resd.sub, resd.rand, alt="g")
  
  pvals = sapply( 1:length(X), function(i) sum ( resd.rand > X[i] )  )  / length(resd.rand)
  pvals.adj = p.adjust( pvals, method="BH")
  print(c(a, d, test$p.value))
  
  results.all = cbind( as.character(attr$name[g][f.ar]),  as.character(attr$entrezID[g][f.ar]), node_degree[f.ar], node_degree_sub, resd.sub, X, pvals, pvals.adj)
  colnames(results.all) = c("Gene", "entrezID", "Node degree full", "Node degree sub", "Residuals", "X", "P-vals", "P-vals adj")
  write.table(results.all, file = paste(dir, label, ".residuals", sep=""))
  write.table(c(N,n), file = paste(dir, label, ".Ns", sep=""))
  print("finale")
  #return(results.all)
  return (list(a, d, test$p.value))
  
}
