#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

d = read.delim(args[1], header=T, sep="\t", as.is=T, row.names=1)

clusters = read.delim("Clusters", header=T, sep="\t", as.is=T)[,-1]
clusters = data.frame(Bait=colnames(clusters), Cluster=as.numeric(clusters[1,]))
nested.clusters = read.delim("NestedClusters", header=F, sep="\t", as.is=T)[1:dim(d)[1],]
nested.phi = read.delim("NestedMu", header=F, sep="\t", as.is=T)[1:dim(d)[1],]
nested.phi2 = read.delim("NestedSigma2", header=F, sep="\t", as.is=T)[1:dim(d)[1],]
mcmc = read.delim("MCMCparameters", header=F, sep="\t", as.is=T)

### distance between bait using phi (also reorder cluster names)
### report nested clusters with positive counts only
### rearrange rows and columns of the raw data matrix according to the back-tracking algorithm

recursivePaste = function(x) {
  n = length(x)
  x = x[order(x)]
  y = x[1]
  if(n > 1) {
    for(i in 2:n) y = paste(y, x[i], sep="-")
  }
  y
}

calcDist = function(x, y) {
  if(length(x) != length(y)) stop("different length\n")
  else res = sum(abs(x-y))
  res
}


#clusters, nested.clusters, nested.phi, d

bcl = clusters
pcl = nested.clusters
phi = nested.phi
phi2 = nested.phi2
dat = d


## bipartite graph
make.graphlet = function(b,p,s) {
  g = NULL
  g$b = b
  g$p = p
  g$s = as.numeric(s)
  g
}

make.hub = function(b,p) {
  g = NULL
  g$b = b
  g$p = p
  g
}

jaccard = function(x,y) {
  j = length(intersect(x,y)) / length(union(x,y))
  j
}

merge.graphlets = function(x, y) {
  g = NULL
  g$b = union(x$b, y$b)
  g$p = union(x$p, y$p)
  g$s1 = rep(0,length(g$p))
  g$s2 = rep(0,length(g$p))
  g$s1[match(x$p, g$p)] = x$s
  g$s2[match(y$p, g$p)] = y$s
  g$s = apply(cbind(g$s1, g$s2), 1, max)
  g
}

summarizeDP = function(bcl, pcl, phi, phi2, dat, hub.size=0.5, ...) {
  pcl = as.matrix(pcl)
  phi = as.matrix(phi)
  phi2 = as.matrix(phi2)
  dat = as.matrix(dat)
  rownames(phi) = rownames(dat)
  rownames(phi2) = rownames(dat)

  ubcl = unique(as.numeric(bcl$Cluster))
  n = length(ubcl)
  pcl = pcl[,ubcl]
  phi = phi[,ubcl]
  phi2 = phi2[,ubcl]
  phi[phi < 0.05] = 0

  bcl$Cluster = match(as.numeric(bcl$Cluster), ubcl)
  colnames(pcl) = colnames(phi) = colnames(phi2) = paste("CL", 1:n, sep="")

  ## remove non-reproducible mean values
  nprey = dim(dat)[1]; nbait = dim(dat)[2]
  preys = rownames(dat); baits = colnames(dat)
  n = length(unique(bcl$Cluster))
  for(j in 1:n) {
    id = c(1:nbait)[bcl$Cluster == j]
    for(k in 1:nprey) {
      do.it = ifelse(mean(as.numeric(dat[k,id]) > 0) <= 0.5,TRUE,FALSE)
      if(do.it) {
        phi[k,j] = 0
      }
    }
  }

  ## create bipartite graphs (graphlets)
  gr = NULL
  for(j in 1:n) {
    id = c(1:nbait)[bcl$Cluster == j]
    id2 = c(1:nprey)[phi[,j] > 0]
    gr[[j]] = make.graphlet(baits[id], preys[id2], phi[id2,j])
  }

  ## intersecting preys between graphlets
  gr2 = NULL
  cur = 1
  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j) {
        combine = jaccard(gr[[i]]$p, gr[[j]]$p) >= 0.75
        if(!is.na(combine) && combine) {
          gr2[[cur]] = merge.graphlets(gr[[i]], gr[[j]])
          cur = cur + 1
        }
      }
    }
  }

  old.phi = phi
  phi = phi[, bcl$Cluster]
  phi2 = phi2[, bcl$Cluster]
  ## find hub preys
  proceed = apply(old.phi, 1, function(x) sum(x>0) >= 2)
  h = NULL
  cur = 1
  for(k in 1:nprey) {
    if(proceed[k]) {
      id = as.numeric(phi[k,]) > 0
      if(mean(id) >= hub.size) {
        h[[cur]] = make.hub(baits[id], preys[k])
        cur = cur + 1
      }
    }
  }
  nhub = cur - 1

  res = list(data=dat, baitCL=bcl, phi=phi, phi2=phi2, gr = gr, gr2 = gr2, hub = h)
  res
}

res = summarizeDP(clusters, nested.clusters, nested.phi, nested.phi2, d)

write.table(res$baitCL[order(res$baitCL$Cluster),], "baitClusters", sep="\t", quote=F, row.names=F)
write.table(res$data, "clusteredData", sep="\t", quote=F)

##### SOFT
library(gplots)
tmpd = res$data
tmpm = res$phi
colnames(tmpm) = paste(colnames(res$data), colnames(tmpm))

pdf("estimated.pdf", height=25, width=8)
tmp.res = heatmap.2(tmpm, Rowv=T, Colv=T, trace="n", col=rev(heat.colors(10)), breaks=seq(0,.5,by=0.05), margins=c(10,10), keysize=0.8, cexRow=0.4)
tmpd = tmpd[rev(tmp.res$rowInd),tmp.res$colInd]
write.table(tmpd, "clustered-matrix.txt", sep="\t", quote=F, col.names=NA)
heatmap.2(tmpd, Rowv=F, Colv=F, trace="n", col=rev(heat.colors(10)), breaks=seq(0,.5,by=0.05), margins=c(10,10), keysize=0.8, cexRow=0.4)
dev.off()

### Statistical Plots
tmp = mcmc[,2]
ymax = max(tmp)
ymin = min(tmp)
pdf("stats.pdf", height=12, width=12)

plot(mcmc[mcmc[,4]=="G",3], type="s", xlab="Iterations", ylab="Number of Clusters", main="")
plot(mcmc[,2], type="l", xlab="Iterations", ylab="Log-Likelihood", main="", ylim=c(ymin,ymax))

dev.off()

createBait2Bait = function() {
  pdf("bait2bait.pdf")
  tmp = res$phi
  colnames(tmp) = paste(colnames(res$phi), res$baitCL$Bait, sep="_")
  dd = cor(tmp)
  heatmap.2(as.matrix(dd), trace="n", breaks=seq(-1,1,by=0.1), col=(greenred(20)), cexRow=0.7, cexCol=0.7)
  dev.off()
}

tryCatch(createBait2Bait(), error=function(e) print(e), warning=function(w) print(w))
