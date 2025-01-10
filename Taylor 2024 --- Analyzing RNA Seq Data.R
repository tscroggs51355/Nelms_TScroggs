
LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]

setwd("C:/Users/taylo/Desktop/Scroggs 2024 - Nelms lab")
load('2022.rda')

## Combine technical replicates 
B[,grepl('02', colnames(B))] = B[,grepl('02', colnames(B))] + B[,grepl('24', colnames(B))]
B[,grepl('23', colnames(B))] = B[,grepl('23', colnames(B))] + B[,grepl('01', colnames(B))]
B = B[,!grepl('24', colnames(B)) & !grepl('01', colnames(B))]


##### 
library(ComplexHeatmap)
library(seriation) 
## seriation flips the dendrogram 
## minmax and minmax2 are sensitive to outiers 
## minmax gets you to an absolute zero 
## minmax 2 does not 
## question whether absolute zero is significant to your interpretation? 
minmax = function(x) {
  sweep(x - log(pseudocount, 10), 1, apply(x - log(pseudocount, 10), 1, max), '/')
}
minmax2 = function(x) {
  x = sweep(x, 1, apply(x, 1, min), '-')
  sweep(x, 1, apply(x, 1, max), '/')
}


SScore = function(A, B, alpha = .2, len = NULL) {
  # xx = match(A, B)  ## For named vectores
  xx = rank(-A)[order(-B)]
  #xx = xx[!is.na(A) & !is.na(B)]
  #if (is.null(leng)) { len = max(c(round(9.21034/alpha + 1), length(xx))) }
  xx = sapply(1:len, function(i) { sum(xx[1:i] <= i) })
  return(sum(exp(-alpha*1:length(xx))*xx))
}

SSmat = function(x, alpha = .2) {
  ncx <- ncol(x)
  len = max(c(round(9.21034/alpha + 1), nrow(x)))
  r <- matrix(0, nrow = ncx, ncol = ncx)
  for (i in seq_len(ncx)) {
    for (j in i:ncx) {
      x2 = x[,i]
      y2 = x[,j]
      ok <- complete.cases(x2, y2)
      r[i, j] <- if (any(ok)) {
        SScore(x2[ok], y2[ok], alpha = alpha, len = len)
      } else NA
      r[j, i] = r[i, j]
    }
  }
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  return(r)
}


controls = c('A02','P02','B03','O03','H04','E06','L06','A09','P09','G10','J10','N12')
controls = c(controls, LayoutPlate[controls])

TFinfo = read.csv('"C:\Users\taylo\OneDrive - University of Georgia\2023 Nelms\gene names.csv\gene names (1).csv')


dim

#########



pseudocount = 100
### good baseline, works a lot of the time, think about when you hit counting error and that is dependent on sequencing depth 
## higher pseudocount gets rid of noise, but less sensitive to low expression of genes 
## low pseudocount for specific genes 

B2 = B[!(rownames(B) %in% unlist(strsplit(TFinfo[,2], ', '))),]
## excluding  164 target TFs 
#heatmap(matrix((colSums(B2) >= 10000)[LayoutPlate], nrow=16)[16:1,]+0, scale='none', Rowv=NA, Colv=NA)
B2 = B2[,colSums(B2) >= 5000]
## need at least 5000 umi's per sample 
TPM = sweep(B2,2,colSums(B2),'/')*10^6
B3 = log(TPM + pseudocount, 10)
## TPM normalize 
## pseudocount and log transform -- this is all standard practice in the lab 
B4 = B3[rowSums(B2 >= 10) >= 2,]
## need at least two samples with 10 UMI's per sample 

#boxplot(by(apply(B3[,colnames(B3) %in% controls], 1, sd), round(rowSums(B2[,colnames(B3) %in% controls])/100)*100, c))


Z = sweep(B4, 1, rowMeans(B4[,colnames(B4) %in% controls]), '-')
# Ztransforming - across the rows, substract rowmeans B4 in the control wells 

Z = Z/mean(apply(Z[,colnames(Z) %in% controls], 2, sd))
# dividing based upon standard deviation of controls 

summary(colSums(B))  # UMI stats; mean UMIs / sample was 32,261

ncol(Z) # 324 out of 352 samples passed QC (92.0%)
ncol(Z)/352


crossval = NULL
for (i in which(controls %in% colnames(B4))) {
  Zcv = sweep(B4, 1, rowMeans(B4[,colnames(B4) %in% controls[-i]]), '-')
  Zcv = Zcv/mean(apply(Zcv[,colnames(Zcv) %in% controls[-i]], 2, sd))
  crossval = cbind(crossval, Zcv[,controls[i]])
}
colnames(crossval) = controls[which(controls %in% colnames(B4))]

Znull = NULL
for (i in 1:(ncol(crossval) - 1)) {
    for (j in (i+1):ncol(crossval)) {
      Znull = c(Znull, rowMeans(crossval[,c(i,j)]))
    }
}
#### How you get the significant genes 
-quantile(-Znull, p = c(.01,.001,.0001,.00001))

numSig = function(well, alpha = .01/100) { sum(rowMeans(Z[,c(well, LayoutPlate[well])]) > -quantile(-Znull, p = alpha)) }
## hist(Znull)
## abline(v=3.456945)
## abline(v=-3.456945)

calcP = function(well) {
  pval = (sapply(rowMeans(Z[,c(well, LayoutPlate[well])]), function(x) { sum(x <= Znull) })+.5)/(length(Znull)+.5)
  padj = p.adjust(pval, method = 'BH') # BH = FDR, holm = p-value
  return(data.frame(p.val = pval, p.adj = padj))
}
#g2pvals = calcP('M04')
#kn1pvals = calcP('J07')
#save(g2pvals,kn1pvals, file = 'pvals.rda')



boxplot(crossval[,order(colSums(B[,colnames(crossval)]))])
abline(h = 0)


######################
Zz = Z[,LayoutPlate[colnames(Z)] %in% colnames(Z)]
Zz = Zz[,!(colnames(Zz) %in% controls)]
cors.wt = cor(sign(Zz)*Zz^2)
diag(cors.wt) = NA
TFcors = diag(cors.wt[,LayoutPlate[colnames(cors.wt)]])
names(TFcors) = rownames(cors.wt)
TFcors = TFcors[order(-TFcors)][seq(1,length(TFcors),2)]

rownames(TFinfo) = TFinfo[,1]

plot(TFcors, pch = 19)
abline(h = quantile(cors.wt, 1 - .05/length(TFcors), na.rm=T), col = 'red')
TFinfo[names(TFcors),][which(TFcors > quantile(cors.wt, 1 - .05/length(TFcors),na.rm=T)),]
######################



par(mfrow = c(2,1))
hist(crossval, breaks = 40, xlim = c(-10,10), ylim = c(0,2000))
hist(Z[,!(colnames(Z) %in% controls)], breaks = 40, xlim = c(-10,10), ylim = c(0,2000))


plotTF = function(well) { plot(Z[,c(well, LayoutPlate[well])], pch = 19, xlim = c(-8,8), ylim = c(-8,8)) }


topDEsTF = function(well, numTops = 172 ) {
Ztf = Z[,c(well, LayoutPlate[well])]
 Ztf = Ztf[order(-rowMeans(Ztf)),]
 return(Ztf[1:numTops,])
}

topDEsTF = function(well, numTops = 22, outputCSV = NULL) {
  Ztf = Z[, c(well, LayoutPlate[well])]
  Ztf = Ztf[order(-rowMeans(Ztf)),]
  topDEs = Ztf[1:numTops,]
  
  if (!is.null(outputCSV)) {
    write.csv(topDEs, file = outputCSV, row.names = TRUE)
  }
  
  return(topDEs)
}
topDEsTF_result <- topDEsTF("J07")
topDEsTF_g2 <- topDEsTF("M04")

controlsB = controls[controls %in% colnames(B2)]

plotTopDEs = function(well, numTops = 20) {
  Heatmap(minmax(B3[rownames(topDEsTF(well, numTops = numTops)), c(well, LayoutPlate[well], controlsB)]))
}


TFplot = function(well) {
  Btf = rowMeans(B2[,c(well, LayoutPlate[well])])
  Bnull = rowMeans(B2[,colnames(B2) %in% controls])
  plot(Bnull + 1, Btf + 1, pch = 19, log = 'xy', xlab = 'Mean mCherry Expression (TPM)', ylab = 'Mean TF Expression (TPM)')
}

topDEsTF('A02')[1:10,]
  
Zx = Z[,!(colnames(Z) %in% controls)]
Zx = Zx[,LayoutPlate[colnames(Zx)] %in% colnames(Zx)]
Zx = sign(Zx)*Zx^2
cors = cor(Zx[rank(-apply(Zx,1,sd)) <= 2000,])
diag(cors) = NA


Zfilt = Z
Zfilt[Zfilt <= 0] = NA
Heatmap(cor(Zfilt, use = "pairwise.complete.obs"))



corRank = apply(-cors, 1, rank)
ReplicateCorRank = cbind(diag(corRank[,LayoutPlate[colnames(corRank)]]), diag(corRank[LayoutPlate[colnames(corRank)],]))
