
files = dir('Mapped_Data/UMIcounts')
A0 = list()
for (f in files) {
	A0[[f]] = read.table(paste('Mapped_Data/UMIcounts/', f, sep = ''), sep = '\t', header=T, row.names=1)
}

gn = unique(unlist(lapply(A0, rownames)))
A = matrix(NA, nrow = length(gn), ncol = length(files))
rownames(A) = gn
colnames(A) = files
for (f in files) {
	A[,f] = A0[[f]][match(gn,rownames(A0[[f]])),1]
}
colnames(A) = sub('TG1-1_','',sub('_S.+L004_dT-', '_', sub('.tsv', '', files)))
A[is.na(A)] = 0
A = A[rowSums(A) > 0,]
#save(A, file = 'TG1-1 repeat.rda')








########################################################################

LayoutBC = rep(rep(paste('Tube',1:4,sep=''), each = 32), 3)
LayoutBC[seq(1,384,2)] = paste(LayoutBC[seq(1,384,2)], '_', c(rep(1:16, 4), rep(1:16 + 16*2, 4), rep(1:16 + 16*4, 4)), 's', sep = '')
LayoutBC[seq(2,384,2)] = paste(LayoutBC[seq(2,384,2)], '_', c(rep(1:16 + 16, 4), rep(1:16 + 16*3, 4), rep(1:16 + 16*5, 4)), 's', sep = '')
LayoutPlate = paste(rep(LETTERS[1:16], 24), 0, rep(1:24, each = 16), sep = '')
LayoutPlate = unlist(lapply(strsplit(LayoutPlate, ''), function(xx) { paste(xx[1], xx[length(xx) - 1], xx[length(xx)], sep = '') }))
names(LayoutPlate) = LayoutPlate[384:1]
load('2022.rda')

LayoutBC67 = rep(NA,384)
LayoutBC67[seq(2,384,2)[rep(rep(c(F,T), each = 8), 8)]] = paste(rep(c('Tube7','Tube6'), each = 48), '_', rep(1:48, 2), 's', sep = '')
LayoutBC67[seq(2,384,2)[c(rep(F, 112), rep(rep(c(T,F), each = 8), 5))]] = paste('Tube6_', 49:88, 's', sep = '')
LayoutBC67[c(72,157,191,313)] = paste('Tube7_', 49:52, 's', sep = '')



load('TG1-1 repeat.rda')
B = A[,grepl('Tube[1-4]', colnames(A))]
colnames(B) = sub('_.+_','_',colnames(B))
B = B[,order(colnames(B))]
B = B[,seq(1,ncol(B),3)] + B[,seq(2,ncol(B),3)] + B[,seq(3,ncol(B),3)]
B = B[,LayoutBC]
colnames(B) = LayoutPlate


rep1 = cbind(B[,grepl('02', colnames(B))], B[,grepl('23', colnames(B))])
rep2 = cbind(B[,grepl('24', colnames(B))], B[,grepl('01', colnames(B))])
keeps = (colSums(rep1) >= 10000) & (colSums(rep2) >= 10000)
diag(cor(rep1,rep2))[keeps]


Bx = A[,grepl('Tube[6-7]', colnames(A))]
colnames(Bx) = sub('_.+_','_',colnames(Bx))
Bx = Bx[,order(colnames(Bx))]
Bx = Bx[,seq(1,ncol(Bx),3)] + Bx[,seq(2,ncol(Bx),3)] + Bx[,seq(3,ncol(Bx),3)]
Bx = Bx[,na.omit(LayoutBC67)]
colnames(Bx) = LayoutPlate[!is.na(LayoutBC67)]
B[,LayoutPlate[!is.na(LayoutBC67)]] = B[,LayoutPlate[!is.na(LayoutBC67)]] + Bx[,LayoutPlate[!is.na(LayoutBC67)]]


#heatmap(matrix(log(colSums(B))[LayoutPlate],nrow=16)[16:1,], scale='none', Rowv=NA, Colv=NA)
#heatmap(matrix((colSums(B)>=10000)[LayoutPlate],nrow=16)[16:1,]+1, scale='none', Rowv=NA, Colv=NA)



### Combine technical replicates
#heatmap(cor(cbind(B[,grepl('02', colnames(B))], B[,grepl('24', colnames(B))],B[,grepl('23', colnames(B))], B[,grepl('01', colnames(B))])), scale='none')
B[,grepl('02', colnames(B))] = B[,grepl('02', colnames(B))] + B[,grepl('24', colnames(B))]
B[,grepl('23', colnames(B))] = B[,grepl('23', colnames(B))] + B[,grepl('01', colnames(B))]
B = B[,!grepl('24', colnames(B)) & !grepl('01', colnames(B))]

##### Make a new script starting with B 

#####################

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

TFinfo = read.csv('gene names (1).csv')


pseudocount = 100
### good baseline, works a lot of the time, think about when you hit counting error and that is dependent on sequencing depth 
## higher pseudocount gets rid of noise, but less sensitive to low expression of genes 
## low pseudocount for specific genes 

B2 = A[!(rownames(B) %in% unlist(strsplit(TFinfo[,2], ', '))),]
## excluding  164 target TFs 
#heatmap(matrix((colSums(B2) >= 10000)[LayoutPlate], nrow=16)[16:1,]+0, scale='none', Rowv=NA, Colv=NA)
B2 = B2[,colSums(B2) >= 5000]
## need at least 5000 umi's per sample 
B3 = log(sweep(B2,2,colSums(B2),'/')*10^6 + pseudocount, 10)
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

## dim(Z) 
##[1] 6652  324

#### EXCLUDE THIS 
if (F) {
Zneg = Z[,colnames(Z) %in% controls]
Zpos = list()
Zpos[['1']] = Z[,(as.numeric(sub('[A-Z]','',colnames(Z))) %in% 2:12) & (colnames(Z) %in% LayoutPlate[colnames(Z)]) & !(colnames(Z) %in% controls)]
Zpos[['2']] = Z[,LayoutPlate[colnames(Zpos[[1]])]]

set.seed(1)
h = 0.5  # h is equivalent to the sd of a random normal added to each individual value
pvals = Zpos[[1]]*0
for (i in 1:2000000) {
	pos = sample(1:2, 2, replace = T)
	pos = (Zpos[[pos[1]]] + Zpos[[pos[2]]])/2 + rnorm(nrow(Zneg), sd = h)
	neg = rowMeans(Zneg[,sample(1:ncol(Zneg), 2, replace = T)])
	pvals = pvals + sweep(pos, 1, neg, '<=')
}
pvals2 = (pvals + .5)/(i+.5)
#save(pvals2, file='pvals_h0_5.rda')
} else { load('pvals_h0_5.rda') }
#### EXCLUDE Above 

table(colSums(pvals2 <= .05/nrow(pvals2)))

#barplot((colSums(apply(pvals2,2,p.adjust, method='holm') <= .05)))

sigs = names(which(colSums(pvals2 <= .05/nrow(pvals2)) > 2))


Zx = Z[,colnames(Z) %in% c(sigs, LayoutPlate[sigs])]
Zx = sign(Zx)*Zx^2
cors = cor(Zx[rank(-apply(Zx,1,sd)) <= 2000,])
diag(cors) = NA

hmcols = colorRampPalette(c('#0571b0','#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020','#ca0020'))(100)
o1 = seriate(dist(cors^3), method = "OLO")
hm=Heatmap(cors, col = hmcols, use_raster = T, raster_device = 'png', show_heatmap_legend = T, show_row_names = F, show_column_names = F, cluster_columns = as.dendrogram(o1[[1]]), cluster_rows = as.dendrogram(o1[[1]]))

#svg('Hm1_75.svg', width = 10, height = 10)
xx=colnames(cors)[row_order(draw(hm))]
#dev.off()

sum(xx[-1] == LayoutPlate[xx[-length(xx)]])
sum(xx %in% LayoutPlate[xx])/2

samps = c(rep(1:140,2),141:163)
out = rep(0,141)
names(out) = 0:140
for (i in 1:10^8) {
	negs = sapply(1:10^5, function(i) { sum(diff(sample(samps)) == 0) })
	negs = table(negs)
	out[names(negs)] = out[names(negs)] + negs
}

ppois(40,280/303, lower.tail=F)


svg('Hm2.svg', width = 4, height = 4)
Heatmap(cors[xx[150:141],xx[150:141]], col = hmcols, show_heatmap_legend = F, cluster_columns = F, cluster_rows = F)
dev.off()

svg('Hm3.svg', width = 4, height = 4)
Heatmap(cors[xx[29:18],xx[29:18]], col = hmcols, show_heatmap_legend = F, cluster_columns = F, cluster_rows = F)
dev.off()





Z2 = Z[,LayoutPlate[colnames(Z)] %in% colnames(Z)]
p0 = Z2[,grepl('0[2-9]', colnames(Z2)) | grepl('1[0-2]', colnames(Z2))]
p0 = p0[,!(colnames(p0) %in% controls)]
p0 = (p0 + Z2[,LayoutPlate[colnames(p0)]])/sqrt(2)
p0 = cbind(p0, Z[,!(colnames(Z) %in% c(colnames(p0), LayoutPlate[colnames(p0)], controls))])



mat = B4[rowSums(Z[,xx[150:141]] >= 3.8) >= 2, colnames(B4) %in% c(xx[150:141], controls)]

o2 = seriate(dist(cor(t(sign(mat[,1:10])*(mat[,1:10]-2)^2))), method = "OLO")
#svg('Hm4.svg', width = 4, height = 5)
Heatmap(minmax2(mat), col = colorRampPalette(c('#0571b0','#4591c0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), show_heatmap_legend = F, show_row_names = F, cluster_rows = as.dendrogram(o2[[1]]))
#dev.off()




mat = B4[rowSums(Z[,c('D21','M04','J13','G12','M20','D05','E11','L14','G06','J19')] >= 3.8) >= 2, colnames(B4) %in% c(c('D21','M04','J13','G12','M20','D05','E11','L14','G06','J19'), controls)]
mat = mat[,c('D21','M04','M20','D05','J13','G12','E11','L14','G06','J19',colnames(mat)[colnames(mat) %in% controls])]

o2 = seriate(dist(cor(t(sign(mat[,1:10])*(mat[,1:10]-2)^3))), method = "OLO")
#svg('Hm4b.svg', width = 4, height = 5)
Heatmap(minmax2(mat), col = colorRampPalette(c('#0571b0','#4591c0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), show_heatmap_legend = F, show_row_names = F, cluster_rows = as.dendrogram(o2[[1]]), cluster_columns = F)
#dev.off()





IDs = read.csv('ID Mapping.txt')
IDs = IDs[!grepl('Old', IDs[,1]),]
IDs[,1] = sub('[.][1-9]', '', IDs[,1])
IDs[,2] = sub('[.][1-9]', '', IDs[,2])
IDs2 = gsub(' ', '', as.character(IDs[,2]))
names(IDs2) = IDs[,1]

IDs = read.csv('ID Mapping2.txt')
IDs = IDs[!grepl('Old', IDs[,1]),]
IDs[,1] = sub('[.][1-9]', '', IDs[,1])
IDs[,2] = sub('[.][1-9]', '', IDs[,2])
IDs3 = gsub(' ', '', as.character(IDs[,2]))
names(IDs3) = IDs[,1]
IDs3 = IDs3[!duplicated(IDs3)]

IDs = read.csv('ID Mapping3.txt')
IDs = IDs[!grepl('Old', IDs[,1]),]
IDs[,1] = sub('[.][1-9]', '', IDs[,1])
IDs[,2] = sub('[.][1-9]', '', IDs[,2])
IDs4 = gsub(' ', '', as.character(IDs[,2]))
names(IDs4) = IDs[,1]
IDs4 = IDs4[IDs3]
names(IDs4) = names(IDs3)
IDs4 = IDs4[!duplicated(IDs4)]
IDs4 = IDs4[!is.na(IDs4)]



testEnr = function(yy, mat = p0, thresh = 4) {
	yy = unique(na.omit(yy))
	out = sapply(1:ncol(mat), function(i) {
		xx = na.omit(mat[,i])
		xx2 = names(xx)[xx >= thresh]
		yy2 = yy[yy %in% names(xx)]
		#-phyper(sum(!(xx2 %in% yy2)), length(xx) - length(yy2), length(yy2), length(xx2), log.p = T)
		sum(xx2 %in% yy2)
	})
	names(out) = colnames(mat)
	out[order(-out)]/log(10)
}


Salvo = read.csv('Embryogenesis/Salvo 2014 Table S4.csv')
Salvo$genev5 = IDs4[Salvo[,2]]
Salvo = Salvo[!is.na(Salvo$genev5),]

SalvoGn = Salvo$genev5[(Salvo[,8] == 'upregulated') & (as.numeric(Salvo[,4]) >= 2)]

embs = testEnr(SalvoGn)

mat = B4[rownames(B4) %in% SalvoGn, colnames(B4) %in% c(names(embs)[embs >= -log(.05/ncol(p0),10)], LayoutPlate[names(embs)[embs >= -log(.05/ncol(p0),10)]], controls)]
mat = mat[rowSums(p0[rownames(mat),colnames(p0) %in% colnames(mat)] >= 3.5) > 0,]

mapping = match(names(embs), TFinfo[,1])
mapping[is.na(mapping)] = match(LayoutPlate[names(embs)[is.na(mapping)]], TFinfo[,1])
mapping = TFinfo[mapping,4]
mapping[mapping == ''] = names(embs)[mapping == '']
names(embs) = mapping

#svg('Embryogenic enrichment.svg', width = 6, height = 4.5)
barplot(embs[order(-embs)], las=2, border = NA, space=0, col = '#777777', ylim = c(0,25), xaxs='i')
abline(h=-log(.05/ncol(p0),10), lty = 3, col = 'red', lwd = 2)
#dev.off()


o2 = seriate(dist(cor(t(mat))), method = "OLO")
#svg('Hm4.svg', width = 4, height = 5)
Heatmap(minmax2(mat), col = colorRampPalette(c('#0571b0','#4591c0','#92c5de','#f7f7f7','#f4a582','#ca0020'))(100), show_heatmap_legend = F, show_row_names = F, cluster_rows = as.dendrogram(o2[[1]]))
#dev.off()



heatmap(p0[rownames(p0) %in% SalvoGn, colnames(p0) %in% c('D21','M04','J13','G12','M20','D05','E11','L14','G06','J19')] >= 3.5)



Pnorm = read.csv('pooledNormalized.csv', row.names=1)
Pnorm2 = apply(as.matrix(Pnorm[,seq(1,ncol(Pnorm),4)]),2,as.numeric)
rownames(Pnorm2) = IDs2[rownames(Pnorm)]
colnames(Pnorm2) = paste('Norm', 1:6, sep = '')

svg('Meiotic enrichment.svg', width = 6, height = 4.5)
barplot(testEnr(rownames(Pnorm2)[Pnorm2[,1] >= 1]), las=2, border = NA, space=0, col = '#777777', ylim = c(0,8), xaxs='i')
abline(h=-log(.05/ncol(p0),10), lty = 3, col = 'red', lwd = 2)
dev.off()

testEnr(rownames(Pnorm2)[which(Pnorm2[,1] >= 1)])[1:20]




