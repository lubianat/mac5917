# This script computes DEGs from the same cells in two conditions in scRNA-Seq using MAST
# Written by Tiago Lubiana, 24/08/2018

library(data.table)
library(MAST)
library(NMF)
library(rsvd)
library(ggplot2)
library(GGally)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(edgeR)

# load data ---------------
load(file = './int/data/train_test_sets_dengue.Rdata')
full_dengue_train <- rbind(expr_control_train, expr_perturbed_train)
annotation <- fread('./int/data/annot_dengue.csv')
annotation <-  annotation[annotation$title %in% rownames(full_dengue_train)]

#scaRaw <- FromMatrix(t(maits$expressionmat), maits$cdat, maits$fdat)
freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)

# remove all zeros
df <- t(full_dengue_train)
df <- df[apply(df[,-1], 1, function(x) !all(x==0)),]
full_dengue_train <- as.data.frame(t(df))

# transforming to log cpm
full_dengue_train <- as.matrix(full_dengue_train)
full_dengue_train_cpm <- as.data.frame(log2(cpm(full_dengue_train) + 1))

#reordering to match annotation
full_dengue_train_cpm <- full_dengue_train_cpm[annotation$title,]

scaRaw <-  FromMatrix(exprsArray = t(full_dengue_train_cpm))

# exploratory data analysis -------
aheatmap(
  assay(scaRaw[1:1000, ]),
  labRow = '',
  distfun = 'spearman'
)

set.seed(123)
plotPCA <- function(sca_obj) {
  projection <- rpca(t(assay(sca_obj)), retx = TRUE, k = 4)$x
  colnames(projection) = c("PC1", "PC2", "PC3", "PC4")
  pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
  print(ggpairs(
    pca,
    columns = c(
      'PC1',
      'PC2',
      'PC3'
    )
  ))
  invisible(pca)
}

plotPCA(scaRaw)

# filtering the data ---------
sca <- scaRaw

sca <- sca[sample(which(freq(sca) > 0),6000),]

cdr2 <- colSums(assay(sca) > 0)

colData(sca)$cngeneson <- scale(cdr2)

plotPCA(sca)



scaSample <- sca[sample(which(freq(sca) > .1), 20),]

flat <- as(scaSample, 'data.table')
ggplot(flat, aes(x = value)) + geom_density() + facet_wrap(~ primerid, scale =
                                                             'free_y')
thres <-
  thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
par(mfrow = c(5, 4))
plot(thres)


assays(sca) <- list(thresh = thres$counts_threshold, tpm = assay(sca))
expressed_genes <- freq(sca) > freq_expressed
sca <- sca[expressed_genes, ]


cond <- factor(colData(sca)$condition)
cond <- relevel(cond, "Unstim")

colData(sca)$condition <- as.factor(ifelse(annotation$`moi:ch1` == 0, 'ctrl', 'infected'))

zlmCond <- zlm( ~ condition + cngeneson, sca)

summaryCond <- summary(zlmCond, doLRT = 'conditioninfected')

print(summaryCond, n = 4)

summaryDt <- summaryCond$datatable

fcHurdle <-
  merge(summaryDt[contrast == 'conditioninfected' &
                    component == 'H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
        summaryDt[contrast == 'conditioninfected' &
                    component == 'logFC', .(primerid, coef, ci.hi, ci.lo)], by = 'primerid') #logFC coefficients

fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <-
  merge(fcHurdle[fdr < .05 &
                   abs(coef) > FCTHRESHOLD], as.data.table(mcols(sca)), by = 'primerid')
setorder(fcHurdleSig, fdr)


entrez_to_plot <- fcHurdleSig[1:20, primerid]
symbols_to_plot <- fcHurdleSig[1:50, symbolid]
flat_dat <- as(sca[entrez_to_plot, ], 'data.table')
ggbase <-
  ggplot(flat_dat, aes(x = condition, y = thresh, color = condition)) + geom_jitter() +
  facet_wrap( ~ primerid, scale = 'free_y') + ggtitle("DE Genes in Activated MAIT Cells")
ggbase + geom_violin()


flat_dat[, lmPred := lm(thresh ~ cngeneson + condition)$fitted, key = symbolid]
ggbase + aes(x = cngeneson) + geom_line(aes(y = lmPred), lty = 1) + xlab('Standardized Cellular Detection Rate')

## This is all rather kludgy at the moment
MM <- model.matrix(~condition,unique(colData(sca)[,c("condition"),drop=FALSE]))
rownames(MM) <- str_extract(rownames(MM), 'Stim|Unstim')
predicted <- predict(zlmCond,modelmatrix=MM)

## Avert your eyes...
predicted[, primerid:=as.character(primerid)]
predicted_sig <- merge(mcols(sca), predicted[primerid%in%entrez_to_plot], by='primerid')
predicted_sig <- as.data.table(predicted_sig)

## plot with inverse logit transformed x-axis
ggplot(predicted_sig)+aes(x=invlogit(etaD),y=muC,xse=seD,yse=seC,col=sample)+
  facet_wrap(~symbolid,scales="free_y")+theme_linedraw()+
  geom_point(size=0.5)+scale_x_continuous("Proportion expression")+
  scale_y_continuous("Estimated Mean")+
  stat_ell(aes(x=etaD,y=muC),level=0.95, invert='x')

mat_to_plot <- assay(sca[entrez_to_plot,])
rownames(mat_to_plot) <- symbols_to_plot
aheatmap(mat_to_plot,annCol=colData(sca)[,"condition"],main="DE genes",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)))



# AUXILAIRY FUNCTIONS -------
write.csv(fcHurdleSig, './int/results/DE_MAST_denge_train.csv')
