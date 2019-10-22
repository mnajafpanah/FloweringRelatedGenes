#### Ball Helix Bioinformatics Scientist Assignmet ###
# Purpose: Identify potential genes that are responsible for a flowering trait 
#          through Differential Gene Expression (DGE) analysis by given a pair-ended 
#          RNA-Seq expresion data using the edgeR, limma, Glimma, etc. 
# Author: Mohammad Najaf-panah
# Date: October 2019

# Loading required packages
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(ggrepel)
library(plyr)

# Set working directory
setwd('/Users/mohammadjavadnajafpanah/Documents/Ball_assignment/')

# Loading expression matrix
data <- read.table("test_gene_expression_matrix.txt", header=TRUE, row.names=1)
data <- cbind(GENENAME=rownames(data), data) # add GENENAME as 1st column
dim(data) 
# So, we have 29333 genes in this plant and 2 samples with 3 biological replicates each.

# Assign conditions 
group <-factor(c(rep("Bud_Rep", 3), rep("Day2_Rep", 3)))

# Stores data in a simple list-based data object called a DGEList 
# NOTE: the '$counts' variable is actually the TPM (Transcript per Million) float values.
deg <- DGEList(data[,-1], group=group, genes = data[,1,drop=FALSE])
options(digits=3)
deg

### Differential expression analysis 
# Creating a design matrix and contrasts
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
deg$design <- design
deg$design
dim(design) # 6 2


contr.matrix <- makeContrasts(
  
  Bud_RepVsDay2_Rep = Bud_Rep - Day2_Rep,
  
  levels = colnames(design))
contr.matrix

## FYI
# A lot of genes are unexpressed throughout all samples. 
table(rowSums(deg$counts==0)==6) 
# FALSE  TRUE 
# 23308  6025  
# In fact, 6025 (21%) of genes (29333) in this dataset have zero counts across all 6 samples.

# Filtering the genes using TPM.filter method by Guo,W. et al. 2019 (ThreeDRNAseq package)
TPM.filter <- function(TPM,sample.n=3,tpm.cut=1){
  rowSums(TPM>=tpm.cut)>=sample.n
}
remove <- TPM.filter(deg$counts, sample.n = 6, tpm.cut = 1)
table(remove)
# FALSE  TRUE 
# 14476 14857

# The DGEList object is then subsetted to retain only the non-filtered genes:
deg <- deg[which(remove), ]
dim(deg) # 14857     6
# Therefore, the number of genes is reduced to 14857, about 51% of the number that we started with.

# Dispersion Estimation
deg <- estimateDisp(deg, design, robust=TRUE)
plotBCV(deg)

fit <- glmQLFit(deg, design, robust=TRUE)
head(fit$coefficients)

plotQLDisp(fit)

summary(fit$df.prior)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.30    3.74    3.74    3.70    3.74    3.74 

## Normalization for composition bias by trimmed mean of M values (TMM) (Robinson and Oshlack 2010).
# calcNormFactors(deg) returns the DGEList argument with only the norm.factors changed.
deg <- calcNormFactors(deg)
deg$samples$norm.factors
# 1.094 1.069 1.077 0.914 0.924 0.940

# Fitting linear models for comparisons of interest
v <- voom(deg, design, plot=TRUE)
v
arrayw <- arrayWeightsSimple(v, design) # Use of array weights increases the significance of the top genes 
                                        # and try to partially recover the library sizes. 
vfit <- lmFit(v, design, weights=arrayw)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, trend = TRUE)
plotSA(efit, main="Final model: Mean-variance trend")

# Examining the number of DE genes
summary(decideTests(efit, p.value = 0.01))
#        Bud_RepVsDay2_Rep p-value=0.01
# Down                4102
# NotSig              6638
# Up                  4117

dt <- decideTests(efit, p.value = 0.01)
dim(dt) # 14857     1

# stor the result tsble
write.fit(efit, dt, file="results.txt")

# Examining individual DE genes from top to bottom
Bud_RepVsDay2_Rep <- topTreat(efit, coef=1, n=Inf)
head(Bud_RepVsDay2_Rep)
dim(Bud_RepVsDay2_Rep) # 14857     7

# Useful graphical representations of differential expression results
plotMD(efit, column=1, status=dt[,1], main=colnames(efit)[1], xlim=c(-8,13))

glMDPlot(efit, coef=1, status=dt, main=colnames(efit)[1],
         side.main="GENENAME", counts=deg$counts, groups=group, launch=TRUE)

# Coolmap:
# coolmap clusters genes by correlation (highly correlated genes are closest) 
# and clusters samples based on Euclidean distance between the expression values. 
ltpm <- log(deg$counts) # log(TPM)
topDEgenes <- Bud_RepVsDay2_Rep$GENENAME[1:30] # top 30 DE genes
i <- which(v$genes$GENENAME %in% topDEgenes)
pdf(file = "coolmap.pdf", width = 9, height = 18)
coolmap(ltpm[i,], margins=c(8,10), lhei=c(1,6), lwid=c(1,3))
dev.off()


## Construct the volcano plot object
Bud_RepVsDay2_Rep[,8] <- ifelse((Bud_RepVsDay2_Rep$P.Value < 0.01 & abs(Bud_RepVsDay2_Rep$logFC) > 0.5), "red", "black")
dim(Bud_RepVsDay2_Rep)
size <- ifelse((Bud_RepVsDay2_Rep$P.Value < 0.01 & abs(Bud_RepVsDay2_Rep$logFC) > 1.5), 4, 2)
vp <- ggplot(data=Bud_RepVsDay2_Rep, aes(x=Bud_RepVsDay2_Rep[,2], y=-log10(Bud_RepVsDay2_Rep$P.Value))) +
  geom_point(size=size, colour=Bud_RepVsDay2_Rep[,8]) +
  xlim(c(-3, 3)) + ylim(c(0,8)) +
  xlab("log2 fold change") + ylab("-log10  p-value") +
  guides(colour = guide_legend(override.aes = list(shape=16)))
vp









