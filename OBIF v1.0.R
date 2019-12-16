# Omics-Based Interaction Framework (OBIF) v1.0 Code ----
# This intended for a self-guided analysis of an Omics dataset with OBIF. For additional inquiries please contact Evans Laboratory in The University of Texas MD Anderson Cancer Center by e-mailing Dr. Scott E. Evans (seevans@mdanderson.org).

# Citing information ----
# Do not forget to cite the use of this code from the following publication:

# Install and/or load the following packages if needed ----
library(readxl)
library(ClassDiscovery)
library(reshape2)
library(ggplot2)
library(devtools)
library(limma)
library(Biobase)
library(lumi)
library(preprocessCore)
library(cowplot)
library(gridGraphics)
library(ggplotify)
library(ggdendro)
library(dendextend)
library(gridExtra)
library(car)
library(stats)
library(rafalib)
library(rmarkdown)
library(knitr)
library(gdata)
library(WriteXLS)
library(readxl)
library(reporttools)
library(survival)
library(DT)
library(colorspace)
library(RColorBrewer)
library(ClassComparison)
library(ClassDiscovery)
library(limma)
library(gplots)
library(kableExtra)
library(GEOquery)
library(Biobase)
library(lumi)
library(ggplot2)
library(limma)
library(heatmap.plus)
library(heatmap3)
library(writexl)
library(devtools)
library(ggfortify)
library(gtools)
library(RColorBrewer)
library(eulerr)
library(ggbiplot)
library(VennDiagram)
library(readxl)
library(ggplot2)
library(devtools)
library(Vennerable)
library(venneuler)
library(rJava)
library(eulerr)


# Load and organize your analysis-ready Omics dataset ----
# Manually import your dataset into your R Environment from any source (i.e., GEO, Excel, CSV, etc.).
HML_AB090342_Metab_Analysis_Ready_ <- read_excel("Documents/R/Metabolomics/HML AB090342 Metab (Analysis-Ready).xlsx")

# Rename your dataset and define factors
data.orig <- HML_AB090342_Metab_Analysis_Ready_
factor.A <- "Colistin"
factor.B <- "Sulbactam"
factor.AB <- "Colistin + Sulbactam"

# Identify the column with the name of the features
feat.ids <- data.orig[,5]
colnames(feat.ids) <- "Feat.ID"

# Identify the columns with the expression values per sample
expr.val <- data.orig[,10:29]

# Identify the columns with additional information about your dataset
add.col <- data.orig[,1:9] 

# Preserve feature names in all subsets of your datasets
rownames(data.orig) <- feat.ids
rownames(feat.ids) <- feat.ids
rownames(expr.val) <- feat.ids
rownames(add.col) <- feat.ids

# Indicate number of replicates per group
n.ctrl <- 5
n.facA <- 5
n.facB <- 5
n.facAB <- 5

# Set color pallete for samples
col.samples <- c(rep("#000000", n.ctrl), rep("#33a02c", n.facA), rep("#ff7f00", n.facB), rep("#1f78b4", n.facAB))

# Create and set column names for samples
cn.ctrl <- paste("ctrl.", seq(1:n.ctrl), sep = "")
cn.facA <- paste("facA.", seq(1:n.facA), sep = "")
cn.facB <- paste("facB.", seq(1:n.facB), sep = "")
cn.facAB <- paste("facAB.", seq(1:n.facAB), sep = "")
cn.expr.val <- c(cn.ctrl, cn.facA, cn.facB, cn.facAB)
colnames(expr.val) <- cn.expr.val

# Quality control review of your dataset ----
## 1. Review data distribution and preliminary model fitting of your original dataset
### Simple dataset transformations (log2 transformation and quantile normalization if needed)
sim.l2.expr.val <- log2(expr.val)
sim.l2qn.expr.val <- normalize.quantiles(data.matrix(sim.l2.expr.val))
sim.l2qn.expr.val <- as.data.frame(sim.l2qn.expr.val)
colnames(sim.l2qn.expr.val) <- cn.expr.val

### Array-based dataset transformations (background correction, log2 transformation and quantile normalization if needed)
arr.eSet <- new("ExpressionSet", exprs=as.matrix(expr.val))
arr.bg.eSet <- lumiB(arr.eSet, method = "bgAdjust.affy")
arr.bgl2.eSet <- lumiT(arr.bg.eSet, "log2")
arr.bgl2qn.eSet <- lumiN(arr.bgl2.eSet, method = "quantile")
arr.bgl2qn.expr.val <- exprs(arr.bgl2qn.eSet)
arr.bgl2qn.expr.val <- as.data.frame(arr.bgl2qn.expr.val)
colnames(arr.bgl2qn.expr.val) <- cn.expr.val

## 1.A Using a violin plot
vp.expr.val.orig <- melt(expr.val, id.vars = NULL)
vp.expr.val.sim <- melt(sim.l2qn.expr.val, id.vars = NULL)
vp.expr.val.arr <- melt(arr.bgl2qn.expr.val, id.vars = NULL)

qc.vp.orig <- ggplot(vp.expr.val.orig, aes(x = variable, y = value, color = variable)) + geom_violin() + geom_boxplot(width = 0.2) + scale_color_manual(values = col.samples) + ggtitle("Original Dataset") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") 
qc.vp.sim <- ggplot(vp.expr.val.sim, aes(x = variable, y = value, color = variable)) + geom_violin() + geom_boxplot(width = 0.2) + scale_color_manual(values = col.samples) + ggtitle("Simple Transformations") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
qc.vp.arr <- ggplot(vp.expr.val.arr, aes(x = variable, y = value, color = variable)) + geom_violin() + geom_boxplot(width = 0.2) + scale_color_manual(values = col.samples) + ggtitle("Array-based Transformations") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")

pdf("HML AB090342 Met - QC 1A.pdf")
qc.vp.orig
dev.off()

pdf("HML AB090342 Met - QC 1B.pdf")
qc.vp.sim
dev.off()

pdf("HML AB090342 Met - QC 1C.pdf")
qc.vp.arr
dev.off()

# 1.B Using a density plot
qc.dp.orig <- ggplot(vp.expr.val.orig, aes(value, color = variable)) + geom_density() + scale_color_manual(values = col.samples) + ggtitle(label="Original Dataset", subtitle = deparse(substitute(vp.expr.val.orig))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top") 
qc.dp.sim <- ggplot(vp.expr.val.sim, aes(value, color = variable)) + geom_density() + scale_color_manual(values = col.samples) + ggtitle(label="Simple Transformations", subtitle = deparse(substitute(vp.expr.val.sim))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top") 
qc.dp.arr <- ggplot(vp.expr.val.arr, aes(value, color = variable)) + geom_density() + scale_color_manual(values = col.samples) + ggtitle(label="Array-based Transformations", subtitle = deparse(substitute(vp.expr.val.arr))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top") 


# 1.C Using a box plot 
pdf("HML AB090342 Met - QC 2A.pdf")
qc.bxpt.orig <- boxplot(expr.val, boxwex=0.5, notch=TRUE, outline=FALSE, las =2, border = col.samples, main = "Original Dataset")
dev.off()

pdf("HML AB090342 Met - QC 2B.pdf")
qc.bxpt.sim <- boxplot(sim.l2qn.expr.val, boxwex=0.5, notch=TRUE, outline=FALSE, las =2, border = col.samples, main = "Simple Transformations")
dev.off()

pdf("HML AB090342 Met - QC 2C.pdf")
qc.bxpt.ar
r <- boxplot(arr.bgl2qn.expr.val, boxwex=0.5, notch=TRUE, outline=FALSE, las =2, border = col.samples, main = "Array-based Transformations")
dev.off()

## 1.D Using hierarchical clutersing
hcl.orig <- hclust(distanceMatrix(
  dataset =  expr.val, 
  metric =  "pearson"), 
  method = "ward.D2")
hcl.sim <- hclust(distanceMatrix(
  dataset =  sim.l2qn.expr.val, 
  metric =  "pearson"), 
  method = "ward.D2")
hcl.arr <- hclust(distanceMatrix(
  dataset =  arr.bgl2qn.expr.val, 
  metric =  "pearson"), 
  method = "ward.D2")

pdf("HML AB090342 Met - QC 3A.pdf")
qc.hcl.orig <- myplclust(hcl.orig, labels = hcl.orig$labels, lab.col = col.samples) 
dev.off()

pdf("HML AB090342 Met - QC 3B.pdf")
qc.hcl.sim <- myplclust(hcl.sim, labels = hcl.orig$labels, lab.col = col.samples)
dev.off()

pdf("HML AB090342 Met - QC 3C.pdf")
qc.hcl.arr <- myplclust(hcl.arr, labels = hcl.orig$labels, lab.col = col.samples)
dev.off()

# Quality control and linear model for interaction statistics ----
### Create factor design for interaction statistics
#### One-way ANOVA dataframe and design
si.orig <- data.frame(Factors=sapply(colnames(expr.val),function(x)substr(x,1,regexpr('\\.',x)[1]-1)))
si.orig$Factors <- factor(si.orig$Factors, levels=c('ctrl','facA','facB','facAB'))
si <- si.orig
design <- model.matrix(~1+si$Factors)
colnames(design) <- levels(si$Factors)

#### Two-way ANOVA dataframe and design
si.orig2x2 <- data.frame("facA" = c(rep("0", n.ctrl), rep("1", n.facA), rep("0", n.facB), rep("1", n.facAB)),
                         "facB" = c(rep("0", n.ctrl), rep("0", n.facA), rep("1", n.facB), rep("1", n.facAB)))
si.orig2x2$facA <- factor(si.orig2x2$facA, levels=c("0", "1"))
si.orig2x2$facB <- factor(si.orig2x2$facB, levels=c("0", "1"))
si2x2 <- si.orig2x2
design2x2 <- model.matrix(~1+si2x2$facA*si2x2$facB)
colnames(design2x2) <- cbind('ctrl','facA','facB','Inter')

### Interaction Statistics
#### Data set up for linear modeling of feature expression
ints.mtx.orig <- as.matrix(expr.val)
ints.mtx.sim <- as.matrix(sim.l2qn.expr.val)
ints.mtx.arr <- as.matrix(arr.bgl2qn.expr.val)

ints.des2x2 <- design2x2

pv.inter <- vector("numeric", length = nrow(ints.mtx.orig))
for (i in 1:nrow(ints.mtx.orig)) {ints.tmp.orig <- data.frame(exp = ints.mtx.orig[i, ], ints.des2x2)}
for (i in 1:nrow(ints.mtx.sim)) {ints.tmp.sim <- data.frame(exp = ints.mtx.sim[i, ], ints.des2x2)}
for (i in 1:nrow(ints.mtx.arr)) {ints.tmp.arr <- data.frame(exp = ints.mtx.arr[i, ], ints.des2x2)}

#### Fitting the expression model for interaction statistics
ints.fit.orig <- lm(exp ~ 0 + facA + facB + facA*facB, data = ints.tmp.orig) 
sum.ints.fit.orig <- summary(ints.fit.orig)
ints.fit.sim <- lm(exp ~ 0 + facA + facB + facA*facB, data = ints.tmp.sim) 
sum.ints.fit.sim <- summary(ints.fit.sim)
ints.fit.arr <- lm(exp ~ 0 + facA + facB + facA*facB, data = ints.tmp.arr) 
sum.ints.fit.arr <- summary(ints.fit.arr)

sum.ints.fit.orig
sum.ints.fit.sim
sum.ints.fit.arr

#### Plotting interaction plots 
par(mfrow = c(2,3))
ints.plot1.orig <- with(ints.tmp.orig, interaction.plot(facA, facB, exp, main="Original Dataset (A vs B)"))
ints.plot1.sim <- with(ints.tmp.sim, interaction.plot(facA, facB, exp, main="Simple Transformations (A vs B)"))
ints.plot1.arr <- with(ints.tmp.arr, interaction.plot(facA, facB, exp, main="Array-based Transformations (A vs B)"))

ints.plot2.orig <- with(ints.tmp.orig, interaction.plot(facB, facA, exp, main="Original Dataset (B vs A)"))
ints.plot2.sim <- with(ints.tmp.sim, interaction.plot(facB, facA, exp, main="Simple Transformations (B vs A)"))
ints.plot2.arr <- with(ints.tmp.arr, interaction.plot(facB, facA, exp, main="Array-based Transformations (B vs A)"))

par(mfrow = c(1,1))

#### Diagnostic plots for linear model fitting of the interaction statistics 
# Supplement the data fitted to a linear model with the model fit statistics
all.for.orig <- fortify(ints.fit.orig)
all.for.sim <- fortify(ints.fit.sim)
all.for.arr <- fortify(ints.fit.arr)

# Retrieve the model fit statistics for the variables without the intercept (facA, facB and facA:facB). 
# In the interaction statistics fit, the ctrl variables are considered the intercept ("exp ~ 0 ...").
# The weights of the ctrl variables in the regresion are considered in the design matrix.
mod.ints.orig <- all.for.orig[-c(1:n.ctrl),]
mod.ints.sim <- all.for.sim[-c(1:n.ctrl),]
mod.ints.arr <- all.for.arr[-c(1:n.ctrl),]

# Modified function of diagPlots from Raju Rimal to represent diagnostic plots of lm models without intercepts
lm0.diagPlot<-function(model){
  p1<-ggplot(model, aes(.fitted, .resid))+geom_point()
  p1<-p1+stat_smooth(method="lm")+geom_hline(yintercept=0, col="red", linetype="dashed")
  p1<-p1+xlab("Fitted values")+ylab("Residuals")
  p1<-p1+ggtitle(label="Residual vs Fitted Plot", subtitle = deparse(substitute(model)))+theme_bw()
  
  p2<-ggplot(model, aes(qqnorm(.stdresid)[[1]], .stdresid))+geom_point(na.rm = TRUE)
  p2<-p2+geom_abline(aes(qqline(.stdresid)))+xlab("Theoretical Quantiles")+ylab("Standardized Residuals")
  p2<-p2+ggtitle(label="Normal Q-Q Plot", subtitle = deparse(substitute(model)))+theme_bw()
  
  p3<-ggplot(model, aes(.fitted, sqrt(abs(.stdresid))))+geom_point(na.rm=TRUE)
  p3<-p3+stat_smooth(method="lm", na.rm = TRUE)+xlab("Fitted Value")
  p3<-p3+ylab(expression(sqrt("|Standardized residuals|")))
  p3<-p3+ggtitle(label="Scale-Location", subtitle = deparse(substitute(model)))+theme_bw()
  
  p4<-ggplot(model, aes(seq_along(.cooksd), .cooksd))+geom_bar(stat="identity", position="identity")
  p4<-p4+xlab("Obs. Number")+ylab("Cook's distance")
  p4<-p4+ggtitle(label="Cook's distance", subtitle = deparse(substitute(model)))+theme_bw()
  
  p5<-ggplot(model, aes(.hat, .stdresid))+geom_point(aes(size=.cooksd), na.rm=TRUE)
  p5<-p5+stat_smooth(method="lm", na.rm=TRUE)
  p5<-p5+xlab("Leverage")+ylab("Standardized Residuals")
  p5<-p5+ggtitle(label="Residual vs Leverage Plot", subtitle = deparse(substitute(model)))
  p5<-p5+scale_size_continuous("Cook's Distance", range=c(1,5))
  p5<-p5+theme_bw()+theme(legend.position="bottom")
  
  p6<-ggplot(model, aes(.hat, .cooksd))+geom_point(na.rm=TRUE)+stat_smooth(method="lm", na.rm=TRUE)
  p6<-p6+xlab("Leverage hii")+ylab("Cook's Distance")
  p6<-p6+ggtitle(label="Cook's dist vs Leverage hii/(1-hii)", subtitle = deparse(substitute(model)))
  p6<-p6+geom_abline(slope=seq(0,3,0.5), color="gray", linetype="dashed")
  p6<-p6+theme_bw()
  
  return(list(ResvsFit=p1, NormQQ=p2, ScLoc=p3, CookD=p4, ResvsLev=p5, CookvsLev=p6))
}

# Diagnostics plots of the interaction statistics fitted linear model with interception samples
qc.diagp.orig.all <- lm0.diagPlot(all.for.orig)
qc.diagp.sim.all <- lm0.diagPlot(all.for.sim)
qc.diagp.arr.all <- lm0.diagPlot(all.for.arr)

qc.diagp.orig.all$ResvsFit
qc.diagp.orig.all$NormQQ
qc.diagp.orig.all$ScLoc
qc.diagp.orig.all$CookD
qc.diagp.orig.all$ResvsLev
qc.diagp.orig.all$CookvsLev

qc.diagp.sim.all$ResvsFit
qc.diagp.sim.all$NormQQ
qc.diagp.sim.all$ScLoc
qc.diagp.sim.all$CookD
qc.diagp.sim.all$ResvsLev
qc.diagp.sim.all$CookvsLev

qc.diagp.arr.all$ResvsFit
qc.diagp.arr.all$NormQQ
qc.diagp.arr.all$ScLoc
qc.diagp.arr.all$CookD
qc.diagp.arr.all$ResvsLev
qc.diagp.arr.all$CookvsLev

# Diagnostics plots of the interaction statistics fitted linear model without interception samples
qc.diagp.orig <- lm0.diagPlot(mod.ints.orig)
qc.diagp.sim <- lm0.diagPlot(mod.ints.sim)
qc.diagp.arr <- lm0.diagPlot(mod.ints.arr)

qc.diagp.orig$ResvsFit
qc.diagp.orig$NormQQ
qc.diagp.orig$ScLoc
qc.diagp.orig$CookD
qc.diagp.orig$ResvsLev
qc.diagp.orig$CookvsLev

qc.diagp.sim$ResvsFit
qc.diagp.sim$NormQQ
qc.diagp.sim$ScLoc
qc.diagp.sim$CookD
qc.diagp.sim$ResvsLev
qc.diagp.sim$CookvsLev

qc.diagp.arr$ResvsFit
qc.diagp.arr$NormQQ
qc.diagp.arr$ScLoc
qc.diagp.arr$CookD
qc.diagp.arr$ResvsLev
qc.diagp.arr$CookvsLev

# Quality control and linear model for expression and contrast analysis ----
# Make contrast comparisons for both One-way and Two-way ANOVAs
cont.mtx.des <- makeContrasts(SME.facA = "facAB-facB", SME.facB = "facAB-facA", Inter.fAxfB = "facAB-facA-facB",  levels=design)
cont.mtx.des2x2 <- makeContrasts(SME.facA = "Inter+facA", SME.facB = "Inter+facB", facAB = "Inter+facA+facB",  levels= design2x2)

# Calculate the statistics for expression and contrast analysis using One-way ANOVA
fit <- lmFit(ints.mtx.arr,design)
fit.ebayes <- eBayes(fit)
adj.Pvalues.lmBH <- apply(fit.ebayes$p.value, 2, p.adjust, method="BH")
adj.Pvalues.lmBonf <- apply(fit.ebayes$p.value, 2, p.adjust, method="bonferroni")

fit2 <- contrasts.fit(fit,cont.mtx.des)
fit2 <- eBayes(fit2)
adj.Pvalues.contrBH <- apply(fit2.eBayes$p.value, 2, p.adjust, method="BH")
adj.Pvalues.contrBonf <- apply(fit2.eBayes$p.value, 2, p.adjust, method="bonferroni")

# Calculate the statistics for expression and contrast analysis using Two-way ANOVA
fit.2x2 <- lmFit(ints.mtx.arr,design2x2)
fit.2x2.ebayes <- eBayes(fit.2x2)
adj.2x2.Pvalues.lmBH <- apply(fit.2x2.ebayes$p.value, 2, p.adjust, method="BH")
adj.2x2.Pvalues.lmBonf <- apply(fit.2x2.ebayes$p.value, 2, p.adjust, method="bonferroni")

fit2.2x2 <- contrasts.fit(fit.2x2,cont.mtx.des2x2)
fit2.2x2 <- eBayes(fit2.2x2)
adj.2x2.Pvalues.contrBH <- apply(fit2.2x2.eBayes$p.value, 2, p.adjust, method="BH")
adj.2x2.Pvalues.contrBonf <- apply(fit2.2x2.eBayes$p.value, 2, p.adjust, method="bonferroni")

# Make sure that all p-values are calculated equally between One-way and Two-way ANOVAs
comp.pvals <- cbind(fit.ebayes$p.value, fit.2x2.ebayes$p.value, fit2.eBayes$p.value, fit2.2x2.eBayes$p.value) 

# Expression value calculations ----
m.l2.vals <- data.frame(m.ctrl = apply(ints.mtx.arr[, si$Factors == "ctrl"], 1, mean), 
                      m.facA = apply(ints.mtx.arr[, si$Factors == "facA"], 1, mean), 
                      m.facB = apply(ints.mtx.arr[, si$Factors == "facB"], 1, mean),
                      m.facAB = apply(ints.mtx.arr[, si$Factors == "facAB"], 1, mean))

d.log2FC <- data.frame(d.facA = m.l2.vals[, "m.facA"] - m.l2.vals[, "m.ctrl"],
                   d.facB = m.l2.vals[, "m.facB"] - m.l2.vals[, "m.ctrl"], 
                   d.facAB = m.l2.vals[, "m.facAB"] - m.l2.vals[, "m.ctrl"])

# Assign Interaction Score ----
# Calculate Interaction Score
IntScore <- cbind.data.frame(RawIScore = (d.log2FC[,3])/(d.log2FC[,1]+d.log2FC[,2]),
                             AbsIScore = abs((d.log2FC[,3])/(d.log2FC[,1]+d.log2FC[,2])))

# Assign double IS category
IntScore$DoubleIS <- ifelse(IntScore$AbsIScore > 1, "Synergistic.Dual", 
                            ifelse(IntScore$AbsIScore < 1, "Antagonistic.Dual", "Error4"))
sum(IntScore$DualIS == "Error4")

# Assign triple IS category
SynCutoff <- 1.1
AntCuoff <- 0.9
IntScore$TripleIS <- ifelse(IntScore$AbsIScore > SynCutoff, "Synergistic", 
                            ifelse(IntScore$AbsIScore < AntCuoff, "Antagonistic",
                                   ifelse(IntScore$AbsIScore < SynCutoff & IntScore$AbsIScore > AntCuoff, "Additive", "Error5")))
sum(IntScore$TripleIS == "Error4")

# Assign Expression Profiles ----
# Downregulated or Upregulated per treatment
ExpProf <- data.frame(UD.A = ifelse (d.log2FC[,1] < 0 , "Down", ifelse ( d.log2FC[,1] > 0, "Up", "Error1")),
                      UD.B = ifelse ( d.log2FC[,2] < 0 , "Down", ifelse ( d.log2FC[,2] > 0, "Up", "Error1")),
                      UD.AB = ifelse ( d.log2FC[,3] < 0 , "Down", ifelse ( d.log2FC[,3] > 0, "Up", "Error1")))
# Verify 0 or NA values
sum(ExpProf$UD.A == "Error1")
sum(ExpProf$UD.B == "Error1")
sum(ExpProf$UD.AB == "Error1")

# Assign & Verify SingleTx group
ExpProf$SingleTx <- ifelse(ExpProf[,1] == ExpProf[,2], "Cooperative", 
                           ifelse(ExpProf[,1] != ExpProf[,2], "Competitive", "Error2"))
sum(ExpProf$SingleTx == "Error2")

# Assign & Verify DualTx group
ExpProf$DualTx <- ifelse(ExpProf$SingleTx == "Cooperative" & ExpProf[,1] == ExpProf[,3], "Concordant", 
                         ifelse(ExpProf$SingleTx == "Cooperative" & ExpProf[,1] != ExpProf[,3], "Discordant", 
                                ifelse(ExpProf$SingleTx == "Competitive" & ExpProf[,1] == ExpProf[,3], "facA-Dominant", 
                                       ifelse(ExpProf$SingleTx == "Competitive" & ExpProf[,2] == ExpProf[,3], "facB-Dominant","Error3"))))
sum(ExpProf$DualTx == "Error3")

# Assign SingleTx & Dual group related to expression in dual exposure
ExpProf$UDSingleTx <- paste(ExpProf$SingleTx, ExpProf$UD.AB, sep = ".")
ExpProf$UDDualTx <- paste(ExpProf$DualTx, ExpProf$UD.AB, sep = ".")

# Final data compilation ----
final.res <- data.frame(feat.ids,
                        add.col,
                        arr.bgl2qn.expr.val,
                        m.l2.vals, 
                        d.log2FC,
                        pv=fit.ebayes$p.value[,2:4],
                        pv=fit2$p.value,
                        adjpvBH=adj.Pvalues.lmBH[,2:4],
                        adjpvBH=adj.Pvalues.contrBH, 
                        adjpvBonf=adj.Pvalues.lmBonf[,2:4],
                        adjpvBonf=adj.Pvalues.contrBonf, 
                        IntScore,
                        ExpProf)

# Venn Diagram of DEMs / iDEMs ----
# Select DEMs per group
A.DEMs <- subset(final.res, 
                 d.facA >= 0.3 | 
                   d.facA <= -0.3)
A.FDR.DEMs <- subset(A.DEMs,
                     adjpvBH.facA <= 0.1)

B.DEMs <- subset(final.res, 
                 d.facB >= 0.3 | 
                   d.facB <= -0.3)
B.FDR.DEMs <- subset(B.DEMs,
                     adjpvBH.facB <= 0.1)

AB.DEMs <- subset(final.res, 
                  d.facAB >= 0.3 | 
                    d.facAB <= -0.3)
AB.FDR.DEMs <- subset(AB.DEMs,
                      adjpvBH.facAB <= 0.1)
AB.iDEMs <- subset(AB.FDR.DEMs, 
                   pv.Inter.fAxfB < 0.05)

# Create list of vectors for DEMs per group for Venn Diagram
venn.list <- list(facA.dems = A.FDR.DEMs$PeakID, 
                  facB.dems = B.FDR.DEMs$PeakID, 
                  facAB.dems =AB.FDR.DEMs$PeakID)

venn.list.iDEMs <- list(facA.dems = A.FDR.DEMs$PeakID, 
                        facB.dems = B.FDR.DEMs$PeakID, 
                        facAB.dems =AB.FDR.DEMs$PeakID,
                        facAB.idems = AB.iDEMs$PeakID)

# Plotting Venn Diagram with Euler 
venn.dems <- venn.diagram(venn.list, filename=NULL)
grid.draw(venn.dems)

plot(euler(venn.list), 
     quantities = TRUE, 
     edges = c("green", "red", "blue"), 
     fills = FALSE, 
     labels= FALSE, 
     legend = TRUE, 
     main = "HML AB090342 Metab Euler-5B1.pdf")

venn.idems <- venn.diagram(venn.list.iDEMs, filename=NULL)
grid.draw(venn.idems)

plot(euler(venn.list.iDEMs), 
     quantities = TRUE, 
     edges = c("green", "red", "blue", "black"), 
     fills = FALSE, 
     labels= FALSE, 
     legend = TRUE, 
     main = "HML AB090342 Metab Euler-5B2.pdf")

# PCA of Expression Profiles ----
# Set colors
col.pal <- c("Concordant.Down"= "#a6cee3",
             "Concordant.Up"= "#1f78b4",
             "Discordant.Down"= "#b2df8a",
             "Discordant.Up"= "#33a02c",
             "facA-Dominant.Down"= "#fb9a99",
             "facA-Dominant.Up"= "#e31a1c",
             "facB-Dominant.Down"= "#fdbf6f", 
             "facB-Dominant.Up"= "#ff7f00")

FC.ExPr <- data.frame(FC.facA = logratio2foldchange(d.log2FC$d.facA, base = 2),
                      FC.facB = logratio2foldchange(d.log2FC$d.facB, base = 2),
                      FC.facC = logratio2foldchange(d.log2FC$d.facAB, base = 2))

ExPr.AB.DEMs <- cbind.data.frame(FC.ExPr,
                                 ExpProf)

# PCA
data.class <- ExPr.AB.DEMs$UDDualTx
data.pca <- prcomp(FC.ExPr, center = TRUE, scale. = TRUE)
g <- ggbiplot(data.pca, 
              obs.scale = 1, 
              var.scale = 1, 
              groups = data.class, 
              ellipse = FALSE, 
              circle = FALSE)
g
fin.g <- g + 
  theme_bw() + 
  xlim(-2.5, 2.5) + 
  ylim(-2.5, 2.5) + 
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  scale_color_manual(name = "Expression Profiles", values = col.pal) +
  scale_fill_manual(values = col.pal)
fin.g

# Summary Interaction Score ----
AB.iDEMs$SummIS <- log(AB.iDEMs$AbsIScore, base = 2)
ISplot <- AB.iDEMs[order(AB.iDEMs$SummIS),]
barplot(ISplot$SummIS, horiz = TRUE)

pdf("JF HCC366 RPPA.pdf")
barplot(ISplot$SummIS, horiz = TRUE)
dev.off()



# Quick R object diagnostics ----
# Use these commands to verify the correct charactericts of any object in you Enviroment if any part of the code is not working 
class(sim.l2qn.expr.val)
typeof(sim.l2qn.expr.val)
attributes(expr.val)
