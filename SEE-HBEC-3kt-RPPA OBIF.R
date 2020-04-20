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
library(circlize)
library(plyr)
library(UpSetR)
library(ape)

# Load and organize your analysis-ready Omics dataset ----
# Manually import your dataset into your R Environment from any source (i.e., GEO, Excel, CSV, etc.).
dataset <- read_excel("~/Documents/R/SEE HBEC-3kt RPPA/Statistical DEGs Analysis/SEE-HBEC3kt-RPPA.xls")

# Rename your dataset and define factors
data.orig <- dataset
factor.A <- "Pam2CSK4"
factor.B <- "ODN2395"
factor.AB <- "Pam2ODN"

# Identify the column with the name of the features
feat.ids <- data.orig[,1] # Select the column numbers 
colnames(feat.ids) <- "Feat.ID"

# Identify the columns with the expression values per sample
expr.val <- data.orig[,3:18] # Select the column numbers

# Identify the columns with additional information about your dataset
add.col <- data.orig[,2]

# Preserve feature names in all subsets of your datasets
rownames(data.orig) <- feat.ids$Feat.ID
rownames(feat.ids) <- feat.ids$Feat.ID
rownames(expr.val) <- feat.ids$Feat.ID
rownames(add.col) <- feat.ids$Feat.ID

# Indicate number of replicates per group
n.ctrl <- 4
n.facA <- 4
n.facB <- 4
n.facAB <- 4

# Set color pallete for samples
col.samples <- c(rep("#000000", n.ctrl), rep("#8dd3c7", n.facA), rep("#fb8072", n.facB), rep("#80b1d3", n.facAB))

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
sim.l2.expr.val <- expr.val # log2(expr.val) # Avoided because data is already in log2
sim.l2qn.expr.val <- normalize.quantiles(data.matrix(sim.l2.expr.val)) # Avoided because this data was already normalize and this step did not improve statistics
sim.l2qn.expr.val <- as.data.frame(sim.l2qn.expr.val)
colnames(sim.l2qn.expr.val) <- cn.expr.val

### Array-based dataset transformations (background correction, log2 transformation and quantile normalization if needed)
arr.eSet <- new("ExpressionSet", exprs=as.matrix(expr.val))
arr.bg.eSet <- arr.eSet # lumiB(arr.eSet, method = "bgAdjust.affy") # Avoided because this data was already pre-processed, and background sustraction greatly reduced the values of the residuals while also reducing the adj.R2
arr.bgl2.eSet <- arr.bg.eSet  # lumiT(arr.bg.eSet, "log2") # Avoided because this data was ready for analysis
arr.bgl2qn.eSet <- lumiN(arr.bgl2.eSet, method = "quantile") #Avoided because this data was ready for analysis
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

qc.vp.orig
qc.vp.sim
qc.vp.arr

# 1.B Using a density plot
qc.dp.orig <- ggplot(vp.expr.val.orig, aes(value, color = variable)) + geom_density() + scale_color_manual(values = col.samples) + ggtitle(label="Original Dataset", subtitle = deparse(substitute(vp.expr.val.orig))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top") 
qc.dp.sim <- ggplot(vp.expr.val.sim, aes(value, color = variable)) + geom_density() + scale_color_manual(values = col.samples) + ggtitle(label="Simple Transformations", subtitle = deparse(substitute(vp.expr.val.sim))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top") 
qc.dp.arr <- ggplot(vp.expr.val.arr, aes(value, color = variable)) + geom_density() + scale_color_manual(values = col.samples) + ggtitle(label="Array-based Transformations", subtitle = deparse(substitute(vp.expr.val.arr))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top") 

qc.dp.orig
qc.dp.sim 
qc.dp.arr

# 1.C Using a box plot 
qc.bxpt.orig <- boxplot(expr.val, boxwex=0.5, notch=TRUE, outline=FALSE, las =2, border = col.samples, main = "Original Dataset")

qc.bxpt.sim <- boxplot(sim.l2qn.expr.val, boxwex=0.5, notch=TRUE, outline=FALSE, las =2, border = col.samples, main = "Simple Transformations")

qc.bxpt.arr <- boxplot(arr.bgl2qn.expr.val, boxwex=0.5, notch=TRUE, outline=FALSE, las =2, border = col.samples, main = "Array-based Transformations")

## 1.D Using hierarchical clutersing
hcl.orig <- hclust(distanceMatrix(
  dataset =  expr.val, 
  metric =  "pearson"), 
  method = "ward.D")
hcl.sim <- hclust(distanceMatrix(
  dataset =  sim.l2qn.expr.val, 
  metric =  "pearson"), 
  method = "ward.D2")
hcl.arr <- hclust(distanceMatrix(
  dataset =  arr.bgl2qn.expr.val, 
  metric =  "pearson"), 
  method = "ward.D")

qc.hcl.orig <- myplclust(hcl.orig, labels = hcl.orig$labels, lab.col = col.samples) 

qc.hcl.sim <- myplclust(hcl.sim, labels = hcl.orig$labels, lab.col = col.samples)

qc.hcl.arr <- myplclust(hcl.arr, labels = hcl.orig$labels, lab.col = col.samples)


# Other dendrograms ----
dend <-hcl.arr %>%
  as.dendrogram %>%
  set("leaves_pch", 19) %>% 
  set("leaves_col",col.samples) %>% 
  set("hang_leaves", -1) 

plot(dend)
ggd1 <- as.ggdend(dend)
ggplot(ggd1) 
ggplot(ggd1, horiz = TRUE, theme = NULL) 
ggplot(ggd1, theme = theme_minimal()) 
ggplot(ggd1, labels = T) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")

plot(as.phylo(dend), type = "fan")


  
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
par(mfrow = c(1,2))
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

#grid.arrange(
 # qc.diagp.arr.all$ResvsFit,
  #qc.diagp.arr.all$NormQQ,
  #qc.diagp.arr.all$ScLoc,
  #qc.diagp.arr.all$CookD,
  #qc.diagp.arr.all$ResvsLev,
  #qc.diagp.arr.all$CookvsLev
#)


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
fit2.ebayes <- eBayes(fit2)
adj.Pvalues.contrBH <- apply(fit2.ebayes$p.value, 2, p.adjust, method="BH")
adj.Pvalues.contrBonf <- apply(fit2.ebayes$p.value, 2, p.adjust, method="bonferroni")

# Calculate the statistics for expression and contrast analysis using Two-way ANOVA
fit.2x2 <- lmFit(ints.mtx.arr,design2x2)
fit.2x2.ebayes <- eBayes(fit.2x2)
adj.2x2.Pvalues.lmBH <- apply(fit.2x2.ebayes$p.value, 2, p.adjust, method="BH")
adj.2x2.Pvalues.lmBonf <- apply(fit.2x2.ebayes$p.value, 2, p.adjust, method="bonferroni")

fit2.2x2 <- contrasts.fit(fit.2x2,cont.mtx.des2x2)
fit2.2x2.ebayes <- eBayes(fit2.2x2)
adj.2x2.Pvalues.contrBH <- apply(fit2.2x2.ebayes$p.value, 2, p.adjust, method="BH")
adj.2x2.Pvalues.contrBonf <- apply(fit2.2x2.ebayes$p.value, 2, p.adjust, method="bonferroni")

# Make sure that all p-values are calculated equally between One-way and Two-way ANOVAs 
# R will compare p-values to an undefined number of significant values, hence we round below to 5 significant digits
comp.pvals <- cbind.data.frame(a = round(fit.ebayes$p.value, digits = 5), 
                               b = round(fit.2x2.ebayes$p.value, digits = 5), 
                               c = round(fit2.ebayes$p.value, digits = 5), 
                               d = round(fit2.2x2.ebayes$p.value, digits = 5))
# Cross comparison between the p-values of facAB are the most useful to use to see difference between one-way and two-way ANOVA calculations, 
# because each one calculates the value differently: One-Way from direct Main Effects in fit.ebayes, and Two-way fron Contrasting in fit2.2x2.ebayes
comp.sum <- ifelse (comp.pvals$a.facAB == comp.pvals$d.facAB, 0, 1)
# The summ of the substraction between the p-values from both methods should be 0
sum(comp.sum)




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
ExpProf <- data.frame(UD.A = ifelse (d.log2FC[,1] < 0 , "Down", ifelse ( d.log2FC[,1] > 0, "Up", "Unchanged")),
                      UD.B = ifelse ( d.log2FC[,2] < 0 , "Down", ifelse ( d.log2FC[,2] > 0, "Up", "Unchanged")),
                      UD.AB = ifelse ( d.log2FC[,3] < 0 , "Down", ifelse ( d.log2FC[,3] > 0, "Up", "Unchanged")))
# Verify 0 or NA values
sum(ExpProf$UD.A == "Unchanged")
sum(ExpProf$UD.B == "Unchanged")
sum(ExpProf$UD.AB == "Unchanged")

# Assign & Verify SingleTx group
ExpProf$SingleTx <- ifelse(ExpProf[,1] == "Up" & ExpProf[,2] == "Up", "Cooperative", 
                           ifelse(ExpProf[,1] == "Down" & ExpProf[,2] == "Down", "Cooperative",
                                  ifelse(ExpProf[,1] == "Unchanged" & ExpProf[,2] == "Unchanged", "Cooperative",
                                         ifelse(ExpProf[,1] == "Up" & ExpProf[,2] == "Down", "Competitive",
                                                ifelse(ExpProf[,1] == "Down" & ExpProf[,2] == "Up", "Competitive", 
                                                       ifelse(ExpProf[,1] == "Up" & ExpProf[,2] == "Unchanged", "Competitive",
                                                              ifelse(ExpProf[,1] == "Down" & ExpProf[,2] == "Unchanged", "Competitive",
                                                                     ifelse(ExpProf[,1] == "Unchanged" & ExpProf[,2] == "Up", "Competitive",
                                                                            ifelse(ExpProf[,1] == "Unchanged" & ExpProf[,2] == "Down", "Competitive",
                                                                                   "Error2")))))))))
sum(ExpProf$SingleTx == "Error2")

# Assign & Verify DualTx group
ExpProf$DualTx <- ifelse(ExpProf$SingleTx == "Cooperative" & ExpProf[,1] == "Up" & ExpProf[,3] == "Up", "Concordant",
                         ifelse(ExpProf$SingleTx == "Cooperative" & ExpProf[,1] == "Down" & ExpProf[,3] == "Down", "Concordant",
                                ifelse(ExpProf$SingleTx == "Cooperative" & ExpProf[,1] == "Up" & ExpProf[,3] == "Down", "Discordant",
                                       ifelse(ExpProf$SingleTx == "Cooperative" & ExpProf[,1] == "Down" & ExpProf[,3] == "Up", "Discordant",
                                              ifelse(ExpProf$SingleTx == "Competitive" & ExpProf[,1] == "Up" & ExpProf[,3] == "Up", "facA-Dominant",
                                                     ifelse(ExpProf$SingleTx == "Competitive" & ExpProf[,1] == "Down" & ExpProf[,3] == "Down", "facA-Dominant",
                                                            ifelse(ExpProf$SingleTx == "Competitive" & ExpProf[,2] == "Up" & ExpProf[,3] == "Up", "facB-Dominant",
                                                                   ifelse(ExpProf$SingleTx == "Competitive" & ExpProf[,2] == "Down" & ExpProf[,3] == "Down", "facB-Dominant",
                                                                          ifelse(ExpProf$SingleTx == "Cooperative" & ExpProf[,3] == "Unchanged", "No.AB.DEM",
                                                                                 ifelse(ExpProf$SingleTx == "Competitive" & ExpProf[,3] == "Unchanged", "No.AB.DEM",
                                                                                        "Error3"))))))))))
sum(ExpProf$DualTx == "No.AB.DEM")
sum(ExpProf$DualTx == "Error3")

# Assign SingleTx & Dual group related to expression in dual exposure
ExpProf$UDSingleTx <- paste(ExpProf$SingleTx, ExpProf$UD.AB, sep = ".")
ExpProf$UDDualTx <- paste(ExpProf$DualTx, ExpProf$UD.AB, sep = ".")

# If too many "Unchanged" features, "No.AB.DEMs", or "Error#" -> Re-consider dataset transformations and undetectable values ----
sum(ExpProf$UD.A == "Unchanged")
sum(ExpProf$UD.B == "Unchanged")
sum(ExpProf$UD.AB == "Unchanged")
sum(ExpProf$DualTx == "No.AB.DEM")
sum(ExpProf$DualTx == "Error3")

# Final data compilation ----
final.res <- data.frame(feat.ids,
                        add.col,
                        arr.bgl2qn.expr.val,
                        m.l2.vals, 
                        d.log2FC,
                        pv=fit.ebayes$p.value[,2:4],
                        pv=fit2.ebayes$p.value,
                        adjpvBH=adj.Pvalues.lmBH[,2:4],
                        adjpvBH=adj.Pvalues.contrBH, 
                        adjpvBonf=adj.Pvalues.lmBonf[,2:4],
                        adjpvBonf=adj.Pvalues.contrBonf, 
                        IntScore,
                        ExpProf)

# Start from final.res file after dataset analysis ----
final.labels <- final.res

# Define cut off values for selection of DEMs
l2fc.co <- 0.3
FDR.co <- 0.01
eff.co <-0.05

# Label your features with by DEM per group
final.labels$A.DEM <- ifelse (abs(final.labels$d.facA) >= l2fc.co & final.labels$adjpvBH.facA <= FDR.co,
                              "Yes", "No")
final.labels$B.DEM <- ifelse (abs(final.labels$d.facB) >= l2fc.co & final.labels$adjpvBH.facB <= FDR.co,
                              "Yes", "No")
final.labels$AB.DEM <- ifelse (abs(final.labels$d.facAB) >= l2fc.co & final.labels$adjpvBH.facAB <= FDR.co,
                               "Yes", "No")
# Label your AB.DEM by factorial effects
final.labels$AB.SME.A <- ifelse (final.labels$AB.DEM == "Yes" & final.labels$pv.SME.facA <= eff.co,
                                 "Yes", "No")
final.labels$AB.SME.B <- ifelse (final.labels$AB.DEM == "Yes" & final.labels$pv.SME.facB <= eff.co,
                                 "Yes", "No")
final.labels$AB.iDEM <- ifelse (final.labels$AB.DEM == "Yes" & final.labels$pv.Inter.fAxfB <= eff.co,
                                "Yes", "No")

# Label your AB.DEM by expression profiles
final.labels$EPs <- ifelse (final.labels$UDDualTx == "Concordant.Up", "I",
                            ifelse (final.labels$UDDualTx == "Concordant.Down", "II",
                                    ifelse (final.labels$UDDualTx == "Discordant.Up", "III",
                                            ifelse (final.labels$UDDualTx == "Discordant.Down", "IV",
                                                    ifelse (final.labels$UDDualTx == "facA-Dominant.Up", "V",
                                                            ifelse (final.labels$UDDualTx == "facA-Dominant.Down", "VI",
                                                                    ifelse (final.labels$UDDualTx == "facB-Dominant.Up", "VII",
                                                                            ifelse (final.labels$UDDualTx == "facB-Dominant.Down", "VIII",
                                                                                    "Error.No.EP"))))))))
final.labels$f.IS <- log2(final.labels$AbsIScore)
final.labels$f.IS.iDEMs <- ifelse (final.labels$AB.iDEM == "Yes", final.labels$f.IS, log2(1))

# Assign a position to each DEM per group to be used as x coordinante for Circos plot in circlize
final.labels <- ddply(final.labels,"A.DEM", transform, x1.A = sample(0:(length(A.DEM)-1), length(A.DEM), replace = FALSE))
final.labels$x2.A <- final.labels$x1.A+1
final.labels <- ddply(final.labels,"B.DEM", transform, x1.B = sample(0:(length(B.DEM)-1), length(B.DEM), replace = FALSE))
final.labels$x2.B <- final.labels$x1.B+1
final.labels <- ddply(final.labels,"AB.DEM", transform, x1.AB = sample(0:(length(AB.DEM)-1), length(AB.DEM), replace = FALSE))
final.labels$x2.AB <- final.labels$x1.AB+1

final.labels$AB.EP <- ifelse (final.labels$AB.DEM == "Yes", final.labels$UDDualTx, 
                              ifelse (final.labels$A.DEM == "Yes", "A.DEMs",
                                      ifelse (final.labels$B.DEM == "Yes", "B.DEMs","No.DEMs")))

final.labels <- ddply(final.labels,"AB.EP", transform, x1.EP = sample(0:(length(AB.EP)-1), length(AB.EP), replace = FALSE))
final.labels$x2.EP <- final.labels$x1.EP+1

final.labels$c.l2fc <- ifelse (final.labels$AB.DEM == "Yes", final.labels$d.facAB, 
                               ifelse (final.labels$A.DEM == "Yes", final.labels$d.facA,
                                       ifelse (final.labels$B.DEM == "Yes", final.labels$d.facB,"0")))

# Assign their positions ordered by log2FC within each group
final.labels <- ddply(final.labels,.(AB.EP), mutate, x1.ord = order(c.l2fc))
final.labels$x2.ord <- final.labels$x1.ord + 1

# Subset the DEMs per group 
A.DEMS <- subset(final.labels, A.DEM == "Yes")
B.DEMS <- subset(final.labels, B.DEM == "Yes")
AB.DEMS <- subset(final.labels, AB.DEM == "Yes")

EP.1 <- subset(AB.DEMS, UDDualTx == "Concordant.Up")
EP.2 <- subset(AB.DEMS, UDDualTx == "Concordant.Down")
EP.3 <- subset(AB.DEMS, UDDualTx == "Discordant.Up")
EP.4 <- subset(AB.DEMS, UDDualTx == "Discordant.Down")
EP.5 <- subset(AB.DEMS, UDDualTx == "facA-Dominant.Up")
EP.6 <- subset(AB.DEMS, UDDualTx == "facA-Dominant.Down")
EP.7 <- subset(AB.DEMS, UDDualTx == "facB-Dominant.Up")
EP.8 <- subset(AB.DEMS, UDDualTx == "facB-Dominant.Down")


No.EP <- subset(AB.DEMS, EPs == "Error.No.EP")
nrow(No.EP) 
# If = 0, continue as usual. Errors in Up and Down labeling does not affect any of the DEMs.
# If > 1, Detect and remove DEMs that could not be classified as EP, or the features that are not DEMs, and consider whether its valuable to maintain or not.
# Verify the No.EP with View(NoEP)

# Assign their unique factors for final dataframe to use in circos
c.A.DEMS <- data.frame(c.factors = rep("A.DEMs", nrow(A.DEMS)))
c.B.DEMS <-data.frame(c.factors = rep("B.DEMs", nrow(B.DEMS)))
c.AB.DEMS <- data.frame(c.factors = rep("AB.DEMs", nrow(AB.DEMS)))

c.EP.1 <- data.frame(c.factors = rep("I", nrow(EP.1)))
c.EP.2 <- data.frame(c.factors = rep("II", nrow(EP.2)))
c.EP.3 <- data.frame(c.factors = rep("III", nrow(EP.3))) 
c.EP.4 <- data.frame(c.factors = rep("IV", nrow(EP.4))) 
c.EP.5 <- data.frame(c.factors = rep("V", nrow(EP.5))) 
c.EP.6 <- data.frame(c.factors = rep("VI", nrow(EP.6))) 
c.EP.7 <- data.frame(c.factors = rep("VII", nrow(EP.7))) 
c.EP.8 <- data.frame(c.factors = rep("VIII", nrow(EP.8))) 

# Assign their unique positions for random order to use in circos
c.A.DEMS$c.x1 <- A.DEMS$x1.A
c.A.DEMS$c.x2 <- A.DEMS$x2.A
c.B.DEMS$c.x1 <- B.DEMS$x1.B
c.B.DEMS$c.x2 <- B.DEMS$x2.B
c.AB.DEMS$c.x1 <- AB.DEMS$x1.AB
c.AB.DEMS$c.x2 <- AB.DEMS$x2.AB

c.EP.1$c.x1 <- EP.1$x1.EP
c.EP.1$c.x2 <- EP.1$x2.EP
c.EP.2$c.x1 <- EP.2$x1.EP
c.EP.2$c.x2 <- EP.2$x2.EP
c.EP.3$c.x1 <- EP.3$x1.EP
c.EP.3$c.x2 <- EP.3$x2.EP
c.EP.4$c.x1 <- EP.4$x1.EP
c.EP.4$c.x2 <- EP.4$x2.EP
c.EP.5$c.x1 <- EP.5$x1.EP
c.EP.5$c.x2 <- EP.5$x2.EP
c.EP.6$c.x1 <- EP.6$x1.EP
c.EP.6$c.x2 <- EP.6$x2.EP
c.EP.7$c.x1 <- EP.7$x1.EP
c.EP.7$c.x2 <- EP.7$x2.EP
c.EP.8$c.x1 <- EP.8$x1.EP
c.EP.8$c.x2 <- EP.8$x2.EP

 # Aqui ----
# Assign their unique log2FC for final dataframe to use in circos
c.A.DEMS$c.lfc <- A.DEMS$d.facA
c.B.DEMS$c.lfc <- B.DEMS$d.facB
c.AB.DEMS$c.lfc <- AB.DEMS$d.facAB

c.EP.1$c.lfc <- EP.1$d.facAB
c.EP.2$c.lfc <- EP.2$d.facAB
c.EP.3$c.lfc <- EP.3$d.facAB
c.EP.4$c.lfc <- EP.4$d.facAB
c.EP.5$c.lfc <- EP.5$d.facAB
c.EP.6$c.lfc <- EP.6$d.facAB
c.EP.7$c.lfc <- EP.7$d.facAB
c.EP.8$c.lfc <- EP.8$d.facAB

# Integrate the rest of dataset with the factors and positions
c.A.DEMS <- cbind.data.frame(c.A.DEMS, A.DEMS)
c.B.DEMS <- cbind.data.frame(c.B.DEMS, B.DEMS)
c.AB.DEMS <- cbind.data.frame(c.AB.DEMS, AB.DEMS)

c.EP.1 <- cbind.data.frame(c.EP.1, EP.1)
c.EP.2 <- cbind.data.frame(c.EP.2, EP.2)
c.EP.3 <- cbind.data.frame(c.EP.3, EP.3)
c.EP.4 <- cbind.data.frame(c.EP.4, EP.4)
c.EP.5 <- cbind.data.frame(c.EP.5, EP.5)
c.EP.6 <- cbind.data.frame(c.EP.6, EP.6)
c.EP.7 <- cbind.data.frame(c.EP.7, EP.7)
c.EP.8 <- cbind.data.frame(c.EP.8, EP.8)

c.AB.iDEMs <- subset(c.AB.DEMS, AB.iDEM == "Yes")

# Combine DEMs dataframes to a final dataframe to be used in circos
c.all.DEMs.df <- rbind.data.frame(c.A.DEMS, c.B.DEMS, c.AB.DEMS)
c.all.EP.df <- rbind.data.frame(c.A.DEMS, c.B.DEMS, c.EP.1,c.EP.2,c.EP.3,c.EP.4,c.EP.5,c.EP.6,c.EP.7,c.EP.8)

#Subset the shared DEMs between groups
links.AnAB <- subset.data.frame(x = c.AB.DEMS,
                                subset = c.AB.DEMS$A.DEM == "Yes" &
                                  c.AB.DEMS$AB.DEM == "Yes")
links.BnAB <- subset.data.frame(x = c.AB.DEMS,
                                subset = c.AB.DEMS$B.DEM == "Yes" &
                                  c.AB.DEMS$AB.DEM == "Yes")
links.AnB <- subset.data.frame(x = c.A.DEMS,
                               subset = c.A.DEMS$A.DEM == "Yes" &
                                 c.A.DEMS$B.DEM == "Yes")

# Create beds of linked DEMs between groups
bed1.AnAB <- data.frame(c.factors = rep("A.DEMs", nrow(links.AnAB)), 
                        start = links.AnAB$x1.A, 
                        end = links.AnAB$x1.A)
bed2.AnAB <- data.frame(c.factors = rep("AB.DEMs", nrow(links.AnAB)), 
                        start = links.AnAB$x1.AB, 
                        end = links.AnAB$x1.AB)

bed1.BnAB <- data.frame(c.factors = rep("B.DEMs", nrow(links.BnAB)), 
                        start = links.BnAB$x1.B, 
                        end = links.BnAB$x1.B)
bed2.BnAB <- data.frame(c.factors = rep("AB.DEMs", nrow(links.BnAB)), 
                        start = links.BnAB$x1.AB, 
                        end = links.BnAB$x1.AB)

bed1.AnB <- data.frame(c.factors = rep("A.DEMs", nrow(links.AnB)), 
                       start = links.AnB$x1.A, 
                       end = links.AnB$x1.A)
bed2.AnB <- data.frame(c.factors = rep("B.DEMs", nrow(links.AnB)), 
                       start = links.AnB$x1.B, 
                       end = links.AnB$x1.B)

# Create beds of linked DEMs between groups and Expression profiles
ep.bed1.AnAB <- data.frame(c.factors = rep("A.DEMs", nrow(links.AnAB)), 
                           start = links.AnAB$x1.A, 
                           end = links.AnAB$x1.A)
ep.bed2.AnAB <- data.frame(c.factors = links.AnAB$EPs, 
                           start = links.AnAB$x1.EP, 
                           end = links.AnAB$x2.EP)

ep.bed1.BnAB <- data.frame(c.factors = rep("B.DEMs", nrow(links.BnAB)), 
                           start = links.BnAB$x1.B, 
                           end = links.BnAB$x1.B)
ep.bed2.BnAB <- data.frame(c.factors = links.BnAB$EPs, 
                           start = links.BnAB$x1.EP, 
                           end = links.BnAB$x2.EP)

ep.bed1.AnB <- data.frame(c.factors = rep("A.DEMs", nrow(links.AnB)), 
                          start = links.AnB$x1.A, 
                          end = links.AnB$x1.A)
ep.bed2.AnB <- data.frame(c.factors = rep("B.DEMs", nrow(links.AnB)), 
                          start = links.AnB$x1.B, 
                          end = links.AnB$x1.B)

# Make the region columns for circos.genomic functions ----
# For random feature order keep c.x1 and c.x2 as $c.x1 and $c.x2, for order by log2 FC change them to $x1.ord and $x2.ord
c.start <- cbind.data.frame(c.factors = c.all.EP.df$c.factors, 
                            c.x1 = c.all.EP.df$c.x1, 
                            c.x2 = c.all.EP.df$c.x2)

# Initialize the circos plot sectors and outer most track
circos.clear()
circos.par(start.degree = 90,
           clock.wise = FALSE,
           gap.degree = 2)
circos.genomicInitialize(data = c.all.EP.df,
                         plotType = NULL)

# Antagonistic and synergistic iDEMs and their Interaction Score ----
# Make the genomic format
c.IS.iDEMs <- cbind.data.frame(c.start, 
                               c.value = ifelse (
                                 c.start$c.factors == "A.DEMs" | 
                                   c.start$c.factors == "B.DEMs",
                                 0, c.all.EP.df$f.IS.iDEMs))

cond_col_is = function (value) {ifelse (value < 0, "orange", 
                                        ifelse (value > 0, "purple",
                                                ifelse (CELL_META$sector.index == "A.DEMs" | CELL_META$sector.index == "B.DEMs", "white", 
                                                        "lightgrey")))}

circos.par(track.margin = c(0.0025,0),
           cell.padding = c(0.0,0,0.0,0))
circos.genomicTrackPlotRegion(data = c.IS.iDEMs,
                              bg.col = NA,
                              bg.border = NA,
                              track.height = 0.25,
                              ylim = range(c.IS.iDEMs$c.value),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, 
                                                    value,
                                                    type = "h",
                                                    baseline = 0,
                                                    lwd = 2,
                                                    col = cond_col_is(value),
                                                    ...)
                              })

c.rect.smeA <- cbind.data.frame(c.start, 
                                c.value = ifelse (c.start$c.factors == "A.DEMs" | c.start$c.factors == "B.DEMs", 0,
                                                  ifelse (c.all.EP.df$AB.SME.A == "Yes", 1, 4)))
c.rect.smeB <- cbind.data.frame(c.start, 
                                c.value = ifelse (c.start$c.factors == "A.DEMs" | c.start$c.factors == "B.DEMs", 0,
                                                  ifelse (c.all.EP.df$AB.SME.B == "Yes", 2, 4)))
c.rect.fAxfB <- cbind.data.frame(c.start, 
                                 c.value = ifelse (c.start$c.factors == "A.DEMs" | c.start$c.factors == "B.DEMs", 0,
                                                   ifelse(c.all.EP.df$AB.iDEM == "Yes", 3, 4)))

c.rect.effs <- cbind.data.frame(c.start,
                                c.v1 = c.rect.smeA$c.value,
                                c.v2 = c.rect.smeB$c.value,
                                c.v3 = c.rect.fAxfB$c.value
)

cond_col_eff = function (value) {ifelse (value == 1, "#1b9e77",
                                         ifelse (value == 2, "#d95f02",
                                                 ifelse (value == 3, "#e7298a", 
                                                         ifelse (value == 4, "#f0f0f0","white"))))}

circos.par(track.margin = c(0, 0.0025),
           cell.padding = c(0,0,0,0))
circos.genomicTrackPlotRegion(data = c.rect.effs,
                              bg.col = "#f0f0f0",
                              stack = TRUE,
                              bg.border = NA,
                              track.height = 0.1,
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region,
                                                   value,
                                                   border = cond_col_eff(value),
                                                   col = cond_col_eff(value),
                                                   ...)})

cond_col_ep = function (value) {ifelse (CELL_META$sector.index == "I", "#1f78b4",
                                        ifelse (CELL_META$sector.index == "II", "#a6cee3",
                                                ifelse (CELL_META$sector.index == "III", "#e31a1c",
                                                        ifelse (CELL_META$sector.index == "IV", "#fb9a99",
                                                                ifelse (CELL_META$sector.index == "V", "#33a02c",
                                                                        ifelse (CELL_META$sector.index == "VI", "#b2df8a",
                                                                                ifelse (CELL_META$sector.index == "VII", "#ff7f00",
                                                                                        ifelse (CELL_META$sector.index == "VIII","#fdbf6f","white"))))))))}
circos.par(track.margin = c(0.0025,0.0025),
           cell.padding = c(0.000,0,0.000,0))

circos.genomicTrackPlotRegion(data = c.all.EP.df,
                              bg.col = NA,
                              bg.border = NA,
                              track.height = 0.05,
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region,
                                                   value,
                                                   border = cond_col_ep(value),
                                                   col = cond_col_ep(Cvalue),
                                                   ...)})
c.l.l2fc <- cbind.data.frame(c.start,
                             c.v1 = c.all.EP.df$c.lfc)
c.l.up <- cbind.data.frame(c.start,
                           c.v1 = ifelse (c.all.EP.df$c.lfc > 0, c.all.EP.df$c.lfc, 0))
c.l.down <- cbind.data.frame(c.start,
                             c.v1 = ifelse (c.all.EP.df$c.lfc < 0, c.all.EP.df$c.lfc, 0))

cond_col_l2fc = function (value) {ifelse (value < 0, "blue",
                                          ifelse (value > 0, "red", "white"))}

circos.par(track.margin = c(0.0025,0.0025),
           cell.padding = c(0.000,0,0.000,0))

circos.genomicTrackPlotRegion(data = c.l.l2fc,
                              bg.col = NA,
                              bg.border = "lightgrey",
                              track.height = 0.20,
                              ylim = range(c.l.l2fc$c.v1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, 
                                                    value,
                                                    type = "h",
                                                    baseline = 0,
                                                    col = cond_col_l2fc(value),
                                                    ...)
                              })

cond_col_dems = function (value) {ifelse (CELL_META$sector.index == "A.DEMs", "#8dd3c7",
                                          ifelse (CELL_META$sector.index == "B.DEMs", "#fb8072",
                                                  "#80b1d3"))}

circos.par(track.margin = c(0.0025,0),
           cell.padding = c(0,0,0,0))

circos.genomicTrackPlotRegion(data = c.l.down,
                              bg.col = NA,
                              bg.border = "lightgrey",
                              track.height = 0.05,
                              ylim = range(c.l.down$c.v1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicRect(region,
                                                   value,
                                                   border = cond_col_dems(value),
                                                   col = cond_col_dems(Cvalue),
                                                   ...)})


circos.par(track.margin = c(0.0025,0.0025),
           cell.padding = c(0.000,0,0.000,0))

circos.genomicLink(region1 = ep.bed1.AnB,
                   region2 = ep.bed2.AnB,
                   border = NA,
                   col = "#d9d9d9",
                   rou1 = circlize:::get_most_inside_radius(),
                   rou2 = circlize:::get_most_inside_radius(),
                   h = 0.5)

circos.genomicLink(region1 = ep.bed1.BnAB,
                   region2 = ep.bed2.BnAB,
                   border = NA,
                   col = "#fb8072",
                   rou1 = circlize:::get_most_inside_radius()*1,
                   rou2 = circlize:::get_most_inside_radius()*0.95)

circos.genomicLink(region1 = ep.bed1.AnAB,
                   region2 = ep.bed2.AnAB,
                   border = NA,
                   col = "#8dd3c7",
                   rou1 = circlize:::get_most_inside_radius()*1,
                   rou2 = circlize:::get_most_inside_radius()*1)

draw.sector(get.cell.meta.data("cell.start.degree", sector.index = "VII") + 1,
            get.cell.meta.data("cell.end.degree", sector.index = "I") - 1,
            rou1 = get.cell.meta.data("cell.top.radius", track.index = 1) + 0.01,
            rou2 = get.cell.meta.data("cell.bottom.radius", track.index = 5) - 0.05,
            col = "NA",
            border = "blue",
            lwd = 1)

# Export as TIFF manually with dimensions 1200 x 1200

# Create list of vectors for DEMs per group for Venn Diagram
venn.list <- list(facA.dems = c.A.DEMS$Feat.ID, 
                  facB.dems = c.B.DEMS$Feat.ID, 
                  facAB.dems =c.AB.DEMS$Feat.ID)

venn.list.iDEMs <- list(facA.dems = c.A.DEMS$Feat.ID, 
                        facB.dems = c.B.DEMS$Feat.ID, 
                        facAB.dems =c.AB.DEMS$Feat.ID,
                        facAB.idems = c.AB.iDEMs$Feat.ID)

# Plotting Venn Diagram with Euler 
venn.dems <- venn.diagram(venn.list, filename=NULL)
grid.draw(venn.dems)

plot(euler(venn.list), 
     quantities = TRUE, 
     edges = c("#8dd3c7","#fb8072","#80b1d3"), 
     fills = FALSE, 
     labels= TRUE)

# Export as PDF manually with dimensions 10 x 10 inches

plot(euler(venn.list.iDEMs), 
     quantities = TRUE, 
     edges = c("#8dd3c7","#fb8072","#80b1d3","#e7298a"), 
     fills = FALSE, 
     labels= FALSE, 
     legend = TRUE)

venn.idems <- venn.diagram(venn.list.iDEMs, filename=NULL)
grid.draw(venn.idems)

# UpSet plot for ALL regions 
intersections <- venn(venn.list)
intersections <- data.frame(n = intersections$original.values)
max.set.val <- sum(intersections)

upset(data = fromList(venn.list),
      matrix.color = "blue",
      main.bar.color = "black",
      sets.bar.color = c("#fb8072","#80b1d3", "#8dd3c7"),
      order.by = "freq",
      decreasing = TRUE,
      empty.intersections = TRUE,
      scale.intersections = "identity",
      scale.sets = "identity",
      set_size.show = TRUE,
      set_size.scale_max = max.set.val)

# UpSet plot for percentages relative to AB.DEMs with empty regions 
AnAB <- subset(c.A.DEMS, AB.DEM == "Yes")
BnAB <- subset(c.B.DEMS, AB.DEM == "Yes")

# Removey any group from comparisons if the intersections (AnAB or BnAB) have 0 elements
venn.ABlist <- list(facA.dems = AnAB$Feat.ID, 
  facB.dems = BnAB$Feat.ID, 
  facAB.dems = c.AB.DEMS$Feat.ID)

int.ABs <- venn(venn.ABlist)
int.per <- int.ABs$original.values
per.AB.input <- round(int.ABs$original.values/sum(int.ABs$original.values)*100, digits = 0)
sum(per.AB.input)

# Match the same number of group colors in sets.bar.col according to the shown intersections. 
upset(fromExpression(per.AB.input),
      matrix.color = "blue",
      main.bar.color = "#80b1d3",
      sets.bar.color = c("#80b1d3", "#fb8072", "#8dd3c7"),
      order.by = "freq",
      decreasing = FALSE,
      empty.intersections = TRUE,
      scale.intersections = "identity",
      scale.sets = "identity",
      set_size.show = TRUE,
      set_size.scale_max = 100,
      mainbar.y.max = 100,
      mainbar.y.label	= "Shared DEMs per intersection (%)",
      sets.x.label = "Shared DEMs per group (%)",
      intersections =  list("facAB.dems",
                            list("facAB.dems","facA.dems"), 
                            list("facAB.dems","facB.dems"),
                            list("facAB.dems","facA.dems", "facB.dems")))

# PCA of Expression Profiles ----
# Set colors
pca.col.pal <- c("I" = "#1f78b4",
                 "II" ="#a6cee3",
                 "III" = "#e31a1c",
                 "IV"="#fb9a99",
                 "V"="#33a02c",
                 "VI"="#b2df8a",
                 "VII"="#ff7f00",
                 "VIII"="#fdbf6f")

FC.ExPr <- data.frame(A = logratio2foldchange(c.AB.DEMS$d.facA, base = 2),
                      B = logratio2foldchange(c.AB.DEMS$d.facB, base = 2),
                      AB  = logratio2foldchange(c.AB.DEMS$d.facAB, base = 2))

# PCA
data.class <- c.AB.DEMS$EPs
data.pca <- prcomp(FC.ExPr, center = T, scale. = T)
g <- ggbiplot(data.pca, 
              obs.scale = 1, 
              var.scale = 1, 
              groups = data.class, 
              ellipse = F, 
              circle = F)
g
fin.g <- g + 
  theme_bw() + 
  xlim(-2, 2) + 
  ylim(-2, 2) + 
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  scale_color_manual(name = "Expression Profiles", values = pca.col.pal, limits = names(pca.col.pal))
fin.g

# Summary Interaction Score ----
c.AB.iDEMs$SummIS <- log(c.AB.iDEMs$AbsIScore, base = 2)
ISplot <- c.AB.iDEMs[order(c.AB.iDEMs$SummIS),]
barplot(ISplot$SummIS, 
        horiz = TRUE,
        col = cond_col_is(ISplot$SummIS),
        border = NA)

write_xlsx(AB.DEMS, "SEE-LungH-RNAseq OBIF AB.DEMs.xlsx")

# Quick R object diagnostics ----
# Use these commands to verify the correct charactericts of any object in you Enviroment if any part of the code is not working 
class(sim.l2qn.expr.val)
typeof(sim.l2qn.expr.val)
attributes(expr.val)

