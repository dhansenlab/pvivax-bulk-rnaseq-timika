# Title: P. vivax bulk RNA-seq 
# Subtitle: Differential epxression analysis
# Authors: Steph Studniberg and Alex Garnham


######## Setup ######
library(janitor)
library(readxl)
library(kableExtra)
library(rtracklayer)
library(edgeR)
library(ggplot2)
library(gridExtra)
library(RUVSeq)
library(ComplexHeatmap)
load("data/Pvcounts.rData")
load("Pvannotation.rData")
targets <- read.csv("data/Pvtargets.csv")


###### Filtering and Normalisation ######
# Create DGEList
dge <- DGEList(counts = counts$counts,
               genes = annot,
               group = targets$Cohort,
               samples = targets)
dge$samples$group <- relevel(dge$samples$group, ref = "HC")
colnames(dge) <- dge$samples$CohortID
rownames(dge$counts) <- dge$genes$Symbol

# Filtering
dge <- dge[!is.na(dge$genes$Symbol),, keep.lib.sizes = FALSE] 
dge <- dge[!duplicated(dge$genes$Symbol),, keep.lib.sizes = FALSE]
dge <- dge[dge$genes$Symbol != "XIST",, keep.lib.sizes = FALSE]
chrY <- dge[dge$genes$Chromosome == "chrY", ]
chrY <- chrY$genes$EnsemblID[rownames(chrY$genes) %in% unlist(strsplit(rownames(chrY$genes), "_"))] 
dge <- dge[!dge$genes$EnsemblID %in% chrY,, keep.lib.sizes = FALSE]
allhb <- dge$genes[grepl("hemoglobin", dge$genes$Description), ]
dge <- dge[!dge$genes$Symbol %in% allhb$Symbol,, keep.lib.sizes = FALSE]
dge <- dge[!dge$genes$Symbol %in% c("HBA2", "HBA1", "HBD", "HBB"),, keep.lib.sizes = FALSE]
keep <- filterByExpr(dge, group = group)
dge <- dge[keep,, keep.lib.sizes = FALSE]

# Normalisation
dge <- calcNormFactors(dge, method = "TMM")


######## Adjusting for variation ######
## Generate negative control genes
design <- model.matrix(~0 + group)
colnames(design) <- c("HC", "AM", "SM")
v <- voom(dge, design = design, plot = FALSE)
fit <- lmFit(v, design)
contrasts <- makeContrasts(SMvHC = SM - HC,
                           SMvAM = SM - AM,
                           AMvHC = AM - HC,
                           levels = design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit, trend = FALSE, robust = TRUE)
tt <- topTable(fit = fit, number = Inf)
genes <- tt[(12359:13358), ]

# Generate residuals
ruv <- estimateGLMCommonDisp(dge, design)
ruv <- estimateGLMTagwiseDisp(ruv, design)
fit <- glmFit(ruv, design)
res <- residuals(fit, type = "deviance")

# Perform RUVr  
ruv_r <- vector(mode = "list", length = 5)
names(ruv_r) <- c(1:5)
for (i in names(ruv_r)){
ruv_r[[i]] <-  RUVr(x = dge$counts,
                    cIdx = dge$genes$Symbol[dge$genes$Symbol %in% rownames(genes)],
                    k = as.numeric(names(ruv_r[i])),
                    residuals = res)
}


# Save log_cpm surrogate variables for data viz
lcpm <- cpm(dge, log = TRUE) 
lcpm_ruv_r <- vector(mode="list", length=5)
names(lcpm_ruv_r) <- c(1:5)
par(mfrow=c(3,2))
for(i in names(lcpm_ruv_r)){
  lcpm_ruv_r[[i]] <- removeBatchEffect(x = lcpm, design = design, covariates = ruv_r[[i]]$W)
  }


######## DE analysis ######
design <- model.matrix(~0+group + ruv_r$`5`$W)
colnames(design) <- c("HC", "AM", "SM", "WW_1", "WW_2", "WW_3","WW_4","WW_5")
rownames(design) <- colnames(dge)
contrasts <- makeContrasts(SMvHC = SM - HC,
                           SMvAM = SM - AM,
                           AMvHC = AM - HC,
                           levels = design)
v <- voom(dge, design = design, plot = FALSE)
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit, robust = TRUE)
dt <- decideTests(fit, p.value = 0.01)

save(v, file = "Pv_expr.rData")
save(fit, file = "Pvfit.rData")
save(dt, file = "Pvdt.rData")
save(tt_SMvHC, file = "topTable_SMvHC.rData")
save(tt_SMvAM, file = "topTable_SMvAM.rData")
save(tt_AMvHC, file = "topTable_AMvHC.rData")