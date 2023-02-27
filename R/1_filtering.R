library(ggplot2)
library(reshape)
library(microbiome)
library(vegan)

### this script does a quality filter (by amplicon) and then some transformations
# TSS --> Total sum scaling (per amplicon) and then merges all amplicon datasets into 1 (PS.TSS)
# T --> Total sum scaling multiplied by sample DNA concentration and then merges all amplicons
# into 1 dataset (PS.T)
# fPS is the filered dataset with no normalization

# We also remove the silva handlers from the taxonomic table (all the "g__", "s__"...), so this
#makes the script run slower than expected.

# finally, we do a bit of cleaning up in the sampla data, by adding inputed immune genes from Fay
# and by inputing the few missing values for BMI and HI


#### preprocessing: filtering and transforming
PS.l <- readRDS(file = "data/PhyloSeqList_HMHZ_All.Rds")
PSLab.l <- readRDS(file="data/PhyloSeqList_All_Tax_New.Rds")

# this is our filtering function
fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.005%
    keepTaxa = (x / sum(x) > 0.00005)
    summary(keepTaxa)
    ps = phyloseq::prune_taxa(keepTaxa, ps)
# plus prevalnce filter at 1%
    KeepTaxap <- microbiome::prevalence(ps)>0.01
    ps <- phyloseq::prune_taxa(KeepTaxap, ps)
# subset samples based on total read count (100 reads)
#ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 100)
    ps <- phyloseq::prune_samples(sample_sums(ps)>100, ps)
    ps
}

## filtering MA by amplicon
fPS.l <- list()
for (i in 1:length(PS.l)) {
    try(fPS.l[[i]] <- fil(PS.l[[i]]), silent=TRUE)
}

# and remove those handlers.
#test <- list()
for (i in 1:length(fPS.l)) {
    try(tax_table(fPS.l[[i]])[, colnames(tax_table(fPS.l[[i]]))] <- gsub(tax_table(fPS.l[[i]])[, colnames(tax_table(fPS.l[[i]]))], pattern="[a-z]__", replacement=""), silent=TRUE)
}

### now pool all amplicons
fPS <- fPS.l[[1]]
for (i in 2:length(fPS.l)){
    fPS <- try(merge_phyloseq(fPS,fPS.l[[i]]))
#    print(fPS)
}

#### let's transform by amplicon
PS.lTSS <- list()
PS.lT <- list()

for (i in 1:length(fPS.l)) {
    try(PS.lTSS[[i]] <- transform_sample_counts(fPS.l[[i]], function(x) x / sum(x)), silent=TRUE)
    try(PS.lT[[i]] <- PS.lTSS[[i]], silent=TRUE)
    try(PS.lT[[i]]@otu_table <- PS.lT[[i]]@otu_table*PS.lT[[i]]@sam_data$Concentration, silent=TRUE)
}

PS.TSS <- PS.lTSS[[1]]
for (i in 2:length(PS.lTSS)){
    PS.TSS <- try(merge_phyloseq(PS.TSS,PS.lTSS[[i]]))
}

PS.T <- PS.lT[[1]]
for (i in 2:length(PS.lT)){
    PS.T <- try(merge_phyloseq(PS.T,PS.lT[[i]]))
}


############# now for lab dataset
## filtering MA by amplicon
fPSLab.l <- list()
for (i in 1:length(PSLab.l)) {
    try(fPSLab.l[[i]] <- fil(PSLab.l[[i]]), silent=TRUE)
}

# and remove those handlers.
for (i in 1:length(fPSLab.l)) {
    try(tax_table(fPSLab.l[[i]])[, colnames(tax_table(fPSLab.l[[i]]))] <- gsub(tax_table(fPSLab.l[[i]])[, colnames(tax_table(fPSLab.l[[i]]))], pattern="[a-z]__", replacement=""), silent=TRUE)
}

# and pooling into one
fPSLab <- fPSLab.l[[1]]
for (i in 2:length(fPSLab.l)){
    fPSLab <- try(merge_phyloseq(fPSLab,fPSLab.l[[i]]))
}

#### let's transform by amplicon
PSLab.lTSS <- list()
PSLab.lT <- list()

for (i in 1:length(fPSLab.l)) {
    try(PSLab.lTSS[[i]] <- transform_sample_counts(fPSLab.l[[i]], function(x) x / sum(x)), silent=TRUE)
    try(PSLab.lT[[i]] <- PSLab.lTSS[[i]], silent=TRUE)
    try(PSLab.lT[[i]]@otu_table <- PSLab.lT[[i]]@otu_table*PSLab.lT[[i]]@sam_data$Concentration, silent=TRUE)
}

PSLab.TSS <- PSLab.lTSS[[1]]
for (i in 2:length(PSLab.lTSS)){
    PSLab.TSS <- try(merge_phyloseq(PSLab.TSS,PSLab.lTSS[[i]]))
}

PSLab.T <- PSLab.lT[[1]]
for (i in 2:length(PSLab.lT)){
    PSLab.T <- try(merge_phyloseq(PSLab.T,PSLab.lT[[i]]))
}

#Wild dataset:
# PS.T Relative abundances rescaled to sample DNA concentration
# PS.TSS Relative abundances - total sum scaling
# fPS filtered dataset
#Lab dataset:
# PSLab.T - Relative abundances rescaled to sample DNA concentration
#PSLab.TSS Relative abundances - total sum scaling
# fPSLab filtered dataset


# Now we want to relabel Eimeria ASV species as in https://github.com/ferreira-scm/Eimeria_AmpSeq/blob/master/R/Wild_8_AssignEim_tree_cor.R

Eim <- readRDS("/SAN/Susanas_den/gitProj/Eimeria_AmpSeq/tmp/Wild/EimeriaSpeciesAssign.RDS")

eim_rn <- which(rownames(PS.T@tax_table) %in% rownames(Eim@tax_table))

# sanity check
rownames(PS.T@tax_table[eim_rn])== rownames(Eim@tax_table)

PS.T@tax_table[eim_rn, 7] <- Eim@tax_table[,7]
fPS@tax_table[eim_rn, 7] <- Eim@tax_table[,7]
PS.TSS@tax_table[eim_rn, 7] <- Eim@tax_table[,7]
PS.CLR@tax_table[eim_rn, 7] <- Eim@tax_table[,7]

## ok now we agglomerate Eimeria species
PS.T1 <- subset_taxa(PS.T, !Genus%in%"Eimeria")
PS.TE <- subset_taxa(PS.T, Genus%in%"Eimeria")
PS.TE <- tax_glom(PS.TE, "Species")
PS.T <- merge_phyloseq(PS.TE, PS.T1)

PS.TSS1 <- subset_taxa(PS.TSS, !Genus%in%"Eimeria")
PS.TSSE <- subset_taxa(PS.TSS, Genus%in%"Eimeria")
PS.TSSE <- tax_glom(PS.TSSE, "Species")
PS.TSS <- merge_phyloseq(PS.TSSE, PS.TSS1)

fPS1 <- subset_taxa(fPS, !Genus%in%"Eimeria")
fPSE <- subset_taxa(fPS, Genus%in%"Eimeria")
fPSE <- tax_glom(fPSE, "Species")
fPS <- merge_phyloseq(fPSE, fPS1)

#### let's adjust the metadata too, to include Fay's immune genes that have been inputed
Immune <- read.csv("https://raw.githubusercontent.com/fayweb/Eimeria_mouse_immunity/main/output_data/2.imputed_MICE_data_set.csv")
#subsetting to filed animals only
Immune <- Immune[Immune$origin=="Field",]
# keeping wanted variables only
keep <- c("Mouse_ID","IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")
Nkeep <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")
Immune <- Immune[,keep]
head(Immune)
# fixing name structure
Immune$Mouse_ID <- gsub("AA", "AA_", Immune$Mouse_ID)
# removing samples that are not in our PS
Immune <- Immune[which(Immune$Mouse_ID%in%PS.T@sam_data$Mouse_ID),]
# sanity checj
Immune[match(PS.T@sam_data$Mouse_ID, Immune$Mouse_ID),"Mouse_ID"]==PS.T@sam_data$Mouse_ID
#and replace
PS.T@sam_data[,Nkeep] <- Immune[match(PS.T@sam_data$Mouse_ID, Immune$Mouse_ID),Nkeep]
fPS@sam_data[,Nkeep] <- Immune[match(PS.T@sam_data$Mouse_ID, Immune$Mouse_ID),Nkeep]
PS.TSS@sam_data[,Nkeep] <- Immune[match(PS.T@sam_data$Mouse_ID, Immune$Mouse_ID),Nkeep]

##### we are going to input that BMI value missing, and the HI values missing too.
library(mice)
df <- data.frame(PS.T@sam_data[,c("Mouse_ID", "HI", "BMI")])
df.t <- mice(df, m=5, maxit=50, meth='pmm', seed=500)
df.t <- complete(df.t,1)
# sanity check
all(PS.T@sam_data$Mouse_ID==df.t$Mouse_ID)
PS.T@sam_data$BMI <- df.t$BMI
PS.TSS@sam_data$BMI <- df.t$BMI
fPS@sam_data$BMI <- df.t$BMI
PS.T@sam_data$HI <- df.t$HI
PS.TSS@sam_data$HI <- df.t$HI
fPS@sam_data$HI <- df.t$HI

