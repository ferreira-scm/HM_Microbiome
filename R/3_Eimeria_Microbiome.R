source("R/2_Parasite_cleaning.R")

library(vegan)

# betadiversity (BC distances), permanova and visualisatio with MDS



#hybridicity
PS.T@sam_data$hi <- abs(PS.T@sam_data$HI-0.5)

# subsetting dataset to not have NA's in permanova
sdata <- sample_data(PS.T)
sdata <- sdata[!is.na(sdata$hi),]
sdata <- sdata[!is.na(sdata$BMI),]
PS.T_sub <- subset_samples(PS.T, sample_names(PS.T)%in%rownames(sdata))
PS.T_sub

dist <- phyloseq::distance(PS.T_sub, method="bray", type="samples")

permaPS <- adonis2(dist~
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality)

permaPS

#### Now I want to know the effect of parasites on the eukaryotic and bacterial communities


PS.T_minusP <- subset_taxa(PS.T_sub, !Genus %in%c("Eimeria", "Cryptosporidium", "Oxyurida", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis"))

rank_names(PS.T)

get_taxa_unique(PS.T, "Kingdom")

Bac <- subset_taxa(PS.T_minusP, Kingdom%in%c("Bacteria", "Archaea"))
Bac <- phyloseq::prune_samples(sample_sums(Bac)>0, Bac)
Euk_minusP <- subset_taxa(PS.T_minusP, !Kingdom%in%c("Bacteria", "Archaea", NA, "unidentified"))
Euk_minusP <- phyloseq::prune_samples(sample_sums(Euk_minusP)>0, Euk_minusP)
get_taxa_unique(Euk_minusP, "Kingdom")

dist_bac <- phyloseq::distance(Bac, method="bray", type="samples")

dist_euk <- phyloseq::distance(Euk_minusP, method="bray", type="samples")

sdata_b <- sample_data(Bac)
sdata_e <- sample_data(Euk_minusP)

permaBac <- adonis2(dist_bac~
                        sdata_b$Eimeria_ferrisi_asv+
                        sdata_b$Eimeria_falciformis_asv+
                        sdata_b$Eimeria_vermiformis_asv+
                        sdata_b$Tritrichomonas_asv+
                        sdata_b$Trichuris_asv+
                        sdata_b$Hymenolepis_asv+
                        sdata_b$Ascaridida_asv+
                        sdata_b$Crypto_asv+
                        sdata_b$Mastophorus_asv+
                        sdata_b$Sex+
                        sdata_b$hi+
                        sdata_b$Year+
                        sdata_b$BMI+
                        sdata_b$Locality)

permaBac

permaEuk <- adonis2(dist_euk~
                        sdata_e$Eimeria_ferrisi_asv+
                        sdata_e$Eimeria_falciformis_asv+
                        sdata_e$Eimeria_vermiformis_asv+
                        sdata_e$Tritrichomonas_asv+
                        sdata_e$Trichuris_asv+
                        sdata_e$Hymenolepis_asv+
                        sdata_e$Ascaridida_asv+
                        sdata_e$Crypto_asv+
                        sdata_e$Mastophorus_asv+
                        sdata_e$Sex+
                        sdata_e$hi+
                        sdata_e$Year+
                        sdata_e$BMI+
                        sdata_e$Locality)

permaEuk

## richness, gotta get the untrimmed!!
erich = estimate_richness(Euk_minusP)

brich = estimate_richness(Bac)

plot_richness(Bac, measures="Shannon")

test <- estimateR(round(Bac@otu_table*100))

test

erich

Richdf <- data.frame(
    brichChao=brich$Observed,
    brichShan=brich$Shannon,
    brichObs=brich$Chao1,
    Year=sdata_b$Year,
    HI=sdata_b$hi,
    Sex=sdata_b$Sex,
    BMI=sdata_b$BMI,
    Locality=sdata_b$Locality,
    Eim_fer=sdata_e$Eimeria_ferrisi_asv,
    Eim_fal=sdata_e$Eimeria_falciformis_asv,
    Eim_ver=sdata_e$Eimeria_vermiformis_asv,
    Tritri=sdata_e$Tritrichomonas_asv,
    Trichuris=sdata_e$Trichuris_asv,
    Hymenolepis=sdata_e$Hymenolepis_asv,
    Ascaridida=sdata_e$Ascaridida_asv,
    Crypto=sdata_e$Crypto_asv,
    Mastophorus=sdata_e$Mastophorus_asv)
                        
