
source("R/2_Parasite_cleaning.R")

library(vegan)

# betadiversity (BC distances), permanova and visualisatio with MDS
dist <- phyloseq::distance(PS.T_sub, method="bray", type="samples")

permaPS <- adonis2(dist~
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality)

permaPS

Parasite_sub <- subset_taxa(PS.T_sub, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

Parasite_sub <- phyloseq::prune_samples(sample_sums(Parasite_sub)>0, Parasite_sub)

dist_para <- phyloseq::distance(Parasite_sub, method="bray", type="samples")

sdata_p <- sample_data(Parasite_sub)

permaPara <- adonis2(dist_para~
                   sdata_p$Sex+
                   sdata_p$hi+
                   sdata_p$BMI+
                   sdata_p$Year+
                   sdata_p$Locality)

permaPara

#### Now I want to know the effect of parasites on the eukaryotic and bacterial communities
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
                        sdata_b$Syphacia_asv+
                        sdata_b$Aspiculuris_asv+
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
                        sdata_b$Syphacia_asv+
                        sdata_b$Aspiculuris_asv+
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


## let's just combine and analysi gut community instead for separate domains
# parasites explain ~ 7%
permaPS_para <- adonis2(dist~
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality+
                        sdata$Eimeria_ferrisi_asv+
                        sdata$Eimeria_falciformis_asv+
                        sdata$Eimeria_vermiformis_asv+
                        sdata$Tritrichomonas_asv+
                        sdata$Syphacia_asv+
                        sdata$Aspiculuris_asv+
                        sdata$Trichuris_asv+
                        sdata$Hymenolepis_asv+
                        sdata$Ascaridida_asv+
                        sdata$Crypto_asv+
                        sdata$Mastophorus_asv)
permaPS_para

PS.ord <- ordinate(PS.T, "MDS", "bray")

p1=plot_ordination(PS.T, PS.ord, type="sample", color="Co_infb")

p1+theme_bw()+stat_ellipse(aes(group=PS.T@sam_data$Co_infb))


PS.df <- as(sample_data(PS.T), "data.frame")

groups <- PS.df$Co_infb

library(vegan)

PS.dis <- phyloseq::distance(PS.T@otu_table, "bray")

mod <- betadisper(PS.dis, groups)



anova(mod) # the dispersion is different between groups, then examine

plot(mod, hull=FALSE, ellipse=TRUE)

boxplot(mod)

mod.HSD <- TukeyHSD(mod )
plot(mod.HSD)

co_adonis <- adonis2(PS.dis~Co_infb, data=PS.df)

co_adonis

permaPS_T <- adonis2(dist~
                   sdata$Co_infb+
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality)

permaPS_T

##### ok, now let's see if co-infection classification is associated with microbiome
permaPS_C <- adonis2(dist~
                   sdata$Co_type+
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality)

permaPS_C

PS.ord_C <- ordinate(PS.T, "MDS", "bray")

p2=plot_ordination(PS.T, PS.ord_C, type="sample", color="Co_type")
p2+theme_bw()+stat_ellipse(aes(group=PS.T@sam_data$Co_type))

PS.df <- as(sample_data(PS.T), "data.frame")

Co_type <- PS.df$Co_type

library(vegan)
PS.dis <- phyloseq::distance(PS.T@otu_table, "bray")

mod_type <- betadisper(PS.dis, Co_type)

anova(mod_type) # the dispersion is different between groups, then examine

plot(mod_type)

boxplot(mod_type)

mod.HSD_type <- TukeyHSD(mod_type)
plot(mod.HSD_type)


### alpha diversity next

get_taxa_unique(Parasite, "Genus")

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
                        
