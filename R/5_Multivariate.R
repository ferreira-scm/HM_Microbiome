source("R/2_Parasite_cleaning.R")

library(mvabund)

PS.T_sub

KeepTaxap <- microbiome::prevalence(PS.T_sub)>0.10
PS10 <- phyloseq::prune_taxa(KeepTaxap, PS.T_sub)

ASVmv <- mvabund(otu_table(PS10))

boxplot(otu_table(PS10), horizontal=TRUE, las=2)

sdata <- sample_data(PS10)

Mod1 <- manyglm(ASVmv~
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
                        sdata$Mastophorus_asv+
                        sdata$Sex+
                        sdata$hi+
                        sdata$Year+
                        sdata$BMI+
                        sdata$Locality)

saveRDS(Mod1, "tmp/Mvabund_PS10.R")

plot(Mod1)

table <- anova(Mod1, p.uni = "adjusted")

which(table$uni.p<0.05)

table

saveRDS(table, "tmp/Mvabund_anovaPS10.R")

#table2 <- anova.manyglm(Mod1, p.uni = "adjusted", show.time="all")

#saveRDS(table2, "tmp/Mvabund_PS10_anova.manyglm.R")
