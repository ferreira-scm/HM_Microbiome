library(phyloseq)
library(vegan)
PS.T <- readRDS("tmp/PS.T.rds")
sdata <- PS.T@sam_data

# betadiversity (BC distances), permanova and visualisatio with MDS
dist_bray <- vegdist(PS.T@otu_table, method="bray")

PS.mP <- subset_taxa(PS.T, !Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))
dist_mP <- vegdist(PS.mP@otu_table, method="bray")

Parasite <- subset_taxa(PS.T, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))
Parasite <- phyloseq::prune_samples(sample_sums(Parasite)>0, Parasite)
dist_para <- phyloseq::distance(Parasite, method="bray", type="samples")

sdata_p <- Parasite@sam_data

permaPara <- adonis2(dist_para~
                   sdata_p$Sex+
                   sdata_p$hi+
                   sdata_p$BMI+
                   sdata_p$Year+
#                   sdata_p$Co_infb+
                   sdata_p$Locality,
                   by="margin")

permaPara

## testing now the Anna Karenina principle: each dysbiotic community is dysbiotic on its own way
PS.mP@sam_data$Co_infb <- as.factor(PS.mP@sam_data$Co_infb)

PS.ord <- ordinate(PS.mP, "MDS", "bray")
p1=plot_ordination(PS.mP, PS.ord, type="sample", color="Co_infb")

p1+
    scale_color_manual(labels=c("uninfected", "single parasite", "multiple parasites"), values = c("wheat", "tomato", "tomato4"))+
    labs(color="Parasite infection")+
    theme_bw()+
    stat_ellipse(aes(group=PS.mP@sam_data$Co_infb)) 

PS.df <- as(sample_data(PS.mP), "data.frame")

groups <- PS.df$Co_infb
#PS.dis <- phyloseq::distance(PS.T@otu_table, "bray")
mod <- betadisper(dist_mP, groups)

anova(mod) # the dispersion is not different between groups
plot(mod, hull=FALSE, ellipse=TRUE)
boxplot(mod) # the dispersion is not different between groups
mod.HSD <- TukeyHSD(mod)
plot(mod.HSD)# the dispersion is not different between groups

#### No evidence for stochastic effect between parasite infection, rather deterministic one

### alpha diversity next
## richness, gotta get the untrimmed!!
shannon <- vegan::diversity(PS.mP@otu_table, index = "shannon")   
simpson <- vegan::diversity(PS.mP@otu_table, index="simpson")

PS.df$shannon <- shannon
PS.df$simpson <- simpson

cor.test(as.numeric(PS.df$Co_inf), PS.df$shannon)# zero!

shan_parasite <- ggplot(data=PS.df, aes(y=shannon, x=Co_infb, fill=Co_infb))+
    geom_boxplot(alpha=0.5, outlier.shape = NA)+
    scale_fill_brewer(palette="Dark2")+
    geom_jitter(position=position_jitter(0.3), alpha=0.5)+
    labs(x="Parasite infection", y="Shannon diversity index")+
    scale_x_discrete(labels=c("uninfected", "single parasite", "multiple parasites"))+
    theme_bw()+
    theme(legend.position="none")

ggplot2::ggsave(file="fig/Shannon_Parasite.pdf", shan_parasite, width = 6, height = 4, dpi = 600)

summary(aov(shannon~Co_infb, PS.df))

library(lme4)
library(lmerTest)
library(effects)

df$Co_inf <- as.numeric(df$Co_inf)

## testing the effect of parasites in shanon diversity

############# This gets complicated, let's do it with BMRS
shaModel <- lmer(shannon~BMI+
                     Eimeria_ferrisi_asv+
                     Eimeria_falciformis_asv+
                     Eimeria_vermiformis_asv+
                     Hymenolepis_asv+
                     Trichuris_asv+
                     Mastophorus_asv+
                     Ascaridida_asv+
                     Crypto_asv+
                     Tritrichomonas_asv+
                     Aspiculuris_asv+
                     Syphacia_asv+
                     hi+Sex+Year+(1|Locality), PS.df)

summary(shaModel)

plot(Effect("Tritrichomonas_asv", shaModel))

0.00107+0.00260+0.00101+0.00129+0.00139+0.00131+0.00147+0.00143+0.00507+0.00257+0.00476+0.00491

sdata$Co_infb <- as.factor(sdata$Co_infb)

permaPS_T_para <- adonis2(dist_mP~
                         sdata$Eimeria_ferrisi_asv+
                         sdata$Eimeria_falciformis_asv+
                         sdata$Eimeria_vermiformis_asv+
                         sdata$Hymenolepis_asv+
                         sdata$Trichuris_asv+
                         sdata$Mastophorus_asv+
                         sdata$Ascaridida_asv+
                         sdata$Crypto_asv+
                         sdata$Tritrichomonas_asv+
                         sdata$Aspiculuris_asv+
                         sdata$Syphacia_asv+
                   sdata$Co_infb+
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality,
                   by="margin")


permaPS_T_para

