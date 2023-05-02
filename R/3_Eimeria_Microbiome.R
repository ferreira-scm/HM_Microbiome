library(phyloseq)
library(vegan)
library(RColorBrewer
library(microshades)

PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")

sdata <- PS.TSS@sam_data
### Let's see what is going on here
Euk <- subset_taxa(PS.TSS, Kingdom%in%"Eukarya")
gen.e <- get_taxa_unique(Euk, "Genus")
gen.e[-grep("Unknown", gen.e)]

Bac <- subset_taxa(PS.TSS, Kingdom%in%"Bacteria")
gen.e <- get_taxa_unique(Bac, "Genus")
gen.e[-grep("Unknown", gen.e)]

Arc <- subset_taxa(PS.TSS, Kingdom%in%"Archaea")
gen.e <- get_taxa_unique(Arc, "Genus")

# Jaccard distances
jac <- vegdist(PS.TSS@otu_table, method="jaccard")

# Bray Curtis dissimilarity matrix
bc <- vegdist(PS.TSS@otu_table, method="bray")

permaJac1 <- adonis2(jac~
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality,
                   by="margin")

permabc1 <- adonis2(bc~
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality,
                   by="margin")


permaJac1
permabc1

### Plotting
PS.TSS@sam_data$t <- "a"

Euk <- subset_taxa(PS.TSS, Kingdom=="Eukarya")
Euk.m <- merge_samples(Euk, "t")
Prok<- subset_taxa(PS.TSS, Kingdom%in%c("Bacteria", "Archaea"))
Prok.m <- merge_samples(Prok, "t")
PS.m <- merge_samples(PS.TSS, "t")

PS.m <- transform_sample_counts(PS.m, function(x) 100*x/sum(x))
Euk.m <- transform_sample_counts(Euk.m, function(x) 100*x/sum(x))
Prok.m <- transform_sample_counts(Prok.m, function(x) 100*x/sum(x))

PS.m <- tax_glom(PS.m, "Kingdom")
ps.m <- psmelt(PS.m)

Euk.m <- tax_glom(Euk.m, "Genus")
Prok.m <- tax_glom(Prok.m, "Genus")

top30 <- names(sort(taxa_sums(Euk.m), TRUE)[1:10])
Euk.30 <- prune_taxa(top30, Euk.m)
euk.m <- psmelt(Euk.30)

euk.m <- euk.m[, c("OTU", "Abundance", "Phylum", "Family", "Genus")]

euk.m <- rbind(data.frame(OTU="other", Abundance=100-sum(euk.m$Abundance), Phylum="other", Family="other", Genus="other"), euk.m)
euk.m <- euk.m[order(euk.m$Abundance),]
euk.m$Genus <-factor(euk.m$Genus)

levels(euk.m$Genus) <- euk.m$Genus

top30 <- names(sort(taxa_sums(Prok.m), TRUE)[1:10])
Prok.30 <- prune_taxa(top30, Prok.m)
prok.m <- psmelt(Prok.30)
prok.m <- prok.m[, c("OTU", "Abundance", "Phylum", "Family", "Genus")]
prok.m <- rbind(data.frame(OTU="other", Abundance=100-sum(prok.m$Abundance), Phylum="other", Family="other", Genus="other"), prok.m)


#levels(euk.m$Genus) <-c("other", as.character(euk.m$Genus[!euk.m$Genus=="other"]))

library(ggpubr)

ggbarplot(ps.m, y="Abundance", fill="Kingdom", position=position_stack())+
        theme_minimal(base_size=12)+
    scale_fill_manual(values=col)+
    ylab("Abundance (%)")+
    xlab("")+
    labs(fill="Domain")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

nb <- length(unique(euk.m$Genus))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb)

ggbarplot(euk.m, y="Abundance", fill="Genus", position=position_stack())+
    theme_minimal(base_size=12)+
    scale_fill_manual(values=mycolors)+
    ylab("Abundance (%)")+
    xlab("")+
    labs(fill="Genus")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

nb <- length(unique(prok.m$Genus))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb)

ggbarplot(prok.m, y="Abundance", fill="Genus", position=position_stack())+
    theme_minimal(base_size=12)+
    scale_fill_manual(values=mycolors)+
    ylab("Abundance (%)")+
    xlab("")+
    labs(fill="Genus")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())



    geom_bar(position="stack")

col <- c(microshades_palette("micro_blue", 1, lightest=FALSE)[1], microshades_palette("micro_orange", 3)[2], microshades_palette("micro_cvd_turquoise", 3)[2], microshades_palette("micro_purple")[1])


a <- plot_bar(PS.m, fill="Kingdom")+
    theme_minimal(base_size=12)+
    scale_fill_manual(values=col)+
    ylab("Abundance (%)")+
    xlab("")+
    labs(fill="Domain")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())




plot_bar(Euk.30, fill="Genus")



ps.19 <- psmelt(PS.19)

names(ps.19)

ps.19 <- ps.19[, c("OTU", "Kingdom", "Family", "Genus")]

oth <- data.frame(OTU="other", Kingdom="other", Family="other", Genus="other")

ps.19 <- rbind(oth, ps.19)

ps.19

plot_bar(PS.19, fill="Genus")

#################### parasite stuff
# removing parasites from gut community
PS.mP <- subset_taxa(PS.TSS, !Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))
# Bray Curtis dissimilarity matrix
dist_mP <- vegdist(PS.mP@otu_table, method="bray")

# Jaccard distances
jac_mP <- vegdist(PS.mP@otu_table, method="jaccard")

# Isolating parasite community
Parasite <- subset_taxa(PS.TSS, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))
Parasite <- phyloseq::prune_samples(sample_sums(Parasite)>0, Parasite)
# BC dissimilary
dist_para <- vegdist(Parasite@otu_table, method="bray")

sdata_p <- Parasite@sam_data

permaPara <- adonis2(dist_para~
                   sdata_p$Sex+
                   sdata_p$hi+
                   sdata_p$BMI+
                   sdata_p$Year+
                   sdata_p$Co_infb,
                   strata=sdata_p$Locality,
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

PS.df

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

permaJac <- adonis2(jac_mP~
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


permaJac

permaPS_T_para
