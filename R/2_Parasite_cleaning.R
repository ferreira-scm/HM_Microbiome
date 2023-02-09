
source("R/1_filtering.R")

source("R/Correlation_net.R")


############# We need to distinguis between Oxyurida ASVs
## create reference sequences for Syhacia and Apiculuris and align them to ASV
library(DECIPHER)
library(phangorn)
library(ShortRead)
library(dada2)
library(dplyr)


Oxy <- subset_taxa(PS.T, Genus%in%"Oxyurida")

seqs <- DNAStringSet(getSequences(colnames(Oxy@otu_table)))

parasite <- "Oxyurida"
amp <-  amp_func(PS.lT, parasite="Oxyurida")
names(seqs) <- amp_func(PS.lT, parasite="Oxyurida")

seqs18 <- seqs[-grep("D3A_", names(seqs))]
seqs28 <- seqs[grep("D3A_", names(seqs))]

#writeFasta(seqs, "tmp/Oxyurida.fa")

## get 18S sequences from Aspiculuris, Syphacia and...
access.aspi18 <- c("EF464551.1", "MH215350.1", "KP338606.1", "KY462827.1", "MT755640.1", "KY462828.1", "MT613322.1")

access.syph18 <- c("EF464553.1", "EF464554.1", "AB629697.1", "KY462829.1", "KY462826.1","OK138900.1", "OK138907.1", "OK138885.1", "OK138886.1", "OK138897.1", "OK138899.1", "OK138904.1", "OK138908.1","EU263105.2", "MT135057.1", "MT135058.1")

# retrieve sequences from NBCI and convert formats and label
aspi18.db <- read.GenBank(access.aspi18)
aspi18.db.ID <- paste(names(aspi18.db), attr(aspi18.db, "species"), sep="_")
aspi18.DB <- aspi18.db %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
syph18.db <- read.GenBank(access.syph18)
syph18.db.ID <- paste(names(syph18.db), attr(syph18.db, "species"), sep="_")
syph18.DB <- syph18.db %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(aspi18.DB) <- aspi18.db.ID
names(syph18.DB) <- syph18.db.ID
outg <- read.GenBank(c("JF934731.1", "FR687850.1", "AB626598.1"))# don't forget an outgroup
outg.ID <- paste(names(outg), attr(outg, "species"), sep="_")
#convert to DNAStringset
outg.db <- outg %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(outg.db) <- outg.ID

# uncomment if rerunning is necessary
#align18 <-  AlignSeqs(c(outg.db, syph18.DB, aspi18.DB, seqs18), anchor=NA, iterations=20, refinements=20, processors=90)
#writeFasta(align18, "tmp/Oxyurida_tree/oxyurida18S.fa")

#run tree outside R
#~/iqtree-2.2.0-Linux/bin/iqtree2 -s tmp/Oxyurida_tree/oxyurida18S.fa -m MFP -B 5000 -T AUTO

# read tree into R
library(ggtree)
phyloOxy <- ggtree::read.tree("tmp/Oxyurida_tree/oxyurida18S.fa.contree")

phyloO <- root(phyloOxy, outg.ID)

clade.nr <- ggtree::MRCA(phyloO, syph18.db.ID)
clade.nr.aspi <- ggtree::MRCA(phyloO, aspi18.db.ID)
syph <- extract.clade(phyloO, clade.nr)
aspi <- extract.clade(phyloO, clade.nr.aspi)

syphASV <- which(names(seqs)%in%syph$tip.label)
aspiASV <- which(names(seqs)%in%aspi$tip.label)

ampdf <-as.data.frame(amp)

ampdf$species <- "name"

ampdf$species[syphASV] <- "Syphacia"
ampdf$species[aspiASV] <- "Aspiculuris"
ampdf$species[ampdf$species=="name"] <- "28S"

ampdf

PS.T@tax_table[which(PS.T@tax_table[,6]=="Oxyurida"),7] <- ampdf$species
fPS@tax_table[which(fPS@tax_table[,6]=="Oxyurida"),7] <- ampdf$species

# now the correlation networks

######## Now we do a correlation network to see the likelihood of ASV within the same genus being from the same species.
# Oxyurida
parasite <- "Oxyurida"
Correlation_net(PS.lT, PS.l, PS.T, "Oxyurida") # several species

PS.T@tax_table[which(PS.T@tax_table[,7]=="28S"),7]<- "Aspiculuris"
PS.T@tax_table[which(PS.T@tax_table[,7]=="Aspiculuris"),6]<- "Aspiculuris"
PS.T@tax_table[which(PS.T@tax_table[,7]=="Syphacia"),6]<- "Syphacia"
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Syphacia")])
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Aspiculuris")])

sample_data(PS.T)$Syphacia_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Syphacia"))
sample_data(PS.T)$Aspiculuris_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Aspiculuris"))

cor.test(log(1+sample_data(PS.T)$Aspiculuris_asv), log(1+sample_data(PS.T)$Aspiculuris_sp))
cor.test(log(1+sample_data(PS.T)$Syphacia_asv), log(1+sample_data(PS.T)$Syphacia_sp))

# Tritrichomonas
parasite <- "Tritrichomonas"
Correlation_net(PS.lT, PS.l, PS.T, "Tritrichomonas")# likely 1 species: merge

PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Tritrichomonas")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Tritrichomonas"),7] <- "Tritrichomonas_sp"
sample_data(PS.T)$Tritrichomonas_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Tritrichomonas"))

# Cryptosporidum
parasite <- "Cryptosporidium"
Correlation_net(PS.lT, PS.l, PS.T, "Cryptosporidium")# too sparce but likely 1 sp, merge?

PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Cryptosporidium")])
sample_data(PS.T)$Crypto_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Cryptosporidium"))


#Ascaridida # heterakis?
parasite <- "Ascaridida"
Correlation_net(PS.lT, PS.l, PS.T, "Ascaridida")# likely 1 sp, merge

PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Ascaridida")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Ascaridida"),7] <- "Ascaridida"
sample_data(PS.T)$Ascaridida_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Ascaridida"))

cor.test(log(1+sample_data(PS.T)$Ascaridida_asv), log(1+sample_data(PS.T)$Heterakis))
#plot(log(1+sample_data(PS.T)$Ascaridida_asv), log(1+sample_data(PS.T)$Heterakis))

# Spirurida, likely Mastophorus (high correlation with worm counts)
parasite <- "Spirurida"
Correlation_net(PS.lT, PS.l, PS.T, "Spirurida")# too sparce, likely 1 sp, merge?

PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Spirurida")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Spirurida"),7] <- "Mastophorus_muris"
PS.T@tax_table[which(PS.T@tax_table[,6]=="Spirurida"),6] <- "Mastophorus"
sample_data(PS.T)$Mastophorus_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Mastophorus"))

cor.test(log(1+sample_data(PS.T)$Mastophorus_asv), log(1+sample_data(PS.T)$Mastophorus_muris))
#plot(sample_data(PS.T)$Mastophorus_asv, sample_data(PS.T)$Mastophorus_muris)

#Trichocephalida, likely Trichuris muris
parasite <- "Trichocephalida"
Correlation_net(PS.lT, PS.l, PS.T, "Trichocephalida")# 1 species, merge

PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Trichocephalida")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Trichocephalida"),7] <- "Trichuris_muri"
PS.T@tax_table[which(PS.T@tax_table[,6]=="Trichocephalida"),6] <- "Trichuris"
sample_data(PS.T)$Trichuris_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Trichuris"))

cor.test(log(1+sample_data(PS.T)$Trichuris_asv), log(1+sample_data(PS.T)$Trichuris_muri))
#plot(log(1+sample_data(PS.T)$Trichuris_asv), log(1+sample_data(PS.T)$Trichuris_muri))

# Cyclophyllidae, likely Hymenolepis, not Catenotaenia_pusilla, nor Taenia
parasite <- "Cyclophyllidea"
Correlation_net(PS.lT, PS.l, PS.T, "Cyclophyllidea")# 1 species, merge

PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Cyclophyllidea")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Cyclophyllidea"),7] <- "Hymenolepis_sp"
PS.T@tax_table[which(PS.T@tax_table[,6]=="Cyclophyllidea"),6] <- "Hymenolepis"

sample_data(PS.T)$Hymenolepis_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Hymenolepis"))

cor.test(log(1+sample_data(PS.T)$Hymenolepis_asv), log(1+sample_data(PS.T)$Hymenolepis_sp))
#plot(log(1+sample_data(PS.T)$Hymenolepis_asv), log(1+sample_data(PS.T)$Hymenolepis_sp))

# adding eimeria intensity to sample data
sample_data(PS.T)$Eimeria_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Eimeria"))
sample_data(PS.T)$Eimeria_ferrisi_asv <- sample_sums(subset_taxa(PS.T, Species %in%"ferrisi"))
sample_data(PS.T)$Eimeria_falciformis_asv <- sample_sums(subset_taxa(PS.T, Species %in%"falciformis"))
sample_data(PS.T)$Eimeria_vermiformis_asv <- sample_sums(subset_taxa(PS.T, Species %in%"vermiformis"))


## for fPS

fPS@tax_table[which(fPS@tax_table[,7]=="28S"),7]<- "Aspiculuris"
fPS@tax_table[which(fPS@tax_table[,7]=="Aspiculuris"),6]<- "Aspiculuris"
fPS@tax_table[which(fPS@tax_table[,7]=="Syphacia"),6]<- "Syphacia"

fPS <- merge_taxa(fPS, taxa_names(fPS)[which(fPS@tax_table[,6]=="Syphacia")])
fPS <- merge_taxa(fPS, taxa_names(fPS)[which(fPS@tax_table[,6]=="Aspiculuris")])

sample_data(fPS)$Syphacia_asv <- sample_sums(subset_taxa(fPS, Genus %in%"Syphacia"))
sample_data(fPS)$Aspiculuris_asv <- sample_sums(subset_taxa(fPS, Genus %in%"Aspiculuris"))

fPS <- merge_taxa(fPS, taxa_names(fPS)[which(fPS@tax_table[,6]=="Tritrichomonas")])
fPS@tax_table[which(fPS@tax_table[,6]=="Tritrichomonas"),7] <- "Tritrichomonas_sp"
sample_data(fPS)$Tritrichomonas_asv <- sample_sums(subset_taxa(fPS, Genus %in%"Tritrichomonas"))

fPS <- merge_taxa(fPS, taxa_names(fPS)[which(fPS@tax_table[,6]=="Cryptosporidium")])
sample_data(fPS)$Crypto_asv <- sample_sums(subset_taxa(fPS, Genus %in%"Cryptosporidium"))

fPS <- merge_taxa(fPS, taxa_names(fPS)[which(fPS@tax_table[,6]=="Ascaridida")])
fPS@tax_table[which(fPS@tax_table[,6]=="Ascaridida"),7] <- "Ascaridida"
sample_data(fPS)$Ascaridida_asv <- sample_sums(subset_taxa(fPS, Genus %in%"Ascaridida"))

fPS <- merge_taxa(fPS, taxa_names(fPS)[which(fPS@tax_table[,6]=="Spirurida")])
fPS@tax_table[which(fPS@tax_table[,6]=="Spirurida"),7] <- "Mastophorus_muris"
fPS@tax_table[which(fPS@tax_table[,6]=="Spirurida"),6] <- "Mastophorus"
sample_data(fPS)$Mastophorus_asv <- sample_sums(subset_taxa(fPS, Genus %in%"Mastophorus"))

fPS <- merge_taxa(fPS, taxa_names(fPS)[which(fPS@tax_table[,6]=="Trichocephalida")])
fPS@tax_table[which(fPS@tax_table[,6]=="Trichocephalida"),7] <- "Trichuris_muri"
fPS@tax_table[which(fPS@tax_table[,6]=="Trichocephalida"),6] <- "Trichuris"

sample_data(fPS)$Trichuris_asv <- sample_sums(subset_taxa(fPS, Genus %in%"Trichuris"))

fPS <- merge_taxa(fPS, taxa_names(fPS)[which(fPS@tax_table[,6]=="Cyclophyllidea")])
fPS@tax_table[which(fPS@tax_table[,6]=="Cyclophyllidea"),7] <- "Hymenolepis_sp"
fPS@tax_table[which(fPS@tax_table[,6]=="Cyclophyllidea"),6] <- "Hymenolepis"

sample_data(fPS)$Hymenolepis_asv <- sample_sums(subset_taxa(fPS, Genus %in%"Hymenolepis"))

# adding eimeria intensity to sample data
sample_data(fPS)$Eimeria_asv <- sample_sums(subset_taxa(fPS, Genus %in%"Eimeria"))
sample_data(fPS)$Eimeria_ferrisi_asv <- sample_sums(subset_taxa(fPS, Species %in%"ferrisi"))
sample_data(fPS)$Eimeria_falciformis_asv <- sample_sums(subset_taxa(fPS, Species %in%"falciformis"))
sample_data(fPS)$Eimeria_vermiformis_asv <- sample_sums(subset_taxa(fPS, Species %in%"vermiformis"))



fPS
PS.T



### Now parasite community
Parasite <- subset_taxa(PS.T, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

PS.T@sam_data$E_fer[PS.T@sam_data$Eimeria_ferrisi_asv>0] <- 1
PS.T@sam_data$E_fer[PS.T@sam_data$Eimeria_ferrisi_asv==0] <- 0

PS.T@sam_data$E_fal[PS.T@sam_data$Eimeria_falciformis_asv>0] <- 1
PS.T@sam_data$E_fal[PS.T@sam_data$Eimeria_falciformis_asv==0] <- 0

PS.T@sam_data$E_ver[PS.T@sam_data$Eimeria_vermiformis_asv>0] <- 1
PS.T@sam_data$E_ver[PS.T@sam_data$Eimeria_vermiformis_asv==0] <- 0

PS.T@sam_data$cryp[PS.T@sam_data$Crypto_asv>0] <- 1
PS.T@sam_data$cryp[PS.T@sam_data$Crypto_asv==0] <- 0

PS.T@sam_data$aspi[PS.T@sam_data$Aspiculuris_asv==0] <- 0
PS.T@sam_data$aspi[PS.T@sam_data$Aspiculuris_asv>0] <- 1

PS.T@sam_data$syph[PS.T@sam_data$Syphacia_asv>0] <- 1
PS.T@sam_data$syph[PS.T@sam_data$Syphacia_asv==0] <- 0

PS.T@sam_data$tri[PS.T@sam_data$Tritrichomonas_asv>0] <- 1
PS.T@sam_data$tri[PS.T@sam_data$Tritrichomonas_asv==0] <- 0

PS.T@sam_data$mas[PS.T@sam_data$Mastophorus_asv>0] <- 1
PS.T@sam_data$mas[PS.T@sam_data$Mastophorus_asv==0] <- 0

PS.T@sam_data$tric[PS.T@sam_data$Trichuris_asv>0] <- 1
PS.T@sam_data$tric[PS.T@sam_data$Trichuris_asv==0] <- 0

PS.T@sam_data$hymn[PS.T@sam_data$Hymenolepis_asv>0] <- 1
PS.T@sam_data$hymn[PS.T@sam_data$Hymenolepis_asv==0] <- 0

PS.T@sam_data$asca[PS.T@sam_data$Ascaridida_asv==0] <- 0
PS.T@sam_data$asca[PS.T@sam_data$Ascaridida_asv>0] <- 1


### now categories of co-infection
PS.T@sam_data$CE_fer[PS.T@sam_data$Eimeria_ferrisi_asv>0] <- "ferrisi"
PS.T@sam_data$CE_fer[PS.T@sam_data$Eimeria_ferrisi_asv==0] <- ""

PS.T@sam_data$CE_fal[PS.T@sam_data$Eimeria_falciformis_asv>0] <- "falciformis"
PS.T@sam_data$CE_fal[PS.T@sam_data$Eimeria_falciformis_asv==0] <- ""

PS.T@sam_data$CE_ver[PS.T@sam_data$Eimeria_vermiformis_asv>0] <- "vermiformis"
PS.T@sam_data$CE_ver[PS.T@sam_data$Eimeria_vermiformis_asv==0] <- ""

PS.T@sam_data$Ccryp[PS.T@sam_data$Crypto_asv>0] <- "Cryptosporidium"
PS.T@sam_data$Ccryp[PS.T@sam_data$Crypto_asv==0] <- ""

PS.T@sam_data$Caspi[PS.T@sam_data$Aspiculuris_asv==0] <- ""
PS.T@sam_data$Caspi[PS.T@sam_data$Aspiculuris_asv>0] <- "Aspiculuris"

PS.T@sam_data$Csyph[PS.T@sam_data$Syphacia_asv>0] <- "Syphacia"
PS.T@sam_data$Csyph[PS.T@sam_data$Syphacia_asv==0] <- ""

PS.T@sam_data$Ctri[PS.T@sam_data$Tritrichomonas_asv>0] <- "Tritrichomonas"
PS.T@sam_data$Ctri[PS.T@sam_data$Tritrichomonas_asv==0] <- ""

PS.T@sam_data$Cmas[PS.T@sam_data$Mastophorus_asv>0] <- "Mastophorus"
PS.T@sam_data$Cmas[PS.T@sam_data$Mastophorus_asv==0] <- ""

PS.T@sam_data$Ctric[PS.T@sam_data$Trichuris_asv>0] <- "Trichuris"
PS.T@sam_data$Ctric[PS.T@sam_data$Trichuris_asv==0] <- ""

PS.T@sam_data$Chymn[PS.T@sam_data$Hymenolepis_asv>0] <- "Hymenolepis"
PS.T@sam_data$Chymn[PS.T@sam_data$Hymenolepis_asv==0] <- ""

PS.T@sam_data$Casca[PS.T@sam_data$Ascaridida_asv==0] <- ""
PS.T@sam_data$Casca[PS.T@sam_data$Ascaridida_asv>0] <- "Ascaridida"


PS.T@sam_data$Co_inf <- PS.T@sam_data$asca+PS.T@sam_data$hymn + PS.T@sam_data$mas + PS.T@sam_data$tri+ PS.T@sam_data$syph + PS.T@sam_data$aspi + PS.T@sam_data$cryp + PS.T@sam_data$E_fer + PS.T@sam_data$E_ver + PS.T@sam_data$E_fal

summary(as.factor(PS.T@sam_data$Co_inf))

PS.T@sam_data$Co_infb <- PS.T@sam_data$Co_inf

PS.T@sam_data$Co_infb[as.numeric(PS.T@sam_data$Co_infb)>3] <- 3

PS.T@sam_data$Co_inf <- as.factor(PS.T@sam_data$Co_inf)

#hybridicity
PS.T@sam_data$hi <- abs(PS.T@sam_data$HI-0.5)


#### let's look at specific co-infections. Which are the most common?
PS.T@sam_data$Co_inf

PS.T@sam_data$CCo_inf <-paste(PS.T@sam_data$CE_fer, PS.T@sam_data$CE_fal, PS.T@sam_data$CE_ver, PS.T@sam_data$Ccryp, PS.T@sam_data$Caspi, PS.T@sam_data$Caspa, PS.T@sam_data$Csyph, PS.T@sam_data$Ctri, PS.T@sam_data$Ctric, PS.T@sam_data$Cmas, PS.T@sam_data$Chymn, PS.T@sam_data$Casca, sep="_")
PS.T@sam_data$CCo_inf <- gsub("_+", "_", PS.T@sam_data$CCo_inf)
PS.T@sam_data$CCo_inf <- gsub("^_", "", PS.T@sam_data$CCo_inf)
PS.T@sam_data$CCo_inf <- gsub("_$", "", PS.T@sam_data$CCo_inf)
PS.T@sam_data$CCo_inf[PS.T@sam_data$CCo_inf==""] <- "Uninfected"
co.df <- table(as.factor(PS.T@sam_data$CCo_inf))
#plot(co.df[order(co.df)])# long tailed!!!

## let's classify between uni and multicellular
unicel <- c("ferrisi|falciformis|vermiformis|Tritrichomonas|Cryptosporidium")
multicel <- c("Aspiculuris|Syphacia|Ascaridida|Mastophorus|Trichuris|Hymenolepis")
PS.T@sam_data$Co_type <- ""
PS.T@sam_data$Co_type[grep(unicel, PS.T@sam_data$CCo_inf)] <- "unicel"
PS.T@sam_data$Co_type[grep(multicel, PS.T@sam_data$CCo_inf)] <- paste(PS.T@sam_data$Co_type[grep(multicel, PS.T@sam_data$CCo_inf)], "multicel", sep="")
PS.T@sam_data$Co_type[PS.T@sam_data$Co_type==""] <- "Uninfected"


# subsetting dataset to not have NA's for the models
sdata <- sample_data(PS.T)
sdata <- sdata[!is.na(sdata$hi),]
sdata <- sdata[!is.na(sdata$BMI),]
PS.T_sub <- subset_samples(PS.T, sample_names(PS.T)%in%rownames(sdata))
PS.T_sub

PS.T_minusP <- subset_taxa(PS.T_sub, !Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))



