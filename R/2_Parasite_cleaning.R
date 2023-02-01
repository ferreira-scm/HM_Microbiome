
source("R/1_filtering.R")

source("R/Correlation_net.R")


############# We need to distinguis between Oxyurida ASVs
## create reference sequences for Syhacia and Apiculuris and align them to ASV
library(DECIPHER)
library(phangorn)
library(ShortRead)
library(dada2)
library(dplyr)

# adding eimeria intensity to sample data
sample_data(PS.T)$Eimeria_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Eimeria"))
sample_data(PS.T)$Eimeria_ferrisi_asv <- sample_sums(subset_taxa(PS.T, Species %in%"ferrisi"))
sample_data(PS.T)$Eimeria_falciformis_asv <- sample_sums(subset_taxa(PS.T, Species %in%"falciformis"))
sample_data(PS.T)$Eimeria_vermiformis_asv <- sample_sums(subset_taxa(PS.T, Species %in%"vermiformis"))

Oxy <- subset_taxa(PS.T, Genus%in%"Oxyurida")

seqs <- DNAStringSet(getSequences(colnames(Oxy@otu_table)))

parasite <- "Oxyurida"
amp <- 

names(seqs) <- amp_func(PS.lT, parasite="Oxyurida")

seqs18 <- seqs[-grep("D3A_", names(seqs))]
seqs28 <- seqs[grep("D3A_", names(seqs))]

writeFasta(seqs, "tmp/Oxyurida.fa")

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

# now the correlation networks

######## Now we do a correlation network to see the likelihood of ASV within the same genus being from the same species.
# Oxyurida
parasite <- "Oxyurida"
Correlation_net(PS.lT, PS.l, PS.T, "Oxyurida") # several species

#sample_data(PS.T)$Syphacia_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Syphacia"))
#sample_data(PS.T)$Aspiculuris_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Aspiculuris"))

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
plot(log(1+sample_data(PS.T)$Hymenolepis_asv), log(1+sample_data(PS.T)$Hymenolepis_sp))
                                                       
### Now parasite community
Parasite <- subset_taxa(PS.T, Genus %in%c("Eimeria", "Cryptosporidium", "Syphachia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis"))


Parasite


