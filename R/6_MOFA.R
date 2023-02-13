source("R/2_Parasite_cleaning.R")

                                        # MOFA+: analysis of multi-modal microbiome data
# remove host, human and plants
PS.Tb <- subset_taxa(PS.TSS, !Phylum%in%c("Anthophyta", "Vertebrata", "Ochrophyta","Phragmoplastophyta"))

## ajusting Domain name
euk <- c("Alveolata", "Fungi", "Viridiplantae" , "Animalia", "Eukaryota", "Stramenopiles", "Chloroplastida", "Metazoa", "Amoebozoa")
pro <- c("Bacteria", "Archaea")
Unk <- c(NA, "unidentified")
PS.Tb@tax_table[which(PS.Tb@tax_table[,1]%in%euk),1] <- "Eukaryote"
PS.Tb@tax_table[which(PS.Tb@tax_table[,1]%in%pro),1] <- "Prokaryote"
PS.Tb@tax_table[which(PS.Tb@tax_table[,1]%in%Unk),1] <- "Unknown"
tax_table(PS.Tb)[PS.Tb@tax_table[,1]%in%Unk] <- "Unknown"
tax_table(PS.Tb)[PS.Tb@tax_table[,1]%in%euk] <- "Eukaryote"
tax_table(PS.Tb)[PS.Tb@tax_table[,1]%in%pro] <- "Prokaryote"

PSmP <- subset_taxa(PS.Tb, !Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

fP <- subset_taxa(PS.Tb, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

KeepTaxap <- microbiome::prevalence(PSmP)>0.10

PS10 <- phyloseq::prune_taxa(KeepTaxap, PSmP)

PS10 <- merge_phyloseq(PS10, fP)

## I think we might have to CLR transform per amplicon, but for now we keep going

PS10 <- microbiome::transform(PS10, transform="rclr")

library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)
#BiocManager::install("MOFA2")
library(MOFA2)

dt <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/microbiome/data.txt.gz")

str(dt)

head(dt)

df <- psmelt(PS10
df1 <- df[,c("Sample", "Abundance", "Kingdom", "OTU")]

df1$feature

names(df1) <- c("sample", "value", "view", "feature")

#df1$OTU <- NULL

ggplot(df1, aes(x=value, colour=view))+
    geom_density()
#    facet_wrap(~Kingdom, nrow=1, scales="free")

metadata <- sample_data(PS10)

names(metadata)

keep <- c("Mouse_ID", "Sex", "Year", "Co_infb", "hi", "Co_type", "IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

metadata <- metadata[,keep]

metadata$sample <- metadata$Mouse_ID

metadata$Mouse_ID <- NULL

mofa <- create_mofa(data=df1, groups=NULL, extract_metadata=FALSE)

plot_data_overview(mofa)

get_default_data_options
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 10

model_opts

mofa <- prepare_mofa(mofa, model_options = model_opts)

mofa <- run_mofa(mofa, use_basilisk = TRUE, outfile="tmp/mofa_training.hdf5", save_data=TRUE)

mofa

samples_metadata(mofa) <- metadata

str(mofa)

plot_data_overview(mofa)

str(mofa)