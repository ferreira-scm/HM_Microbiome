source("R/2_Parasite_cleaning.R")


library(vegan)
library(microbiome)
library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)
#BiocManager::install("MOFA2")
library(MOFA2)


# MOFA+: analysis of multi-modal microbiome data
# remove host, human and plants

PS.Tb <- subset_taxa(PS.T, !Phylum%in%c("Anthophyta", "Vertebrata", "Ochrophyta","Phragmoplastophyta"))
Unk <- c(NA, "unidentified")
PS.Tb <- subset_taxa(PS.Tb, !Kingdom %in%Unk)


## ajusting Domain name
euk <- c("Alveolata", "Fungi", "Viridiplantae" , "Animalia", "Eukaryota", "Stramenopiles", "Chloroplastida", "Metazoa", "Amoebozoa")
pro <- c("Bacteria", "Archaea")
PS.Tb@tax_table[which(PS.Tb@tax_table[,1]%in%euk),1] <- "Eukaryote"
PS.Tb@tax_table[which(PS.Tb@tax_table[,1]%in%pro),1] <- "Prokaryote"

get_taxa_unique(PS.Tb, "Genus")






tax_table(PS.Tb)[PS.Tb@tax_table[,1]%in%Unk] <- "Unknown"
tax_table(PS.Tb)[PS.Tb@tax_table[,1]%in%euk] <- "Eukaryote"
tax_table(PS.Tb)[PS.Tb@tax_table[,1]%in%pro] <- "Prokaryote"

## I think we might have to CLR transform per amplicon, but for now we keep going
PS.Tb <- phyloseq::prune_samples(sample_sums(PS.Tb)>0, PS.Tb)

#PS.Tb@otu_table <- PS.Tb@otu_table+1
#PS.Tb <- microbiome::transform(PS10, transform="clr")

#PS.Tb@otu_table <- decostand(PS.Tb@otu_table, "rclr")

taxa_names(PS.Tb) <- paste("ASV", 1:ntaxa(PS.Tb), PS.Tb@tax_table[,6], sep="_")

df <- psmelt(PS.Tb)

df1 <- df[,c("Sample", "OTU", "Abundance", "Kingdom")]
names(df1) <- c("sample", "feature", "value", "view")
#df1$feature

df1 <- data.table(df1)

df1[,length(unique(feature)),by="view"]

#df1$OTU <- NULL

ggplot(df1, aes(x=value, colour=view))+
    geom_density()
#    facet_wrap(~Kingdom, nrow=1, scales="free")

metadata <- sample_data(PS.Tb)

keep <- c("Mouse_ID", "Sex", "Locality", "Year", "Co_infb", "hi", "Co_type", "BMI", "IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

metadata <- metadata[,keep]

metadata$sample <- metadata$Mouse_ID

metadata$Mouse_ID <- NULL

immune <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

mofa <- create_mofa(data=df1, groups=NULL, extract_metadata=FALSE)

data_opts <- get_default_data_options(mofa)
model_opts <- get_default_model_options(mofa)
train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42


mofa <- prepare_mofa(mofa, model_options = model_opts, data_options=data_opts, training_options=train_opts)
mofa <- run_mofa(mofa, use_basilisk = TRUE)


class(metadata) <- "data.frame"
metadata <- data.table(metadata)
samples_metadata(mofa) <- metadata

plot_variance_explained(mofa, plot_total = T)[[2]]

plot_variance_explained(mofa, max_r2=5)

#pdf("fig/MOFA_immune_hi.pdf")
correlate_factors_with_covariates(mofa,
                     covariates = c(immune, "hi", "BMI"),
                     plot="log_pval"
                     )
#dev.off()

plot_factor(mofa,
               factors = 2,
               color_by = "Factor2"
               )

plot_factors(mofa,
                     factors = c(1,5),
             color_by = "Co_type",
                     dot_size = 4
                     ) #+ guides(fill="none")

plot_factor(mofa,
               factor = 2,
               color_by = "Co_type",
               dot_size = 4,
               dodge = TRUE,
               stroke = 0.4,
               add_violin = T,
               add_boxplot = T
               ) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
              )

plot_factor(mofa,
               factor = 1,
               color_by = "Co_infb",
               dot_size = 4,
               dodge = TRUE,
               stroke = 0.4,
               add_violin = T,
               add_boxplot = T
               ) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
              )


plot_factor(mofa,
               factors = 2,
               color_by = "Sex",
               dodge = TRUE,
               add_violin = TRUE
               )


plot_data_heatmap(mofa,
                        view = "Eukaryote",
                        factor = 5,
                  features = 25,
                  denoise=TRUE,
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        show_rownames = TRUE, show_colnames = FALSE,
                        scale = "row"
                  )

plot_data_heatmap(mofa,
                        view = "Prokaryote",
                        factor = 5,
                  features = 25,
                  denoise=TRUE,
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        show_rownames = TRUE, show_colnames = FALSE,
                        scale = "row"
                        )




plot_top_weights(mofa,
                       view = "Eukaryote",
                       factor = 5,
                       nfeatures = 25,     # Top number of features to highlight
                       scale = T           # Scale weights from -1 to 1
                       )


plot_top_weights(mofa,
                       view = "Prokaryote",
                       factor = 5,
                       nfeatures = 25,     # Top number of features to highlight
                       scale = T           # Scale weights from -1 to 1
                       )

plot_data_heatmap(mofa,
                        factor = 5,
                        view = "Eukaryote",
                        features = 20,
                        denoise = TRUE,
                        cluster_rows = T, cluster_cols = F,
                        show_colnames = F, show_rownames = T,
                        annotation_samples = "Co_type",
 #                       annotation_colors = list("Category"=category.colors),
                        annotation_legend = F,
                        scale = "row"
                        )

plot_data_heatmap(mofa,
                        factor = 5,
                        view = "Prokaryote",
                        features = 20,
                        denoise = TRUE,
                        cluster_rows = T, cluster_cols = F,
                        show_colnames = F, show_rownames = T,
                        annotation_samples = "Co_type",
 #                       annotation_colors = list("Category"=category.colors),
                        annotation_legend = F,
                        scale = "row"
                        )
