library(vegan)
library(microbiome)
library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

#                                        BiocManager::install("MOFA2", force=TRUE)

library(MOFA2)


# MOFA+: analysis of multi-modal microbiome data
# remove host, human and plants

PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")

PS.TSS@otu_table <- decostand(PS.TSS@otu_table, method="clr", pseudocount=1)

tax <- as.data.frame(PS.TSS@tax_table)
tax$Kingdom[tax$Phylum%in%c("Mucoromycota", "Ascomycota", "Basidiomycota")] <- "Fungi"
tax$Kingdom[tax$Phylum%in%c("Anthophyta", "Phragmoplastophyta", "Charophyta", "Ochrophyta")] <- "Plants"
tax$Kingdom[tax$Genus%in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas")] <- "Parasite"
tax$Kingdom[tax$Kingdom%in%"Bacteria"] <- "Bacteria"
tax$Kingdom[!tax$Kingdom%in%c("Bacteria", "Parasite", "Plants", "Fungi")] <- "Other"
tax_table(PS.TSS) <- tax_table(as.matrix(tax))
taxa_names(PS.TSS) <- paste("ASV", 1:ntaxa(PS.TSS), PS.TSS@tax_table[,6], sep="_")

get_taxa_unique(PS.TSS, "Kingdom")
df <- psmelt(PS.TSS)

df1 <- df[,c("Sample", "OTU", "Abundance", "Kingdom")]
names(df1) <- c("sample", "feature", "value", "view")
#df1$feature

df1 <- data.table(df1)

df1[,length(unique(feature)),by="view"]

#df1$OTU <- NULL

ggplot(df1, aes(x=value, colour=view))+
    geom_density()
#    facet_wrap(~Kingdom, nrow=1, scales="free")

metadata <- sample_data(PS.TSS)

keep <- c("Mouse_ID", "Sex", "Locality", "Year", "hi", "HI", "BMI")

metadata <- metadata[,keep]

metadata$sample <- metadata$Mouse_ID

metadata$Mouse_ID <- NULL

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
                     covariates = c("Locality", "hi", "BMI", "Sex", "Co_infb", "Year"),
                     plot="log_pval"
                     )
#dev.off()

plot_factor(mofa,
               factors = 1,
               color_by = "Factor1"
               )

plot_factors(mofa,
                     factors = c(1,2),
             color_by = "hi",
                     dot_size = 4
                     ) #+ guides(fill="none")

plot_factor(mofa,
               factor = 1,
               color_by = "hi",
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

plot_data_heatmap(mofa,
                        view = "Fungi",
                        factor = 1,
                  features = 25,
                  denoise=TRUE,
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        show_rownames = TRUE, show_colnames = FALSE,
                        scale = "row"
                  )

plot_data_heatmap(mofa,
                        view = "Parasite",
                        factor = 1,
                  features = 25,
                  denoise=TRUE,
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        show_rownames = TRUE, show_colnames = FALSE,
                        scale = "row"
                        )




plot_top_weights(mofa,
                       view = "Fungi",
                       factor = 2,
                       nfeatures = 25,     # Top number of features to highlight
                       scale = T           # Scale weights from -1 to 1
                       )


plot_data_heatmap(mofa,
                        factor = 2,
                        view = "Fungi",
                        features = 20,
                        denoise = TRUE,
                        cluster_rows = T, cluster_cols = F,
                        show_colnames = F, show_rownames = T,
                        annotation_samples = "HI",
 #                       annotation_colors = list("Category"=category.colors),
                        annotation_legend = F,
                        scale = "row"
                        )

