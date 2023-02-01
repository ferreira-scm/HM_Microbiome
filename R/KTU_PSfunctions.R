library(KTU)
library(Biostrings)
library(dplyr)


#library(devtools)
#install_github("poyuliu/KTU")
## testing KTU in Eimeria only.

#library(KTU)
#library(Biostrings)
#library(dplyr)

# Extract abundance matrix from the phyloseq object
otu1 <- as(otu_table(fPS), "matrix")
# transpose if necessary

if(!taxa_are_rows(fPS)){otu1 <- t(otu1)}
  # Coerce to data.frame
otu1 <- as.data.frame(otu1)

#feature table (ie OTU table) from phyloseq
ft <- tibble::rownames_to_column(otu1, "asv")

# run clustering
hh <- DNAStringSet(dimnames(otu1)[[1]])
names(hh) <- ft$asv
seqs <- hh
Biostrings::writeXStringSet(seqs, format="fasta", filepath="tmp/testKTU.fasta", width=max(seqs@ranges@width))


ktu_O <- KTU::klustering(repseq = "tmp/testKTU.fasta",
                  feature.table = ft,
                  write.fasta = TRUE,
                  cores = 90)

# merge taxa based on KTU (cluster of ASV)
KTU <- merge_taxa(fPS, taxa_names(fPS)[which(ktu_O$clusters==1)])
for (i in 2:max(ktu_O$clusters)){
    KTU <- try(merge_taxa(KTU, taxa_names(fPS)[which(ktu_O$clusters==i)]), silent=TRUE)
}

KTU

dada2KTU(data=fPS, subset=100, path2fasta="tmp/test.fna", ncores=10)

