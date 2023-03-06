library(phyloseq)
### We shouldn't use spieceasi because of the internal transformations, let's do correlations instead

library(igraph)


library(psych)

PS.T <- readRDS("tmp/PS.T.rds")


# let's subset based on co-infection type
UnI <- subset_samples(PS.T, Co_infb==0)
SiI <- subset_samples(PS.T, Co_infb==1)
MuI <- subset_samples(PS.T, Co_infb==2)

## We will start by doing a network for uninfected mice

Uotu<- UnI@otu_table
Udt <- UnI@sam_data

library(parallel)
library(doParallel)

options('mc.cores' = 4)
registerDoParallel(5)

library(HiClimR)

U.cor <- corr.test(Uotu, use="pairwise", method="spearman", adjust="fdr", alpha=0.05)
saveRDS(U.cor, "tmp/Ucor.rds")




colnames(Uotu)

library(NetCoMi)

UnInf_cor <- netConstruct(UnI,
                          measure="spearman",
                          sparsMethod="bootstrap",
                          seed=1234,
                          verbose=2)
