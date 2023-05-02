library(phyloseq)
library(Hmisc)
library(Matrix)
library(igraph)

### We shouldn't use spieceasi because of the internal transformations, let's do correlations instead

PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")

# let's subset based on co-infection type

UnI <- subset_samples(PS.TSS, Co_infb==0)
SiI <- subset_samples(PS.TSS, Co_infb==1)
MuI <- subset_samples(PS.TSS, Co_infb==2)

UnI <- prune_taxa(taxa_sums(UnI) > 0, UnI)
SiI <- prune_taxa(taxa_sums(SiI) > 0, SiI)
MuI <- prune_taxa(taxa_sums(MuI) > 0, MuI) 

#################################################################
## We will start by doing a network for uninfected mice
net.spear <- function(UnI){
########################## based on spearman ###########
U.cor <- rcorr(as.matrix(UnI@otu_table), type="spearman")
U.pval <- forceSymmetric(U.cor$P) # Self-correlation as NA
U.p <- p.adjust(U.pval, method="BH")# adjust for multiple testing
U.pval@x <- U.p
p.yes <- U.pval<0.01
r.val = U.cor$r # select all the correlation values
p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion
## select asv based on rho
yes.r <- abs(p.yes.r)>0.2 # output is logical vector
p.yes.rr <- p.yes.r*yes.r # use logical vector for subscripting.
adjm <- as.matrix(p.yes.rr)
all(taxa_names(UnI)==colnames(adjm))
adjm[is.na(adjm)] <- 0
net.grph <- graph_from_adjacency_matrix(adjm, mode="undirected", weighted=T)
edgew<-E(net.grph)$weight
E(net.grph)$weight <- abs(E(net.grph)$weight)
V(net.grph)$Genus <- as.vector(UnI@tax_table[,6])
V(net.grph)$Domain <- as.vector(UnI@tax_table[,1])
V(net.grph)$Phylum <- as.vector(UnI@tax_table[,2])
bad.vs<-V(net.grph)[degree(net.grph) == 0] 
net.grph <- delete.vertices(net.grph, bad.vs)

# set edge color postive correlation pink color, negative blue.
E.color.Uni = edgew
E.color.Uni = ifelse(E.color.Uni>0, "pink",ifelse(E.color.Uni<0, "blue","grey"))
E(net.grph)$color = as.character(E.color.Uni)
#change edge width
    E(net.grph)$width = abs(edgew)
    return(net.grph)
}



U.grph <- net.spear(UnI)
S.grph <- net.spear(SiI)
M.grph <- net.spear(MuI)

Parasite <- c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas")

V(S.grph)$Stype <- "shape"
V(S.grph)$Stype[which(V(S.grph)$Domain=="Archaea")] <- "circle"
V(S.grph)$Stype[which(V(S.grph)$Domain=="Bacteria")] <- "circle"
V(S.grph)$Stype[which(V(S.grph)$Domain=="Eukarya")] <- "square"
V(S.grph)$Stype[which(V(S.grph)$Domain=="Unknown_domain")] <- "sphere"
V(S.grph)$Stype[which(V(S.grph)$Genus%in%Parasite)] <- "square"

V(S.grph)$color <- "black"
V(S.grph)$color[which(V(S.grph)$Domain=="Archaea")] <- "deeppink3"
V(S.grph)$color[which(V(S.grph)$Domain=="Bacteria")] <- "deeppink4"
V(S.grph)$color[which(V(S.grph)$Domain=="Eukarya")] <- "darkolivegreen"
V(S.grph)$color[which(V(S.grph)$Domain=="Unknown_domain")] <- "darkorange"
V(S.grph)$color[which(V(S.grph)$Genus%in%Parasite)] <- "darkolivegreen1"

col <- V(S.grph)$color


degS <- degree(S.grph)
plot(U.grph,
     vertex.size=degS/8,
     vertex.shape=V(S.grph)$Stype,
     vertex.frame.color=col,
     vertex.color=col,
     vertex.label="")

V(U.grph)$Stype <- "shape"
V(U.grph)$Stype[which(V(U.grph)$Domain=="Archaea")] <- "circle"
V(U.grph)$Stype[which(V(U.grph)$Domain=="Bacteria")] <- "circle"
V(U.grph)$Stype[which(V(U.grph)$Domain=="Eukarya")] <- "square"
V(U.grph)$Stype[which(V(U.grph)$Domain=="Unknown_domain")] <- "sphere"
V(U.grph)$Stype[which(V(U.grph)$Genus%in%Parasite)] <- "square"

degU <- degree(U.grph)
plot(U.grph,
     vertex.size=degU/8,
     vertex.shape=V(U.grph)$Stype,
     vertex.label="")

V(M.grph)$Stype <- "shape"
V(M.grph)$Stype[which(V(M.grph)$Domain=="Archaea")] <- "circle"
V(M.grph)$Stype[which(V(M.grph)$Domain=="Bacteria")] <- "circle"
V(M.grph)$Stype[which(V(M.grph)$Domain=="Eukarya")] <- "square"
V(M.grph)$Stype[which(V(M.grph)$Domain=="Unknown_domain")] <- "sphere"
V(M.grph)$Stype[which(V(M.grph)$Genus%in%Parasite)] <- "square"
V(M.grph)$color <- "black"
V(M.grph)$color[which(V(S.grph)$Domain=="Archaea")] <- "deeppink3"
V(M.grph)$color[which(V(S.grph)$Domain=="Bacteria")] <- "deeppink4"
V(M.grph)$color[which(V(S.grph)$Domain=="Eukarya")] <- "darkolivegreen"
V(M.grph)$color[which(V(S.grph)$Domain=="Unknown_domain")] <- "darkorange"
V(M.grph)$color[which(V(S.grph)$Genus%in%Parasite)] <- "darkolivegreen1"

degM <- degree(M.grph)
plot(M.grph,
     vertex.size=degM/8,
     vertex.frame.color=col,
     vertex.color=col,
     vertex.shape=V(M.grph)$Stype,
     vertex.label="")

##### betweeness

################# modularity
Umod =cluster_fast_greedy(U.grph)

V(U.grph)$cluster <- Umod$membership

Smod =cluster_fast_greedy(S.grph, weights=E(S.grph)$weight)
S.grph$cluster <- Smod$membership
Mmod =cluster_fast_greedy(M.grph, weights=E(M.grph)$weight)
M.grph$cluster <- Mmod$membership

modularity(Umod)
modularity(Smod)
modularity(Mmod)

sizes(Umod)
sizes(Smod)
sizes(Mmod)

