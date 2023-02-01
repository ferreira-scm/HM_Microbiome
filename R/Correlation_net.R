amp_func <- function(PS.lT, parasite){
nmOxy <- list()
amp <- list()
for (i in 1:length(PS.lT)) {
    try(p <- subset_taxa(PS.lT[[i]],Genus==parasite), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Genus")
        print(paste(i, "- ", names(PS.l[i]), ": ", nrow(p@tax_table), sep=""))
        nmOxy[i] <- names(PS.l[i])
        amp[[i]] <- paste(rep(names(PS.l[i]), nrow(p@tax_table)), "ASV", seq(1, nrow(p@tax_table), 1), sep="_")
    }
    rm(p)
}

nmOxy <- unlist(nmOxy)
amp <- unlist(amp)
}


Correlation_net <- function(PS.lT, PS.l, PS.T, parasite){

print(parasite)
    
#how many primers amplify Oxyurida and which families?
for (i in 1:length(PS.lT)) {
#    print(names(all.PS.l)[[i]])
    try(p <- subset_taxa(PS.lT[[i]],Genus==parasite), silent=TRUE)
#    try(get_taxa_unique(p, "family"), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Species")
        print(paste(i, "- ", names(PS.l[i]), ": ", length(a), " species", sep=""))
        print(a)
    }
    rm(p)
}

nmOxy <- list()
amp <- list()
for (i in 1:length(PS.lT)) {
    try(p <- subset_taxa(PS.lT[[i]],Genus==parasite), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Genus")
        print(paste(i, "- ", names(PS.l[i]), ": ", nrow(p@tax_table), sep=""))
        nmOxy[i] <- names(PS.l[i])
        amp[[i]] <- paste(rep(names(PS.l[i]), nrow(p@tax_table)), "ASV", seq(1, nrow(p@tax_table), 1), sep="_")
    }
    rm(p)
}

nmOxy <- unlist(nmOxy)
amp <- unlist(amp)


    print(amp)
    
############## co-occurrences network
##### Oxy, let's try and disentangle the species here
library(Hmisc)
library(Matrix)
library(igraph)
library(RColorBrewer)

Oxy <-subset_taxa(PS.T, Genus %in%parasite)
Oxy <- prune_samples(sample_sums(Oxy)>0, Oxy)
oxy <- (Oxy@otu_table)
tax <- data.frame(Oxy@tax_table)

otu.cor <- rcorr(as.matrix(oxy), type="spearman")
otu.pval <- forceSymmetric(otu.cor$P)
cor.p <- p.adjust(otu.pval, method="BH")

# consider only significant and strong correlations
#cor.r1[which(!cor.r1 > 0.8 | !cor.r1 < -0.8)]=NA
#cor.p[which(cor.p>0.01)]=NA
otu.pval@x <- cor.p
sel.tax <- tax[rownames(otu.pval),,drop=FALSE]

#sanity check
all.equal(rownames(sel.tax), rownames(otu.pval))
p.yes <- otu.pval<0.05
r.val = otu.cor$r # select all the correlation values
p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion
## select asv based on rho
p.yes.r.str <- abs(p.yes.r)>0.5 # output is logical vector
p.yes.rr <- p.yes.r.str*r.val # use logical vector for subscripting.
adjm <- as.matrix(p.yes.r)

colnames(adjm) <- Oxy@tax_table[,7]
rownames(adjm) <- Oxy@tax_table[,7]

#we also want the node color to code for amplicon
amp_name <- as.factor(gsub("_ASV_[0-9]", "", amp))
#amp_name <- amp_name[-grep("D3A", amp_name)]

nb.col <- length(levels(amp_name))
coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)
mc <- coul[as.numeric(amp_name)]


net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)

V(net.grph)

### negative correlations
summary(adjm<0)
### colour negative edges
E(net.grph)$weight
E(net.grph)$color <- "dodgerblue4"
E(net.grph)$color[which(E(net.grph)$weight<0)] <- "#FF00FF"
E(net.grph)$color


E(net.grph)$weight <- abs(E(net.grph)$weight)


set.seed(1234)
plot(net.grph,
     vertex.label=amp,
     edge.width=E(net.grph)$weight*4,
     vertex.color=adjustcolor(mc, 0.8),
     frame.col="grey")
}
