source("R/2_Parasite_cleaning.R")

library(SpiecEasi)

library(microbiome)

PS.T

# spiec easi takes no normalised taxa
fPS

# remove host reads
# remove diet reads

rank_names(fPS)

get_taxa_unique(fPS, "Phylum")

fPS.b <- subset_taxa(fPS, !Phylum%in%c("Anthophyta", "Vertebrata", "Ochrophyta","Phragmoplastophyta"))

fPS_minusP <- subset_taxa(fPS.b, !Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

fParasite <- subset_taxa(fPS, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

KeepTaxap <- microbiome::prevalence(fPS_minusP)>0.10

PS10 <- phyloseq::prune_taxa(KeepTaxap, fPS_minusP)

PS10 <- merge_phyloseq(PS10, fParasite)

pargs <- list(rep.num=1000, seed=10010, ncores=90, thresh=0.05)

#t1 <- Sys.time()
#se.PS10 <- spiec.easi(PS10, method="mb", pulsar.params=pargs)
#t2 <- Sys.time()
#t2-t1
#saveRDS(se.PS10, file="tmp/se.PS10.rds")

se.PS10 <- readRDS("tmp/se.PS10.rds")

bm=symBeta(getOptBeta(se.PS10), mode="maxabs")
diag(bm) <- 0


#weights <- Matrix::summary(t(bm))[,3] # includes negative weights
weights <- (1-Matrix::summary(t(bm))[,3])/2 # ort

net.ids <- taxa_names(PS10)

net10 <- adj2igraph(Matrix::drop0(getRefit(se.PS10)),
                    edge.attr=list(weight=weights),
                    vertex.attr = list(name=net.ids))

betaMat=as.matrix(bm)

# we want positive edges to be green and negative to be red
edges <- E(net10)
edge.colors=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(net10, edges[e.index])
    xindex=which(net.ids==adj.nodes[1])
    yindex=which(net.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "#1B7837")
    }else if(beta<0){
        edge.colors=append(edge.colors, "#762A83")
    }
}
E(net10)$color=edge.colors

### defining attributes
V(net10)$family=PS10@tax_table[,5]
V(net10)$species=PS10@tax_table[,7]
V(net10)$genus=PS10@tax_table[,6]

#V(net10)$type <- 
#which(V(net10)$species%in%fParasite@tax_table[,7])


PS10@tax_table[,7]

hub.s <- hub_score(net10)$vector

deg <- igraph::degree(net10, mode="all")

# we also want the node color to code for family
nb.col <- length(levels(as.factor(V(net10)$family)))
coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)



set.seed(1002)
plot(net10,
     vertex.label=V(net10)$genus,
     vertex.size=as.integer(cut(hub.s, breaks=10))+2,
     vertex.color=adjustcolor(coul,0.8),
     edge.width=as.integer(cut(E(net10)$weight, breaks=6))/3,
     margin=c(0,1,0,0))
legend(x=-2, y=1, legend=levels(as.factor(V(net10)$family)), col=coul, bty="n", x.intersp=0.25,text.width=0.045, pch=20, pt.cex=1.5)

which(V(net10)$species=="ferrisi")
which(V(net10)$species=="falciformis")

which(V(net10)$genus=="Syphacia")

which(V(net10)$genus=="Aspiculuris")

