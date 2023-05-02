library(Hmsc)
library(snow)
library(corrplot)
library(metagMisc)
library(sjSDM)

#install_sjSDM()


PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")
PS.TSS <- phyloseq_standardize_otu_abundance(PS.TSS, method="pa")


################ todo:
####  I might have to z-transform these variables.


# Environment features (host is the environment)
metadata <- as.data.frame(PS.TSS@sam_data)
class(metadata) <- "data.frame"

X=metadata[,c("hi", "Sex", "BMI")]
XFormula <- ~hi+Sex+BMI

# Species occurrence box
Y <- as.matrix(PS.TSS@otu_table)
class(Y) <- "matrix"

### this is annoying, but apparently we need :
#All spatial locations should be unique. If you have several observations in the same point, they should be identified by the random levels.'
lev <- unique(metadata$Longitude)
for (i in 1:length(lev)){
#length(grep(lev[2], metadata$Longitude))
    metadata$Longitude[which(metadata$Longitude%in%lev[i])] <- paste(metadata$Longitude[which(metadata$Longitude==lev[i])], seq(1:length(grep(lev[i], metadata$Longitude))), sep="")
}

lev <- unique(metadata$Longitude)
for (i in 1:length(lev)){
#length(grep(lev[2], metadata$Longitude))
    metadata$Longitude[which(metadata$Longitude%in%lev[i])] <- paste(metadata$Longitude[which(metadata$Longitude==lev[i])], seq(1:length(grep(lev[i], metadata$Longitude))), sep="")
}

lev <- unique(metadata$Latitude)
for (i in 1:length(lev)){
#length(grep(lev[2], metadata$Longitude))
    metadata$Latitude[which(metadata$Latitude%in%lev[i])] <- paste(metadata$Latitude[which(metadata$Latitude==lev[i])], seq(1:length(grep(lev[i], metadata$Latitude))), sep="")
}

lev <- unique(metadata$Latitude)
for (i in 1:length(lev)){
#length(grep(lev[2], metadata$Longitude))
    metadata$Latitude[which(metadata$Latitude%in%lev[i])] <- paste(metadata$Latitude[which(metadata$Latitude==lev[i])], seq(1:length(grep(lev[i], metadata$Latitude))), sep="")
}

metadata$Longitude <- as.numeric(metadata$Longitude)
metadata$Latitude <- as.numeric(metadata$Latitude)

metadata$Locality <- paste(metadata$Longitude, metadata$Latitude, sep=" ")

xy <- as.matrix(cbind(metadata$Longitude, metadata$Latitude))
colnames(xy) <- c("Latitude", "Longitude")
rownames(xy) <- metadata$Locality

# study design
studyDesign <- metadata[,c("Year", "Locality")]
studyDesign$Year <- as.factor(studyDesign$Year)
studyDesign$Locality <- as.factor(studyDesign$Locality)

rL1=HmscRandomLevel(sData=xy)

rL2=HmscRandomLevel(units=studyDesign$Year)

ranlevels=list(Locality=rL1, Year=rL2)
ranlevels

jModel <- Hmsc(Y=Y, XData=X,
               XFormula=XFormula,
               studyDesign=studyDesign,
               ranLevels=ranlevels,
               distr="probit")

## this takes too long... :(
#thin=50
#samples=1000
#nChains=4
#transient=5000
#
#mod_HMDC <- sampleMcmc(jModel,
#                       samples=samples,
#                       thin=thin,
#                       transient=transient,
#                       nChains=nChains,
#                       nParallel=nChains)
#saveRDS(mod_HMDC, "tmp/JSDM_spatial.R")
#mcoda <- convertToCodaObject(mod_HMDC)
#par(mar=rep(2,4))
#plot(mcoda$Beta[,1:5])
#plot(mcoda$Gamma)
#gelman.diag(mcoda$Beta[,1:50])
#postBeta=getPostEstimate(mod_HMDC, parName="Beta")
#plotBeta(mod_HMDC,
#         post=postBeta,
#         plotTree=F,
#         spNamesNumbers=c(T,F))
#plotBeta(mod_HMDC,
#         post=postBeta,
#         param="Mean",
#         plotTree=F,
#         spNamesNumbers=c(T,F))
#VP=computeVariancePartitioning(mod_HMDC)
#plotVariancePartitioning(mod_HMDC, VP=VP, las=2, horiz=F)


############# Let's try a faster approach

#### to do: let's tune this!!!! 40 random steps with leave-one-out cross validation(LOOCV), 150 iterations and learning rate of 0.01

#sjSDM_cv()

colnames(xy) <- c("XX", "YY")

SPeigen <- generateSpatialEV(xy)

#model1=sjSDM(Y=Y, env=linear(X, ~.), spatial=linear(SPeigen, ~0+.), iter=100L, se=TRUE)

saveRDS(model1, "tmp/sjSDM_model1.rds")

Rsquared(model1)

#getCov(model1)
    
summary(model1)

plot(model1)

#an <- anova(model1)
#saveRDS(an, "tmp/sjSDM_anova.rds")

an <- readRDS("tmp/sjSDM_anova.rds")

print(an)

plot(an)

plot(an, internal=TRUE)


imp <- importance(model1)

plot(imp)


                                        # species-species association matrix
association <- getCor(model1)

#association[abs(association)<0.2] <- NA

plot(association)

adj <- as.matrix(association)

heatmap(adj)

# only strong associations
adj[abs(adj)<0.8] <- 0

library(igraph)

net <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=TRUE, diag=FALSE)

net

bad <- V(net)[degree(net) == 0] 

net <- delete.vertices(net, bad) 

plot(net,
     layout=layout.circle,
     vertex.size=1)

plot(net)
