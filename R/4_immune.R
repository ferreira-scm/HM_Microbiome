library("MCMCglmm")
library(ape)
library(brms)
library(rstan)
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(ggmcmc)
library(ggthemes)
library(ggridges)
library(vegan)

PS.T <- readRDS("tmp/PS.T.rds")

### questions to answer:
# How does immne gene expression interact with host genetics and parasite/microbiome structre

### Approach 1

# subsetting
PS.i <- subset_samples(PS.T, !is.na(PS.T@sam_data$TNF))

# now saving ID names
key <- data.frame(ID=sample_data(PS.i)$Mouse_ID)
# metadata
metadt <- sample_data(PS.i)

# first we need to remove parasites from the microbiome PS object because we will analyse those in separate
PS.iP <- subset_taxa(PS.i, !Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(PS.iP, method="jaccard", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
dimnames(JACM)<- c(key, key)

# distance matrix, we need to control for distance?
distance.df <- metadt[,c("Mouse_ID", "Longitude", "Latitude")]
SPATM <- array(NA, c(length(distance.df$Mouse_ID),length(distance.df$Mouse_ID)))
# derive matrix with spatial distances between each location
for (i in 1:length(distance.df$Mouse_ID)){
    for (j in 1:length(distance.df$Mouse_ID))
    {SPATM[i,j]= sqrt((abs(distance.df$Longitude[i]-distance.df$Longitude[j]))^2+
                      (abs(distance.df$Latitude[i]-distance.df$Latitude[j]))^2)
    }
}
dimnames(SPATM)<- c(key, key)

# 2) pairwise genetic distance based on genetic data
# I actually need to get this from sota again
sota <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv")
sota <- sota[sota$Mouse_ID%in%distance.df$Mouse_ID,c(1,6:25)]
sota <- sota[match(distance.df$Mouse_ID, sota$Mouse_ID),]
all(sota$Mouse_ID==distance.df$Mouse_ID)
rownames(sota) <- sota$Mouse_ID
sota$Mouse_ID <- NULL
gen.dis <- dist.gene(sota, method = "pairwise", pairwise.deletion = TRUE)
gen.dis <- as.matrix(gen.dis)
dimnames(gen.dis) <- c(key, key)

# 3) Making immune gene ex distances
#Create data frame with each sample name (character) and sampling time (numeric)
immune <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")
IMM <- as.matrix(vegdist(metadt[,immune], method="euclidean"))
# transpose disssimilary matrix to similarty matrix
IMM <- 1-IMM
# sanity check
all(rownames(IMM)==key)
dimnames(IMM)<- c(key, key)

# 4) Making BMI distances
#Create data frame with each sample name (character) and sampling time (numeric)
BMI_frame<-metadt[,c("Mouse_ID", "BMI")]
#Create an empty matrix to fill with distances
BMIM<-array(0,c(nrow(BMI_frame),nrow(BMI_frame)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(BMI_frame)){
    for (j in 1:nrow(BMI_frame))
    {BMIM[i,j]=abs(BMI_frame$BMI[i] -BMI_frame$BMI[j])
    }
}
dimnames(BMIM) <- c(key, key)

# 5) Creating parasite matrix with chi-sqr dissimilary
Parasite <- subset_taxa(PS.i, Genus %in%c("Eimeria", "Cryptosporidium", "Tritrichomonas",
                                            "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus",
                                            "Trichuris", "Hymenolepis"))
all(sample_names(Parasite)==key)
PARA <- as.matrix(vegdist(Parasite@otu_table, method="chisq"))
# transpose disssimilary matrix to similarty matrix
PARA <- 1-PARA
# sanity check
all(rownames(PARA)==key)
dimnames(PARA)<- c(key, key)

# 6) Create farm/Locality matrix: 
#Create data frame with each Individual name (character) and their Age (Character)
Loc_frame<-metadt[,c("Mouse_ID","Locality")]
#Create an empty numeric matrix to fill with distances
LocM<-array(0,c(nrow(Loc_frame),nrow(Loc_frame)))
#Derive matrix with binary locality similarity between each sample
for(i in 1:nrow(Loc_frame)){
    for(j in 1:nrow(Loc_frame)){
        if(Loc_frame$Locality[i]==Loc_frame$Locality[j]){
            LocM[i,j]= 1
        } else{
            LocM[i,j]= 0
        }
    }
}
#Note that farm similarity matrix has rownames and colnames in the same order as key
all(rownames(LocM)==key$ID)
dimnames(LocM) <- c(key, key)

# 7) Making HI distances
#Create data frame with each sample name (character) and sampling time (numeric)
HI_frame<-metadt[,c("Mouse_ID", "HI")]
#Create an empty matrix to fill with distances
HIM<-array(0,c(nrow(HI_frame),nrow(HI_frame)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(HI_frame)){
    for (j in 1:nrow(HI_frame))
    {HIM[i,j]=abs(HI_frame$HI[i] -HI_frame$HI[j])
    }
}
dimnames(HIM) <- c(key, key)

# 7) Making hibridicity distances
#Create data frame with each sample name (character) and sampling time (numeric)
hi_frame<-metadt[,c("Mouse_ID", "hi")]
#Create an empty matrix to fill with distances
hiM<-array(0,c(nrow(hi_frame),nrow(hi_frame)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(hi_frame)){
    for (j in 1:nrow(hi_frame))
    {hiM[i,j]=abs(hi_frame$hi[i] -hi_frame$hi[j])
    }
}
dimnames(hiM) <- c(key, key)


jac<-c(as.dist(JACM))
bmi<-c(as.dist(BMIM))
spa<-c(as.dist(SPATM))
imm<-c(as.dist(IMM))
gen<-c(as.dist(gen.dis))
para<-c(as.dist(PARA))
loc <- c(as.dist(LocM))
HIm <- c(as.dist(HIM))
him <- c(as.dist(hiM))


#Combine these vectors into a data frame
data.dyad<-data.frame(Immune=imm, Parasite=para, BMI=bmi,Microbiome_similarity=jac,spatial=spa,
                      genetic_dist=gen, locality=loc, HI_dist=HIm, hi_dist=him)

#Now all we need to do is add the identities of both individuals in each dyad as separate columns into the data frame and exclude self-comparisons (as these are not meaningful).

# extracting Individual-combinations present in the matrices
list<-expand.grid(key$ID, key$ID)

str(list)

# This created individual-to-same-individual pairs as well. Get rid of these:
list<-list[which(list$Var1!=list$Var2),]

# this still has both quantiles in--> add 'unique' key
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))

# sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices
i=nrow(key)
JACM[which(rownames(JACM)==list$Var1[i]),which(colnames(JACM)==list$Var2[i])]==jac[i]

# add the names of both individuals participating in each dyad into the data frame
data.dyad$IDA<-list$Var2
data.dyad$IDB<-list$Var1
# Make sure you have got rid of all self comparisons
data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]

######################### Now we model the data ####################
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","genetic_dist", "Parasite", "BMI", "Immune")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }


preplot1<-ggplot(data = data.dyad, aes(x= Microbiome_similarity, y= Immune))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .8)+
    theme_bw()
preplot1

preplot2<-ggplot(data = data.dyad, aes(x= Parasite, y= Immune))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .8)+
    theme_bw()
preplot2


#### model
modeli<-brm(Microbiome_similarity~1+ spatial+locality+BMI+genetic_dist+Parasite+ Immune*hi_dist+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 50, chains = 10,
                inits=0)
saveRDS(modeli, "tmp/BRMmodeli.rds")
modeli <- readRDS("tmp/BRMmodeli.rds")

conditional_effects(modeli)


modelimm<-brm(Microbiome_similarity~1+ spatial+locality+BMI+genetic_dist+Parasite+ Immune+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 50, chains = 10,
                inits=0)
saveRDS(modelimm, "tmp/BRMmodelimm.rds")
modelimm <- readRDS("tmp/BRMmodelimm.rds")

modelimm

#model_i<-brm(Immune~1+ BMI+genetic_dist+Parasite+Microbiome_similarity+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 50, chains = 10,
#                inits=0)
#saveRDS(model_i, "tmp/BRMmodel_i.rds")
model_i <- readRDS("tmp/BRMmodel_i.rds")

#model_ii<-brm(Immune~1+ BMI+genetic_dist+Parasite+Microbiome_similarity+locality+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 50, chains = 10,
#                inits=0)
#saveRDS(model_ii, "tmp/BRMmodel_ii.rds")
model_ii <- readRDS("tmp/BRMmodel_ii.rds")

modeli

model_ii

resdf1<-summary(modeli)$fixed
resdf1<-as.data.frame(resdf1)
resdf1<-resdf1[c("Estimate","l-95% CI","u-95% CI")]
resdf1<-resdf1[2:nrow(resdf1),]
resdf1$Predictor<-rownames(resdf1)
colnames(resdf1)<-c("Estimate","lCI","uCI","Predictor")
resdf1$Predictor<-factor(resdf1$Predictor)

plot2<-ggplot(resdf1,aes(x=Estimate,y=Predictor,colour=Predictor))+
    geom_linerange(aes(xmin = lCI, xmax = uCI),size=2.5)+
    geom_point(size=4,colour="black", shape=21, fill="white")+
    theme_bw()+
    theme(legend.position='none',text = element_text(size=22))+
    labs(x="Effect on immune structure",y="")+
    geom_vline(xintercept=0, linetype="dashed")

plot2

#################### Eimeria immune

heatmap(cor(metadt[,immune]))

cor(metadt[,immune])


#### testing some stuff



## testing now the Anna Karenina principle: each dysbiotic community is dysbiotic on its own way
PS.mP@sam_data$Co_infb <- as.factor(PS.mP@sam_data$Co_infb)

immuneT <- metadt[,immune]

library(cluster)

model=kmeans(immuneT, 3)

clusplot(immuneT, model$cluster)

dis <- vegdist(immuneT, method="euclidean")

dis2 <- stepacross(dis, path="extended")

I.ord <- metaMDS(immuneT, "euclidean", trymax = 500)

I.rda <- rda(immuneT, data=metadt, scaling="species")

names(metadt)

parasite <- Parasite@otu_table

colnames(parasite) <- Parasite@tax_table[,7]




parasite

parasite_euc <- rda(parasite)

parasite_euc$CA

ordiplot(parasite_euc, display=c("sites", "species"))

orditorp(parasite_euc, display="species")

P.rda <- rda(immuneT, data=metadt, scaling="species")

metadt$Immune_PCA1 <- I.rda$CA$u[,1]
metadt$Immune_PCA2 <- I.rda$CA$u[,2]
metadt$Immune_PCA3 <- I.rda$CA$u[,3]
metadt$para_PCA1 <- parasite_euc$CA$u[,1]
metadt$para_PCA2 <- parasite_euc$CA$u[,2]
metadt$para_PCA3 <- parasite_euc$CA$u[,3]

model_Eimf<-brm(Eimeria_ferrisi_asv~1+ BMI+HI+Immune_PCA1+Immune_PCA2+Immune_PCA3+
                (1|Locality),
                data = metadt,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 50, chains = 4,
                inits=0)
saveRDS(model_Eimf, "tmp/BRMmodel_Eimf.rds")
#
model_Para1<-brm(para_PCA1~1+ BMI+HI+Immune_PCA1+Immune_PCA2+Immune_PCA3+
                (1|Locality),
                data = metadt,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 50, chains = 4,
                inits=0)
saveRDS(model_Para1, "tmp/BRMmodel_Para1.rds")
#
model_Para2<-brm(para_PCA2~1+ BMI+HI+Immune_PCA1+Immune_PCA2+Immune_PCA3+
                (1|Locality),
                data = metadt,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 50, chains = 4,
                inits=0)
saveRDS(model_Para2, "tmp/BRMmodel_Para2.rds")


model_Eimf

model_Para2

head(I.rda)

I.cca <- cca(immuneT~Co_infb, data=metadt)

fit <- envfit(I.cca~metadt$Co_infb, perm=999, display="lc")

plot(fit, p.max=0.05)

summary(I.rda)

plot2 <- ordiplot(I.rda, choices=c(1,2))

plot2

ordiellipse(plot2, groups=Co_infb, display="sites")

biplot(I.rda,
       display=c("sites"))

co.names <- levels(metadt$Co_infb)

ordihull(I.rda,
         group = metadt$Co_infb,
         col = c(1,2,3))

legend("topright",
       col = c(1,2,3),
       lty = 1,
              legend = co.names)

class(metadt) <- "data.frame"

groups <- metadt$Co_infb
mod <- betadisper(dis, groups)

anova(mod) # the dispersion is not different between groups

plot(mod, hull=FALSE, ellipse=TRUE)

boxplot(mod) # the dispersion is not different between groups

mod.HSD <- TukeyHSD(mod)

plot(mod.HSD)# the dispersion is not different between groups


sppscores(I.ord) <- immuneT

plot(I.ord)

p1=plot_ordination(PS.mP, PS.ord, type="sample", color="Co_infb")

p1+
    scale_color_manual(labels=c("uninfected", "single parasite", "multiple parasites"), values = c("wheat", "tomato", "tomato4"))+
    labs(color="Parasite infection")+
    theme_bw()+
    stat_ellipse(aes(group=PS.mP@sam_data$Co_infb)) 

