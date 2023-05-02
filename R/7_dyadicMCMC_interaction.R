#### Testing if there are spatial, host genetics interactions
### we don't include sex because it does not converge
PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")

#### In this script we are going to model Microbiome beta diversity by parasite structure, geographic distance,
#genetic distance, sex combination and body mass index. We do this with dyadic interactions using Bayesian multimodal
#models.

## I follow Aura Raulo approach

### Questions to asnwer:
#1) genetic distance is a better predictor than hybrid index distance? This has
#implications in specific allele combinations.
#2) spatial distance: is it just farm membership, or actual physical distance whithin farms?
# No evidence for sex effect in early models, so we don't include this variable
#

############# First create dyad data#######################
key <- data.frame(ID=sample_data(PS.TSS)$Mouse_ID)
metadt <- sample_data(PS.TSS)


################################# We do not remove parasites
# first we need to remove parasites from the microbiome PS object because we will analyse those in separate
#PS.mP <- subset_taxa(PS.TSS, !Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

####################
## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(PS.mP, method="jaccard", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
dimnames(JACM)<- c(key, key)

## 2) Spatial distance matrix
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

# 3) pairwise genetic distance based on genetic data
# I actually need to get this from sota again
sota <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv")
gen <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
         "Es1","Gpd1","Idh1","Mpi","Np", "Sod1", "Es1C", "Gpd1C", "Idh1C",
         "MpiC","NpC", "Sod1C")
sota <- sota[sota$Mouse_ID%in%distance.df$Mouse_ID,c("Mouse_ID", gen)]
sota <- sota[match(distance.df$Mouse_ID, sota$Mouse_ID),]
all(sota$Mouse_ID==distance.df$Mouse_ID)
rownames(sota) <- sota$Mouse_ID
sota$Mouse_ID <- NULL
gen.dis <- dist.gene(sota, method = "pairwise", pairwise.deletion = TRUE)
gen.dis <- as.matrix(gen.dis)
dimnames(gen.dis) <- c(key, key)

## 3) No need for sex anymore
#Sex_frame<-metadt[,c("Mouse_ID","Sex")]
#Sex_frame$Mouse_ID<-as.character(Sex_frame$Mouse_ID)
#Sex_frame$Sex<-as.character(Sex_frame$Sex)
#Create an empty character matrix to fill with characters
#SEXM<-array(as.character(NA),c(nrow(Sex_frame),nrow(Sex_frame)))
#
#for(i in 1:nrow(Sex_frame)){
#    for(j in 1:nrow(Sex_frame)){
#        if(Sex_frame$Sex[i]=="F" & Sex_frame$Sex[i]==Sex_frame$Sex[j]){
#            SEXM[i,j]= "FF"}
#        if(Sex_frame$Sex[i]=="M" & Sex_frame$Sex[i]==Sex_frame$Sex[j]){
#           SEXM[i,j]= "MM"}
#        if( Sex_frame$Sex[i]!=Sex_frame$Sex[j]){
#            SEXM[i,j]= "FM"}
#    }
#}
#dimnames(SEXM)<-c(key, key)

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
Parasite <- subset_taxa(PS.T, Genus %in%c("Eimeria", "Cryptosporidium", "Tritrichomonas",
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
            LocM[i,j]= "1"
        } else{
            LocM[i,j]= "0"
        }
    }
}
#Note that Locality similarity matrix has rownames and colnames in the same order as key
all(rownames(LocM)==key$ID)
dimnames(LocM) <- c(key, key)


# 6) this matrix will describe the distance in years between samples
#Transform dates into a numeric variable
metadt$Year <- as.numeric(metadt$Year)
#Create data frame with each sample name (character) and sampling time (numeric)
SampleTime_frame<-metadt[,c("Mouse_ID","Year")]
#Create an empty matrix to fill with distances
TEMPM<-array(0,c(nrow(SampleTime_frame),nrow(SampleTime_frame)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(SampleTime_frame)){
 for (j in 1:nrow(SampleTime_frame))
{TEMPM[i,j]=abs(SampleTime_frame$Year[i] -SampleTime_frame$Year[j])
  }
}

dimnames(TEMPM)<-c(key,key)

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


# here are our matrices
str(BMIM)
str(SPATM)
str(JACM)
str(gen.dis)
#str(SEXM)
str(PARA)
str(LocM)
str(HIM)
str(TEMPM)
#First unravel the matrices into vectors matching the lower quantile of each matrix.

#From numeric matrices, this can be done by making a list (c()) of the distance object (dist()) derived from the matrix. as.dist() by default includes only the lower quantile of the matrix and excludes the diagonal.
#From categorical matrices, this can be done by making a list (c()) of the lower quantile of the matrix with lower.tri() -function.

jac<-c(as.dist(JACM))
bmi<-c(as.dist(BMIM))
spa<-c(as.dist(SPATM))
gen<-c(as.dist(gen.dis))
#sex<-c(SEXM[lower.tri(SEXM)])
para<-c(as.dist(PARA))
loc <- c(as.dist(LocM))
HIm <- c(as.dist(HIM))
tempm <- c(as.dist(TEMPM))

#Combine these vectors into a data frame
data.dyad<-data.frame(Parasite=para, BMI=bmi,Microbiome_similarity=jac,spatial=spa, genetic_dist=gen, locality=loc, HI=HIm, year=tempm)

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

## sex combination into a factor
#data.dyad$sex_combination <- factor(data.dyad$sex_combination, levels=c("MM", "FM", "FF"))

#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","genetic_dist", "Parasite", "BMI", "year")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }

#hist(data.dyad$Microbiome_similarity)
# proportional values semi-normally distributed limited between 0 and 1 not including 1 and 0 --> best use betaregression, but gaussian would probably give similar estimates

## now let's take a look at correlation trends
preplot1<-ggplot(data = data.dyad, aes(x= Microbiome_similarity, y= year))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .5, alpha= .8)+
    theme_bw()
preplot1

preplot2<-ggplot(data = data.dyad, aes(x= Microbiome_similarity, y= BMI))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .5, alpha= .8)+
    theme_bw()
preplot2

preplot3<-ggplot(data = data.dyad, aes(x= Microbiome_similarity, y= Parasite))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .5, alpha= .8)+
    theme_bw()
preplot3

preplot4<-ggplot(data = data.dyad, aes(x= Microbiome_similarity, y= genetic_dist))+
    geom_point(size= 1.2, alpha= .8, position= "jitter", aes(fill=sex_combination))+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .8, aes(col=sex_combination))+
    facet_wrap(~sex_combination)+
    theme_bw()
preplot4

preplot5<-ggplot(data = data.dyad, aes(x= genetic_dist, y= HI))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", alpha= .8)+
    theme_bw()

cor.test(data.dyad$HI, data.dyad$genetic_dist)

ggplot(data = data.dyad, aes(x=spatial , y= locality))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .8)+
    theme_bw()

cor.test(data.dyad$spatial, data.dyad$locality)

pdf("fig/genetic_distance_HI.pdf",
    width=6,
    height=4)
preplot5
dev.off()

## Let's model
names(data.dyad)

#### This is our final model
#model1<-brm(Microbiome_similarity~1+ spatial+locality+BMI+Parasite+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 50, chains = 10,
#                inits=0)
#saveRDS(model1, "tmp/BRMmodel1.rds")
#
model1 <- readRDS("tmp/BRMmodel1.rds")

#MCMCglmm gaussian model
mcmcglmm_model<-MCMCglmm(Microbiome_similarity~1+Parasite+locality+BMI+spatial+year,
                            data=data.dyad,
                            family= "gaussian",
                            random =~ mm(IDA+IDB),
                            verbose=FALSE)
saveRDS(mcmcglmm_model, "tmp/mcmcglmm_model.rds")
mcmcglmm_model <- readRDS("tmp/mcmcglmm_model.rds")

summary(mcmcglmm_model)
