# This script follows the super cool approach from Aura Raulo
# https://github.com/nuorenarra/Analysing-dyadic-data-with-brms/blob/main/R_Making_dyadic_data/DYADIC_workshop_data_wrangling.Rmd

library("MCMCglmm")
library(ape)
library(brms)
library(rstan)
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(ggmcmc)
library(ggthemes)
library(ggridges)
library(vegan)
library(phyloseq)
library(ggplot2)
library(bayesplot)
library(bayestestR)
library(brms)
library(MCMCglmm)
library(doParallel)
library(parallel)
library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(posterior)
library(distributional)


PS.T <- readRDS("tmp/PS.T.rds")

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
key <- data.frame(ID=sample_data(PS.T)$Mouse_ID)
metadt <- sample_data(PS.T)

# first we need to remove parasites from the microbiome PS object because we will analyse those in separate
PS.mP <- subset_taxa(PS.T, !Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(PS.mP, method="jaccard", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
dimnames(JACM)<- c(key, key)
# distance matrix
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
            LocM[i,j]= 1
        } else{
            LocM[i,j]= 0
        }
    }
}
#Note that AGE similarity matrix has rownames and colnames in the same order as key
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


# here are our matrices
str(BMIM)
str(SPATM)
str(JACM)
str(gen.dis)
#str(SEXM)
str(PARA)
str(LocM)
str(HIM)
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
         
#Combine these vectors into a data frame
data.dyad<-data.frame(Parasite=para, BMI=bmi,Microbiome_similarity=jac,spatial=spa, genetic_dist=gen, locality=loc, HI=HIm)

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
scalecols<-c("spatial","genetic_dist", "Parasite", "BMI")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }

#hist(data.dyad$Microbiome_similarity)
# proportional values semi-normally distributed limited between 0 and 1 not including 1 and 0 --> best use betaregression, but gaussian would probably give similar estimates

## now let's take a look at correlation trends
preplot1<-ggplot(data = data.dyad, aes(x= Microbiome_similarity, y= genetic_dist))+
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

pdf("fig/genetic_distance_HI.pdf",
    width=6,
    height=4)
preplot5
dev.off()

## Let's model
names(data.dyad)

# some BMI have NAs, this is not a problem for this model, but it will be for MCMCglmm (null model)
#data.dyad <- data.dyad[-which(is.na(data.dyad$BMI)),]  
#
#model1<-brm(Microbiome_similarity~1+ spatial+BMI+genetic_dist+Parasite +
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 50, chains = 10,
#                inits=0)
#saveRDS(model1, "tmp/BRMmodel1.rds")
#
#### This is our final model
#modelLocSpa<-brm(Microbiome_similarity~1+ spatial+locality+BMI+genetic_dist+Parasite+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 50, chains = 10,
#                inits=0)
#saveRDS(modelLocSpa, "tmp/BRMmodelLocSpac.rds")

modelLocSpa <- readRDS("tmp/BRMmodelLocSpac.rds")

#modelHI<-brm(Microbiome_similarity~1+ spatial+locality+BMI+HI+Parasite+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 50, chains = 10,
#                inits=0)
#saveRDS(modelHI, "tmp/BRMmodelHI.rds")
#model<-readRDS("tmp/BRMmodelLocSpa.rds")

#MCMCglmm gaussian model
#mcmcglmm_model<-MCMCglmm(Microbiome_similarity~1+Parasite+BMI+spatial+genetic_dist+locality,
#                            data=data.dyad,
#                            family= "gaussian",
#                            random =~ mm(IDA+IDB),
#                            verbose=FALSE)
#saveRDS(mcmcglmm_model, "tmp/mcmcglmm_model.rds")
mcmcglmm_model <- readRDS("tmp/mcmcglmm_model.rds")

summary(mcmcglmm_model)

#Denisty overlay = # Compare distribution of response variable to distributions of a set of predicted response variable values based on model -- are they a good fit?
#pp1<-pp_check(model)

summary(modelLocSpa)

model1_transformed <- ggs(modelLocSpa)

ggplot(filter(model1_transformed, Parameter %in% c("b_Intercept", "b_spatial", "b_BMI", "b_genetic_dist", "b_sex_combinationFF", "b_sex_combinationFM", "b_Parasite")),
       aes(x   = Iteration,
           y   = value,
           col = as.factor(Chain)))+
    geom_line() +
    geom_vline(xintercept = 1000)+
        facet_grid(Parameter ~ . ,
                      scale  = 'free_y',
                      switch = 'y')+
            labs(title = "Caterpillar Plots",
                        col   = "Chains")

ggplot(filter(model1_transformed, Parameter == "b_Intercept", Iteration > 1000), aes(x = value))+
    geom_density(fill  = "yellow", alpha = .5)+
    geom_vline(xintercept = 0, col  = "red", size = 1)+
#        scale_x_continuous(name= "Value", limits = c(-1, 3)) +
#        geom_vline(xintercept = summary(model1)$fixed[1,3:4],
 #                         col = "blue",
 #                         linetype = 2) +
    theme_light() +
    labs(title = "Posterior Density of Intercept")


summary(modelLocSpa)$fixed


ggplot(filter(model1_transformed, Parameter == "b_spatial", Iteration > 1000), aes(x = value))+
    geom_density(fill = "orange", alpha = .5)+
    geom_vline(xintercept = 0, col = "red", size = 1)+
#        scale_x_continuous(name = "Value", limits = c(-.2, .6))+
#            geom_vline(xintercept = summary(model1)$fixed[6,3:4], col = "blue", linetype = 2)+
                theme_light()+
                      labs(title = "Posterior Density of Regression Coefficient for Extraversion")


ggplot(filter(model1_transformed, Parameter == "b_locality", Iteration > 1000), aes(x = value))+
    geom_density(fill = "orange", alpha = .5)+
    geom_vline(xintercept = 0, col = "red", size = 1)+
      #  scale_x_continuous(name = "Value", limits = c(-.2, .6))+
#            geom_vline(xintercept = summary(model1)$fixed[6,3:4], col = "blue", linetype = 2)+
                theme_light()+
                      labs(title = "Posterior Density of Regression Coefficient for Extraversion")

ggplot(filter(model1_transformed, Parameter == "b_HI", Iteration > 1000), aes(x = value))+
    geom_density(fill = "orange", alpha = .5)+
    geom_vline(xintercept = 0, col = "red", size = 1)+
      #  scale_x_continuous(name = "Value", limits = c(-.2, .6))+
#            geom_vline(xintercept = summary(model1)$fixed[6,3:4], col = "blue", linetype = 2)+
                theme_light()+
                      labs(title = "Posterior Density of Regression Coefficient for Extraversion")

ggplot(filter(model1_transformed, Parameter == "b_genetic_dist", Iteration > 1000), aes(x = value))+
    geom_density(fill = "orange", alpha = .5)+
    geom_vline(xintercept = 0, col = "red", size = 1)+
      #  scale_x_continuous(name = "Value", limits = c(-.2, .6))+
#            geom_vline(xintercept = summary(model1)$fixed[6,3:4], col = "blue", linetype = 2)+
                theme_light()+
                      labs(title = "Posterior Density of Regression Coefficient for Extraversion")



                                        #A. Quick LINEPLOT for non-interaction models
plot1 <- mcmc_plot(modelLocSpa,
                type = "intervals",
                prob = 0.95,
                pars= rownames(fixef(modelLocSpa))[2:nrow(fixef(modelLocSpa))])

plot1

                                        #B. More flexible LINEPLOT in ggplot:
resdf1<-summary(modelLocSpa)$fixed
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
    labs(x="Effect on Microbiota",y="")+
    geom_vline(xintercept=0, linetype="dashed")

plot2

rownames(fixef(modelLocSpa))[2:nrow(fixef(modelLocSpa))]

# Alternative: Construct the model with MCMCglmm package (limitations with response distribution, but quicker)

# cross validation by approximat leave-one-out corss validation (LOO-CV)
#loo_t <- loo(modelLocSpa)


## compare the models
#loo_compare(loo(x), loo(y))


###################

# sanity check
all(sample_names(PS.mP)==key$ID)

# Data wrangling step: Give all unclassified genera a name based on their family or order

tax<-as.data.frame(tax_table(PS.mP))

tax[is.na(tax$Genus),]$Genus <- paste0("Unknown_genus_in_",tax[is.na(tax$Genus),]$Family)

tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Order)

tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Class)

tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Phylum)

tax_G<-as.data.frame(table(tax[,6]))

genuses<- as.character(tax_G$Var1)

genuses <- append('none', genuses)

dropnumber <- length(genuses)#199
model_minutes<-2
cores <- 80
(dropnumber*model_minutes)/cores



#start cluster
dropRes_list<-list()

cl <- parallel::makeCluster(80, type="FORK")

doParallel::registerDoParallel(cl)

## changing names here, to test the loop below,

key<-data.frame(ID=sample_data(PS.mP)$Mouse_ID, Sample_name=sample_data(PS.mP)$Mouse_ID)

data.dyad_REAL <- data.dyad

names(data.dyad_REAL)

dropRes_list<-foreach(i = 1:length(genuses)) %dopar% {
    library(phyloseq)
    #choose the genus to drop and prune the taxa to keep everything else
    gen.i<-genuses[i]
    taxa_tokeep<-rownames(tax[which(tax$Genus!=genuses[i]),])
    mic.i<-prune_taxa(taxa_tokeep, PS.mP)
    #Calculate Jaccard and Bray-Curtis microbiome dissimilarity for each mouse pair
    JACM.i<- as.matrix(phyloseq::distance(mic.i, method="jaccard"))
    BRAY.i<- as.matrix(phyloseq::distance(mic.i, method="bray"))
    #Unravel dissimilarity matrices into vectors
    bray<-c(as.dist(BRAY.i))
    jac<-c(as.dist(JACM.i))

    #Make a new dyadic data frame from these vectors and order it to be in the same order as       the original dyadic data frame
    data.dyad.i<-data.frame(Jaccard=jac,BrayCurtis=bray)
    # extracting Sample_name-combinations of the matrix
    list<-expand.grid(key$Sample_name,key$Sample_name)
    # This created sample-to-same-sample pairs as well. Get rid of these:
    list<-list[which(list$Var1!=list$Var2),]
    # the resulting list still has both quantiles of the original matrix (i.e. all values are doubled) in--> add 'unique' key and subset to one quantile only
    list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
    list<-subset(list, !duplicated(list$key))
    # Sample_name combinations are now in the same order as the lower quantile value vector
    # So we can add dyad name and each participant ID to dyadic dataframe
    data.dyad.i$Sample_A<-list$Var2
    data.dyad.i$Sample_B<-list$Var1

    # extracting combinations of individual IDs for each pair
    keyA<-key[,c("ID","Sample_name")]
    colnames(keyA)<-c("IDA","Sample_A")
    keyB<-key[,c("ID","Sample_name")]
    colnames(keyB)<-c("IDB","Sample_B")

    keyA<-keyA[match(data.dyad.i$Sample_A,keyA$Sample_A),]
    keyB<-keyB[match(data.dyad.i$Sample_B,keyB$Sample_B),]

    data.dyad.i$IDA<-keyA$IDA
    data.dyad.i$IDB<-keyB$IDB

    # Make sure you have no self comparisons in the data (This is the case by default here,       since we are using just one sample per individual)
    data.dyad.i<-data.dyad.i[which(data.dyad.i$IDA!=data.dyad.i$IDB),] #

    ### Combine new Jaccard variable with rest of dyadic data columns
    data.dyad<-data.dyad_REAL
    data.dyad$Jaccard<-data.dyad.i$Jaccard

    #factorize terms used for multimembership random structure and make sure levels are     same and in same order
    data.dyad$IDA<-as.factor(data.dyad$IDA)
    data.dyad$IDB<-as.factor(data.dyad$IDB)
    all(levels(data.dyad$IDA)==levels(data.dyad$IDB))#T

#The MCMCglmm model
dropmodel<-MCMCglmm(Microbiome_similarity~1+Parasite+BMI+spatial+genetic_dist+locality,
                            data=data.dyad,
                            family= "gaussian",
                            random =~ mm(IDA+IDB),
                            verbose=FALSE)
saveRDS(mcmcglmm_model, "tmp/mcmcglmm_model.rds")

   ASVs_dropped.i<-nrow(tax_table(PS.mP))-nrow(tax_table(mic.i))
   resdf.i<-data.frame(Genus_dropped=gen.i,
                      ASVs_dropped=ASVs_dropped.i,
                      Parasite_Estimate=summary(dropmodel)$solutions["Parasite",]["post.mean"],
                      Parasite_lCI=summary(dropmodel)$solutions["Parasite",]["l-95% CI"],
                      Parasite_uCI=summary(dropmodel)$solutions["Parasite",]["u-95% CI"],
                      genetic_Estimate=summary(dropmodel)$solutions["genetic_dist",]["post.mean"],
                      genetic_lCI=summary(dropmodel)$solutions["genetic_dist",]["l-95% CI"],
                      genetic_uCI=summary(dropmodel)$solutions["genetic_dist",]["u-95% CI"]
                      )
    return(resdf.i)
}

parallel::stopCluster(cl)
saveRDS(dropRes_list,"tmp/dropRes_list.rds")

#########################################################################
##rbind the resulting data frames to single master data frame
        dropResults<-data.frame(Genus_dropped=NA,
                                            ASVs_dropped=NA,
                                            Parasite_Estimate=NA,
                                            Parasite_lCI=NA,
                                            Parasite_uCI=NA,
                                            genetic_Estimate=NA,
                                            genetic_lCI=NA,
                                            genetic_uCI=NA)
                                            
  for(j in 1:length(genuses)){
    dropResults<-rbind(dropResults,dropRes_list[[j]])
  }

dropResults<-dropResults[2:nrow(dropResults),]
dropResults$Genus_dropped2<-genuses

dropRes_list[[100]]

saveRDS(dropResults,"dropResults_genus.rds")



#For each microbial genera, calculate their importance on social association effect on microbiome
Parasite_CIbr_baseline<-dropResults[which(dropResults$Genus_dropped=="none"),]$Parasite_CIbr

dropResults$Parasite_CIbr<-abs(dropResults$Parasite_uCI-dropResults$Parasite_lCI)

dropResults$Parasite_CIbr_increase<- dropResults$Parasite_CIbr-Parasite_CIbr_baseline

                                        #Only those genera whose drop increases the uncertainty around an effect will be considered important for that effect:
dropResults[which(dropResults$Parasite_CIbr_increase<0),]$Parasite_CIbr_increase<-0

                                        #Importance value is this increase in uncertainty divided by the square root of how many ASVs were dropped, and multiplied by 100 to increase the scale.
dropResults$IMPORTANCE_PARASITE<-(dropResults$Parasite_CIbr_increase/Parasite_CIbr_baseline)/sqrt(dropResults$ASVs_dropped)*100

dropResults$Parasite_CIbr_increase/Parasite_CIbr_baseline/sqrt(

dropResults$ASVs_dropped

dropResults$Genus_dropped

                                                              
dropResults$IMPORTANCE_PARASITE
