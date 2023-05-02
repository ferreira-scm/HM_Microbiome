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


####################
## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(PS.TSS, method="jaccard", type="samples"))
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

## 3) Sex pairs
Sex_frame<-metadt[,c("Mouse_ID","Sex")]
Sex_frame$Mouse_ID<-as.character(Sex_frame$Mouse_ID)
Sex_frame$Sex<-as.character(Sex_frame$Sex)
#Create an empty character matrix to fill with characters
SEXM<-array(as.character(NA),c(nrow(Sex_frame),nrow(Sex_frame)))

for(i in 1:nrow(Sex_frame)){
    for(j in 1:nrow(Sex_frame)){
        if(Sex_frame$Sex[i]=="F" & Sex_frame$Sex[i]==Sex_frame$Sex[j]){
            SEXM[i,j]= "FF"}
        if(Sex_frame$Sex[i]=="M" & Sex_frame$Sex[i]==Sex_frame$Sex[j]){
           SEXM[i,j]= "MM"}
        if( Sex_frame$Sex[i]!=Sex_frame$Sex[j]){
            SEXM[i,j]= "FM"}
    }
}
dimnames(SEXM)<-c(key, key)

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

# 7) Making hi distances
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


# here are our matrices
str(BMIM)
str(SPATM)
str(JACM)
str(gen.dis)
str(SEXM)
str(LocM)
str(HIM)
str(hiM)
str(TEMPM)
#First unravel the matrices into vectors matching the lower quantile of each matrix.

#From numeric matrices, this can be done by making a list (c()) of the distance object (dist()) derived from the matrix. as.dist() by default includes only the lower quantile of the matrix and excludes the diagonal.
#From categorical matrices, this can be done by making a list (c()) of the lower quantile of the matrix with lower.tri() -function.

jac<-c(as.dist(JACM))
bmi<-c(as.dist(BMIM))
spa<-c(as.dist(SPATM))
gen<-c(as.dist(gen.dis))
sex<-c(SEXM[lower.tri(SEXM)])
loc <- as.character(c(as.dist(LocM)))
HIm <- c(as.dist(HIM))
him <- c(as.dist(hiM))
tempm <- c(as.dist(TEMPM))

#Combine these vectors into a data frame
data.dyad<-data.frame(BMI=bmi,Microbiome_similarity=jac,spatial=spa, genetic_dist=gen, locality=loc, HI=HIm, year=tempm, hi=him, sex=sex)

data.dyad$locality <- as.factor(data.dyad$locality)

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
data.dyad$sex <- factor(data.dyad$sex, levels=c("MM", "FM", "FF"))

#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","genetic_dist", "BMI", "year", "hi")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }

saveRDS(data.dyad, "tmp/data.dyad.RDS")

data.dyad <- readRDS("tmp/data.dyad.RDS")

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

preplot4<-ggplot(data = data.dyad, aes(x= Microbiome_similarity, y= genetic_dist))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .8)+
    theme_bw()
preplot4

ggplot(data = data.dyad, aes(x= genetic_dist, y= hi))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", alpha= .8)+
    theme_bw()

cor.test(data.dyad$HI, data.dyad$genetic_dist)


ggplot(data = data.dyad, aes(x= Microbiome_similarity, y= hi))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .8)+
    theme_bw()

ggplot(data = data.dyad, aes(x=spatial , y= locality))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, se= FALSE, col= "red", size= .8)+
    theme_bw()

## Let's model
names(data.dyad)

#### This is our final model
#model1<-brm(Microbiome_similarity~1+ spatial+locality+genetic_dist*hi+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
 #               inits=0)
#saveRDS(model1, "tmp/BRMmodel1.rds")
model1 <- readRDS("tmp/BRMmodel1.rds")
####
#model2<-brm(Microbiome_similarity~1+ spatial+genetic_dist*hi+year+BMI+sex+
#                (1|mm(IDA,IDB))+(1|locality),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#                inits=0)
#saveRDS(model2, "tmp/BRMmodel2.rds")
model2 <- readRDS("tmp/BRMmodel2.rds")

#model3<-brm(Microbiome_similarity~1+ spatial+locality*genetic_dist+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 50, chains = 10,
#                inits=0)
#saveRDS(model3, "tmp/BRMmodel3.rds")
#model3 <- readRDS("tmp/BRMmodel3.rds")

#model3

#conditional_effects(model3)

library(tidybayes)

#brms::posterior_epred(model1)

head(data.dyad)
library(microshades)

library(cowplot)
    
newdata <- data.frame(spatial=seq_range(data.dyad$spatial, n=51),
                      BMI=rep(median(data.dyad$BMI), n=51),
                      hi=rep(median(data.dyad$hi),n=51),
              year=rep(0, n=51),
              locality=rep(0, n=51),
              IDA=rep("AA_0197", 51),
              sex=rep("MM", 51),
              IDB=rep("AA_0089", 51),
              genetic_dist=rep(median(data.dyad$genetic_dist)))

pred.df <- add_predicted_draws(newdata, model1)

spatial_pred <- ggplot(data.dyad, aes(x=spatial, y=Microbiome_similarity))+
    geom_point(shape=21, size=1, colour="gray", alpha=0.7)+
    stat_lineribbon(data=pred.df, aes(y = .prediction),
                    size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
    scale_fill_manual(values=microshades_palette("micro_cvd_purple"))+
    ylab("Gut community similarity")+
    xlab("Spatial distance")+
            labs(fill="level:")+
    theme_bw(base_size=12)

spatial_pred

newdata <- data.frame(genetic_dist=seq_range(data.dyad$genetic_dist, n=51),
              BMI=rep(median(data.dyad$BMI), n=51),
              year=rep(0, n=51),
                                        #              hi=rep(c(0, 0.5, 1),17),
              hi=rep(median(data.dyad$hi), n=51),
              locality=rep(0, n=51),
              sex=rep("MM", 51),
              IDA=rep("AA_0197", 51),
              IDB=rep("AA_0089", 51),
              spatial=rep(median(data.dyad$spatial)))
pred.df <- add_predicted_draws(newdata, model1)

gen_pred <- ggplot(data.dyad, aes(x=genetic_dist, y=Microbiome_similarity))+
    geom_jitter(width=0.01, shape=21, size=1, colour="gray", alpha=0.7)+
    stat_lineribbon(data=pred.df, aes(y = .prediction),
                    size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
    scale_fill_manual(values=microshades_palette("micro_purple"))+
    ylab("Gut community similarity")+
    xlab("Genetic distance")+
            labs(fill="level:")+
    theme_bw(base_size=12)

gen_pred

fig2 <- plot_grid(spatial_pred, gen_pred, labels="auto", nrow=1)

ggsave("fig/figure2.pdf", fig2, width=170, height=85, units="mm", dpi=300)

names(pred.df)

#MCMCglmm gaussian model
mcmcglmm_model<-MCMCglmm(Microbiome_similarity~1+locality+spatial+year+genetic_dist,
                            data=data.dyad,
                            family= "gaussian",
                            random =~ mm(IDA+IDB),
                            verbose=FALSE)
saveRDS(mcmcglmm_model, "tmp/mcmcglmm_model.rds")
mcmcglmm_model <- readRDS("tmp/mcmcglmm_model.rds")

summary(mcmcglmm_model)

#Denisty overlay = # Compare distribution of response variable to distributions of a set of predicted response variable values based on model -- are they a good fit?
pp_check(model1) # I guess not too bad. Could try gamma distribution I guess.

model1_transformed <- ggs(model1)

conditional_effects(model1)

(unique(model1_transformed$Parameter))[1:10]

ggplot(filter(model1_transformed, Parameter %in% c("b_Intercept", "b_spatial", "b_locality1", "b_genetic_dist", "b_hi", "b_year", "b_BMI", "b_sexFM", "b_sexFF", "b_genetic_dist:hi")),
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


ggplot(filter(model1_transformed, Parameter == "b_spatial", Iteration > 1000), aes(x = value))+
    geom_density(fill = "orange", alpha = .5)+
    geom_vline(xintercept = 0, col = "red", size = 1)+
#        scale_x_continuous(name = "Value", limits = c(-.2, .6))+
#            geom_vline(xintercept = summary(model1)$fixed[6,3:4], col = "blue", linetype = 2)+
                theme_light()+
                      labs(title = "Posterior Density of Regression Coefficient for Extraversion")


ggplot(filter(model1_transformed, Parameter == "b_hi", Iteration > 1000), aes(x = value))+
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
plot1 <- mcmc_plot(model1,
                type = "intervals",
                prob = 0.95,
                pars= rownames(fixef(model1))[2:nrow(fixef(model1))])

plot1

c(microshades_palette("micro_orange", 1, lightest = FALSE))

c(microshades_palette("micro_blue", 1, lightest = FALSE))

c(microshades_palette("micro_purple", 1, lightest = FALSE))

c(microshades_palette("micro_green", 1, lightest = FALSE))

c(microshades_palette("micro_brown", 1, lightest = FALSE))

c(microshades_palette("micro_cvd_orange", 1))

c(microshades_palette("micro_cvd_green", 1))
c(microshades_palette("micro_cvd_purple", 1, lightest = FALSE))
c(microshades_palette("micro_cvd_turquoise", 1, lightest = FALSE))





#B. More flexible LINEPLOT in ggplot:
resdf1<-summary(model1)$fixed
resdf1<-as.data.frame(resdf1)
resdf1<-resdf1[c("Estimate","l-95% CI","u-95% CI")]
resdf1<-resdf1[2:nrow(resdf1),]
resdf1$Predictor<-rownames(resdf1)
colnames(resdf1)<-c("Estimate","lCI","uCI","Predictor")

resdf1$Predictor <- c("Spatial distance", "Locality", "Genetic distance", "Hybridicity", "Year", "BMI", "Female-Male", "Female-Female", "Genetic distance: Hybridiciry")# rename
resdf1$Predictor<-factor(resdf1$Predictor)

plot2<-ggplot(resdf1,aes(x=Estimate,y=Predictor,colour=Predictor))+
    geom_linerange(aes(xmin = lCI, xmax = uCI),size=10)+
    scale_colour_viridis_d()+ 
    geom_point(pch="|", size=7,colour="black")+
    theme_bw()+
    theme(legend.position='none',text = element_text(size=14))+
    labs(x="Effect on the gut community similarity",y="")+
    geom_vline(xintercept=0, linetype="dashed")

plot2

ggsave("fig/figure1.pdf", plot2, width=170, height=75, units="mm", dpi=300)

# Alternative: Construct the model with MCMCglmm package (limitations with response distribution, but quicker)

# cross validation by approximat leave-one-out corss validation (LOO-CV)
#loo_t <- loo(modelLocSpa)


## compare the models
#loo_compare(loo(x), loo(y))


###################

# sanity check
all(sample_names(PS.TSS)==key$ID)

# Data wrangling step: Give all unclassified genera a name based on their family or order
tax <- as.data.frame(PS.TSS@tax_table)
#tax[is.na(tax$Genus),]$Genus
#tax[97,] <- "Unknown_genus_in_Oscillospirales" #  manual change here, need to go back and check
#tax[8,] <- "Unknown_genus_in_Eukarya"

tax_G<-as.data.frame(table(tax[,6]))
genuses<- as.character(tax_G$Var1)

genuses <- c("none", genuses)
#genuses <- append('none', genuses)

dropnumber <- length(genuses)#199
model_minutes<-60
cores <- 70
(dropnumber*model_minutes)/cores

key<-data.frame(ID=sample_data(PS.TSS)$Mouse_ID, Sample_name=sample_data(PS.TSS)$Mouse_ID)
data.dyad_REAL <- readRDS("tmp/data.dyad.RDS")

names(data.dyad_REAL)

#summary(mcmcglmm_model)


#start cluster
dropRes_list<-list()
cl <- parallel::makeCluster(60, type="FORK")
doParallel::registerDoParallel(cl)
dropRes_list<-foreach(i = 1:length(genuses)) %dopar% {
    library(phyloseq)
    #choose the genus to drop and prune the taxa to keep everything else
    gen.i<-genuses[i]
    taxa_tokeep<-rownames(tax[which(tax$Genus!=genuses[i]),])
    mic.i<-prune_taxa(taxa_tokeep, PS.TSS)
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
    data.dyad$Microbiome_similarity<-data.dyad.i$Jaccard

    #factorize terms used for multimembership random structure and make sure levels are     same and in same order
    data.dyad$IDA<-as.factor(data.dyad$IDA)
    data.dyad$IDB<-as.factor(data.dyad$IDB)
    all(levels(data.dyad$IDA)==levels(data.dyad$IDB))#T

# Scale all predictors not between 0-1 already to be between 0-1
    scalecols<-c("spatial","genetic_dist", "Parasite", "BMI", "year")
    range.use <- function(x,min.use,max.use){(x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
    for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
        data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }
#Transpose Jaccard dissimilarity to similarity
        data.dyad$Microbiome_similarity<-1-data.dyad$Microbiome_similarity
    
#The MCMCglmm model
dropmodel<-MCMCglmm(Microbiome_similarity~1+genetic_dist+BMI+spatial+locality,
                            data=data.dyad,
                            family= "gaussian",
                            random =~ mm(IDA+IDB),
                            verbose=FALSE)

   ASVs_dropped.i<-nrow(tax_table(PS.TSS))-nrow(tax_table(mic.i))
   resdf.i<-data.frame(Genus_dropped=gen.i,
                      ASVs_dropped=ASVs_dropped.i,
                      genetic_dist_Estimate=summary(dropmodel)$solutions["genetic_dist",]["post.mean"],
                      genetic_dist_lCI=summary(dropmodel)$solutions["genetic_dist",]["l-95% CI"],
                      genetic_dist_uCI=summary(dropmodel)$solutions["genetic_dist",]["u-95% CI"],
                      spatial_Estimate=summary(dropmodel)$solutions["spatial",]["post.mean"],
                      spatial_lCI=summary(dropmodel)$solutions["spatial",]["l-95% CI"],
                      spatial_uCI=summary(dropmodel)$solutions["spatial",]["u-95% CI"]
                      )
    return(resdf.i)
}


parallel::stopCluster(cl)
saveRDS(dropRes_list,"tmp/dropRes_list.rds")

dropRes_list <- readRDS("tmp/dropRes_list.rds")


    
#########################################################################
##rbind the resulting data frames to single master data frame
dropResults<-data.frame(Genus_dropped=NA,
                                            ASVs_dropped=NA,
                                            genetic_dist_Estimate=NA,
                                            genetic_dist_lCI=NA,
                                            genetic_dist_uCI=NA,
                                            spatial_Estimate=NA,
                                            spatial_lCI=NA,
                                            spatial_uCI=NA)


for(j in 1:length(genuses)){
    dropResults<-rbind(dropResults,dropRes_list[[j]])
  }

dropResults<-dropResults[2:nrow(dropResults),]
dropResults$Genus_dropped2<-genuses

saveRDS(dropResults,"tmp/dropResults_genus.rds")

dropResults <- readRDS("tmp/dropResults_genus.rds")

names(dropResults)

#For each microbial genera, calculate their importance on genetic and spatial effect on microbiome
dropResults$genetic_dist_CIbr<-abs(dropResults$genetic_dist_uCI-dropResults$genetic_dist_lCI)

genetic_dist_CIbr_baseline<-dropResults[which(dropResults$Genus_dropped=="none"),]$genetic_dist_CIbr

dropResults$genetic_dist_CIbr_increase<- dropResults$genetic_dist_CIbr-genetic_dist_CIbr_baseline

#Only those genera whose drop increases the uncertainty around an effect will be considered important for that effect:
dropResults[which(dropResults$genetic_dist_CIbr_increase<0),]$genetic_dist_CIbr_increase<-0

#Importance value is this increase in uncertainty divided by the square root of how many ASVs were dropped, and multiplied by 100 to increase the scale.
dropResults$IMPORTANCE_genetic_dist<-(dropResults$genetic_dist_CIbr_increase/genetic_dist_CIbr_baseline)/sqrt(dropResults$ASVs_dropped)*100

dropResults$IMPORTANCE_genetic_dist



######################## spatial
#For each microbial genera, calculate their importance on social association effect on microbiome
dropResults$spatial_CIbr<-abs(dropResults$spatial_uCI-dropResults$spatial_lCI)
spatial_CIbr_baseline<-dropResults[which(dropResults$Genus_dropped=="none"),]$spatial_CIbr
dropResults$spatial_CIbr_increase<- dropResults$spatial_CIbr-spatial_CIbr_baseline

#Only those genera whose drop increases the uncertainty around an effect will be considered important for that effect:
dropResults[which(dropResults$spatial_CIbr_increase<0),]$spatial_CIbr_increase<-0

#Importance value is this increase in uncertainty divided by the square root of how many ASVs were dropped, and multiplied by 100 to increase the scale.
dropResults$IMPORTANCE_spatial<-(dropResults$spatial_CIbr_increase/spatial_CIbr_baseline)/sqrt(dropResults$ASVs_dropped)*100

dropResults$Genus_dropped[dropResults$IMPORTANCE_spatial>0]

dropResults$ASVs_dropped

head(dropResults)

head(dropResults)

geneticImp <- dropResults[dropResults$IMPORTANCE_genetic_dist>0,]
geneticImp <- geneticImp[order(-geneticImp$IMPORTANCE_genetic_dist),]
geneticImp <- rbind(dropResults[dropResults$Genus_dropped=="none",], geneticImp)
geneticImp <- geneticImp[!is.na(geneticImp$Genus_dropped),]
geneticImp$Genus_dropped <- factor(geneticImp$Genus_dropped, levels=geneticImp$Genus_dropped[order(geneticImp$IMPORTANCE_genetic_dist)])

head(geneticImp)

levels(geneticImp$Genus_dropped)

gen_drop <- ggplot(geneticImp, aes(x=genetic_dist_Estimate, y=Genus_dropped))+
    geom_errorbar(aes(xmin=genetic_dist_lCI, xmax=genetic_dist_uCI), colour="gray", size=1, width=0.4)+
    geom_point()+
    theme_classic()

gen_drop

spatialImp <- dropResults[dropResults$IMPORTANCE_spatial>0,]
spatialImp <- spatialImp[order(-spatialImp$IMPORTANCE_spatial),]
spatialImp <-rbind(dropResults[dropResults$Genus_dropped=="none",], spatialImp)
spatialImp <- spatialImp[!is.na(spatialImp$Genus_dropped),]
spatialImp$Genus_dropped <- factor(spatialImp$Genus_dropped, levels=spatialImp$Genus_dropped[order(spatialImp$IMPORTANCE_spatial)])

spa_drop <- ggplot(spatialImp, aes(x=spatial_Estimate, y=Genus_dropped))+
    geom_errorbar(aes(xmin=spatial_lCI, xmax=spatial_uCI), colour="gray", size=1, width=0.4)+
    geom_point(stat="identity")+
    theme_classic()

spa_drop

figure3 <- plot_grid(gen_drop, spa_drop, labels="auto")

ggsave("fig/figure3.pdf", figure3, width=250, height=350, units="mm", dpi=300)

summary(dropResults$IMPORTANCE_spatial>0)
summary(dropResults$IMPORTANCE_genetic>0)


nrow(spatialImp)

nrow(geneticImp)


summary(spatialImp$Genus_dropped%in%geneticImp$Genus_dropped)

summary(geneticImp$Genus_dropped%in%spatialImp$Genus_dropped)

geneticImp$Genus_dropped[which(!geneticImp$Genus_dropped%in%spatialImp$Genus_dropped)]

tax[tax$Genus=="Anaerostipes",]

tax[tax$Genus=="Enterococcus",]

tax[tax$Genus=="Dorea",]

tax[tax$Genus=="Erwinia",]

Lachnospiraceae

