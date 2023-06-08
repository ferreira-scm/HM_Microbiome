################################### ok now we want to know at a broader sense wether these effects are driven by bacteria or eukarya, 

# sanity check
all(sample_names(PS.TSS)==key$ID)

## changing domain name

Diet <- subset_taxa(PS.TSS, Phylum %in% c("Anthophyta", "Phragmoplastophyta", "Charophyta", "Ochrophyta"))
PS.TSS@tax_table[taxa_names(Diet),1] <- "Diet"


Other <- subset_taxa(PS.TSS, Phylum %in% c("Apicomplexa", "Platyhelminthes", "Incertae_Sedis", "Vertebrata", "Tardigrada", "Parabasalia", "Rhabditida", "Mollusca", "Nematozoa", "Arthropoda"))

Unknown <- subset_taxa(PS.TSS, Phylum %in% c("Unknown_phylum", "Unknown_phylum_in_Eukarya"))

PS.TSS@tax_table[taxa_names(Other),1] <- "Other"

Fungi <- subset_taxa(PS.TSS, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota"))
PS.TSS@tax_table[taxa_names(Fungi),1] <- "Fungi"

Parasite <- subset_taxa(PS.TSS, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))
PS.TSS@tax_table[taxa_names(Parasite),1] <- "Parasite"

get_taxa_unique(PS.TSS, "Kingdom")

tax <- as.data.frame(PS.TSS@tax_table)

tax_D<-as.data.frame(table(tax[,1]))
Domain<- as.character(tax_D$Var1)

Domain <- c("none", Domain)

#start cluster
dropRes_list<-list()
cl <- parallel::makeCluster(30, type="FORK")
doParallel::registerDoParallel(cl)

dropResDOM_list<-foreach(i = 1:length(Domain)) %dopar% {
    library(phyloseq)
    #choose the genus to drop and prune the taxa to keep everything else
    gen.i<-Domain[i]
    taxa_tokeep<-rownames(tax[which(tax$Kingdom!=Domain[i]),])
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
    scalecols<-c("spatial","genetic_dist", "BMI", "year", "hi")
    range.use <- function(x,min.use,max.use){(x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
    for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
        data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }
#Transpose Jaccard dissimilarity to similarity
        data.dyad$Microbiome_similarity<-1-data.dyad$Microbiome_similarity
    
#The MCMCglmm model
dropmodel<-MCMCglmm(Microbiome_similarity~1+spatial+locality+genetic_dist*hi+year+BMI+sex,
                            data=data.dyad,
                            family= "gaussian",
                            random =~ mm(IDA+IDB),
                            verbose=FALSE)

   ASVs_dropped.i<-nrow(tax_table(PS.TSS))-nrow(tax_table(mic.i))
   resdf.i<-data.frame(Domain_dropped=gen.i,
                      ASVs_dropped=ASVs_dropped.i,
                genetic_dist_Estimate=summary(dropmodel)$solutions["genetic_dist",]["post.mean"],
                     genetic_dist_lCI=summary(dropmodel)$solutions["genetic_dist",]["l-95% CI"],
                     genetic_dist_uCI=summary(dropmodel)$solutions["genetic_dist",]["u-95% CI"],
                     spatial_Estimate=summary(dropmodel)$solutions["spatial",]["post.mean"],
                          spatial_lCI=summary(dropmodel)$solutions["spatial",]["l-95% CI"],
                          spatial_uCI=summary(dropmodel)$solutions["spatial",]["u-95% CI"],
                          hi_Estimate=summary(dropmodel)$solutions["hi",]["post.mean"],
                               hi_lCI=summary(dropmodel)$solutions["hi",]["l-95% CI"],
                               hi_uCI=summary(dropmodel)$solutions["hi",]["u-95% CI"],
                    locality_Estimate=summary(dropmodel)$solutions["locality1",]["post.mean"],
                         locality_lCI=summary(dropmodel)$solutions["locality1",]["l-95% CI"],
                         locality_uCI=summary(dropmodel)$solutions["locality1",]["u-95% CI"],
                        year_Estimate=summary(dropmodel)$solutions["year",]["post.mean"],
                             year_lCI=summary(dropmodel)$solutions["year",]["l-95% CI"],
                             year_uCI=summary(dropmodel)$solutions["year",]["u-95% CI"],
                         BMI_Estimate=summary(dropmodel)$solutions["BMI",]["post.mean"],
                              BMI_lCI=summary(dropmodel)$solutions["BMI",]["l-95% CI"],
                              BMI_uCI=summary(dropmodel)$solutions["BMI",]["u-95% CI"],
                gen_hi_Estimate=summary(dropmodel)$solutions["genetic_dist:hi",]["post.mean"],
                         gen_hi_lCI=summary(dropmodel)$solutions["genetic_dist:hi",]["l-95% CI"],
                gen_hi_uCI=summary(dropmodel)$solutions["genetic_dist:hi",]["u-95% CI"]
                      )
    return(resdf.i)
}
parallel::stopCluster(cl)
saveRDS(dropResDOM_list,"tmp/dropResDOMAIN_list.rds")

dropResDOM_list <- readRDS("tmp/dropResDOMAIN_list.rds")
