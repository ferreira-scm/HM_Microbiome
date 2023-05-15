
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


hi_gen <-ggplot(data = data.dyad, aes(x= genetic_dist, y= hi))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, col= "red", size= .8,
                formula=y~x+I(x^2))+
    ylab("Hybridicity distance")+
    xlab("Genetic distance")+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
                    axis.title.y = element_text(vjust = 2, size = 12))


HI_hi <-ggplot(data = data.dyad, aes(x= HI, y= hi))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, col= "red", size= .8,
                formula=y~x+I(x^2))+
    ylab("Hybridicity distance")+
    xlab("Hybrid index distance")+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
                    axis.title.y = element_text(vjust = 2, size = 12))

library("ggpmisc")

library(cowplot)

hy.cor <- cor.test(data.dyad$HI, data.dyad$genetic_dist)

gen_hi <- ggplot(data = data.dyad, aes(x= genetic_dist, y= HI))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
    geom_smooth(method= lm, col= "red", size= .8,)+
    annotate(geom="text", x=0, y=1.05, hjust=0.05, label=paste("Pearson's rho=", round(hy.cor$estimate,2), ",p<0.001", "df= ", hy.cor$parameter, sep=""))+ 
    ylab("Genetic distance")+
    xlab("Hybrid index distance")+
    theme_bw(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 12),
                    axis.title.y = element_text(vjust = 2, size = 12))


gen_hi
    
Figure1 <- cowplot::plot_grid(gen_hi, hi_gen, nrow=1, labels="auto")

ggplot2::ggsave(file="fig/Figure1.pdf", Figure1, width = 170, height = 85, dpi = 200, units="mm")

ggplot2::ggsave(file="fig/hi_HI.pdf", HI_hi, width = 85, height = 85, dpi = 200, units="mm")



ggplot(data = data.dyad, aes(x=spatial , y= locality))+
    geom_point(size= 1.2, alpha= .8, position= "jitter")+
#    geom_smooth(method= lm, se= FALSE, col= "red", size= .8)+
    theme_bw()


## Let's model
names(data.dyad)

newdata <- data.frame(spatial=seq_range(data.dyad$spatial, n=51),
                      BMI=rep(median(data.dyad$BMI), n=51),
                      hi=rep(median(data.dyad$hi),n=51),
              year=rep(0, n=51),
              locality=rep(0, n=51),
              IDA=rep("AA_0197", 51),
              sex=rep("MM", 51),
              IDB=rep("AA_0089", 51),
              genetic_dist=rep(median(data.dyad$genetic_dist)))

pred.df <- add_epred_draws(newdata, model1)

spatial_pred <- ggplot(data.dyad, aes(x=spatial, y=Microbiome_similarity))+
    geom_point(shape=21, size=1, colour="gray", alpha=0.7)+
    stat_lineribbon(data=pred.df, aes(y = .epred),
                    size=0.5, .width=c(.95, .5), alpha=0.5) +
    scale_fill_manual(values=microshades_palette("micro_cvd_purple"))+
    ylab("Gut community similarity")+
    xlab("Spatial distance")+
            labs(fill="level:")+
    theme_bw(base_size=12)

spatial_pred

newdata0 <- data.frame(genetic_dist=seq_range(data.dyad$genetic_dist, n=51),
              BMI=rep(median(data.dyad$BMI), n=51),
              year=rep(0, n=51),
              hi=rep(0, n=51),
 #             hi=rep(median(data.dyad$hi), n=51),
              locality=rep(0, n=51),
              sex=rep("MM", 51),
              IDA=rep("AA_0197", 51),
              IDB=rep("AA_0089", 51),
              spatial=rep(median(data.dyad$spatial)))
pred.df <- add_epred_draws(newdata0, model1)

gen_pred0 <- ggplot(data.dyad, aes(x=genetic_dist, y=Microbiome_similarity))+
    geom_jitter(width=0.01, shape=21, size=1, colour="gray", alpha=0.7)+
    stat_lineribbon(data=pred.df, aes(y = .epred),
                    size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
    scale_fill_manual(values=microshades_palette("micro_purple"))+
    ylab("Gut community similarity")+
    xlab("Genetic distance")+
    labs(fill="level:")+
       ggtitle("hybridicity distance = 0")+
    theme_bw(base_size=12)


newdata1 <- data.frame(genetic_dist=seq_range(data.dyad$genetic_dist, n=51),
              BMI=rep(median(data.dyad$BMI), n=51),
              year=rep(0, n=51),
              hi=rep(1, n=51),
#              hi=rep(median(data.dyad$hi), n=51),
              locality=rep(0, n=51),
              sex=rep("MM", 51),
              IDA=rep("AA_0197", 51),
              IDB=rep("AA_0089", 51),
              spatial=rep(median(data.dyad$spatial)))
pred.df <- add_epred_draws(newdata1, model1)

gen_pred1 <- ggplot(data.dyad, aes(x=genetic_dist, y=Microbiome_similarity))+
    geom_jitter(width=0.01, shape=21, size=1, colour="gray", alpha=0.7)+
    stat_lineribbon(data=pred.df, aes(y = .epred),
                    size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
    scale_fill_manual(values=microshades_palette("micro_purple"))+
    ylab("Gut community similarity")+
    xlab("Genetic distance")+
    labs(fill="level:")+
    ggtitle("hybridicity distance = 1")+
    theme_bw(base_size=12)

newdatam <- data.frame(genetic_dist=seq_range(data.dyad$genetic_dist, n=51),
              BMI=rep(median(data.dyad$BMI), n=51),
              year=rep(0, n=51),
#              hi=rep(0, n=51),
              hi=rep(median(data.dyad$hi), n=51),
              locality=rep(0, n=51),
              sex=rep("MM", 51),
              IDA=rep("AA_0197", 51),
              IDB=rep("AA_0089", 51),
              spatial=rep(median(data.dyad$spatial)))

pred.df <- add_epred_draws(newdatam, model1)
gen_predm <- ggplot(data.dyad, aes(x=genetic_dist, y=Microbiome_similarity))+
    geom_jitter(width=0.01, shape=21, size=1, colour="gray", alpha=0.7)+
    stat_lineribbon(data=pred.df, aes(y = .epred),
                    size=0.5, .width=c(.95, .5), alpha=0.5) +
    scale_fill_manual(values=microshades_palette("micro_purple"))+
    ylab("Gut community similarity")+
    xlab("Genetic distance")+
            labs(fill="level:")+
    theme_bw(base_size=12)


gen_predm

FigureSx <- plot_grid(gen_pred0, gen_pred1, labels="auto", nrow=1)

fig2 <- plot_grid(spatial_pred, gen_predm, labels="auto", nrow=1)

ggsave("fig/figure2.pdf", fig2, width=170, height=85, units="mm", dpi=300)


model1_transformed <- ggs(model1)

#conditional_effects(model1)

(unique(model1_transformed$Parameter))[1:10]

caterpillar <- ggplot(filter(model1_transformed, Parameter %in% c("b_Intercept", "b_spatial", "b_locality1", "b_genetic_dist", "b_hi", "b_year", "b_BMI", "b_sexFM", "b_sexFF", "b_genetic_dist:hi")),
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


ggsave("fig/Sx_caterpillar.pdf", caterpillar, width=170, height=200, units="mm", dpi=300)

#B. More flexible LINEPLOT in ggplot:
resdf1<-summary(model1)$fixed
resdf1<-as.data.frame(resdf1)
resdf1<-resdf1[c("Estimate","l-95% CI","u-95% CI")]
resdf1<-resdf1[2:nrow(resdf1),]
resdf1$Predictor<-rownames(resdf1)
colnames(resdf1)<-c("Estimate","lCI","uCI","Predictor")

resdf1$Predictor <- c("Spatial distance", "Locality", "Genetic distance", "Hybridicity", "Year", "BMI", "Female-Male", "Female-Female", "Genetic distance: Hybridicity")# rename
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

ggsave("fig/figure1.pdf", plot2, width=170, height=90, units="mm", dpi=300)




names(pred.df)
