#### Testing if there are spatial, host genetics interactions

#modelSpaGen_I<-brm(Microbiome_similarity~1+ spatial*genetic_dist+locality+BMI+Parasite+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 50, chains = 4,
#                inits=0)
#saveRDS(modelSpaGen_I, "tmp/BRMmodelSpaGen_I.rds")

modelParasite_BMI_I<-brm(Microbiome_similarity~1+ spatial+locality+BMI*Parasite+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 50, chains = 10,
                inits=0)
saveRDS(modelParasite_BMI_I, "tmp/BRMmodelParasite_BMI_I.rds")
#
modelParasite_BMI_I

readRDS("tmp/BRMmodelParasite_BMI_I.rds")

saveRDS(modelSpaGen_I, "tmp/BRMmodelSpaGen_I.rds")


modelSpaGen_I <- readRDS("tmp/BRMmodelSpaGen_I.rds")

get_variables(modelSpaGen_I)

conditional_effects(modelSpaGen_I)

#pp_check(modelSpaGen_I) # ok ish?

summary(modelSpaGen_I)

library(tidybayes)

(model_fit <- data.dyad %>%
    add_predicted_draws(modelSpaGen_I) %>%
    ggplot(aes(x = genetic_dist, y = Microbiome_similarity)) +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                       alpha = 0.5, colour = "black") +
        geom_point(data = data.dyad, colour = "darkseagreen4", size = 3) +   # raw data
            scale_fill_brewer(palette = "Greys") +
                ylab("Microbiome similarity\n") +
                xlab("\nGenetic distance") +
                theme_bw() +
                    theme(legend.title = element_blank(),
                                    legend.position = c(0.15, 0.85)))





add_predicted_draws(modelSpaGen_I)

add_precited_draws(modelSpaGen_I)
