setwd("G:/Meu Drive/Artigos/Em preparação/Artigo Desmodus x isótopos/mixsiar")

Autor: Marcelo Magioli and Fernando Gonçalves
Project: Foraging choice of common vampire bat on introduced species: a precaution for restoration projects 


library(MixSIAR)

#Load consumer data
mix <- load_mix_data(filename='consumer.csv', 
                     iso_names=c("d13C","d15N"), 
                     factors=NULL, 
                     fac_random=NULL, 
                     fac_nested=NULL, 
                     cont_effects=NULL)

#Load source data
source <- load_source_data(filename='sources.csv',
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="means", 
                           mix)

source <- load_source_data(filename='sources2.csv',
                           source_factors=NULL, 
                           conc_dep=FALSE, 
                           data_type="raw", 
                           mix)

#Load TDF data
discr <- load_discr_data(filename='tdf.csv', mix)

# Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=F, plot_save_png=F, 
          mix,source,discr)

# Calculate the convex hull area, standardized by source variance
calc_area(source=source, mix=mix, discr=discr)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
plot_prior(alpha.prior=1,source, plot_save_pdf = F)
plot_prior(alpha.prior= c(0.11,0.25,4.27,0.37),source, plot_save_pdf = F)


# Write the JAGS model file
model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

jags.1 <- run_model(run="normal", mix, source, discr, model_filename,
                    alpha.prior = c(0.11,0.25,4.27,0.37), resid_err, process_err)

output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)

output_JAGS(jags.1, mix, source, output_options)


###----------- Extracting posterior probabilities
df <- jags.1$BUGSoutput$sims.matrix
write.csv(df, file = "density.model.mixsiar.csv")

library("ggplot2")

df <- read.csv(file = "density.model.mixsiar.csv", header = T)
df1 <- as.data.frame(df)

a <- df1[,7]/0.657677
b <- df1[,8]/0.668018
c <- df1[,9]/1
d <- df1[,10]/0.73936

df2 <- data.frame(a, b, c, d)
summary(df2)

### Plot densities
#plot(density(df2[,1]))

g <- ggplot() +
  theme_bw() +
  theme(legend.position = "bottom", panel.border = element_rect(colour = "black", fill = NA, size = .5),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"))

g <- g + geom_density(data = df2, aes(df2[,1],  y=..scaled.., alpha = 0.9,
                                      fill= "C. penicillata"))
g <- g + geom_density(data = df2, aes(df2[,2], y=..scaled.., fill= "D. aurita",
                                      alpha = 0.9))
g <- g + geom_density(data = df2, aes(df2[,3], y=..scaled.., fill= "H. hydrochaeris",
                                      alpha = 0.9))
g <- g + geom_density(data = df2, aes(df2[,4],  y=..scaled.., fill= "N. nasua",
                                      alpha = 0.9))

g <- g + scale_x_continuous(name = "Proportion of diet") +
  scale_y_continuous(name = "Scaled posterior density")

g + guides(fill=FALSE)

g

p8 = ggplot(airquality_trimmed, aes(x = Ozone, fill = Month.f)) +
  geom_density(position = "identity", alpha = 0.6) +
  scale_x_continuous(name = "Mean ozone in\nparts per billion",
                     breaks = seq(0, 200, 25), limits = c(0, 200)) +
  scale_y_continuous(name = "Density") +
  labs(title = "Frequency histogram of mean ozone",
       subtitle = "Source: New York State Department of Conservation") +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = .5),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"))
