setwd("G:/Meu Drive/Artigos/Em preparação/Artigo Desmodus x isótopos")

library(SIDER)

SIDER_data <- scrumpSider(iso.data = "all")

SIDER_trees <- scrumpSider(tree = "all")

new_data_test <- recipeSider(species = "Leontopithecus_rosalia",
                             habitat = "terrestrial",
                             taxonomic.class = "mammalia",
                             tissue = "hair",
                             diet.type = "herbivore",
                             tree = SIDER_trees)

tdf_data_n <- prepareSider(data.estimate = new_data_test,
                           data.isotope = SIDER_data,
                           tree = SIDER_trees,
                           isotope = "carbon")

formula_n <- delta13C ~ diet.type + habitat
random_terms <- (~ animal + species + tissue)

prior <- list(R = list(V = 1, nu=0.002),
              G = list(G1=list(V = 1, nu=0.002),
                       G2=list(V = 1, nu=0.002),
                       G3=list(V = 1, nu=0.002)))

nitt <- c(1200000)
burnin <- c(200000)
thin <- c(500)
parameters <- c(nitt, thin, burnin)
no.chains <- c(2)

convergence =  c(1.1)
ESS = c(1000)

TDF_est_n <-imputeSider(mulTree.data = tdf_data_n, random_terms,
                        priors = prior,
                        formula = formula_n,
                        output = "test_n_run",
                        parameters = parameters,
                        chains = no.chains,
                        convergence =  convergence, 
                        ESS = ESS)

Credible_intervals <- hdrcde::hdr(TDF_est_n $tdf_global, prob =
                                    c(50, 95, 99))
Credible_intervals

summary(TDF_est_n$tdf_global)

file.remove(list.files(pattern = "test_n_run"))
