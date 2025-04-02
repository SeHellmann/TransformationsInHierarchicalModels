############################################################-
#####          Main script for the analyses in              -
#####      ** HIer awesome paper title einfügen **          -
#__________________________________________________________----

# Sebastian Hellmann, April 2025

REDOALLANALYSIS <- FALSE

## Structure:
# Preamble and imports    
#____________________________________________________________
# A  Refitting Rieskamp-data with original model                   
## 1. Read experimental data                                       
## 2. Fit the hierarchical CPT-model                               
## 3. Compare population means between transformations             
# B  Refitting Rieskamp-data with model with correlated random effs   --> STILL TO BE DONE
## 1. Read experimental data                                       
## 2. Fit the hierarchical CPT-model                               
## 3. Compare population means between transformations             
## 4. Check correlation between parameters                         
# C  Re-do (Extended) original simulation study (unconstrained)    
## 1. Subset data and define JAGS settings                         
## 2. Actual parameter recovery analysis                           
## 3. Visualize original full parameter recovery analysis          
# D  Re-do (Extended) original simulation study (alpha=beta)       
## 1. Subset data and define JAGS settings                         
## 2. Actual parameter recovery analysis                           
## 3. Visualize original full parameter recovery analysis          
# E  Simulation study using wider parameter ranges (alpha=beta)         --> STILL TO BE DONE
# F  Simulation study using wider parameter ranges + correlated pars    --> STILL TO BE DONE


# Preamble and imports                                     ----
rm(list = ls())
# use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
print("Working directory set to:")
print(getwd())

{
  # Tell Rstudio where to find JAGS
  Sys.setenv(JAGS_HOME = "C:/Users/go73jec/AppData/Local/Programs/JAGS/JAGS-4.3.1")
  library(tidyverse)
  library(R2jags)
  library(ggpubr)
  library(viridis)
  library(ggh4x)
  if (grepl("Windows", osVersion)) {
    # To use Times as font on Windows (otherwise ignore the ggplot warnings)
    windowsFonts(Times=windowsFont("Times New Roman"))
  }
  custom_theme <- theme_bw()+
    theme(text = element_text(family = "Times", size=12),
          legend.position = "bottom",
          plot.margin = margin(0, 0, 0, 0, "cm"),
          axis.text = element_text(size=12, family="Times", color="black"),
          legend.spacing = unit(0,"line"),
          legend.key.width=unit(1,"line"),
          panel.spacing=unit(0, "lines"))
  dir.create("figures", showWarnings = FALSE)
  dir.create("saved_details", showWarnings = FALSE)
}
## Import simulation function and define JAGS model file names
source("simulate_CPT.R")
original_full_model = "jags_models/cpt_hierarchical_model.txt"
original_full_model_recovery = "jags_models/cpt_hierarchical_recovery.txt"
original_restricted_model_recovery = "jags_models/cpt_hierarchical_restricted_recovery.txt"


original_full_model = "cpt_hierarchical_model.txt"
original_full_model = "cpt_hierarchical_model.txt"


#___________________________________________________________________----
# A  Refitting Rieskamp-data with original model                    ----
## 1. Read experimental data                                        ----

# Load information about the gamble-pairs used in Rieskamp (2008). 
# GambleA.txt and GambleB.txt are structured as follows: 
# value of outcome 1 (column 1), 
# probability of outcome 1 (column 2), 
# value of outcome 2 (column 3), 
# probability of outcome 2 (column 4) (gambles in rows).
prospects.b.temp <- as.matrix(read.table("GambleB.txt"))
prospects.a.temp <- as.matrix(read.table("GambleA.txt"))

prospects.b <- array(0,dim=c(180,4))
prospects.a <- array(0,dim=c(180,4))

# Arrange so that v and p related to the relatively poor outcome ends 
# up in collumn 1 and 2
for (i in 1:180){
  
  if (prospects.a.temp[i,1] < prospects.a.temp[i,3]){
    prospects.a[i,] <- prospects.a.temp[i,] 
  }else{
    prospects.a[i,1:2] <- prospects.a.temp[i,3:4] 
    prospects.a[i,3:4] <- prospects.a.temp[i,1:2] 
  }
  
  if (prospects.b.temp[i,1] < prospects.b.temp[i,3]){
    prospects.b[i,] <- prospects.b.temp[i,] 
  }else{
    prospects.b[i,1:2] <- prospects.b.temp[i,3:4] 
    prospects.b[i,3:4] <- prospects.b.temp[i,1:2] 
  }
}
prospects.a <- prospects.a[121:180,]
prospects.b <- prospects.b[121:180,]

# Load data (choice made by the first participant when presented the 
# second gamble-pair is saved in column 1 row 2)
rawdata <- as.matrix(read.table("Rieskamp_data.txt"))


## 2. Fit the hierarchical CPT-model                                ----
# Define what information that should be passed on to JAGS
data  = list("prospects.a", "prospects.b", "rawdata") 

# Define initial values for parameters 
inits = function() {
  list(mu.phi.alpha = 0.7, sigma.phi.alpha = 1, 
       mu.phi.beta = 0.7, sigma.phi.beta = 1,
       mu.phi.gamma = 0.7, sigma.phi.gamma = 1, 
       mu.phi.delta = 0.7, sigma.phi.delta = 1,
       lmu.lambda = 0, lsigma.lambda = 0.5, 
       lmu.luce = 0, sigma.phi.luce = 0.5) 
}


# Define the variables of interest. JAGS will return these to R when 
# the analysis is finished (and JAGS is closed).	
parameters = c("alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha", "mu.alpha_sebi",
               "beta", "mu.phi.beta", "mu.beta", "sigma.phi.beta", "mu.beta_sebi",
               "gamma", "mu.phi.gamma", "mu.gamma", "sigma.phi.gamma", "mu.gamma_sebi",
               "delta", "mu.phi.delta", "mu.delta", "sigma.phi.delta", "mu.delta_sebi",
               "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", "mu.luce_sebi",
               "luce", "lmu.luce", "mu.luce", "lsigma.luce", "mu.lambda_sebi")

## To prevent re-fitting when save results are present
if (!file.exists("saved_details/Refitted_Data.RData")) {
  res_rieskamp_1 =  jags.parallel(data,
                                  parameters,  model.file = original_full_model,
                                  inits = inits,
                                  n.chains = 4, n.iter = 20000, n.burnin = 1000, n.thin = 10,
                                  n.cluster = 4, jags.seed = 531)
  save(res_rieskamp_1, file="saved_details/Refitted_Data.RData")
}

load("saved_details/Refitted_Data.RData")


## 3. Compare population means between transformations              ----
temp_summary <- res_rieskamp_1$BUGSoutput$summary
#max(res_rieskamp_1$BUGSoutput$summary[,"Rhat"])
parname <- rownames(temp_summary)
temp_summary <- as_tibble(temp_summary) %>% mutate(parname = parname)
group_pars_summary <- temp_summary %>% 
  filter(grepl(parname, pattern = "mu"))

pd <- position_dodge(width=0.2)

two_colors_transformations <- c("#1b9e77", "#fc8d62")
# color_values <- palette.colors(n = 3, palette = "R4")[-1]
# color_values <- palette.colors(n = 3, palette = "Okab")[-1] %>% as.character()
# color_values <- palette.colors(n = 2, palette = "Okab") %>% as.character()
group_pars_summary %>% 
  filter(!grepl("phi", parname) & !grepl("lmu", parname)) %>%
  mutate(Transformation = ifelse(grepl("sebi", parname), "Correct", "Original"), 
         Parameter = sub("_sebi", "", sub("mu.", "", parname)),
         Parameter = ifelse(Parameter=="luce", "varphi", Parameter)) %>%
  ggplot(aes(x=Parameter, color=Transformation))+
  scale_color_manual(values=two_colors_transformations)+
  geom_point(aes(y=`50%`), size=3, position=pd)+
  scale_x_discrete(labels = scales::parse_format())+
  ylab("Posterior Median (95%CI)")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
  custom_theme
ggsave("figures/Rieskamp_Original.eps",
       width = 17.62, height=9/0.7, units="cm",dpi=600, device = cairo_ps)



# B  Refitting Rieskamp-data with model with correlated random effs ----
## 1. Read experimental data                                        ----
## 2. Fit the hierarchical CPT-model                                ----
## 3. Compare population means between transformations              ----
## 4. Check correlation between parameters                          ----
# C  Re-do (Extended) original simulation study (unconstrained)     ----

## 1. Subset data and define JAGS settings                          ----
mixed_prospects.a <- prospects.a[121:180,]
mixed_prospects.b <- prospects.b[121:180,]
data  = list("mixed_prospects.a", "mixed_prospects.b", "Data", "cur_n") 

# Define initial values for parameter
inits = function() {
  list(mu.phi.alpha = 0.7, sigma.phi.alpha = 1,
       mu.phi.beta = 0.7, sigma.phi.beta = 1,
       mu.phi.gamma = 0.7, sigma.phi.gamma = 1, 
       mu.phi.delta = 0.7, sigma.phi.delta = 1,
       lmu.lambda = 0, lsigma.lambda = 0.5, 
       lmu.luce = 0, sigma.phi.luce = 0.5)
}


# Define the variables of interest. JAGS will return these to R 
# when the analysis is finished (and JAGS is closed).	
parameters = c("alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha", "mu.alpha_sebi",
               "beta", "mu.phi.beta", "mu.beta", "sigma.phi.beta", "mu.beta_sebi",
               "gamma", "mu.phi.gamma", "mu.gamma", "sigma.phi.gamma", "mu.gamma_sebi",
               "delta", "mu.phi.delta", "mu.delta", "sigma.phi.delta", "mu.delta_sebi",
               "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", "mu.luce_sebi",
               "luce", "lmu.luce", "mu.luce", "lsigma.luce", "mu.lambda_sebi"
)

## 2. Actual parameter recovery analysis                            ----

## Set mean parameters for simulation
alpha <- beta <- .88
gamma <- .61 
delta <- .69
lambda <- 2.25

## Define the different settings that should be compared
phis <- c(.04, .14, .40) # choice sensitivity
Nsbjs <- c(20, 50, 90) # number of subjects
variabilities <- c(0.1, 0.5, 1) # btw-sbj variability in parameters

## Actually do the simulation, save simulation, and model fitting
## Only do this, when all analysis should be done again (takes long!)
if (REDOALLANALYSIS) {
  dir.create("saved_details/Recovery_full", showWarnings = FALSE)
  N <- VAR <- PHI <- 1
  for (N  in 1:3) {
    cur_n <- Nsbjs[N]
    Data <- matrix(NA, nrow=60, ncol=cur_n)
    for (VAR in 1:3) {
      cur_var <- variabilities[VAR]
      for ( PHI in 1:3) {
        if (PHI==1 && VAR==1 && N==1) next
        cur_sens <- phis[PHI]
        ## Only if the saved samples do not already exists
        if (!file.exists(paste0("saved_details/Recovery_full/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))) {
          ## Make it reproducible
          seeeed <- 2201 + 100*N + 10*VAR + PHI 
          set.seed(seeeed)
          
          ## Sample from Beta-distribution with mean alpha and scaled variance cur_var (not exactly the variance!)
          Alphas <- rbeta(cur_n, alpha*((alpha*(1-alpha))/cur_var *20 -1), (1-alpha)*((alpha*(1-alpha))/cur_var *20 -1) )
          Betas  <- rbeta(cur_n,  alpha*((alpha*(1-alpha))/cur_var *20 -1), (1-alpha)*((alpha*(1-alpha))/cur_var *20 -1) )
          Gammas <- rbeta(cur_n, gamma*((gamma*(1-gamma))/cur_var *10 -1), (1-gamma)*((gamma*(1-gamma))/cur_var *10 -1) )
          Deltas <- rbeta(cur_n, delta*((delta*(1-delta))/cur_var *10 -1), (1-delta)*((delta*(1-delta))/cur_var *10 -1) )
          # Draw from Gamma distribution with mean lambda and variance cur_var
          Lambdas <- rgamma(cur_n, shape= lambda^2/cur_var , scale=cur_var/lambda)
          
          for (k in 1:cur_n) {
            Data[,k] <- simulate_CPT_individ(Alphas[k], Betas[k], Gammas[k], Deltas[k], Lambdas[k], cur_sens)
          }
          params <- data.frame(alpha=Alphas, beta=Betas, gamma=Gammas, delta=Deltas, lambda=Lambdas, phi=cur_sens)
          simulation_pars <- list(N = cur_n, var=cur_var, sens=cur_sens)
          save(Data, params, simulation_pars,
               file=paste0("saved_details/Recovery_full/SampledData_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
          
          rec_samples =  jags.parallel(data, parameters,
                                       model.file = original_full_model_recovery,
                                       inits = inits,  n.chains = 4,
                                       n.iter = 50000, n.burnin = 1000,
                                       n.thin = 5,  n.cluster = 4, jags.seed = seeeed)
          
          save(Data, params, simulation_pars, rec_samples, 
               file=paste0("saved_details/Recovery_full/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
          
        }
      }
    }
  }
}

## When the fitting is done, we load the results and combine the simulations
if (!file.exists("saved_details/Recovery_full/Collected_recovery_results.RData")) {
  collected_samples <- data.frame()
  collected_summaries <- data.frame()
  collected_true_pop_means <- data.frame()
  for (N  in 1:3) {
    cur_n <- Nsbjs[N]
    for (VAR in 1:3) {
      cur_var <- variabilities[VAR]
      for ( PHI in 1:3) {
        cur_sens <- phis[PHI]
        if (file.exists(paste0("saved_details/Recovery_full/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))) {
          ## Load fit results
          load(paste0("saved_details/Recovery_full/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
          
          ## Combine the whole posterior samples of population parameters
          temp <- rec_samples$BUGSoutput$sims.array[,, (cur_n * 6 + 5 + 1:12)]    
          par_names <- dimnames(temp)[[3]]
          dim(temp) <- c(dim(temp)[1]*dim(temp)[2], dim(temp)[3])
          colnames(temp) <- par_names      
          temp <- as.data.frame(temp) 
          head(temp)
          temp <- cbind(temp, as.data.frame(simulation_pars))
          collected_samples <- rbind(collected_samples, temp)
          
          ## Combine the posterior summaries of population parameters
          temp <- rec_samples$BUGSoutput$summary[(cur_n * 6 + 5 + 1:12),]
          temp <- temp %>% as.data.frame() %>%
            select(c(1,2,3,5,7)) %>% 
            rownames_to_column("parname") 
          temp <- cbind(temp, as.data.frame(simulation_pars))
          collected_summaries <- rbind(collected_summaries, temp)
          
          ## Combine actual sampled population means
          load(paste0("saved_details/Recovery_full/SampledData_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
          temp <- colMeans(params) %>% data.frame()  %>% 
            rownames_to_column("Parameter")
          colnames(temp)[2] <- "value"
          temp <- cbind(temp, as.data.frame(simulation_pars))
          collected_true_pop_means <- rbind(collected_true_pop_means, temp)
          
        }
      }
    }
  }
  ## Clean and Format Parameter Labels
  collected_samples <- collected_samples %>% 
    #filter(!grepl("phi", parname) & !grepl("lmu", parname)) %>%
    pivot_longer(1:12, names_to="parname") %>%
    mutate(Computation = ifelse(grepl("sebi", parname), "Sebis", "Original"), 
           Parameter = sub("_sebi", "", sub("mu.", "", parname)))
  collected_summaries <- collected_summaries %>% 
    mutate(Computation = ifelse(grepl("sebi", parname), "Sebis", "Original"), 
           Parameter = sub("_sebi", "", sub("mu.", "", parname)))
  
  save(collected_samples,collected_summaries, collected_true_pop_means, 
       file="saved_details/Recovery_full/Collected_recovery_results.RData")
} else {
  load("saved_details/Recovery_full/Collected_recovery_results.RData")
}

## 3. Visualize original full parameter recovery analysis           ----
## Reproduce Nilsson et al. (2011), Figure 2:
# Note: variability in Nilsson et al is 0; and N = 30; but the following are
# the values most close to those in Nilsson's paper:
plot_samples <- filter(collected_samples, sens==0.4 & var%in%c(0.1, 1) & N == 20) %>%
  filter(Parameter != "luce")
true_params <- data.frame(Parameter= c("alpha", "beta", "gamma","delta","lambda"), 
                          value    = c(   .88,    .88,       .61,    .69,  2.25))
ggplot(plot_samples, aes(x=value, linetype=as.factor(var), color=Computation))+
  geom_vline(data=true_params, aes(xintercept=value))+
  geom_density(aes(group=interaction(Computation, Parameter, var)), linewidth=1)+
  facet_wrap(.~Parameter, scales = "free")

true_params <- data.frame(Parameter= c("alpha", "beta", "gamma","delta","lambda"), 
                          value    = c(   .88,    .88,       .61,    .69,  2.25))
pd <- position_dodge(width=0.4)
ggplot(filter(collected_summaries,Parameter!="luce"), aes(y=mean, x=interaction(N,var), color=Computation))+
  geom_hline(data=true_params, aes(yintercept=value))+
  geom_point(position=pd)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position=pd, width=0.2)+
  facet_grid(Parameter~sens, scales = "free")


pd <- position_dodge(width=0.4)
# Only take the extreme sampling options for each factor
sub_results <- subset(collected_summaries,Parameter!="luce") %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1)) %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens))
sub_pop_means <- collected_true_pop_means %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1) &
           Parameter != "phi") %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens)) %>%
  merge(data.frame(Computation=c("Original", "Sebis")))
pd <- position_dodge(width=0.2)
ggplot(sub_results,
       aes(y=`50%`, x=as.factor(N), color=Computation))+
  geom_hline(data=true_params, aes(yintercept=value))+
  geom_errorbar(data=sub_pop_means , aes(ymin=value, y=value,ymax=value), linetype="dashed")+
  geom_point(position=pd)+
  geom_line(aes(group=Computation),position=pd)+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
  facet_nested(Parameter~var+sens, scales = "free", labeller = label_parsed )+
  labs(y="Parameter values", x="Simulated Sample Size x Sensitivity")+
  theme_bw()


# 
# sub_results <- subset(collected_summaries,Parameter!="luce") %>%
#   filter(sens %in% c(0.04, 0.4) &
#            var %in% c(0.1, 1)) %>%
#   mutate(var=paste0("Variability: ", var),
#          sens=paste0("Sensitivity: ", sens))
# sub_pop_means <- collected_true_pop_means %>%
#   filter(sens %in% c(0.04, 0.4) &
#            var %in% c(0.1, 1) &
#            Parameter != "phi") %>%
#   mutate(var=paste0("Variability: ", var),
#          sens=paste0("Sensitivity: ", sens)) %>%
#   merge(data.frame(Computation=c("Original", "Sebis")))
# pd <- position_dodge(width=0.2)
# ggplot(sub_results,
#        aes(y=`50%`, x=as.factor(N), color=Computation))+
#   geom_hline(data=true_params, aes(yintercept=value))+
#   geom_errorbar(data=sub_pop_means , aes(ymin=value, y=value,ymax=value), linetype="dashed")+
#   geom_point(position=pd)+
#   geom_line(aes(group=Computation),position=pd)+
#   geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
#   facet_nested(Parameter~sens+var, scales = "free", labeller = label_parsed )+
#   labs(y="Parameter values", x="Simulated Sample Size x Sensitivity")+
#   theme_bw()
# 


# D  Re-do (Extended) original simulation study (alpha=beta)        ----
## 1. Subset data and define JAGS settings                          ----
mixed_prospects.a <- prospects.a[121:180,]
mixed_prospects.b <- prospects.b[121:180,]
data  = list("mixed_prospects.a", "mixed_prospects.b", "Data", "cur_n") 

# Define initial values for parameter
inits = function() {
  list(mu.phi.alpha = 0.7, sigma.phi.alpha = 1,
       mu.phi.gamma = 0.7, sigma.phi.gamma = 1, 
       mu.phi.delta = 0.7, sigma.phi.delta = 1,
       lmu.lambda = 0, lsigma.lambda = 0.5, 
       lmu.luce = 0, sigma.phi.luce = 0.5)
}


# Define the variables of interest. JAGS will return these to R 
# when the analysis is finished (and JAGS is closed).	
parameters = c("alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha", "mu.alpha_sebi",
               "gamma", "mu.phi.gamma", "mu.gamma", "sigma.phi.gamma", "mu.gamma_sebi",
               "delta", "mu.phi.delta", "mu.delta", "sigma.phi.delta", "mu.delta_sebi",
               "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", "mu.luce_sebi",
               "luce", "lmu.luce", "mu.luce", "lsigma.luce", "mu.lambda_sebi"
)


## 2. Actual parameter recovery analysis                            ----
## Set mean parameters for simulation
alpha <- .88
gamma <- .61 
delta <- .69
lambda <- 2.25

## Define the different settings that should be compared
phis <- c(.04, .14, .40) # choice sensitivity
Nsbjs <- c(20, 50, 90) # number of subjects
variabilities <- c(0.1, 0.5, 1) # btw-sbj variability in parameters

## Actually do the simulation, save simulation, and model fitting
## Only do this, when all analysis should be done again (takes long!)
if (REDOALLANALYSIS) {
  dir.create("saved_details/Recovery_restricted", showWarnings = FALSE)
  N <- VAR <- PHI <- 1
  for (N  in 1:3) {
    cur_n <- Nsbjs[N]
    Data <- matrix(NA, nrow=60, ncol=cur_n)
    for (VAR in 1:3) {
      cur_var <- variabilities[VAR]
      for ( PHI in 1:3) {
        cur_sens <- phis[PHI]
        seeeed <- 2201 + 100*N + 10*VAR + PHI 
        set.seed(seeeed)
        ## Sample from Beta-distribution with mean alpha and scaled variance cur_var (not exactly the variance!)
        Alphas <- rbeta(cur_n, alpha*((alpha*(1-alpha))/cur_var *20 -1), (1-alpha)*((alpha*(1-alpha))/cur_var *20 -1) )
        Gammas <- rbeta(cur_n, gamma*((gamma*(1-gamma))/cur_var *10 -1), (1-gamma)*((gamma*(1-gamma))/cur_var *10 -1) )
        Deltas <- rbeta(cur_n, delta*((delta*(1-delta))/cur_var *10 -1), (1-delta)*((delta*(1-delta))/cur_var *10 -1) )
        # Draw from Gamma distribution with mean lambda and variance cur_var
        Lambdas <- rgamma(cur_n, shape= lambda^2/cur_var , scale=cur_var/lambda)
        
        for (k in 1:cur_n) {
          Data[,k] <- simulate_CPT_individ(Alphas[k], Alphas[k], Gammas[k], Deltas[k], Lambdas[k], cur_sens)
        }
        params <- data.frame(alpha=Alphas, gamma=Gammas, delta=Deltas, lambda=Lambdas, phi=cur_sens)
        simulation_pars <- list(N = cur_n, var=cur_var, sens=cur_sens)
        save(Data, params, simulation_pars,
             file=paste0("saved_details/Recovery_restricted/SampledData_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
        
        rec_samples =  jags.parallel(data, parameters,
                                     model.file = original_restricted_model_recovery,
                                     inits = inits,  n.chains = 4,
                                     n.iter = 50000, n.burnin = 1000,
                                     n.thin = 5,  n.cluster = 4, jags.seed = seeeed)
        
        
        save(Data, params, simulation_pars, rec_samples, 
             file=paste0("saved_details/Recovery_restricted/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
        
      }
    }
  }
}

## When the fitting is done, we load the results and combine the simulations
if (!file.exists("saved_details/Recovery_restricted/Collected_recovery_results_restricted.RData")) {
  collected_samples_restricted <- data.frame()
  collected_summaries_restricted <- data.frame()
  collected_true_pop_means_restricted <- data.frame()
  for (N  in 1:3) {
    cur_n <- Nsbjs[N]
    for (VAR in 1:3) {
      cur_var <- variabilities[VAR]
      for ( PHI in 1:3) {
        cur_sens <- phis[PHI]
        if (file.exists(paste0("saved_details/Recovery_restricted/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))) {
          ## Load fit results
          load(paste0("saved_details/Recovery_restricted/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
          
          ## Combine the whole posterior samples of population parameters
          temp <- rec_samples$BUGSoutput$sims.array[,, (cur_n * 5 + 5 + 1:10)]    
          par_names <- dimnames(temp)[[3]]
          dim(temp) <- c(dim(temp)[1]*dim(temp)[2], dim(temp)[3])
          colnames(temp) <- par_names      
          temp <- as.data.frame(temp) 
          head(temp)
          temp <- cbind(temp, as.data.frame(simulation_pars))
          collected_samples_restricted <- rbind(collected_samples_restricted, temp)
          
          ## Combine the posterior summaries of population parameters
          temp <- rec_samples$BUGSoutput$summary[(cur_n * 5 + 5 + 1:10),]
          temp <- temp %>% as.data.frame() %>%
            select(c(1,2,3,5,7)) %>% 
            rownames_to_column("parname") 
          temp <- cbind(temp, as.data.frame(simulation_pars))
          collected_summaries_restricted <- rbind(collected_summaries_restricted, temp)
          
          ## Combine actual sampled population means
          load(paste0("saved_details/Recovery_restricted/SampledData_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
          temp <- colMeans(params) %>% data.frame()  %>% 
            rownames_to_column("Parameter")
          colnames(temp)[2] <- "value"
          temp <- cbind(temp, as.data.frame(simulation_pars))
          collected_true_pop_means_restricted <- rbind(collected_true_pop_means_restricted, temp)
          
        }
      }
    }
  }
  ## Clean and Format Parameter Labels
  collected_samples_restricted <- collected_samples_restricted %>% 
    #filter(!grepl("phi", parname) & !grepl("lmu", parname)) %>%
    pivot_longer(1:10, names_to="parname") %>%
    mutate(Computation = ifelse(grepl("sebi", parname), "Sebis", "Original"), 
           Parameter = sub("_sebi", "", sub("mu.", "", parname)))
  collected_summaries_restricted <- collected_summaries_restricted %>% 
    mutate(Computation = ifelse(grepl("sebi", parname), "Sebis", "Original"), 
           Parameter = sub("_sebi", "", sub("mu.", "", parname)))
  
  save(collected_samples_restricted,collected_summaries_restricted, collected_true_pop_means_restricted, 
       file="saved_details/Recovery_restricted/Collected_recovery_results_restricted.RData")
} else {
  load("saved_details/Recovery_restricted/Collected_recovery_results_restricted.RData")
}

## 3. Visualize original full parameter recovery analysis           ----
## Reproduce Nilsson et al. (2011), Figure 2:
# Note: variability in Nilsson et al is 0; and N = 30; but the following are
# the values most close to those in Nilsson's paper:
plot_samples <- filter(collected_samples_restricted, sens==0.4 & var==0.1 & N == 20) %>%
  filter(Parameter != "luce")
true_params <- data.frame(Parameter= c("alpha", "gamma","delta","lambda"), 
                          value    = c(   .88,     .61,    .69,  2.25))
ggplot(plot_samples, aes(x=value, linetype=Computation, color=Computation))+
  geom_vline(data=true_params, aes(xintercept=value))+
  geom_density(aes(group=interaction(Computation, Parameter)))+
  facet_wrap(.~Parameter, scales = "free")

## Extend Nilsson et al. (2011), Figure 2:
# Note: variability in Nilsson et al is 0; and N = 30; but the following are
# the values most close to those in Nilsson's paper:
plot_samples <- filter(collected_samples_restricted, sens==0.4 & var%in%c(0.1, 1) & N == 20) %>%
  filter(Parameter != "luce") %>%
  mutate(var=paste0("Variability: ", var))
true_params <- data.frame(Parameter= c("alpha","gamma","delta","lambda"), 
                          value    = c(   .88,     .61,    .69,  2.25))
p1<- ggplot(subset(plot_samples, Parameter!="lambda"),
            aes(x=value, color=Computation))+
  geom_vline(data=subset(true_params, Parameter!="lambda"), aes(xintercept=value))+
  geom_density(aes(group=interaction(Computation, Parameter, var)), linewidth=1)+
  facet_grid2(var~Parameter, scales = "free", labeller = label_parsed)+
  theme(strip.text.y = element_blank(), strip.background.y = element_blank())
p2 <- ggplot(subset(plot_samples, Parameter=="lambda"),
             aes(x=value, color=Computation))+
  geom_vline(data=subset(true_params, Parameter=="lambda"), aes(xintercept=value))+
  geom_density(aes(group=interaction(Computation, Parameter, var)), linewidth=1)+
  xlim(c(1, 3.7))+theme(axis.title.y = element_blank())+
  facet_grid2(var~Parameter, scales = "free", labeller = label_parsed)
ggarrange(p1, p2, 
          widths = c(0.7, 0.3), nrow=1, common.legend = TRUE, legend="bottom")


true_params <- data.frame(Parameter= c("alpha","gamma","delta","lambda"), 
                          value    = c(   .88,    .61,    .69,  2.25))
pd <- position_dodge(width=0.4)
ggplot(filter(collected_summaries_restricted,Parameter!="luce"), aes(y=mean, x=interaction(N,var), color=Computation))+
  geom_hline(data=true_params, aes(yintercept=value))+
  geom_point(position=pd)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position=pd, width=0.2)+
  facet_grid(Parameter~sens, scales = "free")


pd <- position_dodge(width=0.4)
# Only take the extreme sampling options for each factor
sub_results <- subset(collected_summaries_restricted,Parameter!="luce") %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1)) %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens))
sub_pop_means <- collected_true_pop_means_restricted %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1) &
           Parameter != "phi") %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens)) %>%
  merge(data.frame(Computation=c("Original", "Sebis")))
pd <- position_dodge(width=0.2)
ggplot(sub_results,
       aes(y=`50%`, x=as.factor(N), color=Computation))+
  geom_hline(data=true_params, aes(yintercept=value))+
  geom_errorbar(data=sub_pop_means , aes(ymin=value, y=value,ymax=value), linetype="dashed")+
  geom_point(position=pd)+
  geom_line(aes(group=Computation),position=pd)+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
  facet_nested(Parameter~var+sens, scales = "free", labeller = label_parsed )+
  labs(y="Parameter values", x="Simulated Sample Size x Sensitivity")+
  theme_bw()





sub_results <- subset(collected_summaries_restricted,Parameter!="luce") %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1)) %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens))
sub_pop_means <- collected_true_pop_means_restricted %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1) &
           Parameter != "phi") %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens)) %>%
  merge(data.frame(Computation=c("Original", "Sebis")))
pd <- position_dodge(width=0.2)
ggplot(sub_results,
       aes(y=`50%`, x=as.factor(N), color=Computation))+
  geom_hline(data=true_params, aes(yintercept=value))+
  geom_errorbar(data=sub_pop_means , aes(ymin=value, y=value,ymax=value), linetype="dashed")+
  geom_point(position=pd)+
  geom_line(aes(group=Computation),position=pd)+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
  facet_nested(Parameter~sens+var, scales = "free", labeller = label_parsed )+
  labs(y="Parameter values", x="Simulated Sample Size")+
  custom_theme





differences_df <- collected_samples_restricted %>% 
  mutate(iter=1:n(),.by=c(N, var, sens, parname, Parameter, Computation)) %>% select(-parname) %>% 
  pivot_wider(names_from="Computation", values_from = value) %>%
  mutate(difference = Sebis-Original) %>%
  group_by(N, sens, var, Parameter) %>%
  reframe(Med=median(difference), 
          lower = quantile(difference, probs = 0.025),
          upper = quantile(difference, probs = 0.975)) 
plot_differences_df <- differences_df %>% filter(sens==0.4) %>%
  mutate(Parameter = ifelse(Parameter=="luce", "varphi", Parameter))
ggplot(plot_differences_df, aes(x=as.factor(var), y=Med))+
  geom_point()+geom_line(aes(group=1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  facet_nested(Parameter~"Sample~Size"+N, scales = "free_y", labeller=label_parsed)+ 
  labs(x="Variability between individuals", y="Differences in the mean estimates")+
  custom_theme

# E  Simulation study using wider parameter ranges (alpha=beta)  -----
# F  Simulation study using wider parameter ranges + correlated pars ----