############################################################-
#####          Main script for the analyses in              -
#####                  Manuscript:                          -
#####  Bias by Variance: How Commonly Used Computations To  -
#####     Constrain Parameters in Hierarchical Modeling     -
#####                 Distort Estimation                    -


# Use Alt+o in RStudio to collapse all folds!

# Sebastian Hellmann, June 2025
rm(list = ls())
REDOALLANALYSIS <- FALSE

#________       Structure of the script         ____________----
# Preamble and imports    
#___________________________________________________________________
# A  Read and prepare experimental data define JAGS inputs          
## 1. Rieskamp (2008) data and gambles                              
## 2. Pachur et al (2017) age data                                  
#___________________________________________________________________
#______    Re-doing analysis of Nilsson et al (2011)       _________
#___________________________________________________________________
# B  Refitting Rieskamp-data with original model (alpha=beta)       
## 1. Fit the hierarchical CPT-model                                
## 2. Compare population means between transformations              
#___________________________________________________________________
# C  Re-do (Extended) original simulation study (alpha=beta)        
## 1. Actual parameter recovery analysis                            
## 2. Visualize original restricted parameter recovery analysis     
#___________________________________________________________________
# D  Re-doing age difference analysis in Pachur et al (2017)  ______
#___________________________________________________________________
## 1. Fit the hierarchical CPT-model                                
## 2. Compare means between young and old                           
#___________________________________________________________________
#_______                 For Supplement                     ________
#___________________________________________________________________
# E  Refitting Rieskamp-data with original model                    
## 1. Fit the hierarchical CPT-model                                
## 2. Compare population means between transformations              
#___________________________________________________________________
# F  Re-do (Extended) original simulation study (unconstrained)     
## 1. Actual parameter recovery analysis                            
## 2. Visualize original full parameter recovery analysis           










# Preamble and imports                                     ----

# use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
print("Working directory set to:")
print(getwd())

{
  # Tell Rstudio where to find JAGS
  #Sys.setenv(JAGS_HOME = "C:/Users/go73jec/AppData/Local/Programs/JAGS/JAGS-4.3.1")
  library(tidyverse)
  library(R2jags)
  library(ggpubr)
  library(viridis)
  library(ggh4x)
  library(cowplot)
  library(tensr)
  library(readxl)
  library(kableExtra)# For the first table of posteriors
  library(xtable)    # For the second table of posteriors
  source('helper_fcts/custom_theme.R') # import custom ggplot theme
  dir.create("figures", showWarnings = FALSE)
  dir.create("saved_details", showWarnings = FALSE)
  
  par_names <- c("alpha", "beta", "gamma.loss", "gamma.gain", "lambda", "sens")
  par_labels <- c("alpha","beta",  "gamma^'-'", "gamma^'+'","lambda", "phi" )
}
## Import simulation function and define JAGS model file names
source("helper_fcts/simulate_CPT.R")
original_restricted_model <- "jags_models/cpt_hierarchical_restricted.txt"
original_full_model <- "jags_models/cpt_hierarchical_model.txt"

original_full_model_recovery <- "jags_models/cpt_hierarchical_recovery.txt"
original_restricted_model_recovery <- "jags_models/cpt_hierarchical_restricted_recovery.txt"

Pachur_age_model <- "jags_models/cpt_hierarchical_age_model.txt"


## Generate Figure 1 and 2 for the theoretical part
source("helper_fcts/Fig1_transformation_viz.R")
source("helper_fcts/Fig2_theoretical_bias_viz.R")

#___________________________________________________________________----
# A  Read and prepare experimental data define JAGS inputs          ----

## 1. Rieskamp (2008) data and gambles                              ----

# Load information about the gamble-pairs used in Rieskamp (2008). 
# GambleA.txt and GambleB.txt are structured as follows: 
# value of outcome 1 (column 1), 
# probability of outcome 1 (column 2), 
# value of outcome 2 (column 3), 
# probability of outcome 2 (column 4) (gambles in rows).
prospects.b.temp <- as.matrix(read.table("data/Rieskamp_2008_data/GambleB.txt"))
prospects.a.temp <- as.matrix(read.table("data/Rieskamp_2008_data/GambleA.txt"))

prospects.b <- array(0,dim=c(180,4))
prospects.a <- array(0,dim=c(180,4))

# Arrange so that v and p related to the relatively poor outcome ends 
# up in column 1 and 2
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

# Load data (choice made by the first participant when presented the 
# second gamble-pair is saved in column 1 row 2; 180 problems x 30 participants)
rawdata <- as.matrix(read.table("data/Rieskamp_2008_data/Rieskamp_data.txt"))


# Define what information that should be passed on to JAGS for the empirical data analysis
data <- list("prospects.a", "prospects.b", "rawdata") 


# Subset mixed gambles for the recovery study and define JAGS-relevant objects
mixed_prospects.a <- prospects.a[121:180,]
mixed_prospects.b <- prospects.b[121:180,]
simu_data  <- list("mixed_prospects.a", "mixed_prospects.b", "Data", "cur_n") 


## 2. Pachur et al (2017) age data                                  ----
# Read the data 
choice_data <- read_xlsx("data/PachurEtAl_Who errs, who dares_Data.xlsx",
                              sheet = "Choice task", range="B2:EB107")
## Bring data in correct format
gambles <- as.matrix(choice_data[,124:ncol(choice_data)])
lotteries_a <- gambles[,1:4]
lotteries_b <- gambles[,5:8]

## Ensure that the smaller outcome (and corresponding probability) is always left
## (lotteries_a[73:80,] entail only outcome 0 (o1=0,p1=.5,o2=0,p2=.5) )
for (i in 1:nrow(lotteries_a)) {
  if (lotteries_a[i,1] > lotteries_a[i,3]) {
    lotteries_a[i,] <- lotteries_a[i,c(3,4,1,2)]
  }
  if (lotteries_b[i,1] > lotteries_b[i,3]) {
    lotteries_b[i,] <- lotteries_b[i,c(3,4,1,2)]
  }
}
## Check order of the positive, negative, and mixed gambles
all(lotteries_a[1:41,] >= 0 & lotteries_b[1:41,] >=0)
all(lotteries_a[42:72, c(1, 3)] <= 0 & lotteries_b[42:72, c(1, 3)] <=0)
all(lotteries_a[73:105, 1] <= 0 & lotteries_a[73:105, 3] >= 0) #  a includes 0-0 outcomes
all(lotteries_b[73:105, 1] < 0 & lotteries_b[73:105, 3] > 0)

all_choices <- choice_data[,1:122] %>% as.matrix()

age_data <- read_xlsx("data/PachurEtAl_Who errs, who dares_Data.xlsx",
                      sheet = "Data", 
                      range="A1:B12811" # include to omit the NA warnings for Speed
                      )
age_data <- age_data %>% 
  select(sbj=Subject, group=Age_group) %>%
  distinct()
young_choices <- all_choices[,{
  age_data %>% filter(group=="younger") %>% pull("sbj")
  }]
older_choices <- all_choices[,{
  age_data %>% filter(group=="older") %>% pull("sbj")
}]

age_data <- list("lotteries_a", "lotteries_b", "age_choices", "N_parts")




#___________________________________________________________________----
#______    Re-doing analysis of Nilsson et al (2011)       _________----
#___________________________________________________________________----
# B  Refitting Rieskamp-data with original model (alpha=beta)       ----

## 1. Fit the hierarchical CPT-model                                ----

# Define initial values for parameters 
inits = function() {
  list(mu.phi.alpha = 0, sigma.phi.alpha = 1, 
       mu.phi.gamma.gain = 0, sigma.phi.gamma.gain = 1, 
       mu.phi.gamma.loss = 0, sigma.phi.gamma.loss = 1,       
       lmu.lambda = 0, lsigma.lambda = 0.5, 
       lmu.sens = 0, sigma.phi.sens = 0.5)
}


# Define the variables of interest. JAGS will return these to R when 
# the analysis is finished (and JAGS is closed).	
parameters = c("alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha", "mu.alpha_sebi",
               "gamma.gain", "mu.phi.gamma.gain", "mu.gamma.gain", "sigma.phi.gamma.gain", "mu.gamma.gain_sebi",
               "gamma.loss", "mu.phi.gamma.loss", "mu.gamma.loss", "sigma.phi.gamma.loss", "mu.gamma.loss_sebi",
               "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", "mu.sens_sebi",
               "sens", "lmu.sens", "mu.sens", "lsigma.sens", "mu.lambda_sebi")

## To prevent re-fitting when save results are present
if (!file.exists("saved_details/Refitted_Data.RData")) {
  res_rieskamp_restricted =  jags.parallel(data, parameters,
     model.file = original_restricted_model,
     inits = inits,
     n.chains = 5, n.iter = 70000,
     n.burnin = 6000, n.thin = 20,
     n.cluster = 5, jags.seed = 10042025)
  res_rieskamp_restricted <- list(samples=res_rieskamp_restricted$BUGSoutput$sims.array,
                                  summaries = res_rieskamp_restricted$BUGSoutput$summary)
  save(res_rieskamp_restricted, file="saved_details/Refitted_Data_restricted.RData")
}

load("saved_details/Refitted_Data_restricted.RData")

## 2. Compare population means between transformations              ----
temp_summary <- res_rieskamp_restricted$summaries
#max(res_rieskamp_restricted$BUGSoutput$summary[,"Rhat"])
parname <- rownames(temp_summary)
temp_summary <- as_tibble(temp_summary) %>% mutate(parname = parname)
group_pars_summary <- temp_summary %>% 
  filter(grepl(parname, pattern = "mu"))

pd <- position_dodge(width=0.2)
plt_group_pars_summary <- group_pars_summary %>% 
  filter(!grepl("phi", parname) & !grepl("lmu", parname)) %>%
  mutate(Computation = ifelse(grepl("sebi", parname), "Correct", "Incorrect"),
         Computation = factor(Computation, levels=c("Incorrect", "Correct")),
         Parameter = sub("_sebi", "", sub("mu.", "", parname)),
         Parameter = factor(Parameter, levels=par_names, labels=par_labels))

ggplot(plt_group_pars_summary, aes(x=Parameter, color=Computation))+
  scale_color_manual(values=two_colors_transformations)+
  geom_point(aes(y=`50%`), size=3, position=pd)+
  scale_x_discrete(labels = scales::parse_format())+
  ylab("Posterior Median (95%CI)")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
  custom_theme
ggsave("figures/Rieskamp_restricted.eps",
       width = 17.62, height=9/0.7, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/Rieskamp_restricted.png",
       width = 17.62, height=9/0.7, units="cm",dpi=900)

#___________________________________________________________________----
# C  Re-do (Extended) original simulation study (alpha=beta)        ----
## 1. Actual parameter recovery analysis                            ----

# Define initial values for parameter
inits = function() {
  list(mu.phi.alpha = 0.7, sigma.phi.alpha = 1,
       mu.phi.gamma.gain = 0.7, sigma.phi.gamma.gain = 1, 
       mu.phi.gamma.loss = 0.7, sigma.phi.gamma.loss = 1,
       lmu.lambda = 0, lsigma.lambda = 0.5, 
       lmu.sens = 0, sigma.phi.sens = 0.5)
}


# Define the variables of interest. JAGS will return these to R 
# when the analysis is finished (and JAGS is closed).	
parameters = c("alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha", "mu.alpha_sebi",
               "gamma.gain", "mu.phi.gamma.gain", "mu.gamma.gain", "sigma.phi.gamma.gain", "mu.gamma.gain_sebi",
               "gamma.loss", "mu.phi.gamma.loss", "mu.gamma.loss", "sigma.phi.gamma.loss", "mu.gamma.loss_sebi",
               "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", "mu.sens_sebi",
               "sens", "lmu.sens", "mu.sens", "lsigma.sens", "mu.lambda_sebi"
)



## Set mean parameters for simulation
alpha <- .88
gamma.gain <- .61 
gamma.loss <- .69
lambda <- 2.25

## Define the different settings that should be compared
phis <- c(.04, .14, .40) # choice sensitivity
Nsbjs <- c(20, 50, 90) # number of subjects
variabilities <- c(0.1, 0.5, 1) # btw-sbj variability in parameters

## Actually do the simulation, save simulation, and model fitting
## Only do this, when all analysis should be done again (takes long!)
if (REDOALLANALYSIS) {
  collected_samples_restricted <- data.frame()
  collected_summaries_restricted <- data.frame()
  collected_true_pop_means_restricted <- data.frame()
  getpars <-c("alpha", "gamma.loss", "gamma.gain", "lambda", "sens")
  getpars <- paste0("mu.", getpars, rep(c("", "_sebi"), each=length(getpars)))
  
  dir.create("saved_details/Recovery_restricted", showWarnings = FALSE)
  N <- VAR <- PHI <- 1
  for (N  in 1:3) { # for each subject ... 
    cur_n <- Nsbjs[N]
    Data <- matrix(NA, nrow=60, ncol=cur_n)
    for (VAR in 1:3) { # ... loop over all levels of variability ... 
      cur_var <- variabilities[VAR]
      for ( PHI in 1:3) { # ... and over all levels of choice sensitivity 
        cur_sens <- phis[PHI]
        seeeed <- 2201 + 100*N + 10*VAR + PHI 
        set.seed(seeeed)
        ## Sample from Beta-distribution with mean alpha and scaled variance cur_var (not exactly the variance!)
        Alphas <- rbeta(cur_n, alpha*((alpha*(1-alpha))/cur_var *20 -1), (1-alpha)*((alpha*(1-alpha))/cur_var *20 -1) )
        # check means and variance of the beta distribution (statements below should be true) 
        # shape_1 <- alpha*((alpha*(1-alpha))/cur_var *20 -1)
        # shape_2 <- (1-alpha)*((alpha*(1-alpha))/cur_var *20 -1)
        # mean <- shape_1/(shape_1+shape_2)
        # variance <- shape_1*shape_2/((shape_1+shape_2)^2 * (shape_1+shape_2+1))
        # mean == alpha
        # round(variance,2) == round(cur_var/20,2)
        Gammas.gain <- rbeta(cur_n, 
                             gamma.gain*((gamma.gain*(1-gamma.gain))/cur_var *10 -1), 
                             (1-gamma.gain)*((gamma.gain*(1-gamma.gain))/cur_var *10 -1) )
        Gammas.loss <- rbeta(cur_n, 
                             gamma.loss*((gamma.loss*(1-gamma.loss))/cur_var *10 -1), 
                             (1-gamma.loss)*((gamma.loss*(1-gamma.loss))/cur_var *10 -1) )
        # Draw from Gamma distribution with mean lambda and variance cur_var
        Lambdas <- rgamma(cur_n, shape= lambda^2/cur_var , scale=cur_var/lambda)
        
        for (k in 1:cur_n) {
          Data[,k] <- simulate_CPT_individ(Alphas[k], Alphas[k], Gammas.gain[k], Gammas.loss[k], Lambdas[k], cur_sens)
        }
        params <- data.frame(alpha=Alphas, gamma.gain=Gammas.gain, gamma.loss=Gammas.loss, lambda=Lambdas, phi=cur_sens)
        simulation_pars <- list(N = cur_n, var=cur_var, sens=cur_sens)
        save(Data, params, simulation_pars,
             file=paste0("saved_details/Recovery_restricted/SampledData_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
        
        rec_samples =  jags.parallel(simu_data, parameters,
                                     model.file = original_restricted_model_recovery,
                                     inits = inits,  n.chains = 4,
                                     n.iter = 50000, n.burnin = 1000,
                                     n.thin = 5,  n.cluster = 4, jags.seed = seeeed)
        
        rec_summary <- rec_samples$BUGSoutput$summary
        rec_samples <- rec_samples$BUGSoutput$sims.array
        save(Data, params, simulation_pars, rec_samples, rec_summary,
             file=paste0("saved_details/Recovery_restricted/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
        
        ### Put everything in common data frames
        
        ## Combine the whole posterior samples of population parameters
        temp <- rec_samples[,, getpars]    
        dim(temp) <- c(dim(temp)[1]*dim(temp)[2], dim(temp)[3])
        colnames(temp) <- getpars     
        temp <- as.data.frame(temp) 
        #head(temp)
        temp <- cbind(temp, as.data.frame(simulation_pars))
        collected_samples_restricted <- rbind(collected_samples_restricted, temp)
        
        ## Combine the posterior summaries of population parameters
        temp <- rec_summary[getpars,]
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
  ## Clean and Format Parameter Labels
  collected_samples_restricted <- collected_samples_restricted %>% 
    #filter(!grepl("phi", parname) & !grepl("lmu", parname)) %>%
    pivot_longer(1:10, names_to="parname") %>%
    mutate(Computation = ifelse(grepl("sebi", parname), "Correct", "Incorrect"), 
           Computation = factor(Computation, levels=c("Incorrect", "Correct")),
           Parameter = sub("_sebi", "", sub("mu.", "", parname)))
  collected_summaries_restricted <- collected_summaries_restricted %>% 
    mutate(Computation = ifelse(grepl("sebi", parname), "Correct", "Incorrect"), 
           Computation = factor(Computation, levels=c("Incorrect", "Correct")),
           Parameter = sub("_sebi", "", sub("mu.", "", parname)))
  
  save(collected_samples_restricted,collected_summaries_restricted, collected_true_pop_means_restricted, 
       file="saved_details/Collected_recovery_results_restricted.RData")
}

load("saved_details/Collected_recovery_results_restricted.RData")

## 2. Visualize original restricted parameter recovery analysis     ----

# Only take the extreme sampling options for each factor (lowest/highest variance and sensitivity)
true_params <- data.frame(Parameter= c("alpha","gamma.gain","gamma.loss","lambda"),
                          value    = c(   .88,     .61,    .69,  2.25)) %>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))

sub_results <- subset(collected_summaries_restricted,Parameter!="sens") %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1)) %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens),
         Parameter = factor(Parameter, levels=par_names, labels=par_labels))

sub_pop_means <- collected_true_pop_means_restricted %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1) &
           Parameter != "phi") %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens),
         Parameter = factor(Parameter, levels=par_names, labels=par_labels)) %>%
  merge(data.frame(Computation=c("Incorrect", "Correct")))

pd <- position_dodge(width=0.2)
ggplot(sub_results,
       aes(y=`50%`, x=as.factor(N), color=Computation))+
  geom_hline(data=subset(true_params,Parameter!="beta"), aes(yintercept=value))+
  geom_errorbar(data=sub_pop_means , aes(ymin=value, y=value,ymax=value), linetype="dashed", color="gray20")+
  geom_point(position=pd)+
  geom_line(aes(group=Computation),position=pd)+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
  scale_color_manual(values=two_colors_transformations)+
  facet_nested(Parameter~var+sens, scales = "free", labeller = label_parsed , drop = TRUE)+
  labs(y="Parameter values", x="Simulated sample size")+
  custom_theme+
  theme(panel.spacing= unit(0.1, "cm"))
ggsave("figures/Recovery_restricted_posteriorCIs.eps",
       width = 17.62, height=22.62, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/Recovery_restricted_posteriorCIs.png",
       width = 17.62, height=22.62, units="cm",dpi=900)


## Extend Nilsson et al. (2011), Figure 2:
# Note: variability in Nilsson et al is 0; and N = 30; but the following are
# the values most close to those in Nilsson's paper:
plot_samples <- filter(collected_samples_restricted, sens==0.4 & var%in%c(0.1, 1) & N == 50) %>%
  #mutate(Parameter =ifelse(Parameter!= "sens", Parameter, "phi")) %>%
  filter(Parameter!= "sens") %>%
  mutate(var=factor(var))%>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))
true_params <- data.frame(Parameter= c("alpha","gamma.gain","gamma.loss","lambda"),
                          value    = c(   .88,     .61,    .69,  2.25)) %>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))

sub_pop_means <- collected_true_pop_means_restricted %>%
  filter(sens == c(0.4) & N==50 &
           var %in% c(0.1, 1) &
           Parameter != "phi") %>% mutate(var=factor(var))%>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))

p1<- ggplot(subset(plot_samples),
            aes(x=value, color=Computation, linetype=var))+
  #geom_vline(data=subset(true_params), aes(xintercept=value))+
  geom_density(aes(group=interaction(Computation, Parameter, var)), linewidth=1)+
  geom_vline(data=subset(sub_pop_means, sens==0.4 & var%in% c(0.1, 1) & N==50), 
             aes(xintercept=value, linetype=var))+
  scale_color_manual(values=two_colors_transformations)+
  scale_linetype_manual(values=c("dashed", "solid"))+
  labs(linetype="Variability", y="Posterior density", x="Value")+
  facet_wrap(.~Parameter, scales = "free", labeller = label_parsed, drop=TRUE)+
  custom_theme+
  theme(legend.direction = "vertical", legend.box = "vertical",
        legend.position = "right")
p1
ggsave("figures/Recovery_restricted_distributions.eps",
       width = 17.62, height=8, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/Recovery_restricted_distributions.png",
       width = 17.62, height=8, units="cm",dpi=900)



differences_df <- collected_samples_restricted %>% 
  mutate(iter=row_number(),.by=c(N, var, sens, parname, Parameter, Computation)) %>% 
  select(-parname) %>% 
  pivot_wider(names_from="Computation", values_from = value) %>%
  mutate(difference = Incorrect-Correct) %>%
  group_by(N, sens, var, Parameter) %>%
  reframe(Med=median(difference), 
          lower = quantile(difference, probs = 0.025),
          upper = quantile(difference, probs = 0.975)) 

plot_differences_df <- differences_df %>% 
  filter(sens==0.4) %>%
#  mutate(Parameter = ifelse(Parameter=="sens", "phi", Parameter))
  filter(Parameter!="sens")%>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))

ggplot(plot_differences_df, aes(x=as.factor(var), y=Med))+
  geom_point()+geom_line(aes(group=1))+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  facet_nested(Parameter~"Sample~Size"+N, scales = "free_y", labeller=label_parsed)+ 
  labs(x="Variability between individuals", 
       y="Induced bias in group-level means (Incorrect-Correct)")+
  custom_theme+
  theme(panel.spacing= unit(0.1, "cm"))
ggsave("figures/Recovery_restricted_trafodifferences_SUPPLEMENT.eps",
       width = 12, height=12, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/Recovery_restricted_trafodifferences_SUPPLEMENT.png",
       width = 12, height=12, units="cm",dpi=900)




#___________________________________________________________________----
# D  Re-doing age difference analysis in Pachur et al (2017)  ______----
#___________________________________________________________________----
## 1. Fit the hierarchical CPT-model                                ----

# Define initial values for parameters 
inits = function() {
  list(mu.phi.alpha = 0, sigma.phi.alpha = 1, 
       mu.phi.gamma = 0, sigma.phi.gamma = 1, 
       mu.phi.delta_p = 0, sigma.phi.delta_p = 1,
       mu.phi.delta_m = 0, sigma.phi.delta_m = 1,
       mu.phi.lambda = 0, sigma.phi.lambda = 1, 
       mu.phi.sens = 0, sigma.phi.sens = 1) 
}


# Define the variables of interest. JAGS will return these to R when 
# the analysis is finished (and JAGS is closed).	
parameters = c("alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha", "mu.alpha_sebi",
               "gamma", "mu.phi.gamma", "mu.gamma", "sigma.phi.gamma", "mu.gamma_sebi",
               "delta_p", "mu.phi.delta_p", "mu.delta_p", "sigma.phi.delta_p", "mu.delta_p_sebi",
               "delta_m", "mu.phi.delta_m", "mu.delta_m", "sigma.phi.delta_m", "mu.delta_m_sebi",
               "lambda", "mu.phi.lambda", "mu.lambda", "sigma.phi.lambda", "mu.sens_sebi",
               "sens", "mu.phi.sens", "mu.sens", "sigma.phi.sens", "mu.lambda_sebi")

## To prevent re-fitting when save results are present
if (!file.exists("saved_details/Refitted_Age_Data.RData")) {
  ## Fit younger group
  age_choices = young_choices
  N_parts <- ncol(young_choices)
  res_younger =  jags.parallel(age_data,
                               parameters,  model.file = Pachur_age_model,
                               inits = inits,
                               n.chains = 6, n.iter = 40000, n.burnin = 3000, n.thin = 20,
                               n.cluster = 6, jags.seed = 771)
  res_younger <- list(samples=res_younger$BUGSoutput$sims.array,
                      summaries = res_younger$BUGSoutput$summary)
  
  ## Fit older group
  age_choices = older_choices
  N_parts <- ncol(older_choices)
  res_older =  jags.parallel(age_data,
                             parameters,  model.file = Pachur_age_model,
                             inits = inits,
                             n.chains = 6, n.iter = 40000, n.burnin = 3000, n.thin = 20,
                             n.cluster = 6, jags.seed = 188)
  res_older <- list(samples=res_older$BUGSoutput$sims.array,
                    summaries = res_older$BUGSoutput$summary)
  
  save(res_older, res_younger, file="saved_details/Refitted_Age_Data.RData")
}

load("saved_details/Refitted_Age_Data.RData")

#res_younger$summaries[order(-res_younger$summaries[,"Rhat"]),] # Rhat(sigma.phi.alpha) = 1.007
#res_older$summaries[order(-res_older$summaries[,"Rhat"]),]  # Rhat(sigma.phi.sens) = 1.007

## 2. Compare means between young and old                           ----
plot_parameters <- c("mu.alpha",   "mu.alpha_sebi",
                    "mu.gamma",   "mu.gamma_sebi",
                    "mu.delta_p", "mu.delta_p_sebi",
                    "mu.delta_m", "mu.delta_m_sebi",
                    "mu.lambda", "mu.sens_sebi",
                    "mu.sens",   "mu.lambda_sebi",
                    "sigma.phi.alpha", "sigma.phi.delta_m", 
                    "sigma.phi.delta_p", "sigma.phi.gamma", "sigma.phi.lambda", "sigma.phi.sens")
collected_age_summaries <- rbind(
  cbind(as_tibble(res_older$summaries[plot_parameters,c("mean", "50%", "2.5%", "97.5%")]), `Age group`="Older", parname = plot_parameters),
  cbind(as_tibble(res_younger$summaries[plot_parameters,c("mean", "50%", "2.5%", "97.5%")]),`Age group`="Younger", parname = plot_parameters)
  ) %>%
  mutate(Statistic = ifelse(grepl("sigma", parname), "Variability","Mean"),
         Computation = ifelse(Statistic=="Variability", "Variability", 
                                 ifelse(grepl("sebi", parname), "Correct Mean", "Incorrect Mean")), 
         Computation = factor(Computation, levels=c("Incorrect Mean", "Correct Mean", "Variability")),
         Parameter = sub("sigma.phi.", "", sub("_sebi", "", sub("mu.", "", parname))),
         Parameter = factor(Parameter, 
                            levels= c("alpha", "delta_m", "delta_p", "gamma", "lambda", "sens"),
                            labels= c("alpha", "delta^'-'", "delta^'+'", "gamma", "lambda", "phi"))) 

pd <- position_dodge(width=0.4)
collected_age_summaries %>% 
  ggplot(aes(x=`Age group`, color=Computation, shape=`Age group`))+
  scale_color_manual(name="",values=three_colors_trafovar)+
  geom_point(aes(y=`50%`), size=3, position=pd)+
  scale_x_discrete(labels = scales::parse_format())+
  ylab("Posterior median (95%CI)")+guides(shape="none")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
  facet_wrap(~Parameter, scales = "free_y",
             labeller = label_parsed, nrow=2)+
  custom_theme

ggsave("figures/Age_Comparison.eps",
       width = 17.62, height=9/0.6, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/Age_Comparison.png",
       width = 17.62, height=9/0.6, units="cm",dpi=900)


## Re-produce Table 5 in Pachur et al. (2017)
collected_age_samples <- rbind(
  cbind(as.data.frame(apply(res_older$samples[,,plot_parameters], 3, c)), group="Older"),
  cbind(as.data.frame(apply(res_younger$samples[,,plot_parameters], 3, c)), group="Younger"))

age_samples_long <- collected_age_samples %>% 
  pivot_longer(cols = -group, names_to = "parname", values_to = "samples") %>%
  mutate(Computation = ifelse(grepl("sebi", parname), "Correct", "Incorrect"), 
         Parameter = sub("_sebi", "", sub("mu.", "", parname)),
         Parameter = factor(Parameter, 
                            levels= c("alpha", "delta_m", "delta_p", "gamma", "lambda", "sens"),
                            labels= c("alpha", "delta^'-'", "delta^'+'", "gamma", "lambda", "phi"))) %>%
  filter(!grepl("sigma", parname))

summary_differences <- age_samples_long %>% 
  group_by(Parameter, Computation, group) %>% 
  mutate(N=1:n()) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=c(Parameter, Computation, N), 
              values_from = samples, names_from = group) %>%
  mutate(diff=Older-Younger) %>% 
  group_by(Parameter, Computation) %>% 
  reframe(Lower = quantile(diff, 0.025), 
          Upper = quantile(diff, 0.975),
          value=paste0(format(round(mean(diff), 2), nsmall=2), " [", 
                       format(round(Lower, 2), nsmall=2), ", ",
                       format(round(Upper, 2), nsmall=2), "]")) %>%
  ## Make credible differences get printed in bold in Latex:
  mutate(bold = Lower>0 | Upper < 0,
         value = ifelse(bold, paste0("\\textbf{", value,"}"), value)) %>% 
  select(-Lower, -Upper, -bold) %>% # remove unncesseray columns
  mutate(`Age group`= "Difference\n(older-younger)")

## Kept, in case we decide for a different format of the table
# 
# table_comparison <- collected_age_summaries %>%
#   mutate(value= paste0(format(round(mean, 2), nsmall=2), " [", 
#                        format(round(`2.5%`, 2), nsmall=2), ", ",
#                        format(round(`97.5%`, 2), nsmall=2), "]")) %>%
#   select( `Age group`, Parameter, Computation, value) %>%
#   rbind(summary_differences) %>%
#   mutate(Parameter = factor(Parameter, 
#                             labels=c("$\\alpha$", "$\\delta^-$", "$\\delta^+$", "$\\gamma$", "$\\lambda$", "$\\phi$")))%>%
#   pivot_wider(names_from = Parameter)
# post_table <- kable(table_comparison, format = "latex", escape = FALSE,
#       caption="Posterior means (and 95\\% CIs) for the fitted parameters in younger and older individuals and their difference.")
# writeLines(post_table, 'figures/TableAgeComparison.tex')




table_comparison2 <- collected_age_summaries %>% filter(!grepl("sigma", parname)) %>%
  mutate(Computation=sub(" Mean", "", as.character(Computation))) %>%
  mutate(value= paste0(format(round(mean, 2), nsmall=2), " [", 
                       format(round(`2.5%`, 2), nsmall=2), ", ",
                       format(round(`97.5%`, 2), nsmall=2), "]")) %>%
  select( `Age group`, Parameter, Computation, value) %>%
  rbind(summary_differences) %>%
  mutate(Parameter = factor(Parameter, 
                            labels=c("$\\alpha$", "$\\delta^-$", "$\\delta^+$", "$\\gamma$", "$\\lambda$", "$\\phi$")) ,
         `Age group`= ifelse(grepl("Diff", `Age group`), "(Older-Younger)",# "\\baselineskip=15pt Difference\\newline (Older-Younger)",
                             `Age group`)
         ) %>% 
  select(Parameter, Computation, `Age group`, value) %>%
  pivot_wider(names_from = c(`Age group`)) %>% 
  arrange(Parameter, Computation) %>%
  group_by(Parameter) %>% 
  mutate(Parameter=c(paste0("\\multirow{ 2}{*}{", Parameter[1],"}"), "")) %>%
  ungroup() %>% 
  mutate(Computation = as.character(Computation)) %>% 
  rename(Computation=Computation)
#table_comparison2 <- rbind(names(table_comparison2), table_comparison2)
#names(table_comparison2) <- c("Parameter", "Computation", "\\multicolumn{2}{c}{Age Group}", "Difference")
table_comparison2 <- xtable(table_comparison2, align = c("l", "l", "l","|", "c", "c","|", "c"),
                            caption="\\raggedright Mean and 95\\% CI for posterior distribution of parameters.")
addtorow <- list()
addtorow$pos <- list(c(-1),c(2, 4, 6, 8, 10, 12)) 
addtorow$command <- c('&&\\multicolumn{2}{|c|}{Age Group}&Difference \\\\', '\\midrule')
#addtorow$pos <- list(c(2, 3, 5, 7, 9, 11, 13)) 
#addtorow$command <- c('\\midrule')

print(table_comparison2, type="latex",sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(table_comparison2)), booktabs = TRUE,
      caption.placement="top")
dir.create("figures", showWarnings = FALSE)
print(table_comparison2, type="latex",
      file="figures/TableAgeComparison2.tex",
      sanitize.text.function=function(x){x},
      include.rownames=FALSE,
      add.to.row=addtorow,
      hline.after = c(0, nrow(table_comparison2)), booktabs = TRUE,
      caption.placement="top",label = "tab:age",
      table.placement="hp")
#print.xtable()

# ### Overloaded plot
# 
# library(ggpattern)
# is_in_range <- function(x, range) return(x > min(range) & x < max(range))
# 
# plot_parameters = c("mu.phi.alpha", "mu.alpha", "sigma.phi.alpha", "mu.alpha_sebi",
#                     "mu.phi.gamma", "mu.gamma", "sigma.phi.gamma", "mu.gamma_sebi",
#                     "mu.phi.delta_p", "mu.delta_p", "sigma.phi.delta_p", "mu.delta_p_sebi",
#                     "mu.phi.delta_m", "mu.delta_m", "sigma.phi.delta_m", "mu.delta_m_sebi",
#                     "mu.phi.lambda", "mu.lambda", "sigma.phi.lambda", "mu.sens_sebi",
#                     "mu.phi.sens", "mu.sens", "sigma.phi.sens", "mu.lambda_sebi")
# 
# collected_age_samples <- rbind(
#   cbind(as.data.frame(apply(res_older$samples[,,plot_parameters], 3, c)), group="older"),
#   cbind(as.data.frame(apply(res_younger$samples[,,plot_parameters], 3, c)), group="younger"))
# age_samples_long <- collected_age_samples %>%
#   pivot_longer(cols = -group, names_to = "parameter", values_to = "samples")
# 
# quantiles_age_samples <- age_samples_long %>%
#   group_by(parameter, group) %>% 
#   reframe(quantiles=quantile(samples, probs=c(0.025, 0.975)))
# densities_age_samples <- age_samples_long %>%
#   group_by(parameter, group) %>%
#   reframe(densx = density(samples)$x,
#           densy = density(samples)$y)
# 
# ## Clean and Format Parameter Labels
# densities_age_samples <- densities_age_samples %>% 
#   #filter(!grepl("phi", parname) & !grepl("lmu", parname)) %>%
#   mutate(Statistic = ifelse(grepl("mu", parameter), "Mean", "SD"),
#          Scale = ifelse(grepl("phi", parameter), "Real", "Parameter"),
#          Computation = ifelse(grepl("sebi", parameter), "Correct", "Incorrect"), 
#          Parameter = sub("_sebi", "", sub("mu.", "", parameter)))
# 
# densities_age_samples_HDI <- densities_age_samples %>%
#   group_by(parameter, group) %>% 
#   filter(is_in_range(densx, subset(quantiles_age_samples, parameter==cur_group()$parameter & group==cur_group()$group)$quantiles)) %>%
#   ungroup()
# 
# 
# 
# 
# p_group_comparison <-ggplot(subset(densities_age_samples, Scale=="Parameter" & Statistic=="Mean"), 
#                             aes(x=densx, y=densy))+
#   geom_line(aes(color=Computation, linetype=group))+
#   geom_area_pattern(data =subset(densities_age_samples_HDI, Scale=="Parameter" & Statistic=="Mean"),
#                     mapping=aes(pattern_density=group, pattern_spacing=group,
#                                 color=Computation,fill= Computation,
#                                 group=interaction(Computation, group)),
#                     alpha=0.5, position="identity",
#                     pattern_fill="gray20", pattern_spacing=0.06,
#                     show.legend=c(pattern_density=TRUE, color=FALSE, fill=TRUE))+
#   scale_pattern_density_manual(name="",values = c(`older` = 0, `younger`=0.004))+
#   scale_discrete_manual(aesthetics = c("color", "fill"), name="", values = two_colors_transformations)+
#   scale_y_continuous(name="Posterior density",
#                      expand = expansion(mult=c(0.01, 0.05)))+# c(0.01))+
#   #    expand_limits(y=c(0.01, 23))+
#   facet_nested(Parameter~., scales = "free", independent = "x")+
#   xlab("Parameter value") +
#   custom_theme#+ #ylab("Posterior density")+
# #theme_bw()+theme(legend.position = "bottom")+
# # ggtitle(paste0("Posterior distributions of mean coefficients (shaded area represents 95%-HDI)",
# #                "\nStudy", study, "; Model: ", model))
# p_group_comparison

#___________________________________________________________________----
#_______                 For Supplement                     ________----
#___________________________________________________________________----
# E  Refitting Rieskamp-data with original model                    ----
## 1. Fit the hierarchical CPT-model                                ----

# Define initial values for parameters 
inits = function() {
  list(mu.phi.alpha = 0.7, sigma.phi.alpha = 1, 
       mu.phi.beta = 0.7, sigma.phi.beta = 1,
       mu.phi.gamma.gain = 0.7, sigma.phi.gamma.gain = 1, 
       mu.phi.gamma.loss = 0.7, sigma.phi.gamma.loss = 1,
       lmu.lambda = 0, lsigma.lambda = 0.5, 
       lmu.sens = 0, sigma.phi.sens = 0.5) 
}


# Define the variables of interest. JAGS will return these to R when 
# the analysis is finished (and JAGS is closed).	
parameters = c("alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha", "mu.alpha_sebi",
               "beta", "mu.phi.beta", "mu.beta", "sigma.phi.beta", "mu.beta_sebi",
               "gamma.gain", "mu.phi.gamma.gain", "mu.gamma.gain", "sigma.phi.gamma.gain", "mu.gamma.gain_sebi",
               "gamma.loss", "mu.phi.gamma.loss", "mu.gamma.loss", "sigma.phi.gamma.loss", "mu.gamma.loss_sebi",
               "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", "mu.sens_sebi",
               "sens", "lmu.sens", "mu.sens", "lsigma.sens", "mu.lambda_sebi")

## To prevent re-fitting when save results are present
if (!file.exists("saved_details/Refitted_Data.RData")) {
  res_rieskamp_1 =  jags.parallel(data,
                                  parameters,  model.file = original_full_model,
                                  inits = inits,
                                  n.chains = 4, n.iter = 20000, n.burnin = 1000, n.thin = 10,
                                  n.cluster = 4, jags.seed = 531)
  res_rieskamp_1 <- list(samples=res_rieskamp_1$BUGSoutput$sims.array,
                         summaries = res_rieskamp_1$BUGSoutput$summary)
  save(res_rieskamp_1, file="saved_details/Refitted_Data.RData")
}

load("saved_details/Refitted_Data.RData")


## 2. Compare population means between transformations              ----
temp_summary <- res_rieskamp_1$summaries
#max(res_rieskamp_1$BUGSoutput$summary[,"Rhat"])
parname <- rownames(temp_summary)
temp_summary <- as_tibble(temp_summary) %>% mutate(parname = parname)
group_pars_summary <- temp_summary %>% 
  filter(grepl(parname, pattern = "mu"))

pd <- position_dodge(width=0.2)
plt_group_pars_summary <- group_pars_summary %>% 
  filter(!grepl("phi", parname) & !grepl("lmu", parname)) %>%
  mutate(Computation = ifelse(grepl("sebi", parname), "Correct", "Incorrect"), 
         Computation = factor(Computation, levels=c("Incorrect", "Correct")),
         Parameter = sub("_sebi", "", sub("mu.", "", parname)))%>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))
ggplot(plt_group_pars_summary , aes(x=Parameter, color=Computation))+
  scale_color_manual(values=two_colors_transformations)+
  geom_point(aes(y=`50%`), size=3, position=pd)+
  scale_x_discrete(labels = scales::parse_format())+
  ylab("Posterior median (95%CI)")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
  custom_theme
ggsave("figures/Rieskamp_Original.eps",
       width = 17.62, height=9/0.7, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/Rieskamp_Original.png",
       width = 17.62, height=9/0.7, units="cm",dpi=900)

#___________________________________________________________________----
# F  Re-do (Extended) original simulation study (unconstrained)     ----

## 1. Actual parameter recovery analysis                            ----

# Define initial values for parameter
inits = function() {
  list(mu.phi.alpha = 0.7, sigma.phi.alpha = 1,
       mu.phi.beta = 0.7, sigma.phi.beta = 1,
       mu.phi.gamma.gain = 0.7, sigma.phi.gamma.gain = 1, 
       mu.phi.gamma.loss = 0.7, sigma.phi.gamma.loss = 1,
       lmu.lambda = 0, lsigma.lambda = 0.5, 
       lmu.sens = 0, sigma.phi.sens = 0.5)
}


# Define the variables of interest. JAGS will return these to R 
# when the analysis is finished (and JAGS is closed).	
parameters = c("alpha", "mu.phi.alpha", "mu.alpha", "sigma.phi.alpha", "mu.alpha_sebi",
               "beta", "mu.phi.beta", "mu.beta", "sigma.phi.beta", "mu.beta_sebi",
               "gamma.gain", "mu.phi.gamma.gain", "mu.gamma.gain", "sigma.phi.gamma.gain", "mu.gamma.gain_sebi",
               "gamma.loss", "mu.phi.gamma.loss", "mu.gamma.loss", "sigma.phi.gamma.loss", "mu.gamma.loss_sebi",
               "lambda", "lmu.lambda", "mu.lambda", "lsigma.lambda", "mu.sens_sebi",
               "sens", "lmu.sens", "mu.sens", "lsigma.sens", "mu.lambda_sebi"
)


## Set mean parameters for simulation
alpha <- beta <- .88
gamma.gain <- .61 
gamma.loss <- .69
lambda <- 2.25

## Define the different settings that should be compared
phis <- c(.04, .14, .40) # choice sensitivity
Nsbjs <- c(20, 50, 90) # number of subjects
variabilities <- c(0.1, 0.5, 1) # btw-sbj variability in parameters

## Actually do the simulation, save simulation, and model fitting
## Only do this, when all analysis should be done again (takes long!)
if (REDOALLANALYSIS) {
  collected_samples <- data.frame()
  collected_summaries <- data.frame()
  collected_true_pop_means <- data.frame()
  getpars <-c("alpha", "beta", "gamma.loss", "gamma.gain", "lambda", "sens")
  getpars <- paste0("mu.", getpars, rep(c("", "_sebi"), each=length(getpars)))
  
  dir.create("saved_details/Recovery_full", showWarnings = FALSE)
  N <- VAR <- PHI <- 1
  for (N  in 1:3) {
    cur_n <- Nsbjs[N]
    Data <- matrix(NA, nrow=60, ncol=cur_n)
    for (VAR in 1:3) {
      cur_var <- variabilities[VAR]
      for ( PHI in 1:3) {
        cur_sens <- phis[PHI]
        ## Only if the saved samples do not already exists
        if (!file.exists(paste0("saved_details/Recovery_full/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))) {
          ## Make it reproducible
          seeeed <- 2201 + 100*N + 10*VAR + PHI 
          set.seed(seeeed)
          
          ## Sample from Beta-distribution with mean alpha and scaled variance cur_var (not exactly the variance!)
          Alphas <- rbeta(cur_n, alpha*((alpha*(1-alpha))/cur_var *20 -1), (1-alpha)*((alpha*(1-alpha))/cur_var *20 -1) )
          Betas  <- rbeta(cur_n,  alpha*((alpha*(1-alpha))/cur_var *20 -1), (1-alpha)*((alpha*(1-alpha))/cur_var *20 -1) )
          Gammas.gain <- rbeta(cur_n, gamma.gain*((gamma.gain*(1-gamma.gain))/cur_var *10 -1), 
                               (1-gamma.gain)*((gamma.gain*(1-gamma.gain))/cur_var *10 -1) )
          Gammas.loss <- rbeta(cur_n, gamma.loss*((gamma.loss*(1-gamma.loss))/cur_var *10 -1), 
                               (1-gamma.loss)*((gamma.loss*(1-gamma.loss))/cur_var *10 -1) )
          # Draw from Gamma distribution with mean lambda and variance cur_var
          Lambdas <- rgamma(cur_n, shape= lambda^2/cur_var , scale=cur_var/lambda)
          
          for (k in 1:cur_n) {
            Data[,k] <- simulate_CPT_individ(Alphas[k], Betas[k], Gammas.gain[k], Gammas.loss[k], Lambdas[k], cur_sens)
          }
          params <- data.frame(alpha=Alphas, beta=Betas, gamma.gain=Gammas.gain, gamma.loss=Gammas.loss, lambda=Lambdas, phi=cur_sens)
          simulation_pars <- list(N = cur_n, var=cur_var, sens=cur_sens)
          save(Data, params, simulation_pars,
               file=paste0("saved_details/Recovery_full/SampledData_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
          
          rec_samples =  jags.parallel(simu_data, parameters,
                                       model.file = original_full_model_recovery,
                                       inits = inits,  n.chains = 4,
                                       n.iter = 50000, n.burnin = 1000,
                                       n.thin = 5,  n.cluster = 4, jags.seed = seeeed)
          rec_summary <- rec_samples$BUGSoutput$summary
          rec_samples <- rec_samples$BUGSoutput$sims.array
          save(Data, params, simulation_pars, rec_summary, rec_samples, 
               file=paste0("saved_details/Recovery_full/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
        } else {
          load(paste0("saved_details/Recovery_full/RecoveryResult_N_", cur_n,"_var_", cur_var, "_phi_", cur_sens,".RData"))
        }
        
        ## Combine the whole posterior samples of population parameters
        temp <- rec_samples[,, getpars]    
        dim(temp) <- c(dim(temp)[1]*dim(temp)[2], dim(temp)[3])
        colnames(temp) <- getpars     
        temp <- as.data.frame(temp) 
        #head(temp)
        temp <- cbind(temp, as.data.frame(simulation_pars))
        collected_samples <- rbind(collected_samples, temp)
        
        ## Combine the posterior summaries of population parameters
        temp <- rec_summary[getpars,] 
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
  ## Clean and Format Parameter Labels
  collected_samples <- collected_samples %>% 
    #filter(!grepl("phi", parname) & !grepl("lmu", parname)) %>%
    pivot_longer(1:12, names_to="parname") %>%
    mutate(Computation = ifelse(grepl("sebi", parname), "Correct", "Incorrect"), 
           Computation = factor(Computation, levels=c("Incorrect", "Correct")),
           Parameter = sub("_sebi", "", sub("mu.", "", parname)))
  collected_summaries <- collected_summaries %>% 
    mutate(Computation = ifelse(grepl("sebi", parname), "Correct", "Incorrect"), 
           Computation = factor(Computation, levels=c("Incorrect", "Correct")),
           Parameter = sub("_sebi", "", sub("mu.", "", parname)))
  
  save(collected_samples,collected_summaries, collected_true_pop_means, 
       file="saved_details/Collected_recovery_results.RData")
}


load("saved_details/Collected_recovery_results.RData")

## 2. Visualize original full parameter recovery analysis           ----
## Reproduce Nilsson et al. (2011), Figure 2:
# Note: variability in Nilsson et al is 0; and N = 30; but the following are
# the values most close to those in Nilsson's paper:
plot_samples <- filter(collected_samples, sens==0.4 & var%in%c(0.1, 1) & N == 20) %>%
  filter(Parameter != "sens")%>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))
true_params <- data.frame(Parameter= c("alpha", "beta", "gamma.gain","gamma.loss","lambda"), 
                          value    = c(   .88,    .88,       .61,    .69,  2.25))%>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))

ggplot(plot_samples, aes(x=value, linetype=as.factor(var), color=Computation))+
  geom_vline(data=true_params, aes(xintercept=value))+
  geom_density(aes(group=interaction(Computation, Parameter, var)), linewidth=1)+
  scale_color_manual(values=two_colors_transformations)+
  facet_wrap(.~Parameter, scales = "free", labeller=label_parsed)+
  labs(y="Posterior density", x="Parameter value", linetype="Variability")+
  custom_theme+
  theme(plot.margin = margin(0, 0.3, 0, 0, "cm"))
ggsave("figures/Recovery_full_posteriordists.eps",
       width = 23, height=9/0.6, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/Recovery_full_posteriordists.png",
       width = 23, height=9/0.6, units="cm",dpi=900)

# true_params <- data.frame(Parameter= c("alpha", "beta", "gamma.gain","gamma.loss","lambda"), 
#                           value    = c(   .88,    .88,       .61,    .69,  2.25))
# pd <- position_dodge(width=0.4)
# ggplot(filter(collected_summaries,Parameter!="sens"), aes(y=mean, x=interaction(N,var), color=Computation))+
#   geom_hline(data=true_params, aes(yintercept=value))+
#   geom_point(position=pd)+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position=pd, width=0.2)+
#   facet_grid(Parameter~sens, scales = "free")



# Only take the extreme sampling options for each factor
sub_results <- subset(collected_summaries,Parameter!="sens") %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1)) %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens))%>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))
sub_pop_means <- collected_true_pop_means %>%
  filter(sens %in% c(0.04, 0.4) &
           var %in% c(0.1, 1) &
           Parameter != "phi") %>%
  mutate(var=paste0("Variability: ", var),
         sens=paste0("Sensitivity: ", sens)) %>%
  merge(data.frame(Computation=c("Incorrect", "Correct")))%>%
  mutate(Parameter = factor(Parameter, levels=par_names, labels=par_labels))
pd <- position_dodge(width=0.2)
ggplot(sub_results,
       aes(y=`50%`, x=as.factor(N), color=Computation))+
  geom_hline(data=true_params, aes(yintercept=value))+
  geom_errorbar(data=sub_pop_means , aes(ymin=value, y=value,ymax=value), linetype="dashed", color="gray20")+
  geom_line(aes(group=Computation),position=pd)+
  geom_point(position=pd)+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), position=pd, width=0.2)+
  scale_color_manual(values=two_colors_transformations)+
  facet_nested(Parameter~var+sens, scales = "free", labeller = label_parsed )+
  labs(y="Parameter values", x="Simulated sample size")+
  custom_theme+
  theme(panel.spacing= unit(0.1, "cm"))
ggsave("figures/Recovery_full_posteriorCIs_SUPPLEMENT.eps",
       width = 17.62, height=17.62, units="cm",dpi=600, device = cairo_ps)

ggsave("figures/Recovery_full_posteriorCIs_SUPPLEMENT.png",
       width = 17.62, height=17.62, units="cm",dpi=900)


# 
# sub_results <- subset(collected_summaries,Parameter!="sens") %>%
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
#   merge(data.frame(Computation=c("Incorrect", "Correct")))
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