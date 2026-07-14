library(ggh4x)
library(tidyverse)
par1 <- expand.grid(mu = seq(-5, 5, by = 0.1), 
                   sigma = c(0.1, 1, 2)) %>% 
  mutate(alpha_Incorrect= pnorm(mu), 
         alpha_Correct = pnorm(mu/sqrt(1+sigma^2)),
         alpha_Bias = alpha_Incorrect-alpha_Correct ) %>%
  pivot_longer(cols = 3:5, names_sep = "_", names_to = c("Parameter", "Computation"))  %>%
  mutate(Transformation = ifelse(Parameter=="alpha", "Phi", "exp"),
         Computation = factor(Computation, levels=c("Incorrect", "Correct", "Bias")))
par2 <- expand.grid(mu = seq(-2, 1, by = 0.1), 
                   sigma = c(0.1, 1, 2)) %>% 
  mutate(lambda_Incorrect = exp(mu),
         lambda_Correct= exp(mu + sigma^2/2),
         lambda_Bias = lambda_Incorrect-lambda_Correct) %>%
  pivot_longer(cols = 3:5, names_sep = "_", names_to = c("Parameter", "Computation"))  %>%
  mutate(Transformation = ifelse(Parameter=="alpha", "Phi", "exp"),
         Computation = factor(Computation, levels=c("Incorrect", "Correct", "Bias")))
par <- rbind(par1, par2)

ggplot(subset(par,Computation!="Bias"), aes(x=mu, y=value))+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  #geom_errorbar(data=sub_pop_means , aes(ymin=value, y=value,ymax=value), linetype="dashed", color="gray20")+
  #geom_point(position=pd)+
  geom_line(aes( color=Computation, group=Computation))+
  geom_line(data=subset(par, Computation=="Bias"),aes(group=Computation, linetype="Bias"))+
  scale_color_manual(values=c(two_colors_transformations))+
  scale_linetype_discrete(name="")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.05)))+
  facet_nested("Standard~deviation"+sigma~"Transformation"+Transformation, scales = "free", independent="y", labeller = label_parsed , drop = TRUE)+
  labs(y="Parameter scale mean", x="Real scale mean")+
  custom_theme+
  theme(panel.spacing.y= unit(0.1, "cm"), legend.box = "vertical", legend.direction = "horizontal",
        legend.spacing.y = unit(0, "cm"), legend.box.margin = margin(0, 0, 0, 0, "pt"),
        legend.key.spacing.y = unit(0.1, "cm"),legend.key.height = unit(0.1, "cm"), legend.box.spacing = unit(0.1, "cm"))

ggsave("figures/biasillustration.eps", 
       width = 14, height=13, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/biasillustration.png", 
       width = 14, height=13, units="cm",dpi=600, device=png)
