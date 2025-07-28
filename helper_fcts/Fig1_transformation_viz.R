library(tidyverse)
library(cowplot)

# data -------------------------------------------------------------
N <- 1e5
conditions <- expand_grid(tibble(sign = c("positive", "negative") , mu = c(1, -1)) , 
                          tibble(vars = c("small", "large") , sigma = c(0.1, 1)) 
                          )
# simulate data
dat <- pmap(conditions, function(sign, mu, vars, sigma){
  tibble(sign, vars , real = rnorm(N,mu,sigma) , parameter = pnorm(real))
  }) |>
  list_rbind()

# compute means ...
means <- dat |> 
  group_by(sign, vars) |> 
  summarise(Correct = mean(parameter) , 
            Original = pnorm(mean(real))) |> 
  ungroup() |> 
  pivot_longer(cols=Correct:Original, values_to = "mean", names_to = "Computation") |> 
  mutate(dens = NA)

# ... and their densities
for(i in 1:nrow(means)) { 
  dat_temp <- dat |> filter(vars==means[[i,'vars']] & sign==means[[i,'sign']])
  den <- density(dat_temp$parameter)
  den <- data.frame(x = den$x, y = den$y)
  means[[i, 'dens']] <- approx(den$x, den$y, xout = means[[i, 'mean']])$y
  }

# plots --------------------------------------------------------------------
## theme ----------------------------------------------------------
if (!exists("custom_theme")) {
  source('helper_fcts/custom_theme.R') # get Sebi's custom theme
}
expand <- c(0.01,0.01) # remove space around axis limits
two_lines_variance <- c('dotdash', 'solid')
linewidths <- c(1,1.5)

# to grob legend: https://github.com/wilkelab/cowplot/issues/202#issuecomment-1981765769
get_legend2 <- function(plot, legend_number = 1) {
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  idx <- which(vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE))
  if (length(idx) >= legend_number) {
    return(legends[[idx[legend_number]]])
  } else if (length(idx) >= 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}

## base --------------------------------------------------------------------

og_hline_p <- data.frame(xbegin=1, xend=5, ybegin=pnorm(1), yend=pnorm(1))
og_vline_p <- data.frame(xbegin=1, xend=1, ybegin=pnorm(1), yend=1)
og_hline_n <- data.frame(xbegin=-1, xend=5, ybegin=pnorm(-1), yend = pnorm(-1))
og_vline_n <- data.frame(xbegin=-1, xend=-1, ybegin=pnorm(-1), yend=1)

baseplot <- 
  tibble(real=seq(-5,5,.01) , parameter=pnorm(real)) |> 
  ggplot(aes(real, parameter)) + 
  geom_line(linewidth=linewidths[1]) +
  geom_segment(data=og_vline_p, aes(x=xbegin, xend=xend, y=ybegin, yend=yend), linetype=two_lines_variance[2], linewidth=linewidths[1], color=two_colors_transformations[2]) +
  geom_segment(data=og_hline_p, aes(x=xbegin, xend=xend, y=ybegin, yend=yend), linetype=two_lines_variance[2], linewidth=linewidths[1], color=two_colors_transformations[2], arrow = arrow(length = unit(0.06, "npc"), angle=20, type="open")) +
  geom_segment(data=og_vline_n, aes(x=xbegin, xend=xend, y=ybegin, yend=yend), linetype=two_lines_variance[1], linewidth=linewidths[1], color=two_colors_transformations[2]) +
  geom_segment(data=og_hline_n, aes(x=xbegin, xend=xend, y=ybegin, yend=yend), linetype=two_lines_variance[1], linewidth=linewidths[1], color=two_colors_transformations[2], arrow = arrow(length = unit(0.06, "npc"), type="open")) +
  scale_x_continuous(breaks = seq(-5,5,1),  expand=expand) +
  scale_y_continuous(breaks = seq(0,1,.25),  expand=expand) +
  labs(y = "Parameter Scale" , 
       x = "Real Scale") +
  custom_theme

baseplot_l <- 
  baseplot + 
  annotate("label", label = expression(Phi[mu * "," *sigma^2](X)) , 
           x = 2, y = .5, family="Times", label.padding = unit(0.4, "lines"), size = 5 ) 

## left (small variance) --------------------------------------------------------------------
### x marginal --------------------------------------------------------
xma_l <- 
  dat |> 
  filter(vars=='small') |> 
  ggplot(aes(x=real, linetype=sign)) +
  geom_density(alpha=.5, linewidth=linewidths[2]) +
  #geom_segment(dat=means |> filter(vars=='small' & Computation=="Original"), aes(x = mean, xend = mean, y = 0, yend = dens, color=Computation), linewidth=linewidths[2]) +
  scale_x_continuous(limits=c(-5,5),  expand=expand) +
  scale_y_continuous(expand=expand) + 
  scale_linetype_manual(values=two_lines_variance) +
  theme_void() +
  theme(legend.position = "none") +
  annotate("text", x = -3, y = 2, label = expression(italic(N)(-1, .01)),  family="Times") + 
  annotate("text", x = 3, y = 2, label = expression(italic(N)(1, .01)),  family="Times") 

### y marginal ---------------------------------------------------------
yma_l <- 
  dat |> 
  filter(vars=='small') |> 
  ggplot(aes(x=parameter, linetype=sign)) +
  geom_density(alpha = .5, linewidth=linewidths[2]) +
  geom_segment(dat=means |> filter(vars=='small' & Computation=="Correct"), aes(x = mean, xend = mean, y = 0, yend = dens, color=Computation), linewidth=linewidths[2]) +
  scale_color_manual(values=two_colors_transformations)+
  scale_x_continuous(limits = c(0,1), expand=expand) +
  scale_y_continuous(expand=expand) +
  scale_linetype_manual(values=two_lines_variance) +
  theme_void() +
  theme(legend.position = "none", 
        axis.title.x = element_blank())
yma_l90 <- yma_l + coord_flip()

## right (large variance) ----------------------------------------------------------
### x marginal ----------------------------------------------------------
xma_r <- 
  dat |> 
  filter(vars=='large') |> 
  ggplot(aes(x=real, linetype=sign)) +
  geom_density(aes(fill=sign), alpha=.5, linewidth=linewidths[2]) +
  scale_x_continuous(limits = c(-5,5), expand=expand) +
  scale_y_continuous(expand=expand) + 
  scale_linetype_manual(values=two_lines_variance) +
  theme_void() +
  theme(legend.position = "none") +
  annotate("text", x = -3.8, y = .2, label = expression(italic(N)(-1, 1)), family="Times") + 
  annotate("text", x = 3.8, y = .2, label = expression(italic(N)(1, 1)),  family="Times") +
  scale_fill_manual(values = c("#FFFFFF00", "#fc8d62")) 

### y marginal -------------------------------------------------------

yma_r <-  
  dat |> 
  filter(vars=='large') |> 
  ggplot(aes(x=parameter, linetype=sign)) +
  geom_density(aes(fill=sign), alpha = .3, linewidth=linewidths[2]) +
  geom_segment(dat=means |> filter(vars=='large' & Computation=="Correct"), aes(x = mean, xend = mean, y = 0, yend = dens, color=Computation), linewidth=linewidths[2]) +
  scale_color_manual(values=two_colors_transformations)+
  scale_fill_manual(values = c("#FFFFFF00", "#1b9e77")) +
  scale_linetype_manual(name="Mean Original Scale", values=two_lines_variance, 
                        labels = c("p" = expression(mu = 1), "n" = expression(mu = -1))) +
  scale_x_continuous(limits = c(0,1), expand=expand) +
  scale_y_continuous(expand=expand) +
  theme_void() + 
  guides(color = guide_legend(nrow = 1), linetype = "none", fill='none') +
  theme(legend.position = 'none')
yma_r90 <- yma_r + coord_flip()

# collect plots -------------------------------------------------------------------
# small variance
g_xma_l <- ggplotGrob(xma_l)
g_yma_l <- ggplotGrob(yma_l90)
base_xma_l <- insert_xaxis_grob(baseplot_l, g_xma_l, position = "top", grid::unit(0.2, "null"))
plot_l <- insert_yaxis_grob(base_xma_l, g_yma_l, position = "right", grid::unit(0.2, "null"))

# large variance
g_xma_r <- ggplotGrob(xma_r)
g_yma_r <- ggplotGrob(yma_r90)
base_xma_r <- insert_xaxis_grob(baseplot, g_xma_r, position = "top", grid::unit(0.2, "null"))
plot_r <- insert_yaxis_grob(base_xma_r, g_yma_r, position = "right", grid::unit(0.2, "null"))

legend <- 
  ggplot(dat=means) +
  geom_segment(aes(x = mean, xend = mean, y = 0, yend = dens, color=Computation), linewidth=linewidths[2]) +
  scale_color_manual(name="Computation", values=two_colors_transformations) +
  theme_void()  
legend <- get_legend2(legend + theme(legend.position = "right", legend.spacing.x = unit(2,'cm')))


p_combined <- plot_grid(plot_l, NULL, plot_r, NULL, legend, 
                        nrow=1, ncol=5, rel_widths = c(1,.1,1,.01,.3) , 
                        labels = c("A","","B","","")) 

# safe
ggsave("figures/transformations.eps", width=10, height=3.9, unit='in', dpi=900, device = cairo_ps)
ggsave("figures/transformations.jpg", width=10, height=3.9, unit='in', dpi=900)
