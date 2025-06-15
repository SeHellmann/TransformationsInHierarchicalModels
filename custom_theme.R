if (grepl("Windows", osVersion)) {
  # To use Times as font on Windows (otherwise ignore the ggplot warnings)
  windowsFonts(Times=windowsFont("Times New Roman"))
}
two_colors_transformations <- c("#1b9e77", "#fc8d62")
# color_values <- palette.colors(n = 3, palette = "R4")[-1]
# color_values <- palette.colors(n = 3, palette = "Okab")[-1] %>% as.character()
# color_values <- palette.colors(n = 2, palette = "Okab") %>% as.character()
custom_theme <- theme_bw()+
  theme(text = element_text(family = "Times", size=12),
        strip.text = element_text(family="Times", size=12),
        legend.position = "bottom",
        plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text = element_text(size=12, family="Times", color="black"),
        legend.spacing = unit(0,"line"),
        legend.key.width=unit(1,"line"),
        panel.spacing=unit(1, "lines"))
