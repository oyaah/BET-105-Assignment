# alternative R plotting script (need ggplot2, dplyr, readr)
library(ggplot2)
library(dplyr)
library(readr)

df <- read_tsv("final/angles.tsv", show_col_types = FALSE)

df$size_class <- factor(df$size_class,
                        levels = c("Tiny", "Small", "Intermediate", "Large", "Bulky"))

df <- df %>% mutate(angle = ((angle + 180) %% 360) - 180)

pal <- c(
  "Tiny"         = "#c8c4b5",
  "Small"        = "#e8b44c",
  "Intermediate" = "#e07840",
  "Large"        = "#c93d1d",
  "Bulky"        = "#8b0a1a"
)

n <- nrow(df)

p <- ggplot(df, aes(x = angle, color = size_class)) +
  geom_density(linewidth = 1.2, adjust = 1.2) +
  scale_color_manual(values = pal) +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-150, 150, 50)) +
  labs(
    title = paste0("Tripeptide (XRX) in Helix (n = ", n, ")"),
    x = "Angle between adjacent CA-Centroid vectors (degrees)",
    y = "Norm. Freq. [A.U.]",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background  = element_rect(fill = "#bdbdbd", color = NA),
    plot.background   = element_rect(fill = "#bdbdbd", color = NA),
    panel.grid.major.x = element_line(color = "blue", linetype = "dotted"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = c(0.15, 0.85),
    legend.background  = element_rect(fill = "white", color = "black"),
    plot.title         = element_text(hjust = 0.5)
  )

ggsave("final/angle_plot.png", plot = p, width = 10, height = 6, dpi = 300)
