# Plot restitution portrait of simulation

library(ggplot2)
library(ggthemes)
library(plyr)
library(extrafont)

data.input <-
  read.table('../APD_BCL.dat', header = TRUE, sep = ",")

data.plot <- data.frame(
  apd = data.input$APD,
  bcl = data.input$BCL
)

p <- ggplot(data.plot, aes(bcl, apd)) +
  geom_point(size = 0.5) +
  scale_x_continuous(
    breaks = seq(0, max(data.plot$bcl), by = 200)) +
  scale_y_continuous(
    breaks = round_any(
      seq(min(data.plot$apd), max(data.plot$apd), by = 20), 10)) +
  xlab("Basic Cycle Length (ms)") + 
  ylab("Action Potential Duration (ms)") +
  # Theme
  theme_hc() +
  theme(text = element_text(size = 10),
        axis.line.y = element_line(size = 0.5, color = "black"),
        axis.line.x = element_line(size = 0.5, color = "black"),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        axis.ticks.y = element_line(size = 0.5, color = "black"),
        axis.ticks.x = element_line(size = 0.5, color = "black"),
        axis.ticks.length = unit(5, "point"),
        panel.grid.major = element_line(size = 0.25, linetype = "dashed",
                                        color = "black"),
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical"
  )

ggsave(p, filename = "Restitution.pdf", 
       family = "Arial", 
       width = 7.5, height = 5, 
       units = "in", useDingbats = FALSE)