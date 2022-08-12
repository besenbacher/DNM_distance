#!/usr/bin/env Rscript
library(tidyverse)

infile <- commandArgs(TRUE)[1]
outfile <- commandArgs(TRUE)[2]

d <- read.table(infile, head=T)

plot_dist_btwn_mut <- function(d){
  p <- ggplot(d, aes(y=obs,x=median.random, color=type)) + 
    geom_abline(slope=1, linetype="dashed") +
    geom_point() + 
    scale_x_log10("Expected distance between mutations") + 
    scale_y_log10("Observed distance between mutations") +
    theme_bw() + theme(legend.title=element_text(size=9)) +
    scale_colour_discrete(name="Distances between\nmutations in:",
                          breaks=c("any", "same", "notsame"),
                          labels=c("Any individual", "Same individual", "Different individual")) 
  p + theme(legend.justification=c(0,0), legend.position=c(0.1,0.7)) +
    theme(legend.background = element_rect(size=0.2, linetype="solid", color="black"))
}

p <- plot_dist_btwn_mut(d)

ggsave(outfile, p, width=5, height=5)
