#!/usr/bin/env Rscript
library(Hmisc)
library(tidyverse)

infile <- commandArgs(TRUE)[1]
n_samples <- as.integer(commandArgs(TRUE)[2])
outfile <- commandArgs(TRUE)[3]

d <- read_table(infile)

dist.breaks <- c(0,1,10,100,1000,5000,10000,25000,50000,100000,200000,300000,Inf)

p <- d %>% 
  mutate(dist.bin=cut(dist,breaks=dist.breaks,include.lowest=T,dig.lab=5)) %>%
  group_by(dist.bin, insame) %>%
  summarise(count=sum(count)) %>%
  group_by(dist.bin) %>%
  pivot_wider(names_from = insame, values_from = count, values_fill = 0) %>%
  rename(n_insame = `1`, n_notinsame = `0`) %>%
  mutate(fraction = n_insame / (n_insame + n_notinsame)) %>%
  mutate(ymin=binconf(n_insame,n_insame+n_notinsame)[,2],
         ymax=binconf(n_insame,n_insame+n_notinsame)[,3]) %>%
  ggplot(aes(x=dist.bin,y=fraction, ymin=ymin, ymax=ymax)) + 
  geom_point() + 
  geom_errorbar(width=0.5) +
  geom_hline(yintercept=1/(n_samples-1), linetype="dashed") +
  scale_y_log10("Fraction of pairs of mutations\nwhere both mutations occured\nin the same individual") +
  scale_x_discrete("Distance between mutations") +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust = 1))
    
ggsave(outfile, p, width=10, height=5)
