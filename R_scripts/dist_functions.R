library(Hmisc)
library(ggplot2)
library(dplyr)
library(reshape2)
data.dir <- '/Volumes/GenomeDK/Div/MutationAnalysisPipeline/data/'
#data.dir <- '/Users/besen/Projects/MutationAnalysisPipeline/data/'

read_dist_data <- function(dataset, n_random){
  filename <- paste('dist_btwn_mut_concatted_chrom_callability_', n_random,'.dat', sep="")
  distfile <- file.path(data.dir,dataset,'dist',filename)  
  read.table(distfile, head=T)
}

read_pair_dist_data <- function(dataset, n_random){
  filename <- 'all_pairwise_w_anno.txt'
  distfile <- file.path(data.dir,dataset,'dist',filename)
  read.table(distfile, head=T)#head=F, col.name=c("dist","in.same"))
}

plot_dist_btwn_mut <- function(datasets, n_random=500, labels=NULL){
  d <- NULL
  for (dataset in datasets){
    d.next <- read_dist_data(dataset, n_random)
    d.next$dataset <- dataset
    d <- rbind(d, d.next)
  }
  p <- ggplot(d, aes(y=obs,x=median.random, color=type)) + geom_point() + 
    scale_x_log10("Expected distance between mutations") + 
    scale_y_log10("Observed distance between mutations") +
    theme_bw() + theme(legend.title=element_text(size=9)) +
    scale_colour_discrete(name="Distances between\nmutations in:",
                          breaks=c("any", "same", "notsame"),
                          labels=c("Any individual", "Same individual", "Different individual")) 
  if(length(datasets)==1){
    p #+ theme(legend.justification=c(0,0), legend.position=c(0.5,0.5))
  } else {
    p + facet_wrap(~dataset, scales = "free") 
  }
}

plot_fraction_in_same <- function(datasets, dist.breaks){
  d <- NULL
  for (dataset in datasets){
    d.next <- read_pair_dist_data(dataset)
    d.next$dataset <- dataset
    d <- rbind(d, d.next)
  }
  p <- d %>% 
    mutate(dist.bin=cut(dist,breaks=dist.breaks,include.lowest=T,dig.lab=5)) %>%
    group_by(dataset, dist.bin) %>%
    summarise(n_insame = sum(in.same==1),
              n_notinsame = sum(in.same==0),
              fraction = mean(in.same==1)) %>%
    mutate(ymin=binconf(n_insame,n_insame+n_notinsame)[,2],
           ymax=binconf(n_insame,n_insame+n_notinsame)[,3]) %>%
    ggplot(aes(x=dist.bin,y=fraction, ymin=ymin, ymax=ymax)) + 
    #geom_bar(fill="lightblue",stat="identity",position="dodge") + 
    geom_point() + 
    geom_errorbar(width=0.5) +
    geom_hline(yintercept=1/49, linetype="dashed") +
    #scale_y_log10("Fraction of pairs of mutations\nwhere both mutations occured\nin the same individual",breaks=c(0.01,0.02,0.05,0.1,0.2,0.5,1.0)) +
    scale_y_log10("Fraction of mutations pairs",breaks=c(0.01,0.02,0.05,0.1,0.2,0.5,1.0)) +
    scale_x_discrete("Distance between mutations") +
    theme_bw() + theme(axis.text.x=element_text(angle=45, hjust = 1))
  if(length(datasets)==1){
    p
  } else {
    p+facet_grid(~dataset)
  }
}

#d <- read_pair_dist_data("decode_september")
# dist.breaks <- c(0,1,10,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000,16000,20000,25000,30000,40000,50000,100000,300000,400000,500000,1e9)
#dist.breaks <- c(0,1,10,100,1000,5000,10000,25000,50000,100000,200000,300000,Inf)
# 
#   
# #facet_grid(~dataset) +
# #d.next <- read_dist_data("kong_2012", 500)
# #plot_dist_btwn_mut(c("kong_2012","michaelson_2012"), 500)
#dist.breaks <- c(0,1,10,100,1000,5000,10000,25000,50000,100000,200000,300000,Inf)
# plot_fraction_in_same(c("decode_2014","kong_2012"), dist.breaks)
#plot_fraction_in_same(c("decode_september"), dist.breaks)
#plot_dist_btwn_mut(c("decode_september")) + theme(legend.justification=c(1,0), legend.position=c(0.5,0.5))

#do_plot <- function(d, dist.breaks){

# d2<- d %>% 
#   mutate(dist.bin=cut(dist,breaks=dist.breaks,include.lowest=T,dig.lab=5)) %>%
#   group_by(dataset, dist.bin) %>%
#   summarise(n_insame = sum(in.same==1),
#             n_notinsame = sum(in.same==0),
#             fraction = mean(in.same==1)) #%>%
# levels(d2$dist.bin)[12] <- "(3e+05,Inf]"
# d2 %>% mutate(ymin=binconf(n_insame,n_insame+n_notinsame)[,2],
#               ymax=binconf(n_insame,n_insame+n_notinsame)[,3]) %>%
#     ggplot(aes(x=dist.bin,y=fraction, ymin=ymin, ymax=ymax)) + 
#     geom_point() + 
#     geom_errorbar(width=0.5) +
#     geom_hline(yintercept=1/282, linetype="dashed") +
#     scale_y_log10("Fraction of mutation pairs",breaks=c(0.01,0.02,0.05,0.1,0.2,0.5,1.0)) +
#     scale_x_discrete("Distance between mutations (bp)") +
#     theme_bw() + theme(axis.text.x=element_text(angle=45, hjust = 1))
# 
# 
# do_plot(d, dist.breaks)
