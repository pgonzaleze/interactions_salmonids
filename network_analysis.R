
############ =========================== ############
########### Network analysis on salmonids ###########
############ =========================== ############

#####
# The University of British Columbia
# Institute for the Oceans and Fisheries
# Author: Pedro G. González-Espinosa
# Date: 20/ July /2022

#####
#Load libraries
library(tidyverse) 
library(readxl) 
library(igraph) 
library(dplyr)  
library(ggplot2)  
library(esc) # Hedges´d
library(forcats) # reorder axis of plots

# load database
int_db <- read_excel('Hedges.xlsx') # Former "interaction_database"
int_db$Paired <- paste(int_db$`Stressor A`, ' - ', int_db$`Stressor B`) # combine columns

# Count the co-occurrences, interactions and species
cooc <-  int_db %>% 
  count(Paired)
species <-  int_db %>% 
  count(Organism)
interaction <- int_db %>% 
  count(Interaction)
year <- int_db %>% 
  count(Year)
life_stage <- int_db %>% 
  count(Stage)

# Create a new column combining stressors as a string 
int_db$PairedJoin <- paste(int_db$`Stressor A`, ' ',int_db$`Stressor B`)
# count the co-occurences from a string (from previous step)
coocJoin <-  int_db %>% 
  count(PairedJoin)
# separate the columns to analyse the data
coocount <- coocJoin %>% 
  separate(PairedJoin, c('Stressor A', 'Stressor B'))

# create a igraph object
g <- graph.data.frame(coocount)
# vcount(g) # count vertices
# ecount(g) # count edges
# # some meassurements
# edge.betweenness(g) # measure of how important the node is to the flow of information through a network
# edge_density(g)
# graph.density(g)

as_data_frame(g, what="edges")  # co-occurrence of paired variables
as_data_frame(g, what="vertices")
# account for the # of co-occurrences
#E(g)$width <- coocount$n/1.5 # Num of cooc relates to the width of the line


####################### ================== ###################
#######################       PLOTS        ###################
####################### ================== ###################

# tempcenter <- layout_as_star(g, center=V(g)[[10]]) # centre temperature 
# set.seed(202)
# plot.igraph(g, edge.arrow.size=0, vertex.color="gold", vertex.size=20, 
#             vertex.frame.color="blue", vertex.label.color="black", 
#             vertex.label.cex=0.4, vertex.label.dist=0, 
#             layout=tempcenter, edge.curved=0.2) #layout=layout_in_circle

set.seed(105)
id <- plot(g, edge.arrow.size=0, vertex.color="gold", vertex.size=20, 
             vertex.frame.color="gray", vertex.label.color="black", 
             vertex.label.cex=1, vertex.label.dist=0, 
             layout=layout_with_gem, edge.curved=0.2)

# Basic barplot ordered from lowest frequency to highest frequency
g <- ggplot(int_db) +
  theme_minimal()

# Number of species:
# use "fct_rev" to reverse the order to from low to high
g + geom_bar(aes(y = fct_rev(fct_infreq(Organism)), fill = Organism)) +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = 'none') +
  labs(y = "Species")

# Number of interaction type:
g + geom_bar(aes(x = fct_rev(fct_infreq(Interaction)), fill = Interaction)) +
  scale_fill_brewer(palette = 'Blues') +
  #scale_fill_grey(start = 0.25, end = 0.75) +
  labs(x = "Interaction")

# Number of paired variables:
g + geom_bar(aes(y = fct_rev(fct_infreq(Paired))), fill='gray50') +
  theme(legend.position = 'none') +
  labs(y = "Paired stressors")

# by year
g + geom_bar(aes(Year))

# by life stage
g <- int_db %>%
  ggplot(aes(x = fct_rev(fct_infreq(Stage)), fill=Stage)) + 
  theme_minimal()
g +  geom_bar() +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = 'none') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "Life stage")


################################################################
#############       Stressor - stressor        #################
################################################################

# subset to stressor-stressor relation dataset
int_db_SS <-  int_db %>% 
  subset(Effect=='Stressor-stressor')

# count the co-occurences from a string (from previous step)
coocJoin_SS <-  int_db_SS %>% 
  count(PairedJoin)
# separate the columns to analyse the data
coocount_SS <- coocJoin_SS %>% 
  separate(PairedJoin, c('Stressor A', 'Stressor B'))

# create a igraph object for stressor-stressor dataseset
ss <- graph.data.frame(coocount_SS)
vcount(ss) # count vertices
ecount(ss) # count edges
# some measurements
betweenness(ss) # measure of how important the node is to the flow of information through a network
edge_density(ss)
graph.density(ss)

as_data_frame(ss, what="edges")  # co-occurrence of paired variables
as_data_frame(ss, what="vertices")
# account for the # of co-occurrences
E(ss)$width <- coocount_SS$n # Num of cooc relates to the width of the line

# The function degree() has a mode of in for in-degree, out for out-degree, and all or total for total degree.
degin <- degree(ss, mode="in") # nodes that influence a particular stressor (most influenced)
degout <- degree(ss, mode="out") # nodes that a particular stressor can influence (most influential)
degall <- degree(ss, mode="all")

####################### ================== ###################
#######################       PLOTS        ###################
####################### ================== ###################

set.seed(45)
plot(ss, vertex.size=degin*2, edge.arrow.size=.75, edge.curved=0.25,
     vertex.label.cex=1, vertex.label.dist=1, layout=layout.gem)

set.seed(45)
plot(ss, vertex.size=degout*2, edge.arrow.size=.75, edge.curved=0.25,
     vertex.label.cex=1, vertex.label.dist=1, layout=layout.gem)

# set.seed(1984)
# plot(ss, vertex.size=degall*3, edge.arrow.size=.7, edge.curved=0.25,
#       vertex.label.cex=0.8, vertex.label.dist=1, layout=layout_on_grid)
