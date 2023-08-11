
############ =========================== ############
########### Network analysis on salmonids ###########
############ =========================== ############

#####
# The University of British Columbia
# Institute for the Oceans and Fisheries
# Author: Pedro G. Gonz?lez-Espinosa
# Date: 13/ AUG /2022

#####
#Load libraries
library(tidyverse) 
library(readxl) 
library(igraph) 
library(dplyr)  
library(ggplot2)  
library(esc) # Hedges?d
library(tidyr)
library(forcats) # reorder axis of plots
library(imputeTS) # to change NAs to any other number
library(FSA) # Dunn Test
library(metafor) # for mixed-effects meta-regression model & Rank Correlation Test
library(ggwordcloud) # make a word cloud image
library(compute.es)

# load database
int_db <- read_excel('Hedges.xlsx')

# Subset the data frame to compute the effect size 
int_db <- int_db[1:59, ] # adjust according to the length of those records with full data 

# save lists of medias, SD and sample sizes, use "unlist" to  convert a list to vector
media1 <- as.numeric(unlist(int_db[,12]))   
media2 <- as.numeric(unlist(int_db[,13]))
sd1 <- as.numeric(unlist(int_db[,14]))
sd2 <- as.numeric(unlist(int_db[,15]))
replicates <- as.numeric(unlist(int_db[,21]))
replicates <- na_replace(replicates, 1) # change NA?s to 1
sample_ctrl <- as.numeric(unlist(int_db[,19])) #* replicates 
sample_n <- as.numeric(unlist(int_db[,20])) #* replicates

# compute the Hedges?d
g=(esc_mean_sd(grp1m=media1,grp1sd=sd1,grp1n=sample_n,
               grp2m=media2,grp2sd=sd2,grp2n=sample_n, es.type = 'g'))

effect_size = as.data.frame(cbind(g$es, g$ci.lo, g$ci.hi)) 
colnames(effect_size) <- c('ES', 'CI_lo', 'CI_hi')  # rename columns

# paste both dataframes original dataframe with the effects size and CI dataframe
hedges <- cbind(int_db, effect_size)

# add data of two cases from Jackson et al. 2015 
hedges[58,'ES'] = -0.10 
#hedges[59,'ES'] = -0.43
hedges[58,'CI_lo'] = -0.25  
#hedges[59,'CI_lo'] = -0.74
hedges[58,'CI_hi'] =  0.04
#hedges[59,'CI_hi'] = -0.13

### Compute the overall effect sizes
AveInt = hedges %>% 
  summarise(ES = mean(ES),
            CI_lo = mean(CI_lo),
            CI_hi = mean(CI_hi))
AveInt$Author <- "Overall"
AveInt$study <- "Overall"
hedges <- full_join(hedges, AveInt)

# Append a column of the estimated "interaction". A condition is applied.
hedges <- hedges %>%
  add_column(EstimInter = if_else(.$'CI_lo' < 0 & .$'CI_hi' > 0, 'Additive', 
                                 if_else(.$'CI_lo' < 0 & .$'CI_hi' < 0, 'Antagonistic', 'Synergistic')))
hedges$index <- 1:nrow(hedges) # create an index column
hedges$study <- str_c(hedges$Author, ' S' ,hedges$index) # create a column with author and index ID

###########################################
### create general forest plot by study ###
###########################################
ggplot(data=hedges, aes(y=reorder(as.factor(study),-ES), x=ES, xmin=CI_lo, xmax=CI_hi, 
                        colour=EstimInter)) +
  geom_point() + 
  geom_errorbarh(height=.25) +
  labs(x='Effect Size', y = 'Study') + #title='Effect Size by Study'
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic() +
  theme(legend.position = 'top')

### Average the effects sizes by grouping the interactions ###
SumInt = hedges %>% group_by(EstimInter) %>%
  summarise(ES_avg = mean(ES),
            CI_lo_avg = mean(CI_lo),
            CI_hi_avg = mean(CI_hi),
            .groups = 'drop')
#View(SumInt)

ggplot(data=SumInt, aes(y=reorder(as.factor(EstimInter),-ES_avg), x=ES_avg, xmin=CI_lo_avg, xmax=CI_hi_avg, 
                        colour=EstimInter)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  #scale_x_continuous(limits=c(-10,10)) +
  #scale_y_discrete(breaks=1:nrow(hedges), labels=hedges$study) +
  labs(x='Effect Size', y = 'Type') + #title='Pooled Effect Size',
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic() +
  theme(legend.position = 'none')

#### Average the effects sizes by grouping the interactions by species ###
SumIntOrg = hedges %>% group_by(EstimInter, Organism) %>%
  summarise(ES_avg = mean(ES),
            CI_lo_avg = mean(CI_lo),
            CI_hi_avg = mean(CI_hi),
            .groups = 'drop')
SumIntOrg$index <- 1:nrow(SumIntOrg) # create a column based on the index number
SumIntOrg$OrgID <- str_c(SumIntOrg$Organism, ' S' ,SumIntOrg$index)
#View(SumIntOrg)

ggplot(data=SumIntOrg, aes(y=reorder(as.factor(OrgID),-ES_avg), x=ES_avg, xmin=CI_lo_avg, xmax=CI_hi_avg, 
                        colour=EstimInter)) +
  geom_point() + 
  geom_errorbarh(height=.1) +
  #scale_x_continuous(limits=c(-10,10)) +
  #scale_y_discrete(breaks=1:nrow(hedges), labels=hedges$study) +
  labs(title='Pooled Effect Size Species', x='Effect Size', y = 'Species') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_classic()

# Basic barplot of estimated interactions ordered from lowest frequency to highest frequency
g <- ggplot(hedges) +
  theme_minimal()
# Number of interaction type:
g + geom_bar(aes(x = fct_rev(fct_infreq(EstimInter)), 
                 colour=EstimInter, fill = EstimInter)) +
  theme(legend.position = "none") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  #scale_fill_brewer(palette = 'Blues') +
  #scale_fill_grey(start = 0.25, end = 0.75) +
  labs(x = "Interaction")


##################### ====================== ###################
#####################       Map minder       ###################
##################### ====================== ###################
# Subset the data frame to compute the effect size 
hedges <- hedges[1:59, ] # adjust according to the length of those records with full data 
# Create a new column combining stressors as a string 
hedges$PairedJoin <- paste(hedges$`Stressor A`, "-",hedges$`Stressor B`, sep="")
hedges$InterTypePaired <- paste(hedges$EstimInter, ' ',hedges$PairedJoin)
hedges[1,'InterTypePaired'] = "Antagonistic   Temperature-DO" # match names and order of stressors
hedges[2,'InterTypePaired'] = "Synergistic   Temperature-DO" # match names and order of stressors
hedges[3,'InterTypePaired'] = "Additive   Temperature-DO" # match names and order of stressors
hedges[11,'InterTypePaired'] = "Synergistic   Temperature-Metal" # match names and order of stressors
hedges[12,'InterTypePaired'] = "Synergistic   Temperature-Metal" # match names and order of stressors
hedges[13,'InterTypePaired'] = "Synergistic   Temperature-Pesticides" # match names and order of stressors
hedges[25,'InterTypePaired'] = "Synergistic   Temperature-Nitrate" # match names and order of stressors
hedges[51,'InterTypePaired'] = "Synergistic   Temperature-pH" # match names and order of stressors
hedges[56,'InterTypePaired'] = "Synergistic   Temperature-pH" # match names and order of stressors

coocInterType <-  hedges %>% 
  count(InterTypePaired)
# separate the columns to analyse the data
coocount <- coocInterType %>% 
  separate(InterTypePaired, c('interaction type', "pair"), sep = "^\\S*\\K\\s+")

# subset to interaction type
additive_df <-  coocount %>% 
  subset(`interaction type`=='Additive')
synergistic_df <-  coocount %>% 
  subset(`interaction type`=='Synergistic')
antagonistic_df <-  coocount %>% 
  subset(`interaction type`=='Antagonistic')

### plot paired stressors by interaction type
g <- graph.data.frame(additive_df)
E(g)$width <- additive_df$n*2 # Num of cooc relates to the width of the line
set.seed(105)
id <- plot(g, edge.arrow.size=0, vertex.color="#FF6666", 
           vertex.frame.color="gray", vertex.label.color="black", 
           vertex.label.cex=1, vertex.label.dist=0,
           layout=layout_with_dh)

g <- graph.data.frame(antagonistic_df)
E(g)$width <- antagonistic_df$n*2 # Num of cooc relates to the width of the line
set.seed(105)
id <- plot(g, edge.arrow.size=0, vertex.color="#33FF00", 
           vertex.frame.color="gray", vertex.label.color="black", 
           vertex.label.cex=1, vertex.label.dist=0,
           layout=layout_with_dh)

g <- graph.data.frame(synergistic_df)
E(g)$width <- synergistic_df$n*2 # Num of cooc relates to the width of the line
set.seed(105)
id <- plot(g, edge.arrow.size=0, vertex.color="#3399CC", 
           vertex.frame.color="gray", vertex.label.color="black", 
           vertex.label.cex=1, vertex.label.dist=0,
           layout=layout_with_dh)


############################################################
##############  Kruskal-Wallis and Dunn test   #############
############################################################

 
kruskal.test(n ~ as.factor(`interaction type`), data=coocount)
dunnTest(n ~ as.factor(`interaction type`), data=coocount, method = "bonferroni")
cor.test(x=hedges$`Sample size`, y=hedges$ES, method = 'spearman', exact = FALSE) # test correlation between effects size and sample size


# Subset the data frame to compute the effect size 
hedges <- hedges[1:59, ] # adjust according to the length of those records with full data 

# comparing multiple columns
as.data.frame(table(hedges$`Stressor A`))
as.data.frame(table(hedges$`Stressor B`))

(total_obs <- as.matrix(table(hedges$`Stressor B`, hedges$`Stressor A`)))
#write.csv(total_obs, "observation_matrix.csv")
#sum(total_obs)

### most studiesd pair of stressors
hedges[1,'PairedJoin'] = "Temperature-DO" # match names and order of stressors
hedges[2,'PairedJoin'] = "Temperature-DO" # match names and order of stressors
hedges[3,'PairedJoin'] = "Temperature-DO" # match names and order of stressors
hedges[11,'PairedJoin'] = "Temperature-Metal" # match names and order of stressors
hedges[12,'PairedJoin'] = "Temperature-Metal" # match names and order of stressors
hedges[13, 'PairedJoin'] = "Temperature-Pesticides" # match names and order of stressors
hedges[25, 'PairedJoin'] = "Temperature-Nitrate" # match names and order of stressors
hedges[51,'PairedJoin'] = "Temperature-pH" # match names and order of stressors
hedges[56,'PairedJoin'] = "Temperature-pH" # match names and order of stressors

coocPair <-  hedges %>% 
  count(PairedJoin)
#write.csv(coocPair, "paireedCount.csv")


##### =========== #####
##### Funnel plot #####
##### =========== #####

# compute the SE of the sample
funnel_plot_data <- hedges
funnel_plot_data <- funnel_plot_data[1:57, ]
funnel_plot_data$SE_funnel <- (as.numeric(funnel_plot_data$SD2) / sqrt(funnel_plot_data$`Sample size`))
funnel_plot_data$variance <- as.numeric(funnel_plot_data$SD2)^2
colnames(funnelDF) <- c("SE", "ES", "SS")  # Standard error SE; Effect size ES
funnel_plot_data$ES <- hedges_g(funnelDF$ES, funnelDF$SS)
funnelDF <- as.data.frame(cbind(funnel_plot_data$SE_funnel, funnel_plot_data$ES, 
                                funnel_plot_data$`Sample size`))

funnelDF[28, 'SE'] = 0


funnel(g$es, g$se, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), legend=TRUE)
regtest(g$es, sei=g$se)
ranktest(g$es, sei=g$se)


######### ================= #########
######### word cloud image  #########
######### ================= #########

drivers <-c("Temperature", "Dissolved Oxygen", "Pesticides", "Metal", 
            "Salinity", "Nitrate", "Aerial Expossure", "chemical", "Density",
            "pH", "Radiation", "Suspended Sediments", "UVB", "Ammonia",
            "Dissease", "Water Quality", "Invasive Species", "Water Flow", 
            "Cyanobacteria")
weight <- c(61, 23, 17, 16, 15, 10, 10, 10, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8)

driverDf <- cbind(drivers, weight)

driverDf <- driverDf %>%
  mutate(angle = 15 * sample(-2:2, n(), replace = TRUE, prob = c(1, 1, 4, 1, 1)))

set.seed(1984)
ggplot(drivers, aes(label = drivers, size = weight,
                    color = factor(sample.int(10, nrow(driverDf), replace = TRUE)),
                    angle = driverDf$angle)) +
  geom_text_wordcloud(eccentricity = 0.5) +
  scale_size_area(max_size = 10) +
  theme_minimal()


######### =================================================== ###########
#########                  Descriptive plots                  ###########
######### =================================================== ###########

# Number of species:
# use "fct_rev" to reverse the order to from low to high
# Basic barplot ordered from lowest frequency to highest frequency
plt <- ggplot(hedges) +
  theme_minimal()

plt + geom_bar(aes(y = fct_rev(fct_infreq(Organism)), fill = Organism)) +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = 'none') +
  labs(y = "Species")

# by life stage ()
plt <- hedges %>%
  ggplot(aes(x = fct_rev(fct_infreq(Stage)), fill=Stage)) + 
  theme_minimal()
plt +  geom_bar() +
  scale_fill_brewer(palette = "Spectral") +
  theme(legend.position = 'none') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "Life stage")


stage_count <- hedges %>% 
  count(Stage)

species_count <-  hedges %>% 
  count(Organism)



