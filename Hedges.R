
############ =========================== ############
########### Network analysis on salmonids ###########
############ =========================== ############

#####
# The University of British Columbia
# Institute for the Oceans and Fisheries
# Author: Pedro G. Gonzalez-Espinosa
# Created: 13/ AUG /2022
# Last update: 21/ AUG /2023

#####
#Load libraries
library(tidyverse) 
library(readxl) 
library(igraph) 
library(dplyr)  
library(ggplot2)  
library(esc) # To compute Hedges d. It wont work using R Version 4.3.1 but the issue has been fixed in this script
library(tidyr)
library(forcats) # reorder axis of plots
library(imputeTS) # to change NAs to any other number
library(FSA) # Dunn Test
library(metafor) # for mixed-effects meta-regression model & Rank Correlation Test
#library(ggwordcloud) # make a word cloud image
#library(compute.es) 

# load database
int_db <- read_excel('Hedges.xlsx')

# Subset the data frame to compute the effect size 
int_db <- int_db[1:59, ] # adjust according to the length of those records with full data 

# save lists of medias, SD and sample sizes, use "unlist" to  convert a list to vector
media1 <- as.numeric(unlist(int_db[,13]))   
media2 <- as.numeric(unlist(int_db[,14]))
sd1 <- as.numeric(unlist(int_db[,15]))
sd2 <- as.numeric(unlist(int_db[,16]))
replicates <- as.numeric(unlist(int_db[,19]))
replicates <- na_replace(replicates, 1) # change NA?s to 1
sample_ctrl <- as.numeric(unlist(int_db[,17])) #* replicates 
sample_n <- as.numeric(unlist(int_db[,18])) #* replicates


##### Hot fix to avoid error using R Version 4.3.1 ####
# Functions extracted from https://github.com/strengejacke/esc

# Compute variance of d-type effect size
esc.vd <- function(d, grp1n, grp2n) {
  (grp1n + grp2n) / (grp1n * grp2n) + (d * d) / (2 * (grp1n + grp2n))
}

# 95% confidence interval
#' @importFrom stats qnorm
lower_d <- function(d, v) d - stats::qnorm(.975) * sqrt(v)
upper_d <- function(d, v) d + stats::qnorm(.975) * sqrt(v)

# small sample size bias correction
sssbc <- function(totaln) return(1 - (3 / (4 * totaln - 9)))

# generic conversion function
esc_generic <- function(es, v, grp1n, grp2n, es.type, info, study) {
  # compute total n
  totaln <- grp1n + grp2n
  
  
  # return effect size as odds ratio
  
  if (es.type == "or")
    return(convert_d2or(
      d = es,
      v = v,
      totaln = totaln,
      es.type = "logit",
      info = paste0(info, " to effect size odds ratios"),
      study = study
    ))
  
  
  # return effect size as cox odds ratio
  
  if (es.type == "cox.or")
    return(convert_d2or(
      d = es,
      v = v,
      totaln = totaln,
      es.type = "cox",
      info = paste0(info, " to effect size Cox odds ratios"),
      study = study
    ))
  
  
  # return effect size as f
  
  if (es.type == "f")
    return(convert_d2f(
      d = es,
      v = v,
      totaln = totaln,
      info = paste0(info, " to effect size f"),
      study = study
    ))
  
  
  # return effect size as eta squared
  
  if (es.type == "eta")
    return(convert_d2etasq(
      d = es,
      v = v,
      grp1n = grp1n,
      grp2n = grp2n,
      info = paste0(info, " to effect size eta squared"),
      study = study
    ))
  
  
  # return effect size as log odds
  
  if (es.type == "logit")
    return(convert_d2logit(
      d = es,
      v = v,
      totaln = totaln,
      es.type = "logit",
      info = paste0(info, " to effect size logits"),
      study = study
    ))
  
  
  # return effect size as cox log odds
  
  if (es.type == "cox.log")
    return(convert_d2logit(
      d = es,
      v = v,
      totaln = totaln,
      es.type = "cox",
      info = paste0(info, " to effect size Cox logits"),
      study = study
    ))
  
  
  # return effect size as correlation
  
  if (es.type == "r")
    return(convert_d2r(
      d = es,
      v = v,
      grp1n = grp1n,
      grp2n = grp2n,
      info = paste0(info, " to effect size correlation"),
      study = study
    ))
  
  
  # return d or Hedges' g
  
  if (es.type == "d") {
    info.suffix <- " to effect size d"
  } else if (es.type == "g") {
    info.suffix <- " to effect size Hedges' g"
    es <- hedges_g(d = es, totaln = grp1n + grp2n)
  }
  
  
  # return effect size as standardized mean difference d or Hedges' g
  
  structure(
    class = c("esc", "esc_d"),
    list(es = es, se = sqrt(v), var = v, ci.lo = lower_d(es, v), ci.hi = upper_d(es, v),
         w = 1 / v, totaln = totaln, measure = es.type,
         info = paste0(info, info.suffix), study = study
    ))
}

#' @title Compute effect size from Mean and Standard Deviation
#' @name esc_mean_sd
#'
#' @description Compute effect size from mean and either group-based standard
#'              deviations or full sample standard deviation.
#'
#' @param grp1m The mean of the first group.
#' @param grp1sd The standard deviation of the first group.
#' @param grp1n The sample size of the first group.
#' @param grp2m The mean of the second group.
#' @param grp2sd The standard deviation of the second group.
#' @param grp2n The sample size of the second group.
#' @param r Correlation for within-subject designs (paired samples, repeated measures).
#' @param totalsd The full sample standard deviation. Either \code{grp1sd} and
#'        \code{grp2sd}, or \code{totalsd} must be specified.
#'
#' @inheritParams esc_beta
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @note If \code{es.type = "r"}, Fisher's transformation for the effect size
#'       \code{r} and their confidence intervals are also returned.
#'
#' @examples
#' # with standard deviations for each group
#' esc_mean_sd(
#'   grp1m = 7, grp1sd = 2, grp1n = 50,
#'   grp2m = 9, grp2sd = 3, grp2n = 60,
#'   es.type = "logit"
#' )
#'
#' # effect-size d, within-subjects design
#' esc_mean_sd(
#'   grp1m = 7, grp1sd = 2, grp1n = 50,
#'   grp2m = 9, grp2sd = 3, grp2n = 60, r = .7
#' )
#'
#' # with full sample standard deviations
#' esc_mean_sd(grp1m = 7, grp1n = 50, grp2m = 9, grp2n = 60, totalsd = 4)
#'
#' @export
esc_mean_sd <- function(grp1m, grp1sd, grp1n, grp2m, grp2sd, grp2n, totalsd, r,
                        es.type = c("d", "g", "or", "logit", "r", "cox.or", "cox.log"), study = NULL) {
  es.type <- match.arg(es.type)
  
  # check if parameter are complete
#  if ((missing(totalsd) || is.null(totalsd) || any(is.na(totalsd))) &&
#      ((missing(grp1sd) || is.null(grp1sd) || any(is.na(grp1sd))) ||
#       (missing(grp2sd) || is.null(grp2sd) || any(is.na(grp2sd))))) {
#    warning("Either `totalsd` or both `grp1sd` and `grp2sd` must be specified.", call. = F)
#    return(esc_generic(es = NA, v = NA, es.type = es.type, grp1n = NA, grp2n = NA, info = NA, study = NA))
#  }
  
  # compute totaln, better overview
  totaln <- grp1n + grp2n
  
  # compute mean difference
  dm <- grp1m - grp2m
  info <- "mean and sd"
  
  # compute pooled standard deviation.
  if (!missing(totalsd) && !is.null(totalsd)) {
    # pooled sd from full sample sd, formula from book
    sdp <- ((totalsd^2 * (totaln - 1) - ((dm^2 * grp1n * grp2n) / totaln)) / (totaln - 1))
    
    # pooled sd from full sample sd, formula from unpublished manuscript. formulas vary,
    # email-correspondence with author suggests that book-formula should be correct
    # however, in some case value might be negative, so sqrt is not possible, use
    # alternative formula then
    if (sdp < 0)
      sdp <- (totalsd^2 * (totaln - 1) - ((grp1m^2 + grp2m^2 - 2 * grp1m * grp2m) / totaln)) / totaln
    
    sd_pooled <- sqrt(sdp)
  } else if (!missing(r)) {
    # pooled sd, within-subject
    sd_pooled <- sqrt(grp1sd^2 + grp2sd^2 - (2 * r * grp1sd * grp2sd))
    info <- "mean and sd (within-subject)"
  } else {
    # pooled sd from group sd's
    sd_pooled <- sqrt((grp1sd^2 * (grp1n - 1) + grp2sd^2 * (grp2n - 1)) / (grp1n + grp2n - 2))
  }
  
  # compute effect size
  es <- (grp1m - grp2m) / sd_pooled
  # compute variance
  v <- esc.vd(es, grp1n, grp2n)
  
  # return effect size
  esc_generic(
    es = es,
    v = v,
    es.type = es.type,
    grp1n = grp1n,
    grp2n = grp2n,
    info = info,
    study = study
  )
}


#' @title Compute effect size from Mean and Standard Error
#' @name esc_mean_se
#'
#' @description Compute effect size from Mean and Standard Error.
#'
#' @param grp1se The standard error of the first group.
#' @param grp2se The standard error of the second group.
#' @inheritParams esc_mean_sd
#'
#' @return The effect size \code{es}, the standard error \code{se}, the variance
#'         of the effect size \code{var}, the lower and upper confidence limits
#'         \code{ci.lo} and \code{ci.hi}, the weight factor \code{w} and the
#'         total sample size \code{totaln}.
#'
#' @note If \code{es.type = "r"}, Fisher's transformation for the effect size
#'       \code{r} and their confidence intervals are also returned.
#'
#' @references Lipsey MW, Wilson DB. 2001. Practical meta-analysis. Thousand Oaks, Calif: Sage Publications
#'             \cr \cr
#'             Wilson DB. 2016. Formulas Used by the "Practical Meta-Analysis Effect Size Calculator". Unpublished manuscript: George Mason University
#'
#' @examples
#' esc_mean_se(grp1m = 7, grp1se = 1.5, grp1n = 50,
#'             grp2m = 9, grp2se = 1.8, grp2n = 60, es.type = "or")
#'
#' @export
esc_mean_se <- function(grp1m, grp1se, grp1n, grp2m, grp2se, grp2n, r,
                        es.type = c("d", "g", "or", "logit", "r", "f", "eta", "cox.or", "cox.log"), study = NULL) {
  es.type <- match.arg(es.type)
  
  grp1sd <- grp1se * sqrt(grp1n - 1)
  grp2sd <- grp2se * sqrt(grp2n - 1)
  
  info <- "mean and se"
  
  if (!missing(r)) {
    # pooled sd, within-subject
    sd_pooled <- sqrt(grp1sd^2 + grp2sd^2 - (2 * r * grp1sd * grp2sd))
    info <- "mean and se (within-subject)"
  } else {
    sd_pooled <- sqrt((grp1sd^2 * (grp1n - 1) + grp2sd^2 * (grp2n - 1)) / (grp1n + grp2n - 2))
  }
  
  es <- (grp1m - grp2m) / sd_pooled
  v <- esc.vd(es, grp1n, grp2n)
  
  # return effect size
  esc_generic(
    es = es,
    v = v,
    es.type = es.type,
    grp1n = grp1n,
    grp2n = grp2n,
    info = info,
    study = study
  )
}


#####

# compute the Hedges?d
g=(esc_mean_sd(grp1m=media1,grp1sd=sd1,grp1n=sample_n,
               grp2m=media2,grp2sd=as.numeric(sd2),grp2n=sample_n, es.type = 'g'))

effect_size = as.data.frame(cbind(g$es, g$ci.lo, g$ci.hi)) 
colnames(effect_size) <- c('ES', 'CI_lo', 'CI_hi')  # rename columns

# paste both dataframes original dataframe with the effects size and CI dataframe
hedges <- cbind(int_db, effect_size)

# add data a single case from Jackson et al. 2015 
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

row_to_drop <- nrow(hedges)
hedges <- hedges[-row_to_drop, ] # drop row with overall effect to avoid wrong estimates

# This step is only as an explatory analysis. The meta-regression is below.
#funnel(g$es, g$se, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), legend=TRUE)
regtest(g$es, sei=g$se)
ranktest(g$es, sei=g$se)


##### ========================================================== ##### 
##### Meta-Regression to show the similarity between a subgroup  #####
##### ========================================================== #####

# Perform meta-regression analysis
# adjust the mods to analyse a desired subgroup using the "$" sign
meta_regression_model <- rma(yi = ES, vi = (CI_hi - CI_lo)^2 / (4 * log(2)),
                             mods = ~ hedges$EstimInter,
                             data = hedges)

# Print the summary of the meta-regression analysis
summary(meta_regression_model)


# generate a funnel plot for our meta-regression analysis output
funnel(meta_regression_model, xlab = "Hedges' g")
regtest(meta_regression_model)
ranktest(meta_regression_model)
funnel(meta_regression_model, level=c(90, 95, 99), 
       shade=c("white", "yellow", "orange"),
       xlab = "Hedges' g", legend=F)
