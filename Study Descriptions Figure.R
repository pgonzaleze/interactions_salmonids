
############ =========================== ############
########### Network analysis on salmonids ###########
############ =========================== ############

#####
# The University of British Columbia
# Institute for the Oceans and Fisheries
# Author: Sara Cannon
# Created: 31 August, 2023
# Last update: 31 August, 2023

#load libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(readxl) 
library(RColorBrewer)
library(forcats)

#load and explore data
int_db <- read_excel('Hedges.xlsx')
head(int_db)

int_db$`Type of study` <- as.factor(int_db$`Type of study`)
summary(int_db$`Type of study`)

int_db$Organism <- as.factor(int_db$Organism)
summary(int_db$Organism)

int_db$Stage <- as.factor(int_db$Stage)
summary(int_db$Stage)

#create new column for common names (to use for plot instead of full organism name)
int_db[c('Common Name', 'Scientific Name')] <- str_split_fixed(int_db$Organism, pattern = regex("\\([^)]*\\)"), 2)
head(int_db)
int_db$`Common Name`
int_db$`Scientific Name` #somethign wrong with the regex code above, but we will only use common name for now anyway
#would need to correct if using scientific names on the plot instead

#reorder "stage" for correct position in stacked bars (e.g. want to move "Not specified" to the end)
new_stage_levels <- c(levels(int_db$Stage)[levels(int_db$Stage) != "Not specified"], "Not specified")
int_db$Stage <- factor(int_db$Stage, levels = new_stage_levels)
head(int_db)

#plot
ggplot(int_db) +
  aes(x = forcats::fct_infreq(`Common Name`), fill = Stage) +
  geom_bar(color = "black", linewidth = 0.25) +
  scale_fill_brewer(palette = "Spectral") +
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 3
  ) +
  labs(
    x = "Species (Common Names)",
    y = "Studies (n)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(1,1), 
        legend.justification = c(1,1)
        )

