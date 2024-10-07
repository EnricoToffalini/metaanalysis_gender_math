
####################################################

#### INITIAL SETTINGS

# Empty workspace
rm(list=ls())

# Load custom functions
source("_Custom functions.R")

# Load necessary packages
install_and_load("readxl")

####################################################

#### IMPORT AND FILTER DATA

# Import and prepare data; calculate Cohen's d
d = data.frame(read_excel("math_gender_full_dataset_effects_to_be_computed.xlsx"))
d$N = as.numeric(d$N)
d$GGGI = as.numeric(d$GGGI)
d$GradeM = as.numeric(d$GradeM)
d$id_effect = 1:nrow(d)

# compute Cohen's d for all rows and determine effect size sign
vars = c("male_mean", "male_sd", "fem_mean", "fem_sd")
for (v in vars) d[, v] = as.numeric(d[, v])
x = cd(d$male_mean, d$male_sd, d$male_N, d$fem_mean, d$fem_sd, d$fem_N)
d$eff = x$eff
d$vi = x$vd
d$vi_0 = x$vd.0
d$eff[d$Sign_score == -1] = d$eff[d$Sign_score == -1] * -1 # revert sign where necessary

# approximate confidence intervals
d$lb = d$eff + sqrt(d$vi) * qt(.025, d$N)
d$ub = d$eff + sqrt(d$vi) * qt(.975, d$N)

# see number of studies, samples, and rows
length(unique(d$ID_ARTICLE)); length(unique(d$ID_SAMPLE)); nrow(d)

####################################################

#### EXPORT DATA
write.table(d, file="math_gender_full_dataset.csv",sep=",",row.names=F)

####################################################

#### DESCRIPTIVE STATISTICS

# Aggregate effects
sum(aggregate(d$fem_N, by = list(d$ID_SAMPLE), FUN = max)$x)
sum(aggregate(d$male_N, by = list(d$ID_SAMPLE), FUN = max)$x)
median(aggregate(rowSums(d[,c("fem_N","male_N")]), by = list(d$ID_SAMPLE), FUN = max)$x)
median(aggregate(rowSums(d[,c("fem_N","male_N")]), by = list(d$ID_ARTICLE), FUN = max)$x)

# Mean and standard deviation of mean age
x = d; x$MeanAge_years = as.numeric(d$MeanAge_months)/12
tb = aggregate(x$MeanAge_years, by = list(x$ID_SAMPLE), FUN = max)
(mean_age = mean(tb$x, na.rm = TRUE))
(sd_age = sd(tb$x, na.rm = TRUE))

# Descriptive statistics for moderators
(country_code_count = colSums(table(d$ID_SAMPLE, d$GeoArea) != 0) )

gradeM_count = colSums(table((d$ID_SAMPLE), d$GradeM) != 0)
sum(gradeM_count[names(gradeM_count) == 0])
sum(gradeM_count[names(gradeM_count) %in% seq(1,5,.5)])
sum(gradeM_count[names(gradeM_count) %in% seq(6,8,.5)])
sum(gradeM_count[names(gradeM_count) %in% seq(9,13,.5)])
sum(gradeM_count[names(gradeM_count) == 14 ])

(math_content_count = colSums(table((d$ID_SAMPLE), d$MathContent) != 0) )


####################################################


