##load libraries
library(slider)
library(stats)
library(coda)
library(odin2)
library(dust2)
library(monty)
library(RPMG)
library(deSolve)
library(DescTools)
library(ggpmisc)
library(mgcv)
library(lubridate)
library(neonUtilities)
library(neonOS)
library(dplyr)
library(tidyr)
library(purrr)
library(extraDistr)
library(BayesianTools)
library(tidyverse)
library(broom)

n_posterior_samples = 120
##load in parameter estimates
someyrs_10k_out<-readRDS("/model_estimates.rds")
someyrs_10k_pars <- getSample(someyrs_10k_out, parametersOnly = TRUE)

posterior_pars <- someyrs_10k_pars[sample(nrow(someyrs_10k_pars), n_posterior_samples), ]

log_alpha0_B_L_2016<-posterior_pars[,"log_alpha0_L_BLAN_2016"]
log_alpha0_B_L_2017<-posterior_pars[,"log_alpha0_L_BLAN_2017"]
log_alpha0_B_L_2018<-posterior_pars[,"log_alpha0_L_BLAN_2018"]
log_alpha0_B_L_2019<-posterior_pars[,"log_alpha0_L_BLAN_2019"]

log_alpha0_B_N_2016<-posterior_pars[,"log_alpha0_N_BLAN_2016"]
log_alpha0_B_N_2017<-posterior_pars[,"log_alpha0_N_BLAN_2017"]
log_alpha0_B_N_2018<-posterior_pars[,"log_alpha0_N_BLAN_2018"]
log_alpha0_B_N_2019<-posterior_pars[,"log_alpha0_N_BLAN_2019"]


log_alpha0_H_L_2017<-posterior_pars[,"log_alpha0_L_HARV_2017"]
log_alpha0_H_L_2018<-posterior_pars[,"log_alpha0_L_HARV_2017"]
log_alpha0_H_L_2021<-posterior_pars[,"log_alpha0_L_HARV_2021"]
log_alpha0_H_L_2023<-posterior_pars[,"log_alpha0_L_HARV_2023"]

log_alpha0_H_N_2017<-posterior_pars[,"log_alpha0_N_HARV_2017"]
log_alpha0_H_N_2018<-posterior_pars[,"log_alpha0_N_HARV_2017"]
log_alpha0_H_N_2021<-posterior_pars[,"log_alpha0_N_HARV_2021"]
log_alpha0_H_N_2023<-posterior_pars[,"log_alpha0_N_HARV_2023"]


log_alpha0_S_L_2016<-posterior_pars[,"log_alpha0_L_SCBI_2016"]
log_alpha0_S_L_2017<-posterior_pars[,"log_alpha0_L_SCBI_2017"]
log_alpha0_S_L_2018<-posterior_pars[,"log_alpha0_L_SCBI_2018"]
log_alpha0_S_L_2019<-posterior_pars[,"log_alpha0_L_SCBI_2019"]
log_alpha0_S_L_2023<-posterior_pars[,"log_alpha0_L_SCBI_2023"]

log_alpha0_S_N_2016<-posterior_pars[,"log_alpha0_N_SCBI_2016"]
log_alpha0_S_N_2017<-posterior_pars[,"log_alpha0_N_SCBI_2017"]
log_alpha0_S_N_2018<-posterior_pars[,"log_alpha0_N_SCBI_2018"]
log_alpha0_S_N_2019<-posterior_pars[,"log_alpha0_N_SCBI_2019"]
log_alpha0_S_N_2023<-posterior_pars[,"log_alpha0_N_SCBI_2023"]


##load functions
# helper function
impute_with_overall_mean <- function(year_list, time_col = "time", value_col = "temperature",
                                     time_range = 1:366) {
  stopifnot(is.list(year_list), length(year_list) > 0)
  
  # 1) Build a full grid (all years x all times), then merge each year's data onto it
  years <- names(year_list)
  if (is.null(years) || any(years == "")) years <- seq_along(year_list)
  
  full <- do.call(rbind, lapply(seq_along(year_list), function(i) {
    df <- year_list[[i]]
    df <- df[, c(time_col, value_col)]
    names(df) <- c("time", "value")
    
    out <- merge(data.frame(time = time_range), df, by = "time", all.x = TRUE, sort = TRUE)
    out$year <- years[i]
    out
  }))
  
  # 2) Overall average at each time across years (ignoring NAs)
  overall_mean <- aggregate(value ~ time, data = full, FUN = function(x) mean(x, na.rm = TRUE))
  names(overall_mean)[2] <- "overall_mean"
  
  # 3) Impute missing values in each year with the overall mean for that time
  full <- merge(full, overall_mean, by = "time", all.x = TRUE, sort = TRUE)
  full$value_imputed <- ifelse(is.na(full$value), full$overall_mean, full$value)
  
  # 4) Return:
  #    - overall mean by time
  #    - imputed data split back into the original list structure
  imputed_list <- split(full[, c("time", "value_imputed")], full$year)
  imputed_list <- lapply(imputed_list, function(d) {
    names(d) <- c(time_col, value_col)
    d[order(d[[time_col]]), ]
  })
  
  list(
    overall_mean = overall_mean[order(overall_mean$time), ],
    imputed_years = imputed_list,
    long = full[order(full$year, full$time), c("year", "time", "value", "overall_mean", "value_imputed")]
  )
}


impute_from_3yravg <- function(year_df, avg_df) {
  tibble(time = 1:365) %>%
    left_join(year_df, by = "time") %>%
    left_join(avg_df, by = "time", suffix = c("_year", "_avg")) %>%
    mutate(RH = coalesce(RH_year, RH_avg)) %>%
    dplyr::select(time, RH)
}


##input year to run to initialize
total=365*10-1

##and to plot
yr=7
i=yr*365

# #load in data for fitting 

##LOAD TICK ABUNDANCE DATA
stackByTable("/NEON_count-ticks.zip")

tckabun <- readTableNEON(
  dataFile="/NEON_count-ticks/stackedFiles/tck_fielddata.csv",
  varFile="/NEON_count-ticks/stackedFiles/variables_10093.csv")

ticktaxon<-read.csv("/NEON_count-ticks/stackedFiles/tck_taxonomyProcessed.csv")

field <- tckabun
tax   <- ticktaxon

# 1) Build lab-derived life stage counts per sampleID
lab_counts <- tax %>%
  mutate(
    lifeStage = case_when(
      sexOrAge %in% c("Male", "Female", "Adult") ~ "Adult",
      sexOrAge %in% c("Nymph", "Larva")          ~ sexOrAge,
      TRUE                                       ~ NA_character_
    )
  ) %>%
  filter(!is.na(sampleID), !is.na(lifeStage)) %>%
  group_by(sampleID, lifeStage) %>%
  summarise(n = sum(individualCount, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = lifeStage, values_from = n, values_fill = 0) %>%
  rename(
    adultCount_lab = Adult,
    nymphCount_lab = Nymph,
    larvaCount_lab = Larva
  )

# 2) Join to fielddata and fill only where field counts are NA
field_filled <- field %>%
  left_join(lab_counts, by = "sampleID") %>%
  mutate(
    adultCount = coalesce(adultCount, adultCount_lab),
    nymphCount = coalesce(nymphCount, nymphCount_lab),
    larvaCount = coalesce(larvaCount, larvaCount_lab)
  ) %>%
  select(-adultCount_lab, -nymphCount_lab, -larvaCount_lab)




BLAN_a<-subset(field_filled,field_filled$siteID == "BLAN")
HARV_a<-subset(field_filled,field_filled$siteID == "HARV")
SCBI_a<-subset(field_filled,field_filled$siteID == "SCBI")

##BLAN
#Format data
BLAN_a$DOY<-as.numeric(strftime(BLAN_a$collectDate, format = "%j"))
BLAN_a$year<-as.numeric(format(as.POSIXct(BLAN_a$collectDate), "%Y"))
BLAN_a$DOY_year<-paste(BLAN_a$DOY,"_",BLAN_a$year,sep="")

#remove rows w NAs in larval, nymphal, and adult count columes
BLAN_a<-BLAN_a[!is.na(BLAN_a$larvaCount ),]
BLAN_a<-BLAN_a[!is.na(BLAN_a$nymphCount ),]
BLAN_a<-BLAN_a[!is.na(BLAN_a$adultCount ),]

#calculate the total area sampled on each day
BLAN_area<-aggregate(BLAN_a$totalSampledArea, by=list(Category=BLAN_a$DOY_year), FUN=sum)

#calculate number of each life stage collected on each day
BLAN_larvae<-aggregate(BLAN_a$larvaCount, by=list(Category=BLAN_a$DOY_year), FUN=sum)
BLAN_nymphs<-aggregate(BLAN_a$nymphCount, by=list(Category=BLAN_a$DOY_year), FUN=sum)
BLAN_adults<-aggregate(BLAN_a$adultCount, by=list(Category=BLAN_a$DOY_year), FUN=sum)

##calculate per capita based on area sampled
BLAN_larvae_pc<-BLAN_larvae[,2]/BLAN_area[,2]
BLAN_nymphs_pc<-BLAN_nymphs[,2]/BLAN_area[,2]
BLAN_adults_pc<-BLAN_adults[,2]/BLAN_area[,2]

##splot DOY_year into 2 strings
split_strings <- sapply(BLAN_area[,1], function(x) strsplit(x, "_")[[1]])

# Convert the result to a data frame for better readability
split_df <- data.frame(
  first_part = sapply(split_strings, `[`, 1)
)

BLAN_DOY<- as.numeric(split_df$first_part[seq(1, nrow(split_df), by = 2)])
BLAN_year <- as.numeric(split_df$first_part[seq(2, nrow(split_df), by = 2)])

#combine columns of interest into a dataframe
BLAN_a<-cbind(BLAN_DOY, BLAN_year, BLAN_area[,2],BLAN_larvae_pc,BLAN_nymphs_pc,BLAN_adults_pc, BLAN_larvae[,2], BLAN_nymphs[,2])
colnames(BLAN_a)<-c("DOY","year","area","larvae_perarea","nymphs_perarea","adults_perarea","larvae", "nymphs")
BLAN_a<-as.data.frame(BLAN_a)

#Cchange format into long form
BLAN_a_larv<-BLAN_a[,c(1,2,3,7,4)]
BLAN_a_larv$stage<-rep("L",nrow(BLAN_a_larv))
colnames(BLAN_a_larv)<-c("DOY","year","area","ticks","ticks_pc","stage")

BLAN_a_nymphs<-BLAN_a[,c(1,2,3,8,5)]
BLAN_a_nymphs$stage<-rep("N",nrow(BLAN_a_nymphs))
colnames(BLAN_a_nymphs)<-c("DOY","year","area","ticks","ticks_pc","stage")

BLAN_a_adult<-BLAN_a[,c(1,2,3,6)]
BLAN_a_adult$stage<-rep("A",nrow(BLAN_a_adult))
colnames(BLAN_a_adult)<-c("DOY","year","area","ticks_pc","stage")

#get rid of adult rows
BLAN_a<-as.data.frame(rbind(BLAN_a_larv,BLAN_a_nymphs))


##HARV
#Format data
HARV_a$DOY<-as.numeric(strftime(HARV_a$collectDate, format = "%j"))
HARV_a$year<-as.numeric(format(as.POSIXct(HARV_a$collectDate), "%Y"))
HARV_a$DOY_year<-paste(HARV_a$DOY,"_",HARV_a$year,sep="")

#remove rows w NAs in larval, nymphal, and adult count columes
HARV_a<-HARV_a[!is.na(HARV_a$larvaCount ),]
HARV_a<-HARV_a[!is.na(HARV_a$nymphCount ),]
HARV_a<-HARV_a[!is.na(HARV_a$adultCount ),]
HARV_area<-aggregate(HARV_a$totalSampledArea, by=list(Category=HARV_a$DOY_year), FUN=sum)
HARV_larvae<-aggregate(HARV_a$larvaCount, by=list(Category=HARV_a$DOY_year), FUN=sum)
HARV_nymphs<-aggregate(HARV_a$nymphCount, by=list(Category=HARV_a$DOY_year), FUN=sum)
HARV_adults<-aggregate(HARV_a$adultCount, by=list(Category=HARV_a$DOY_year), FUN=sum)
HARV_larvae_pc<-HARV_larvae[,2]/HARV_area[,2]
HARV_nymphs_pc<-HARV_nymphs[,2]/HARV_area[,2]
HARV_adults_pc<-HARV_adults[,2]/HARV_area[,2]

##splot DOY_year into 2 strings
split_strings <- sapply(HARV_area[,1], function(x) strsplit(x, "_")[[1]])

# Convert the result to a data frame for better readability
split_df <- data.frame(
  first_part = sapply(split_strings, `[`, 1)
)

HARV_DOY<- as.numeric(split_df$first_part[seq(1, nrow(split_df), by = 2)])
HARV_year <- as.numeric(split_df$first_part[seq(2, nrow(split_df), by = 2)])


HARV_a<-cbind(HARV_DOY, HARV_year, HARV_area[,2],HARV_larvae_pc,HARV_nymphs_pc,HARV_adults_pc, HARV_larvae[,2], HARV_nymphs[,2])
colnames(HARV_a)<-c("DOY","year","area","larvae_pc","nymphs_pc","adults_pc", "larvae", "nymphs")



HARV_a<-as.data.frame(HARV_a)
HARV_a_larv<-HARV_a[,c(1,2,3,7,4)]
HARV_a_larv$stage<-rep("L",nrow(HARV_a_larv))
colnames(HARV_a_larv)<-c("DOY","year","area", "ticks", "ticks_pc","stage")
HARV_a_nymphs<-HARV_a[,c(1,2,3,8,5)]
HARV_a_nymphs$stage<-rep("N",nrow(HARV_a_nymphs))
colnames(HARV_a_nymphs)<-c("DOY","year","area", "ticks", "ticks_pc","stage")
HARV_a_adult<-HARV_a[,c(1,2,3,6)]
HARV_a_adult$stage<-rep("A",nrow(HARV_a_adult))
colnames(HARV_a_adult)<-c("DOY","year","area", "ticks_pc","stage")
HARV_a<-as.data.frame(rbind(HARV_a_larv,HARV_a_nymphs))



##SCBI
#For each day calculate: total nymphs, total larvae, total adults, total area flagged

#Format data
SCBI_a$DOY<-as.numeric(strftime(SCBI_a$collectDate, format = "%j"))
SCBI_a$year<-as.numeric(format(as.POSIXct(SCBI_a$collectDate), "%Y"))
SCBI_a$DOY_year<-paste(SCBI_a$DOY,"_",SCBI_a$year,sep="")

#remove rows w NAs in larval, nymphal, and adult count columes
SCBI_a<-SCBI_a[!is.na(SCBI_a$larvaCount ),]
SCBI_a<-SCBI_a[!is.na(SCBI_a$nymphCount ),]
SCBI_a<-SCBI_a[!is.na(SCBI_a$adultCount ),]
SCBI_area<-aggregate(SCBI_a$totalSampledArea, by=list(Category=SCBI_a$DOY_year), FUN=sum)
SCBI_larvae<-aggregate(SCBI_a$larvaCount, by=list(Category=SCBI_a$DOY_year), FUN=sum)
SCBI_nymphs<-aggregate(SCBI_a$nymphCount, by=list(Category=SCBI_a$DOY_year), FUN=sum)
SCBI_adults<-aggregate(SCBI_a$adultCount, by=list(Category=SCBI_a$DOY_year), FUN=sum)
SCBI_larvae_pc<-SCBI_larvae[,2]/SCBI_area[,2]
SCBI_nymphs_pc<-SCBI_nymphs[,2]/SCBI_area[,2]
SCBI_adults_pc<-SCBI_adults[,2]/SCBI_area[,2]

##splot DOY_year into 2 strings
split_strings <- sapply(SCBI_area[,1], function(x) strsplit(x, "_")[[1]])

# Convert the result to a data frame for better readability
split_df <- data.frame(
  first_part = sapply(split_strings, `[`, 1)
)

SCBI_DOY<- as.numeric(split_df$first_part[seq(1, nrow(split_df), by = 2)])
SCBI_year <- as.numeric(split_df$first_part[seq(2, nrow(split_df), by = 2)])
SCBI_a<-cbind(SCBI_DOY, SCBI_year, SCBI_area[,2],SCBI_larvae_pc,SCBI_nymphs_pc,SCBI_adults_pc, SCBI_larvae[,2], SCBI_nymphs[,2])
colnames(SCBI_a)<-c("DOY","year","area","larvae_pc","nymphs_pc","adults_pc", "larvae", "nymphs")
SCBI_a<-as.data.frame(SCBI_a)


SCBI_a_larv<-SCBI_a[,c(1,2,3,7,4)]
SCBI_a_larv$stage<-rep("L",nrow(SCBI_a_larv))
colnames(SCBI_a_larv)<-c("DOY","year","area", "ticks", "ticks_pc","stage")
SCBI_a_nymphs<-SCBI_a[,c(1,2,3,8,5)]
SCBI_a_nymphs$stage<-rep("N",nrow(SCBI_a_nymphs))
colnames(SCBI_a_nymphs)<-c("DOY","year","area", "ticks", "ticks_pc","stage")
SCBI_a_adult<-SCBI_a[,c(1,2,3,6)]
SCBI_a_adult$stage<-rep("A",nrow(SCBI_a_adult))
colnames(SCBI_a_adult)<-c("DOY","year","area", "ticks_pc","stage")
SCBI_a<-as.data.frame(rbind(SCBI_a_larv,SCBI_a_nymphs))


##Subset for each year

BLAN_a_2016<-subset(BLAN_a,BLAN_a$year == 2016)
BLAN_a_2017<-subset(BLAN_a,BLAN_a$year == 2017)
BLAN_a_2018<-subset(BLAN_a,BLAN_a$year == 2018)
BLAN_a_2019<-subset(BLAN_a,BLAN_a$year == 2019)
BLAN_a_2020<-subset(BLAN_a,BLAN_a$year == 2020)
BLAN_a_2021<-subset(BLAN_a,BLAN_a$year == 2021)
BLAN_a_2022<-subset(BLAN_a,BLAN_a$year == 2022)
BLAN_a_2023<-subset(BLAN_a,BLAN_a$year == 2023)
BLAN_a_2024<-subset(BLAN_a,BLAN_a$year == 2024)


HARV_a_2016<-subset(HARV_a,HARV_a$year == 2016)
HARV_a_2017<-subset(HARV_a,HARV_a$year == 2017)
HARV_a_2018<-subset(HARV_a,HARV_a$year == 2018)
HARV_a_2019<-subset(HARV_a,HARV_a$year == 2019)
HARV_a_2020<-subset(HARV_a,HARV_a$year == 2020)
HARV_a_2021<-subset(HARV_a,HARV_a$year == 2021)
HARV_a_2022<-subset(HARV_a,HARV_a$year == 2022)
HARV_a_2023<-subset(HARV_a,HARV_a$year == 2023)
HARV_a_2024<-subset(HARV_a,HARV_a$year == 2024)

SCBI_a_2016<-subset(SCBI_a,SCBI_a$year == 2016)
SCBI_a_2017<-subset(SCBI_a,SCBI_a$year == 2017)
SCBI_a_2018<-subset(SCBI_a,SCBI_a$year == 2018)
SCBI_a_2019<-subset(SCBI_a,SCBI_a$year == 2019)
SCBI_a_2020<-subset(SCBI_a,SCBI_a$year == 2020)
SCBI_a_2021<-subset(SCBI_a,SCBI_a$year == 2021)
SCBI_a_2022<-subset(SCBI_a,SCBI_a$year == 2022)
SCBI_a_2023<-subset(SCBI_a,SCBI_a$year == 2023)
SCBI_a_2024<-subset(SCBI_a,SCBI_a$year == 2024)



##format data for fitting
BLAN_a_2016$time<- ((365*7) + BLAN_a_2016$DOY - 1)

BLAN_a_2017$time<- ((365*7) + BLAN_a_2017$DOY - 1)

BLAN_a_2018$time<- ((365*7) + BLAN_a_2018$DOY - 1)

BLAN_a_2019$time<- ((365*7) + BLAN_a_2019$DOY - 1)

BLAN_a_2020$time<- ((365*7) + BLAN_a_2020$DOY - 1)

BLAN_a_2021$time<- ((365*7) + BLAN_a_2021$DOY - 1)

BLAN_a_2022$time<- ((365*7) + BLAN_a_2022$DOY - 1)

BLAN_a_2023$time<- ((365*7) + BLAN_a_2023$DOY - 1)

BLAN_a_2024$time<- ((365*7) + BLAN_a_2024$DOY - 1)

HARV_a_2016$time<- ((365*7) + HARV_a_2016$DOY - 1)

HARV_a_2017$time<- ((365*7) + HARV_a_2017$DOY - 1)

HARV_a_2018$time<- ((365*7) + HARV_a_2018$DOY - 1)

HARV_a_2019$time<- ((365*7) + HARV_a_2019$DOY - 1)

HARV_a_2020$time<- ((365*7) + HARV_a_2020$DOY - 1)

HARV_a_2021$time<- ((365*7) + HARV_a_2021$DOY - 1)

HARV_a_2022$time<- ((365*7) + HARV_a_2022$DOY - 1)

HARV_a_2023$time<- ((365*7) + HARV_a_2023$DOY - 1)

HARV_a_2024$time<- ((365*7) + HARV_a_2024$DOY - 1)









###SCBI
SCBI_a_2016$time<- ((365*7) + SCBI_a_2016$DOY - 1)

SCBI_a_2017$time<- ((365*7) + SCBI_a_2017$DOY - 1)

SCBI_a_2018$time<- ((365*7) + SCBI_a_2018$DOY - 1)

SCBI_a_2019$time<- ((365*7) + SCBI_a_2019$DOY - 1)

SCBI_a_2020$time<- ((365*7) + SCBI_a_2020$DOY - 1)

SCBI_a_2021$time<- ((365*7) + SCBI_a_2021$DOY - 1)

SCBI_a_2022$time<- ((365*7) + SCBI_a_2022$DOY - 1)

SCBI_a_2023$time<- ((365*7) + SCBI_a_2023$DOY - 1)

SCBI_a_2024$time<- ((365*7) + SCBI_a_2024$DOY - 1)


##load in functions to create temperature splines
# Wrap the spline function to handle periodicity
periodic_spline_BLAN_2016 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_BLAN_2016(t_mod))
}

periodic_spline_BLAN_2017 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_BLAN_2017(t_mod))
}

periodic_spline_BLAN_2018 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_BLAN_2018(t_mod))
}

periodic_spline_BLAN_2019 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_BLAN_2019(t_mod))
}

periodic_spline_BLAN_2020 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_BLAN_2020(t_mod))
}

periodic_spline_BLAN_2021 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_BLAN_2021(t_mod))
}

periodic_spline_BLAN_2022 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_BLAN_2022(t_mod))
}

periodic_spline_BLAN_2023 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_BLAN_2023(t_mod))
}

periodic_spline_BLAN_2024 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_BLAN_2024(t_mod))
}

periodic_spline_HARV_2016 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2016(t_mod))
}

periodic_spline_HARV_2017 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2017(t_mod))
}

periodic_spline_HARV_2018 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2018(t_mod))
}

periodic_spline_HARV_2018 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2018(t_mod))
}

periodic_spline_HARV_2019 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2019(t_mod))
}

periodic_spline_HARV_2020 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2020(t_mod))
}

periodic_spline_HARV_2021 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2021(t_mod))
}

periodic_spline_HARV_2022 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2022(t_mod))
}

periodic_spline_HARV_2023 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2023(t_mod))
}

periodic_spline_HARV_2024 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_HARV_2024(t_mod))
}

periodic_spline_SCBI_2016 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_SCBI_2016(t_mod))
}

periodic_spline_SCBI_2017 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_SCBI_2017(t_mod))
}

periodic_spline_SCBI_2018 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_SCBI_2018(t_mod))
}

periodic_spline_SCBI_2019 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_SCBI_2019(t_mod))
}

periodic_spline_SCBI_2020 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_SCBI_2020(t_mod))
}

periodic_spline_SCBI_2021 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_SCBI_2021(t_mod))
}

periodic_spline_SCBI_2022 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_SCBI_2022(t_mod))
}

periodic_spline_SCBI_2023 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_SCBI_2023(t_mod))
}

periodic_spline_SCBI_2024 <- function(t) {
  t_mod <- ((t - 1) %% 365) + 1  # Map t to the range 1 to 365
  return(spline_fit_SCBI_2024(t_mod))
}

#Temperature simulation
f_temp_BLAN_2016<-function(t)
{
  periodic_spline_BLAN_2016(t)
}

f_temp_BLAN_2017<-function(t)
{
  periodic_spline_BLAN_2017(t)
}

f_temp_BLAN_2018<-function(t)
{
  periodic_spline_BLAN_2018(t)
}

f_temp_BLAN_2019<-function(t)
{
  periodic_spline_BLAN_2019(t)
}

f_temp_BLAN_2020<-function(t)
{
  periodic_spline_BLAN_2020(t)
}

f_temp_BLAN_2021<-function(t)
{
  periodic_spline_BLAN_2021(t)
}

f_temp_BLAN_2022<-function(t)
{
  periodic_spline_BLAN_2022(t)
}

f_temp_BLAN_2023<-function(t)
{
  periodic_spline_BLAN_2023(t)
}

f_temp_BLAN_2024<-function(t)
{
  periodic_spline_BLAN_2024(t)
}

f_temp_HARV_2016<-function(t)
{
  periodic_spline_HARV_2016(t)
}

f_temp_HARV_2017<-function(t)
{
  periodic_spline_HARV_2017(t)
}

f_temp_HARV_2018<-function(t)
{
  periodic_spline_HARV_2018(t)
}

f_temp_HARV_2019<-function(t)
{
  periodic_spline_HARV_2019(t)
}

f_temp_HARV_2020<-function(t)
{
  periodic_spline_HARV_2020(t)
}

f_temp_HARV_2021<-function(t)
{
  periodic_spline_HARV_2021(t)
}

f_temp_HARV_2022<-function(t)
{
  periodic_spline_HARV_2022(t)
}

f_temp_HARV_2023<-function(t)
{
  periodic_spline_HARV_2023(t)
}

f_temp_HARV_2024<-function(t)
{
  periodic_spline_HARV_2024(t)
}



f_temp_SCBI_2016<-function(t)
{
  periodic_spline_SCBI_2016(t)
}

f_temp_SCBI_2017<-function(t)
{
  periodic_spline_SCBI_2017(t)
}

f_temp_SCBI_2018<-function(t)
{
  periodic_spline_SCBI_2018(t)
}

f_temp_SCBI_2019<-function(t)
{
  periodic_spline_SCBI_2019(t)
}

f_temp_SCBI_2020<-function(t)
{
  periodic_spline_SCBI_2020(t)
}

f_temp_SCBI_2021<-function(t)
{
  periodic_spline_SCBI_2021(t)
}

f_temp_SCBI_2022<-function(t)
{
  periodic_spline_SCBI_2022(t)
}

f_temp_SCBI_2023<-function(t)
{
  periodic_spline_SCBI_2023(t)
}

f_temp_SCBI_2024<-function(t)
{
  periodic_spline_SCBI_2024(t)
}

##BLAN

##load in temp data and make it periodic 
BLAN<-read.csv("/BLAN.csv")

#Format data
BLAN$DOY<-as.numeric(strftime(BLAN$startDateTime, format = "%j"))
BLAN<-na.omit(BLAN)
BLAN$year <- year(ymd_hms(BLAN$startDateTime))

##subset by year
BLAN_2016<-subset(BLAN, BLAN$year == 2016)
BLAN_2017<-subset(BLAN, BLAN$year == 2017)
BLAN_2018<-subset(BLAN, BLAN$year == 2018)
BLAN_2019<-subset(BLAN, BLAN$year == 2019)
BLAN_2020<-subset(BLAN, BLAN$year == 2020)
BLAN_2021<-subset(BLAN, BLAN$year == 2021)
BLAN_2022<-subset(BLAN, BLAN$year == 2022)
BLAN_2023<-subset(BLAN, BLAN$year == 2023)
BLAN_2024<-subset(BLAN, BLAN$year == 2024)




BLAN_2016<-subset(BLAN_2016,BLAN_2016$finalQF == 0)
BLAN_2017<-subset(BLAN_2017,BLAN_2017$finalQF == 0)
BLAN_2018<-subset(BLAN_2018,BLAN_2018$finalQF == 0)
BLAN_2019<-subset(BLAN_2019,BLAN_2019$finalQF == 0)
BLAN_2020<-subset(BLAN_2020,BLAN_2020$finalQF == 0)
BLAN_2021<-subset(BLAN_2021,BLAN_2021$finalQF == 0)
BLAN_2022<-subset(BLAN_2022,BLAN_2022$finalQF == 0)
BLAN_2023<-subset(BLAN_2023,BLAN_2023$finalQF == 0)


BLAN_mean_2016<-aggregate(BLAN_2016$tempSingleMaximum, by=list(Category=BLAN_2016$DOY), FUN=mean)
BLAN_mean_2017<-aggregate(BLAN_2017$tempSingleMaximum, by=list(Category=BLAN_2017$DOY), FUN=mean)
BLAN_mean_2018<-aggregate(BLAN_2018$tempSingleMaximum, by=list(Category=BLAN_2018$DOY), FUN=mean)
BLAN_mean_2019<-aggregate(BLAN_2019$tempSingleMaximum, by=list(Category=BLAN_2019$DOY), FUN=mean)
BLAN_mean_2020<-aggregate(BLAN_2020$tempSingleMaximum, by=list(Category=BLAN_2020$DOY), FUN=mean)
BLAN_mean_2021<-aggregate(BLAN_2021$tempSingleMaximum, by=list(Category=BLAN_2021$DOY), FUN=mean)
BLAN_mean_2022<-aggregate(BLAN_2022$tempSingleMaximum, by=list(Category=BLAN_2022$DOY), FUN=mean)
BLAN_mean_2023<-aggregate(BLAN_2023$tempSingleMaximum, by=list(Category=BLAN_2023$DOY), FUN=mean)



colnames(BLAN_mean_2016)<-c("time","temperature")
colnames(BLAN_mean_2017)<-c("time","temperature")
colnames(BLAN_mean_2018)<-c("time","temperature")
colnames(BLAN_mean_2019)<-c("time","temperature")
colnames(BLAN_mean_2020)<-c("time","temperature")
colnames(BLAN_mean_2021)<-c("time","temperature")
colnames(BLAN_mean_2022)<-c("time","temperature")
colnames(BLAN_mean_2023)<-c("time","temperature")

##create df to compute overall daily average for imputations

temps_BLAN<-list(
  '2016' = BLAN_mean_2016,
  '2017' = BLAN_mean_2017,
  '2018' = BLAN_mean_2018,
  '2019' = BLAN_mean_2019,
  '2020' = BLAN_mean_2020,
  '2021' = BLAN_mean_2021,
  '2022' = BLAN_mean_2022,
  '2023' = BLAN_mean_2023
)

BLAN_res <- impute_with_overall_mean(temps_BLAN)  # temps is your list of yearly dfs
BLAN_overall <- BLAN_res$overall_mean             # average temp for each time (1:366)
BLAN_temps_imputed <- BLAN_res$imputed_years      # list of yearly dfs with missing times filled

BLAN_mean_2016<-BLAN_temps_imputed[["2016"]]
BLAN_mean_2017<-BLAN_temps_imputed[["2017"]]
BLAN_mean_2018<-BLAN_temps_imputed[["2018"]]
BLAN_mean_2019<-BLAN_temps_imputed[["2019"]]
BLAN_mean_2020<-BLAN_temps_imputed[["2020"]]
BLAN_mean_2021<-BLAN_temps_imputed[["2021"]]
BLAN_mean_2022<-BLAN_temps_imputed[["2022"]]
BLAN_mean_2023<-BLAN_temps_imputed[["2023"]]

BLAN_data<- list(BLAN_mean_2016, BLAN_mean_2017, BLAN_mean_2018,
                 BLAN_mean_2019, BLAN_mean_2020, BLAN_mean_2021,
                 BLAN_mean_2022, BLAN_mean_2023)


##loop through years
yr_i = 1

days_BLAN_2016 <- BLAN_data[[yr_i]]$time - 1
temp_BLAN_2016 <- BLAN_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_BLAN_2016 <- c(days_BLAN_2016, (max(days_BLAN_2016) + (days_BLAN_2016[2] - days_BLAN_2016[1])))
temp_extended_BLAN_2016 <- c(temp_BLAN_2016, temp_BLAN_2016[1])

# Fit the spline
spline_fit_BLAN_2016 <- splinefun(days_extended_BLAN_2016, temp_extended_BLAN_2016)

##create the amount of temp data that I need for odin interpolation
DOY_BLAN_2016<-c(0:total)
temperature_BLAN_2016<-f_temp_BLAN_2016(DOY_BLAN_2016)

##loop through years
yr_i = 2

days_BLAN_2017 <- BLAN_data[[yr_i]]$time - 1
temp_BLAN_2017 <- BLAN_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_BLAN_2017 <- c(days_BLAN_2017, (max(days_BLAN_2017) + (days_BLAN_2017[2] - days_BLAN_2017[1])))
temp_extended_BLAN_2017 <- c(temp_BLAN_2017, temp_BLAN_2017[1])

# Fit the spline
spline_fit_BLAN_2017 <- splinefun(days_extended_BLAN_2017, temp_extended_BLAN_2017)

##create the amount of temp data that I need for odin interpolation
DOY_BLAN_2017<-c(0:total)
temperature_BLAN_2017<-f_temp_BLAN_2017(DOY_BLAN_2017)

##loop through years
yr_i = 3

days_BLAN_2018 <- BLAN_data[[yr_i]]$time - 1
temp_BLAN_2018 <- BLAN_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_BLAN_2018 <- c(days_BLAN_2018, (max(days_BLAN_2018) + (days_BLAN_2018[2] - days_BLAN_2018[1])))
temp_extended_BLAN_2018 <- c(temp_BLAN_2018, temp_BLAN_2018[1])

# Fit the spline
spline_fit_BLAN_2018 <- splinefun(days_extended_BLAN_2018, temp_extended_BLAN_2018)

##create the amount of temp data that I need for odin interpolation
DOY_BLAN_2018<-c(0:total)
temperature_BLAN_2018<-f_temp_BLAN_2018(DOY_BLAN_2018)


yr_i = 4

days_BLAN_2019 <- BLAN_data[[yr_i]]$time - 1
temp_BLAN_2019 <- BLAN_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_BLAN_2019 <- c(days_BLAN_2019, (max(days_BLAN_2019) + (days_BLAN_2019[2] - days_BLAN_2019[1])))
temp_extended_BLAN_2019 <- c(temp_BLAN_2019, temp_BLAN_2019[1])

# Fit the spline
spline_fit_BLAN_2019 <- splinefun(days_extended_BLAN_2019, temp_extended_BLAN_2019)

##create the amount of temp data that I need for odin interpolation
DOY_BLAN_2019<-c(0:total)
temperature_BLAN_2019<-f_temp_BLAN_2019(DOY_BLAN_2019)

yr_i = 5

days_BLAN_2020 <- BLAN_data[[yr_i]]$time - 1
temp_BLAN_2020 <- BLAN_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_BLAN_2020 <- c(days_BLAN_2020, (max(days_BLAN_2020) + (days_BLAN_2020[2] - days_BLAN_2020[1])))
temp_extended_BLAN_2020 <- c(temp_BLAN_2020, temp_BLAN_2020[1])

# Fit the spline
spline_fit_BLAN_2020 <- splinefun(days_extended_BLAN_2020, temp_extended_BLAN_2020)

##create the amount of temp data that I need for odin interpolation
DOY_BLAN_2020<-c(0:total)
temperature_BLAN_2020<-f_temp_BLAN_2020(DOY_BLAN_2020)

yr_i = 6

days_BLAN_2021 <- BLAN_data[[yr_i]]$time - 1
temp_BLAN_2021 <- BLAN_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_BLAN_2021 <- c(days_BLAN_2021, (max(days_BLAN_2021) + (days_BLAN_2021[2] - days_BLAN_2021[1])))
temp_extended_BLAN_2021 <- c(temp_BLAN_2021, temp_BLAN_2021[1])

# Fit the spline
spline_fit_BLAN_2021 <- splinefun(days_extended_BLAN_2021, temp_extended_BLAN_2021)

##create the amount of temp data that I need for odin interpolation
DOY_BLAN_2021<-c(0:total)
temperature_BLAN_2021<-f_temp_BLAN_2021(DOY_BLAN_2021)

yr_i = 7

days_BLAN_2022 <- BLAN_data[[yr_i]]$time - 1
temp_BLAN_2022 <- BLAN_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_BLAN_2022 <- c(days_BLAN_2022, (max(days_BLAN_2022) + (days_BLAN_2022[2] - days_BLAN_2022[1])))
temp_extended_BLAN_2022 <- c(temp_BLAN_2022, temp_BLAN_2022[1])

# Fit the spline
spline_fit_BLAN_2022 <- splinefun(days_extended_BLAN_2022, temp_extended_BLAN_2022)

##create the amount of temp data that I need for odin interpolation
DOY_BLAN_2022<-c(0:total)
temperature_BLAN_2022<-f_temp_BLAN_2022(DOY_BLAN_2022)

yr_i = 8

days_BLAN_2023 <- BLAN_data[[yr_i]]$time - 1
temp_BLAN_2023 <- BLAN_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_BLAN_2023 <- c(days_BLAN_2023, (max(days_BLAN_2023) + (days_BLAN_2023[2] - days_BLAN_2023[1])))
temp_extended_BLAN_2023 <- c(temp_BLAN_2023, temp_BLAN_2023[1])

# Fit the spline
spline_fit_BLAN_2023 <- splinefun(days_extended_BLAN_2023, temp_extended_BLAN_2023)

##create the amount of temp data that I need for odin interpolation
DOY_BLAN_2023<-c(0:total)
temperature_BLAN_2023<-f_temp_BLAN_2023(DOY_BLAN_2023)


##HARV
##load in temp data and make it periodic 
HARV<-read.csv("/HARV.csv")

#Format data
HARV$DOY<-as.numeric(strftime(HARV$startDateTime, format = "%j"))
HARV<-na.omit(HARV)
HARV$year <- year(ymd_hms(HARV$startDateTime))

##subset by year
HARV_2016<-subset(HARV, HARV$year == 2016)
HARV_2017<-subset(HARV, HARV$year == 2017)
HARV_2018<-subset(HARV, HARV$year == 2018)
HARV_2019<-subset(HARV, HARV$year == 2019)
HARV_2020<-subset(HARV, HARV$year == 2020)
HARV_2021<-subset(HARV, HARV$year == 2021)
HARV_2022<-subset(HARV, HARV$year == 2022)
HARV_2023<-subset(HARV, HARV$year == 2023)
HARV_2024<-subset(HARV, HARV$year == 2024)




HARV_2016<-subset(HARV_2016,HARV_2016$finalQF == 0)
HARV_2017<-subset(HARV_2017,HARV_2017$finalQF == 0)
HARV_2018<-subset(HARV_2018,HARV_2018$finalQF == 0)
HARV_2019<-subset(HARV_2019,HARV_2019$finalQF == 0)
HARV_2020<-subset(HARV_2020,HARV_2020$finalQF == 0)
HARV_2021<-subset(HARV_2021,HARV_2021$finalQF == 0)
HARV_2022<-subset(HARV_2022,HARV_2022$finalQF == 0)
HARV_2023<-subset(HARV_2023,HARV_2023$finalQF == 0)


HARV_mean_2016<-aggregate(HARV_2016$tempSingleMaximum, by=list(Category=HARV_2016$DOY), FUN=mean)
HARV_mean_2017<-aggregate(HARV_2017$tempSingleMaximum, by=list(Category=HARV_2017$DOY), FUN=mean)
HARV_mean_2018<-aggregate(HARV_2018$tempSingleMaximum, by=list(Category=HARV_2018$DOY), FUN=mean)
HARV_mean_2019<-aggregate(HARV_2019$tempSingleMaximum, by=list(Category=HARV_2019$DOY), FUN=mean)
HARV_mean_2020<-aggregate(HARV_2020$tempSingleMaximum, by=list(Category=HARV_2020$DOY), FUN=mean)
HARV_mean_2021<-aggregate(HARV_2021$tempSingleMaximum, by=list(Category=HARV_2021$DOY), FUN=mean)
HARV_mean_2022<-aggregate(HARV_2022$tempSingleMaximum, by=list(Category=HARV_2022$DOY), FUN=mean)
HARV_mean_2023<-aggregate(HARV_2023$tempSingleMaximum, by=list(Category=HARV_2023$DOY), FUN=mean)



colnames(HARV_mean_2016)<-c("time","temperature")
colnames(HARV_mean_2017)<-c("time","temperature")
colnames(HARV_mean_2018)<-c("time","temperature")
colnames(HARV_mean_2019)<-c("time","temperature")
colnames(HARV_mean_2020)<-c("time","temperature")
colnames(HARV_mean_2021)<-c("time","temperature")
colnames(HARV_mean_2022)<-c("time","temperature")
colnames(HARV_mean_2023)<-c("time","temperature")

##create df to compute overall daily average for imputations

temps_HARV<-list(
  '2016' = HARV_mean_2016,
  '2017' = HARV_mean_2017,
  '2018' = HARV_mean_2018,
  '2019' = HARV_mean_2019,
  '2020' = HARV_mean_2020,
  '2021' = HARV_mean_2021,
  '2022' = HARV_mean_2022,
  '2023' = HARV_mean_2023
)

HARV_res <- impute_with_overall_mean(temps_HARV)  # temps is your list of yearly dfs
HARV_overall <- HARV_res$overall_mean             # average temp for each time (1:366)
HARV_temps_imputed <- HARV_res$imputed_years      # list of yearly dfs with missing times filled

HARV_mean_2016<-HARV_temps_imputed[["2016"]]
HARV_mean_2017<-HARV_temps_imputed[["2017"]]
HARV_mean_2018<-HARV_temps_imputed[["2018"]]
HARV_mean_2019<-HARV_temps_imputed[["2019"]]
HARV_mean_2020<-HARV_temps_imputed[["2020"]]
HARV_mean_2021<-HARV_temps_imputed[["2021"]]
HARV_mean_2022<-HARV_temps_imputed[["2022"]]
HARV_mean_2023<-HARV_temps_imputed[["2023"]]

HARV_data<- list(HARV_mean_2016, HARV_mean_2017, HARV_mean_2018,
                 HARV_mean_2019, HARV_mean_2020, HARV_mean_2021,
                 HARV_mean_2022, HARV_mean_2023)


##loop through years
yr_i = 1

days_HARV_2016 <- HARV_data[[yr_i]]$time - 1
temp_HARV_2016 <- HARV_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_HARV_2016 <- c(days_HARV_2016, (max(days_HARV_2016) + (days_HARV_2016[2] - days_HARV_2016[1])))
temp_extended_HARV_2016 <- c(temp_HARV_2016, temp_HARV_2016[1])

# Fit the spline
spline_fit_HARV_2016 <- splinefun(days_extended_HARV_2016, temp_extended_HARV_2016)

##create the amount of temp data that I need for odin interpolation
DOY_HARV_2016<-c(0:total)
temperature_HARV_2016<-f_temp_HARV_2016(DOY_HARV_2016)

##loop through years
yr_i = 2

days_HARV_2017 <- HARV_data[[yr_i]]$time - 1
temp_HARV_2017 <- HARV_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_HARV_2017 <- c(days_HARV_2017, (max(days_HARV_2017) + (days_HARV_2017[2] - days_HARV_2017[1])))
temp_extended_HARV_2017 <- c(temp_HARV_2017, temp_HARV_2017[1])

# Fit the spline
spline_fit_HARV_2017 <- splinefun(days_extended_HARV_2017, temp_extended_HARV_2017)

##create the amount of temp data that I need for odin interpolation
DOY_HARV_2017<-c(0:total)
temperature_HARV_2017<-f_temp_HARV_2017(DOY_HARV_2017)

##loop through years
yr_i = 3

days_HARV_2018 <- HARV_data[[yr_i]]$time - 1
temp_HARV_2018 <- HARV_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_HARV_2018 <- c(days_HARV_2018, (max(days_HARV_2018) + (days_HARV_2018[2] - days_HARV_2018[1])))
temp_extended_HARV_2018 <- c(temp_HARV_2018, temp_HARV_2018[1])

# Fit the spline
spline_fit_HARV_2018 <- splinefun(days_extended_HARV_2018, temp_extended_HARV_2018)

##create the amount of temp data that I need for odin interpolation
DOY_HARV_2018<-c(0:total)
temperature_HARV_2018<-f_temp_HARV_2018(DOY_HARV_2018)


yr_i = 4

days_HARV_2019 <- HARV_data[[yr_i]]$time - 1
temp_HARV_2019 <- HARV_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_HARV_2019 <- c(days_HARV_2019, (max(days_HARV_2019) + (days_HARV_2019[2] - days_HARV_2019[1])))
temp_extended_HARV_2019 <- c(temp_HARV_2019, temp_HARV_2019[1])

# Fit the spline
spline_fit_HARV_2019 <- splinefun(days_extended_HARV_2019, temp_extended_HARV_2019)

##create the amount of temp data that I need for odin interpolation
DOY_HARV_2019<-c(0:total)
temperature_HARV_2019<-f_temp_HARV_2019(DOY_HARV_2019)

yr_i = 5

days_HARV_2020 <- HARV_data[[yr_i]]$time - 1
temp_HARV_2020 <- HARV_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_HARV_2020 <- c(days_HARV_2020, (max(days_HARV_2020) + (days_HARV_2020[2] - days_HARV_2020[1])))
temp_extended_HARV_2020 <- c(temp_HARV_2020, temp_HARV_2020[1])

# Fit the spline
spline_fit_HARV_2020 <- splinefun(days_extended_HARV_2020, temp_extended_HARV_2020)

##create the amount of temp data that I need for odin interpolation
DOY_HARV_2020<-c(0:total)
temperature_HARV_2020<-f_temp_HARV_2020(DOY_HARV_2020)

yr_i = 6

days_HARV_2021 <- HARV_data[[yr_i]]$time - 1
temp_HARV_2021 <- HARV_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_HARV_2021 <- c(days_HARV_2021, (max(days_HARV_2021) + (days_HARV_2021[2] - days_HARV_2021[1])))
temp_extended_HARV_2021 <- c(temp_HARV_2021, temp_HARV_2021[1])

# Fit the spline
spline_fit_HARV_2021 <- splinefun(days_extended_HARV_2021, temp_extended_HARV_2021)

##create the amount of temp data that I need for odin interpolation
DOY_HARV_2021<-c(0:total)
temperature_HARV_2021<-f_temp_HARV_2021(DOY_HARV_2021)

yr_i = 7

days_HARV_2022 <- HARV_data[[yr_i]]$time - 1
temp_HARV_2022 <- HARV_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_HARV_2022 <- c(days_HARV_2022, (max(days_HARV_2022) + (days_HARV_2022[2] - days_HARV_2022[1])))
temp_extended_HARV_2022 <- c(temp_HARV_2022, temp_HARV_2022[1])

# Fit the spline
spline_fit_HARV_2022 <- splinefun(days_extended_HARV_2022, temp_extended_HARV_2022)

##create the amount of temp data that I need for odin interpolation
DOY_HARV_2022<-c(0:total)
temperature_HARV_2022<-f_temp_HARV_2022(DOY_HARV_2022)

yr_i = 8

days_HARV_2023 <- HARV_data[[yr_i]]$time - 1
temp_HARV_2023 <- HARV_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_HARV_2023 <- c(days_HARV_2023, (max(days_HARV_2023) + (days_HARV_2023[2] - days_HARV_2023[1])))
temp_extended_HARV_2023 <- c(temp_HARV_2023, temp_HARV_2023[1])

# Fit the spline
spline_fit_HARV_2023 <- splinefun(days_extended_HARV_2023, temp_extended_HARV_2023)

##create the amount of temp data that I need for odin interpolation
DOY_HARV_2023<-c(0:total)
temperature_HARV_2023<-f_temp_HARV_2023(DOY_HARV_2023)


##SCBI
##load in temp data and make it periodic 
SCBI<-read.csv("/SCBI.csv")

#Format data
SCBI$DOY<-as.numeric(strftime(SCBI$startDateTime, format = "%j"))
SCBI<-na.omit(SCBI)
SCBI$year <- year(ymd_hms(SCBI$startDateTime))

##subset by year
SCBI_2016<-subset(SCBI, SCBI$year == 2016)
SCBI_2017<-subset(SCBI, SCBI$year == 2017)
SCBI_2018<-subset(SCBI, SCBI$year == 2018)
SCBI_2019<-subset(SCBI, SCBI$year == 2019)
SCBI_2020<-subset(SCBI, SCBI$year == 2020)
SCBI_2021<-subset(SCBI, SCBI$year == 2021)
SCBI_2022<-subset(SCBI, SCBI$year == 2022)
SCBI_2023<-subset(SCBI, SCBI$year == 2023)
SCBI_2024<-subset(SCBI, SCBI$year == 2024)




SCBI_2016<-subset(SCBI_2016,SCBI_2016$finalQF == 0)
SCBI_2017<-subset(SCBI_2017,SCBI_2017$finalQF == 0)
SCBI_2018<-subset(SCBI_2018,SCBI_2018$finalQF == 0)
SCBI_2019<-subset(SCBI_2019,SCBI_2019$finalQF == 0)
SCBI_2020<-subset(SCBI_2020,SCBI_2020$finalQF == 0)
SCBI_2021<-subset(SCBI_2021,SCBI_2021$finalQF == 0)
SCBI_2022<-subset(SCBI_2022,SCBI_2022$finalQF == 0)
SCBI_2023<-subset(SCBI_2023,SCBI_2023$finalQF == 0)


SCBI_mean_2016<-aggregate(SCBI_2016$tempSingleMaximum, by=list(Category=SCBI_2016$DOY), FUN=mean)
SCBI_mean_2017<-aggregate(SCBI_2017$tempSingleMaximum, by=list(Category=SCBI_2017$DOY), FUN=mean)
SCBI_mean_2018<-aggregate(SCBI_2018$tempSingleMaximum, by=list(Category=SCBI_2018$DOY), FUN=mean)
SCBI_mean_2019<-aggregate(SCBI_2019$tempSingleMaximum, by=list(Category=SCBI_2019$DOY), FUN=mean)
SCBI_mean_2020<-aggregate(SCBI_2020$tempSingleMaximum, by=list(Category=SCBI_2020$DOY), FUN=mean)
SCBI_mean_2021<-aggregate(SCBI_2021$tempSingleMaximum, by=list(Category=SCBI_2021$DOY), FUN=mean)
SCBI_mean_2022<-aggregate(SCBI_2022$tempSingleMaximum, by=list(Category=SCBI_2022$DOY), FUN=mean)
SCBI_mean_2023<-aggregate(SCBI_2023$tempSingleMaximum, by=list(Category=SCBI_2023$DOY), FUN=mean)



colnames(SCBI_mean_2016)<-c("time","temperature")
colnames(SCBI_mean_2017)<-c("time","temperature")
colnames(SCBI_mean_2018)<-c("time","temperature")
colnames(SCBI_mean_2019)<-c("time","temperature")
colnames(SCBI_mean_2020)<-c("time","temperature")
colnames(SCBI_mean_2021)<-c("time","temperature")
colnames(SCBI_mean_2022)<-c("time","temperature")
colnames(SCBI_mean_2023)<-c("time","temperature")

##create df to compute overall daily average for imputations

temps_SCBI<-list(
  '2016' = SCBI_mean_2016,
  '2017' = SCBI_mean_2017,
  '2018' = SCBI_mean_2018,
  '2019' = SCBI_mean_2019,
  '2020' = SCBI_mean_2020,
  '2021' = SCBI_mean_2021,
  '2022' = SCBI_mean_2022,
  '2023' = SCBI_mean_2023
)

SCBI_res <- impute_with_overall_mean(temps_SCBI)  # temps is your list of yearly dfs
SCBI_overall <- SCBI_res$overall_mean             # average temp for each time (1:366)
SCBI_temps_imputed <- SCBI_res$imputed_years      # list of yearly dfs with missing times filled

SCBI_mean_2016<-SCBI_temps_imputed[["2016"]]
SCBI_mean_2017<-SCBI_temps_imputed[["2017"]]
SCBI_mean_2018<-SCBI_temps_imputed[["2018"]]
SCBI_mean_2019<-SCBI_temps_imputed[["2019"]]
SCBI_mean_2020<-SCBI_temps_imputed[["2020"]]
SCBI_mean_2021<-SCBI_temps_imputed[["2021"]]
SCBI_mean_2022<-SCBI_temps_imputed[["2022"]]
SCBI_mean_2023<-SCBI_temps_imputed[["2023"]]

SCBI_data<- list(SCBI_mean_2016, SCBI_mean_2017, SCBI_mean_2018,
                 SCBI_mean_2019, SCBI_mean_2020, SCBI_mean_2021,
                 SCBI_mean_2022, SCBI_mean_2023)


##loop through years
yr_i = 1

days_SCBI_2016 <- SCBI_data[[yr_i]]$time - 1
temp_SCBI_2016 <- SCBI_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_SCBI_2016 <- c(days_SCBI_2016, (max(days_SCBI_2016) + (days_SCBI_2016[2] - days_SCBI_2016[1])))
temp_extended_SCBI_2016 <- c(temp_SCBI_2016, temp_SCBI_2016[1])

# Fit the spline
spline_fit_SCBI_2016 <- splinefun(days_extended_SCBI_2016, temp_extended_SCBI_2016)

##create the amount of temp data that I need for odin interpolation
DOY_SCBI_2016<-c(0:total)
temperature_SCBI_2016<-f_temp_SCBI_2016(DOY_SCBI_2016)

##loop through years
yr_i = 2

days_SCBI_2017 <- SCBI_data[[yr_i]]$time - 1
temp_SCBI_2017 <- SCBI_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_SCBI_2017 <- c(days_SCBI_2017, (max(days_SCBI_2017) + (days_SCBI_2017[2] - days_SCBI_2017[1])))
temp_extended_SCBI_2017 <- c(temp_SCBI_2017, temp_SCBI_2017[1])

# Fit the spline
spline_fit_SCBI_2017 <- splinefun(days_extended_SCBI_2017, temp_extended_SCBI_2017)

##create the amount of temp data that I need for odin interpolation
DOY_SCBI_2017<-c(0:total)
temperature_SCBI_2017<-f_temp_SCBI_2017(DOY_SCBI_2017)

##loop through years
yr_i = 3

days_SCBI_2018 <- SCBI_data[[yr_i]]$time - 1
temp_SCBI_2018 <- SCBI_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_SCBI_2018 <- c(days_SCBI_2018, (max(days_SCBI_2018) + (days_SCBI_2018[2] - days_SCBI_2018[1])))
temp_extended_SCBI_2018 <- c(temp_SCBI_2018, temp_SCBI_2018[1])

# Fit the spline
spline_fit_SCBI_2018 <- splinefun(days_extended_SCBI_2018, temp_extended_SCBI_2018)

##create the amount of temp data that I need for odin interpolation
DOY_SCBI_2018<-c(0:total)
temperature_SCBI_2018<-f_temp_SCBI_2018(DOY_SCBI_2018)


yr_i = 4

days_SCBI_2019 <- SCBI_data[[yr_i]]$time - 1
temp_SCBI_2019 <- SCBI_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_SCBI_2019 <- c(days_SCBI_2019, (max(days_SCBI_2019) + (days_SCBI_2019[2] - days_SCBI_2019[1])))
temp_extended_SCBI_2019 <- c(temp_SCBI_2019, temp_SCBI_2019[1])

# Fit the spline
spline_fit_SCBI_2019 <- splinefun(days_extended_SCBI_2019, temp_extended_SCBI_2019)

##create the amount of temp data that I need for odin interpolation
DOY_SCBI_2019<-c(0:total)
temperature_SCBI_2019<-f_temp_SCBI_2019(DOY_SCBI_2019)

yr_i = 5

days_SCBI_2020 <- SCBI_data[[yr_i]]$time - 1
temp_SCBI_2020 <- SCBI_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_SCBI_2020 <- c(days_SCBI_2020, (max(days_SCBI_2020) + (days_SCBI_2020[2] - days_SCBI_2020[1])))
temp_extended_SCBI_2020 <- c(temp_SCBI_2020, temp_SCBI_2020[1])

# Fit the spline
spline_fit_SCBI_2020 <- splinefun(days_extended_SCBI_2020, temp_extended_SCBI_2020)

##create the amount of temp data that I need for odin interpolation
DOY_SCBI_2020<-c(0:total)
temperature_SCBI_2020<-f_temp_SCBI_2020(DOY_SCBI_2020)

yr_i = 6

days_SCBI_2021 <- SCBI_data[[yr_i]]$time - 1
temp_SCBI_2021 <- SCBI_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_SCBI_2021 <- c(days_SCBI_2021, (max(days_SCBI_2021) + (days_SCBI_2021[2] - days_SCBI_2021[1])))
temp_extended_SCBI_2021 <- c(temp_SCBI_2021, temp_SCBI_2021[1])

# Fit the spline
spline_fit_SCBI_2021 <- splinefun(days_extended_SCBI_2021, temp_extended_SCBI_2021)

##create the amount of temp data that I need for odin interpolation
DOY_SCBI_2021<-c(0:total)
temperature_SCBI_2021<-f_temp_SCBI_2021(DOY_SCBI_2021)

yr_i = 7

days_SCBI_2022 <- SCBI_data[[yr_i]]$time - 1
temp_SCBI_2022 <- SCBI_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_SCBI_2022 <- c(days_SCBI_2022, (max(days_SCBI_2022) + (days_SCBI_2022[2] - days_SCBI_2022[1])))
temp_extended_SCBI_2022 <- c(temp_SCBI_2022, temp_SCBI_2022[1])

# Fit the spline
spline_fit_SCBI_2022 <- splinefun(days_extended_SCBI_2022, temp_extended_SCBI_2022)

##create the amount of temp data that I need for odin interpolation
DOY_SCBI_2022<-c(0:total)
temperature_SCBI_2022<-f_temp_SCBI_2022(DOY_SCBI_2022)

yr_i = 8

days_SCBI_2023 <- SCBI_data[[yr_i]]$time - 1
temp_SCBI_2023 <- SCBI_data[[yr_i]]$temperature

# Make it periodic by duplicating the boundary
days_extended_SCBI_2023 <- c(days_SCBI_2023, (max(days_SCBI_2023) + (days_SCBI_2023[2] - days_SCBI_2023[1])))
temp_extended_SCBI_2023 <- c(temp_SCBI_2023, temp_SCBI_2023[1])

# Fit the spline
spline_fit_SCBI_2023 <- splinefun(days_extended_SCBI_2023, temp_extended_SCBI_2023)

##create the amount of temp data that I need for odin interpolation
DOY_SCBI_2023<-c(0:total)
temperature_SCBI_2023<-f_temp_SCBI_2023(DOY_SCBI_2023)


#Format data
BLAN_humidity<-read.csv("/BLAN_RH.csv")

#Format data
BLAN_humidity$DOY<-as.numeric(strftime(BLAN_humidity$startDateTime, format = "%j"))
BLAN_humidity<-na.omit(BLAN_humidity)

BLAN_humidity_clean <- BLAN_humidity %>%
  filter(RHMaximum <= 100) %>%
  mutate(year = year(as.POSIXct(startDateTime, tz = "UTC")))

# Split by year

BLAN_humidity_avg <- BLAN_humidity_clean %>% filter(year %in% c(2016:2023))
BLAN_humidity_2016 <- BLAN_humidity_clean %>% filter(year == 2016)
BLAN_humidity_2017 <- BLAN_humidity_clean %>% filter(year == 2017)
BLAN_humidity_2018 <- BLAN_humidity_clean %>% filter(year == 2018)
BLAN_humidity_2019 <- BLAN_humidity_clean %>% filter(year == 2019)
BLAN_humidity_2020 <- BLAN_humidity_clean %>% filter(year == 2020)
BLAN_humidity_2021 <- BLAN_humidity_clean %>% filter(year == 2021)
BLAN_humidity_2022 <- BLAN_humidity_clean %>% filter(year == 2022)
BLAN_humidity_2023 <- BLAN_humidity_clean %>% filter(year == 2023)



BLAN_humidity_mean<-aggregate(BLAN_humidity_avg$RHMaximum, by=list(Category=BLAN_humidity_avg$DOY), FUN=mean)
BLAN_humidity_2016_mean<-aggregate(BLAN_humidity_2016$RHMaximum, by=list(Category=BLAN_humidity_2016$DOY), FUN=mean)
BLAN_humidity_2017_mean<-aggregate(BLAN_humidity_2017$RHMaximum, by=list(Category=BLAN_humidity_2017$DOY), FUN=mean)
BLAN_humidity_2018_mean<-aggregate(BLAN_humidity_2018$RHMaximum, by=list(Category=BLAN_humidity_2018$DOY), FUN=mean)
BLAN_humidity_2019_mean<-aggregate(BLAN_humidity_2019$RHMaximum, by=list(Category=BLAN_humidity_2019$DOY), FUN=mean)
BLAN_humidity_2020_mean<-aggregate(BLAN_humidity_2020$RHMaximum, by=list(Category=BLAN_humidity_2020$DOY), FUN=mean)
BLAN_humidity_2021_mean<-aggregate(BLAN_humidity_2021$RHMaximum, by=list(Category=BLAN_humidity_2021$DOY), FUN=mean)
BLAN_humidity_2022_mean<-aggregate(BLAN_humidity_2022$RHMaximum, by=list(Category=BLAN_humidity_2022$DOY), FUN=mean)
BLAN_humidity_2023_mean<-aggregate(BLAN_humidity_2023$RHMaximum, by=list(Category=BLAN_humidity_2023$DOY), FUN=mean)


colnames(BLAN_humidity_mean)<-c("time","RH")
colnames(BLAN_humidity_2016_mean)<-c("time","RH")
colnames(BLAN_humidity_2017_mean)<-c("time","RH")
colnames(BLAN_humidity_2018_mean)<-c("time","RH")
colnames(BLAN_humidity_2019_mean)<-c("time","RH")
colnames(BLAN_humidity_2020_mean)<-c("time","RH")
colnames(BLAN_humidity_2021_mean)<-c("time","RH")
colnames(BLAN_humidity_2022_mean)<-c("time","RH")
colnames(BLAN_humidity_2023_mean)<-c("time","RH")


##add in missing values

BLAN_humidity_2016_mean <- impute_from_3yravg(
  BLAN_humidity_2016_mean,
  BLAN_humidity_mean
)

BLAN_humidity_2017_mean <- impute_from_3yravg(
  BLAN_humidity_2017_mean,
  BLAN_humidity_mean
)

BLAN_humidity_2018_mean <- impute_from_3yravg(
  BLAN_humidity_2018_mean,
  BLAN_humidity_mean
)

BLAN_humidity_2019_mean <- impute_from_3yravg(
  BLAN_humidity_2019_mean,
  BLAN_humidity_mean
)
BLAN_humidity_2020_mean <- impute_from_3yravg(
  BLAN_humidity_2020_mean,
  BLAN_humidity_mean
)

BLAN_humidity_2021_mean <- impute_from_3yravg(
  BLAN_humidity_2021_mean,
  BLAN_humidity_mean
)

BLAN_humidity_2022_mean <- impute_from_3yravg(
  BLAN_humidity_2022_mean,
  BLAN_humidity_mean
)

BLAN_humidity_2023_mean <- impute_from_3yravg(
  BLAN_humidity_2023_mean,
  BLAN_humidity_mean
)

##Humididty sine curves 

BLAN_humidity_nls_2016 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = BLAN_humidity_2016_mean,
  start = list(
    a = mean(BLAN_humidity_2016_mean$RH),
    b = sd(BLAN_humidity_2016_mean$RH),
    phi = 155
  )
)

BLAN_humidity_nls_2017 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = BLAN_humidity_2017_mean,
  start = list(
    a = mean(BLAN_humidity_2017_mean$RH),
    b = sd(BLAN_humidity_2017_mean$RH),
    phi = 155
  )
)

BLAN_humidity_nls_2018 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = BLAN_humidity_2018_mean,
  start = list(
    a = mean(BLAN_humidity_2018_mean$RH),
    b = sd(BLAN_humidity_2018_mean$RH),
    phi = 155
  )
)

BLAN_humidity_nls_2019 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = BLAN_humidity_2019_mean,
  start = list(
    a = mean(BLAN_humidity_2019_mean$RH),
    b = sd(BLAN_humidity_2019_mean$RH),
    phi = 155
  )
)

BLAN_humidity_nls_2020 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = BLAN_humidity_2020_mean,
  start = list(
    a = mean(BLAN_humidity_2020_mean$RH),
    b = sd(BLAN_humidity_2020_mean$RH),
    phi = 155
  )
)

BLAN_humidity_nls_2021 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = BLAN_humidity_2021_mean,
  start = list(
    a = mean(BLAN_humidity_2021_mean$RH),
    b = sd(BLAN_humidity_2021_mean$RH),
    phi = 155
  )
)

BLAN_humidity_nls_2022 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = BLAN_humidity_2022_mean,
  start = list(
    a = mean(BLAN_humidity_2022_mean$RH),
    b = sd(BLAN_humidity_2022_mean$RH),
    phi = 155
  )
)

BLAN_humidity_nls_2023 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = BLAN_humidity_2023_mean,
  start = list(
    a = mean(BLAN_humidity_2023_mean$RH),
    b = sd(BLAN_humidity_2023_mean$RH),
    phi = 155
  )
)

BLAN_humid_int_2016 = coef(BLAN_humidity_nls_2016)["a"]
BLAN_humid_amp_2016 = coef(BLAN_humidity_nls_2016)["b"]
BLAN_phaseshift_2016 = coef(BLAN_humidity_nls_2016)["phi"]

BLAN_humid_int_2017 = coef(BLAN_humidity_nls_2017)["a"]
BLAN_humid_amp_2017 = coef(BLAN_humidity_nls_2017)["b"]
BLAN_phaseshift_2017 = coef(BLAN_humidity_nls_2017)["phi"]

BLAN_humid_int_2018 = coef(BLAN_humidity_nls_2018)["a"]
BLAN_humid_amp_2018 = coef(BLAN_humidity_nls_2018)["b"]
BLAN_phaseshift_2018 = coef(BLAN_humidity_nls_2018)["phi"]

BLAN_humid_int_2019 = coef(BLAN_humidity_nls_2019)["a"]
BLAN_humid_amp_2019 = coef(BLAN_humidity_nls_2019)["b"]
BLAN_phaseshift_2019 = coef(BLAN_humidity_nls_2019)["phi"]

BLAN_humid_int_2020 = coef(BLAN_humidity_nls_2020)["a"]
BLAN_humid_amp_2020 = coef(BLAN_humidity_nls_2020)["b"]
BLAN_phaseshift_2020 = coef(BLAN_humidity_nls_2020)["phi"]

BLAN_humid_int_2021 = coef(BLAN_humidity_nls_2021)["a"]
BLAN_humid_amp_2021 = coef(BLAN_humidity_nls_2021)["b"]
BLAN_phaseshift_2021 = coef(BLAN_humidity_nls_2021)["phi"]

BLAN_humid_int_2022 = coef(BLAN_humidity_nls_2022)["a"]
BLAN_humid_amp_2022 = coef(BLAN_humidity_nls_2022)["b"]
BLAN_phaseshift_2022 = coef(BLAN_humidity_nls_2022)["phi"]

BLAN_humid_int_2023 = coef(BLAN_humidity_nls_2023)["a"]
BLAN_humid_amp_2023 = coef(BLAN_humidity_nls_2023)["b"]
BLAN_phaseshift_2023 = coef(BLAN_humidity_nls_2023)["phi"]






#HARV
HARV_humidity<-read.csv("/HARV_RH.csv")

#Format data
HARV_humidity$DOY<-as.numeric(strftime(HARV_humidity$startDateTime, format = "%j"))
HARV_humidity<-na.omit(HARV_humidity)

HARV_humidity_clean <- HARV_humidity %>%
  filter(RHMaximum <= 100) %>%
  mutate(year = year(as.POSIXct(startDateTime, tz = "UTC")))

# Split by year

HARV_humidity_avg <- HARV_humidity_clean %>% filter(year %in% c(2016:2023))
HARV_humidity_2016 <- HARV_humidity_clean %>% filter(year == 2016)
HARV_humidity_2017 <- HARV_humidity_clean %>% filter(year == 2017)
HARV_humidity_2018 <- HARV_humidity_clean %>% filter(year == 2018)
HARV_humidity_2019 <- HARV_humidity_clean %>% filter(year == 2019)
HARV_humidity_2020 <- HARV_humidity_clean %>% filter(year == 2020)
HARV_humidity_2021 <- HARV_humidity_clean %>% filter(year == 2021)
HARV_humidity_2022 <- HARV_humidity_clean %>% filter(year == 2022)
HARV_humidity_2023 <- HARV_humidity_clean %>% filter(year == 2023)



HARV_humidity_mean<-aggregate(HARV_humidity_avg$RHMaximum, by=list(Category=HARV_humidity_avg$DOY), FUN=mean)
HARV_humidity_2016_mean<-aggregate(HARV_humidity_2016$RHMaximum, by=list(Category=HARV_humidity_2016$DOY), FUN=mean)
HARV_humidity_2017_mean<-aggregate(HARV_humidity_2017$RHMaximum, by=list(Category=HARV_humidity_2017$DOY), FUN=mean)
HARV_humidity_2018_mean<-aggregate(HARV_humidity_2018$RHMaximum, by=list(Category=HARV_humidity_2018$DOY), FUN=mean)
HARV_humidity_2019_mean<-aggregate(HARV_humidity_2019$RHMaximum, by=list(Category=HARV_humidity_2019$DOY), FUN=mean)
HARV_humidity_2020_mean<-aggregate(HARV_humidity_2020$RHMaximum, by=list(Category=HARV_humidity_2020$DOY), FUN=mean)
HARV_humidity_2021_mean<-aggregate(HARV_humidity_2021$RHMaximum, by=list(Category=HARV_humidity_2021$DOY), FUN=mean)
HARV_humidity_2022_mean<-aggregate(HARV_humidity_2022$RHMaximum, by=list(Category=HARV_humidity_2022$DOY), FUN=mean)
HARV_humidity_2023_mean<-aggregate(HARV_humidity_2023$RHMaximum, by=list(Category=HARV_humidity_2023$DOY), FUN=mean)


colnames(HARV_humidity_mean)<-c("time","RH")
colnames(HARV_humidity_2016_mean)<-c("time","RH")
colnames(HARV_humidity_2017_mean)<-c("time","RH")
colnames(HARV_humidity_2018_mean)<-c("time","RH")
colnames(HARV_humidity_2019_mean)<-c("time","RH")
colnames(HARV_humidity_2020_mean)<-c("time","RH")
colnames(HARV_humidity_2021_mean)<-c("time","RH")
colnames(HARV_humidity_2022_mean)<-c("time","RH")
colnames(HARV_humidity_2023_mean)<-c("time","RH")


##add in missing values

HARV_humidity_2016_mean <- impute_from_3yravg(
  HARV_humidity_2016_mean,
  HARV_humidity_mean
)

HARV_humidity_2017_mean <- impute_from_3yravg(
  HARV_humidity_2017_mean,
  HARV_humidity_mean
)

HARV_humidity_2018_mean <- impute_from_3yravg(
  HARV_humidity_2018_mean,
  HARV_humidity_mean
)

HARV_humidity_2019_mean <- impute_from_3yravg(
  HARV_humidity_2019_mean,
  HARV_humidity_mean
)
HARV_humidity_2020_mean <- impute_from_3yravg(
  HARV_humidity_2020_mean,
  HARV_humidity_mean
)

HARV_humidity_2021_mean <- impute_from_3yravg(
  HARV_humidity_2021_mean,
  HARV_humidity_mean
)

HARV_humidity_2022_mean <- impute_from_3yravg(
  HARV_humidity_2022_mean,
  HARV_humidity_mean
)

HARV_humidity_2023_mean <- impute_from_3yravg(
  HARV_humidity_2023_mean,
  HARV_humidity_mean
)

##Humididty sine curves 

HARV_humidity_nls_2016 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = HARV_humidity_2016_mean,
  start = list(
    a = mean(HARV_humidity_2016_mean$RH),
    b = sd(HARV_humidity_2016_mean$RH),
    phi = 155
  )
)

HARV_humidity_nls_2017 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = HARV_humidity_2017_mean,
  start = list(
    a = mean(HARV_humidity_2017_mean$RH),
    b = sd(HARV_humidity_2017_mean$RH),
    phi = 155
  )
)

HARV_humidity_nls_2018 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = HARV_humidity_2018_mean,
  start = list(
    a = mean(HARV_humidity_2018_mean$RH),
    b = sd(HARV_humidity_2018_mean$RH),
    phi = 155
  )
)

HARV_humidity_nls_2019 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = HARV_humidity_2019_mean,
  start = list(
    a = mean(HARV_humidity_2019_mean$RH),
    b = sd(HARV_humidity_2019_mean$RH),
    phi = 155
  )
)

HARV_humidity_nls_2020 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = HARV_humidity_2020_mean,
  start = list(
    a = mean(HARV_humidity_2020_mean$RH),
    b = sd(HARV_humidity_2020_mean$RH),
    phi = 155
  )
)

HARV_humidity_nls_2021 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = HARV_humidity_2021_mean,
  start = list(
    a = mean(HARV_humidity_2021_mean$RH),
    b = sd(HARV_humidity_2021_mean$RH),
    phi = 155
  )
)

HARV_humidity_nls_2022 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = HARV_humidity_2022_mean,
  start = list(
    a = mean(HARV_humidity_2022_mean$RH),
    b = sd(HARV_humidity_2022_mean$RH),
    phi = 155
  )
)

HARV_humidity_nls_2023 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = HARV_humidity_2023_mean,
  start = list(
    a = mean(HARV_humidity_2023_mean$RH),
    b = sd(HARV_humidity_2023_mean$RH),
    phi = 155
  )
)

HARV_humid_int_2016 = coef(HARV_humidity_nls_2016)["a"]
HARV_humid_amp_2016 = coef(HARV_humidity_nls_2016)["b"]
HARV_phaseshift_2016 = coef(HARV_humidity_nls_2016)["phi"]

HARV_humid_int_2017 = coef(HARV_humidity_nls_2017)["a"]
HARV_humid_amp_2017 = coef(HARV_humidity_nls_2017)["b"]
HARV_phaseshift_2017 = coef(HARV_humidity_nls_2017)["phi"]

HARV_humid_int_2018 = coef(HARV_humidity_nls_2018)["a"]
HARV_humid_amp_2018 = coef(HARV_humidity_nls_2018)["b"]
HARV_phaseshift_2018 = coef(HARV_humidity_nls_2018)["phi"]

HARV_humid_int_2019 = coef(HARV_humidity_nls_2019)["a"]
HARV_humid_amp_2019 = coef(HARV_humidity_nls_2019)["b"]
HARV_phaseshift_2019 = coef(HARV_humidity_nls_2019)["phi"]

HARV_humid_int_2020 = coef(HARV_humidity_nls_2020)["a"]
HARV_humid_amp_2020 = coef(HARV_humidity_nls_2020)["b"]
HARV_phaseshift_2020 = coef(HARV_humidity_nls_2020)["phi"]

HARV_humid_int_2021 = coef(HARV_humidity_nls_2021)["a"]
HARV_humid_amp_2021 = coef(HARV_humidity_nls_2021)["b"]
HARV_phaseshift_2021 = coef(HARV_humidity_nls_2021)["phi"]

HARV_humid_int_2022 = coef(HARV_humidity_nls_2022)["a"]
HARV_humid_amp_2022 = coef(HARV_humidity_nls_2022)["b"]
HARV_phaseshift_2022 = coef(HARV_humidity_nls_2022)["phi"]

HARV_humid_int_2023 = coef(HARV_humidity_nls_2023)["a"]
HARV_humid_amp_2023 = coef(HARV_humidity_nls_2023)["b"]
HARV_phaseshift_2023 = coef(HARV_humidity_nls_2023)["phi"]



#SCBI
SCBI_humidity<-read.csv("~/SCBI_RH.csv")

#Format data
SCBI_humidity$DOY<-as.numeric(strftime(SCBI_humidity$startDateTime, format = "%j"))
SCBI_humidity<-na.omit(SCBI_humidity)

SCBI_humidity_clean <- SCBI_humidity %>%
  filter(RHMaximum <= 100) %>%
  mutate(year = year(as.POSIXct(startDateTime, tz = "UTC")))

# Split by year

SCBI_humidity_avg <- SCBI_humidity_clean %>% filter(year %in% c(2016:2023))
SCBI_humidity_2016 <- SCBI_humidity_clean %>% filter(year == 2016)
SCBI_humidity_2017 <- SCBI_humidity_clean %>% filter(year == 2017)
SCBI_humidity_2018 <- SCBI_humidity_clean %>% filter(year == 2018)
SCBI_humidity_2019 <- SCBI_humidity_clean %>% filter(year == 2019)
SCBI_humidity_2020 <- SCBI_humidity_clean %>% filter(year == 2020)
SCBI_humidity_2021 <- SCBI_humidity_clean %>% filter(year == 2021)
SCBI_humidity_2022 <- SCBI_humidity_clean %>% filter(year == 2022)
SCBI_humidity_2023 <- SCBI_humidity_clean %>% filter(year == 2023)



SCBI_humidity_mean<-aggregate(SCBI_humidity_avg$RHMaximum, by=list(Category=SCBI_humidity_avg$DOY), FUN=mean)
SCBI_humidity_2016_mean<-aggregate(SCBI_humidity_2016$RHMaximum, by=list(Category=SCBI_humidity_2016$DOY), FUN=mean)
SCBI_humidity_2017_mean<-aggregate(SCBI_humidity_2017$RHMaximum, by=list(Category=SCBI_humidity_2017$DOY), FUN=mean)
SCBI_humidity_2018_mean<-aggregate(SCBI_humidity_2018$RHMaximum, by=list(Category=SCBI_humidity_2018$DOY), FUN=mean)
SCBI_humidity_2019_mean<-aggregate(SCBI_humidity_2019$RHMaximum, by=list(Category=SCBI_humidity_2019$DOY), FUN=mean)
SCBI_humidity_2020_mean<-aggregate(SCBI_humidity_2020$RHMaximum, by=list(Category=SCBI_humidity_2020$DOY), FUN=mean)
SCBI_humidity_2021_mean<-aggregate(SCBI_humidity_2021$RHMaximum, by=list(Category=SCBI_humidity_2021$DOY), FUN=mean)
SCBI_humidity_2022_mean<-aggregate(SCBI_humidity_2022$RHMaximum, by=list(Category=SCBI_humidity_2022$DOY), FUN=mean)
SCBI_humidity_2023_mean<-aggregate(SCBI_humidity_2023$RHMaximum, by=list(Category=SCBI_humidity_2023$DOY), FUN=mean)


colnames(SCBI_humidity_mean)<-c("time","RH")
colnames(SCBI_humidity_2016_mean)<-c("time","RH")
colnames(SCBI_humidity_2017_mean)<-c("time","RH")
colnames(SCBI_humidity_2018_mean)<-c("time","RH")
colnames(SCBI_humidity_2019_mean)<-c("time","RH")
colnames(SCBI_humidity_2020_mean)<-c("time","RH")
colnames(SCBI_humidity_2021_mean)<-c("time","RH")
colnames(SCBI_humidity_2022_mean)<-c("time","RH")
colnames(SCBI_humidity_2023_mean)<-c("time","RH")


##add in missing values

SCBI_humidity_2016_mean <- impute_from_3yravg(
  SCBI_humidity_2016_mean,
  SCBI_humidity_mean
)

SCBI_humidity_2017_mean <- impute_from_3yravg(
  SCBI_humidity_2017_mean,
  SCBI_humidity_mean
)

SCBI_humidity_2018_mean <- impute_from_3yravg(
  SCBI_humidity_2018_mean,
  SCBI_humidity_mean
)

SCBI_humidity_2019_mean <- impute_from_3yravg(
  SCBI_humidity_2019_mean,
  SCBI_humidity_mean
)
SCBI_humidity_2020_mean <- impute_from_3yravg(
  SCBI_humidity_2020_mean,
  SCBI_humidity_mean
)

SCBI_humidity_2021_mean <- impute_from_3yravg(
  SCBI_humidity_2021_mean,
  SCBI_humidity_mean
)

SCBI_humidity_2022_mean <- impute_from_3yravg(
  SCBI_humidity_2022_mean,
  SCBI_humidity_mean
)

SCBI_humidity_2023_mean <- impute_from_3yravg(
  SCBI_humidity_2023_mean,
  SCBI_humidity_mean
)

##Humididty sine curves 

SCBI_humidity_nls_2016 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = SCBI_humidity_2016_mean,
  start = list(
    a = mean(SCBI_humidity_2016_mean$RH),
    b = sd(SCBI_humidity_2016_mean$RH),
    phi = 155
  )
)

SCBI_humidity_nls_2017 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = SCBI_humidity_2017_mean,
  start = list(
    a = mean(SCBI_humidity_2017_mean$RH),
    b = sd(SCBI_humidity_2017_mean$RH),
    phi = 155
  )
)

SCBI_humidity_nls_2018 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = SCBI_humidity_2018_mean,
  start = list(
    a = mean(SCBI_humidity_2018_mean$RH),
    b = sd(SCBI_humidity_2018_mean$RH),
    phi = 155
  )
)

SCBI_humidity_nls_2019 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = SCBI_humidity_2019_mean,
  start = list(
    a = mean(SCBI_humidity_2019_mean$RH),
    b = sd(SCBI_humidity_2019_mean$RH),
    phi = 155
  )
)

SCBI_humidity_nls_2020 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = SCBI_humidity_2020_mean,
  start = list(
    a = mean(SCBI_humidity_2020_mean$RH),
    b = sd(SCBI_humidity_2020_mean$RH),
    phi = 155
  )
)

SCBI_humidity_nls_2021 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = SCBI_humidity_2021_mean,
  start = list(
    a = mean(SCBI_humidity_2021_mean$RH),
    b = sd(SCBI_humidity_2021_mean$RH),
    phi = 155
  )
)

SCBI_humidity_nls_2022 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = SCBI_humidity_2022_mean,
  start = list(
    a = mean(SCBI_humidity_2022_mean$RH),
    b = sd(SCBI_humidity_2022_mean$RH),
    phi = 155
  )
)

SCBI_humidity_nls_2023 <- nls(
  RH ~ a + b * sin((2 * pi / 365) * (time - phi)),
  data = SCBI_humidity_2023_mean,
  start = list(
    a = mean(SCBI_humidity_2023_mean$RH),
    b = sd(SCBI_humidity_2023_mean$RH),
    phi = 155
  )
)

SCBI_humid_int_2016 = coef(SCBI_humidity_nls_2016)["a"]
SCBI_humid_amp_2016 = coef(SCBI_humidity_nls_2016)["b"]
SCBI_phaseshift_2016 = coef(SCBI_humidity_nls_2016)["phi"]

SCBI_humid_int_2017 = coef(SCBI_humidity_nls_2017)["a"]
SCBI_humid_amp_2017 = coef(SCBI_humidity_nls_2017)["b"]
SCBI_phaseshift_2017 = coef(SCBI_humidity_nls_2017)["phi"]

SCBI_humid_int_2018 = coef(SCBI_humidity_nls_2018)["a"]
SCBI_humid_amp_2018 = coef(SCBI_humidity_nls_2018)["b"]
SCBI_phaseshift_2018 = coef(SCBI_humidity_nls_2018)["phi"]

SCBI_humid_int_2019 = coef(SCBI_humidity_nls_2019)["a"]
SCBI_humid_amp_2019 = coef(SCBI_humidity_nls_2019)["b"]
SCBI_phaseshift_2019 = coef(SCBI_humidity_nls_2019)["phi"]

SCBI_humid_int_2020 = coef(SCBI_humidity_nls_2020)["a"]
SCBI_humid_amp_2020 = coef(SCBI_humidity_nls_2020)["b"]
SCBI_phaseshift_2020 = coef(SCBI_humidity_nls_2020)["phi"]

SCBI_humid_int_2021 = coef(SCBI_humidity_nls_2021)["a"]
SCBI_humid_amp_2021 = coef(SCBI_humidity_nls_2021)["b"]
SCBI_phaseshift_2021 = coef(SCBI_humidity_nls_2021)["phi"]

SCBI_humid_int_2022 = coef(SCBI_humidity_nls_2022)["a"]
SCBI_humid_amp_2022 = coef(SCBI_humidity_nls_2022)["b"]
SCBI_phaseshift_2022 = coef(SCBI_humidity_nls_2022)["phi"]

SCBI_humid_int_2023 = coef(SCBI_humidity_nls_2023)["a"]
SCBI_humid_amp_2023 = coef(SCBI_humidity_nls_2023)["b"]
SCBI_phaseshift_2023 = coef(SCBI_humidity_nls_2023)["phi"]





##Probability of questing larvae (get max)
theta_i_max<-function(x)
{
  ifelse( x < 0, 0,
          ifelse(x < 26.6,
                 ifelse(-0.2912015 + 0.001641181 * x^2 < 0, 0, 
                        -0.2912015 + 0.001641181 * x^2), 
                 ifelse(3.111262 - 0.08418802 * x < 0, 0,
                        3.111262 - 0.08418802 * x)))
}


##Probability of nymphal questing (from Gilbert) (get max)
theta_n_max<-function(x)
{ 
  ifelse((-0.06003227 + 0.1140132 * x - 0.003410052 * x^2)> 0,
         (-0.06003227 + 0.1140132 * x - 0.003410052 * x^2), 0)
}

##Probability of adult questing (get max)
theta_a_max<-function(x)
{
  
  ifelse((-0.6370415 + 0.3792334 * x - 0.02169659 * x^2) < 0, 0,
         (-0.6370415 + 0.3792334 * x - 0.02169659 * x^2))
}


##solve for theta_max to normalize temperature curve
t_theta<-seq(-30,50,1)
theta_i_norm<-max(theta_i_max(t_theta))
theta_n_norm<-max(theta_n_max(t_theta))
theta_a_norm<-max(theta_a_max(t_theta))

##run to reach equilibrium
t <- seq(0, total-1, by = 1)


BLAN_a_2016_L_new_data <- data.frame(x = c(min(BLAN_a_2016$DOY):max(BLAN_a_2016$DOY)))
BLAN_a_2017_L_new_data <- data.frame(x = c(min(BLAN_a_2017$DOY):max(BLAN_a_2017$DOY)))
BLAN_a_2018_L_new_data <- data.frame(x = c(min(BLAN_a_2018$DOY):max(BLAN_a_2018$DOY)))
BLAN_a_2019_L_new_data <- data.frame(x = c(min(BLAN_a_2019$DOY):max(BLAN_a_2019$DOY)))
BLAN_a_2020_L_new_data <- data.frame(x = c(min(BLAN_a_2020$DOY):max(BLAN_a_2020$DOY)))
BLAN_a_2021_L_new_data <- data.frame(x = c(min(BLAN_a_2021$DOY):max(BLAN_a_2021$DOY)))
BLAN_a_2022_L_new_data <- data.frame(x = c(min(BLAN_a_2022$DOY):max(BLAN_a_2022$DOY)))
BLAN_a_2023_L_new_data <- data.frame(x = c(min(BLAN_a_2023$DOY):max(BLAN_a_2023$DOY)))


HARV_a_2016_L_new_data <- data.frame(x = c(min(HARV_a_2016$DOY):max(HARV_a_2016$DOY)))
HARV_a_2017_L_new_data <- data.frame(x = c(min(HARV_a_2017$DOY):max(HARV_a_2017$DOY)))
HARV_a_2018_L_new_data <- data.frame(x = c(min(HARV_a_2018$DOY):max(HARV_a_2018$DOY)))
HARV_a_2019_L_new_data <- data.frame(x = c(min(HARV_a_2019$DOY):max(HARV_a_2019$DOY)))
HARV_a_2020_L_new_data <- data.frame(x = c(min(HARV_a_2020$DOY):max(HARV_a_2020$DOY)))
HARV_a_2021_L_new_data <- data.frame(x = c(min(HARV_a_2021$DOY):max(HARV_a_2021$DOY)))
HARV_a_2022_L_new_data <- data.frame(x = c(min(HARV_a_2022$DOY):max(HARV_a_2022$DOY)))
HARV_a_2023_L_new_data <- data.frame(x = c(min(HARV_a_2023$DOY):max(HARV_a_2023$DOY)))


SCBI_a_2016_L_new_data <- data.frame(x = c(min(SCBI_a_2016$DOY):max(SCBI_a_2016$DOY)))
SCBI_a_2017_L_new_data <- data.frame(x = c(min(SCBI_a_2017$DOY):max(SCBI_a_2017$DOY)))
SCBI_a_2018_L_new_data <- data.frame(x = c(min(SCBI_a_2018$DOY):max(SCBI_a_2018$DOY)))
SCBI_a_2019_L_new_data <- data.frame(x = c(min(SCBI_a_2019$DOY):max(SCBI_a_2019$DOY)))
SCBI_a_2020_L_new_data <- data.frame(x = c(min(SCBI_a_2020$DOY):max(SCBI_a_2020$DOY)))
SCBI_a_2021_L_new_data <- data.frame(x = c(min(SCBI_a_2021$DOY):max(SCBI_a_2021$DOY)))
SCBI_a_2022_L_new_data <- data.frame(x = c(min(SCBI_a_2022$DOY):max(SCBI_a_2022$DOY)))
SCBI_a_2023_L_new_data <- data.frame(x = c(min(SCBI_a_2023$DOY):max(SCBI_a_2023$DOY)))


sample_window_BLAN_2016 <- integer(365)
sample_window_BLAN_2017 <- integer(365)
sample_window_BLAN_2018 <- integer(365)
sample_window_BLAN_2019 <- integer(365)
sample_window_BLAN_2020 <- integer(365)
sample_window_BLAN_2021 <- integer(365)
sample_window_BLAN_2022 <- integer(365)
sample_window_BLAN_2023 <- integer(365)

sample_window_HARV_2016 <- integer(365)
sample_window_HARV_2017 <- integer(365)
sample_window_HARV_2018 <- integer(365)
sample_window_HARV_2019 <- integer(365)
sample_window_HARV_2020 <- integer(365)
sample_window_HARV_2021 <- integer(365)
sample_window_HARV_2022 <- integer(365)
sample_window_HARV_2023 <- integer(365)

sample_window_SCBI_2016 <- integer(365)
sample_window_SCBI_2017 <- integer(365)
sample_window_SCBI_2018 <- integer(365)
sample_window_SCBI_2019 <- integer(365)
sample_window_SCBI_2020 <- integer(365)
sample_window_SCBI_2021 <- integer(365)
sample_window_SCBI_2022 <- integer(365)
sample_window_SCBI_2023 <- integer(365)




# Set positions corresponding to sample_dates to 1
sample_window_BLAN_2016[BLAN_a_2016_L_new_data[,1]] <- 1
sample_window_BLAN_2017[BLAN_a_2017_L_new_data[,1]] <- 1
sample_window_BLAN_2018[BLAN_a_2018_L_new_data[,1]] <- 1
sample_window_BLAN_2019[BLAN_a_2019_L_new_data[,1]] <- 1
sample_window_BLAN_2020[BLAN_a_2020_L_new_data[,1]] <- 1
sample_window_BLAN_2021[BLAN_a_2021_L_new_data[,1]] <- 1
sample_window_BLAN_2022[BLAN_a_2022_L_new_data[,1]] <- 1
sample_window_BLAN_2023[BLAN_a_2023_L_new_data[,1]] <- 1


sample_window_HARV_2016[HARV_a_2016_L_new_data[,1]] <- 1
sample_window_HARV_2017[HARV_a_2017_L_new_data[,1]] <- 1
sample_window_HARV_2018[HARV_a_2018_L_new_data[,1]] <- 1
sample_window_HARV_2019[HARV_a_2019_L_new_data[,1]] <- 1
sample_window_HARV_2020[HARV_a_2020_L_new_data[,1]] <- 1
sample_window_HARV_2021[HARV_a_2021_L_new_data[,1]] <- 1
sample_window_HARV_2022[HARV_a_2022_L_new_data[,1]] <- 1
sample_window_HARV_2023[HARV_a_2023_L_new_data[,1]] <- 1

sample_window_SCBI_2016[SCBI_a_2016_L_new_data[,1]] <- 1
sample_window_SCBI_2017[SCBI_a_2017_L_new_data[,1]] <- 1
sample_window_SCBI_2018[SCBI_a_2018_L_new_data[,1]] <- 1
sample_window_SCBI_2019[SCBI_a_2019_L_new_data[,1]] <- 1
sample_window_SCBI_2020[SCBI_a_2020_L_new_data[,1]] <- 1
sample_window_SCBI_2021[SCBI_a_2021_L_new_data[,1]] <- 1
sample_window_SCBI_2022[SCBI_a_2022_L_new_data[,1]] <- 1
sample_window_SCBI_2023[SCBI_a_2023_L_new_data[,1]] <- 1


##sort data by DOY
BLAN_a_2016 <- BLAN_a_2016[order(BLAN_a_2016$DOY), ]
BLAN_a_2017 <- BLAN_a_2017[order(BLAN_a_2017$DOY), ]
BLAN_a_2018 <- BLAN_a_2018[order(BLAN_a_2018$DOY), ]
BLAN_a_2019 <- BLAN_a_2019[order(BLAN_a_2019$DOY), ]
BLAN_a_2020 <- BLAN_a_2020[order(BLAN_a_2020$DOY), ]
BLAN_a_2021 <- BLAN_a_2021[order(BLAN_a_2021$DOY), ]
BLAN_a_2022 <- BLAN_a_2022[order(BLAN_a_2022$DOY), ]
BLAN_a_2023 <- BLAN_a_2023[order(BLAN_a_2023$DOY), ]

HARV_a_2016 <- HARV_a_2016[order(HARV_a_2016$DOY), ]
HARV_a_2017 <- HARV_a_2017[order(HARV_a_2017$DOY), ]
HARV_a_2018 <- HARV_a_2018[order(HARV_a_2018$DOY), ]
HARV_a_2019 <- HARV_a_2019[order(HARV_a_2019$DOY), ]
HARV_a_2020 <- HARV_a_2020[order(HARV_a_2020$DOY), ]
HARV_a_2021 <- HARV_a_2021[order(HARV_a_2021$DOY), ]
HARV_a_2022 <- HARV_a_2022[order(HARV_a_2022$DOY), ]
HARV_a_2023 <- HARV_a_2023[order(HARV_a_2023$DOY), ]


SCBI_a_2016 <- SCBI_a_2016[order(SCBI_a_2016$DOY), ]
SCBI_a_2017 <- SCBI_a_2017[order(SCBI_a_2017$DOY), ]
SCBI_a_2018 <- SCBI_a_2018[order(SCBI_a_2018$DOY), ]
SCBI_a_2019 <- SCBI_a_2019[order(SCBI_a_2019$DOY), ]
SCBI_a_2020 <- SCBI_a_2020[order(SCBI_a_2020$DOY), ]
SCBI_a_2021 <- SCBI_a_2021[order(SCBI_a_2021$DOY), ]
SCBI_a_2022 <- SCBI_a_2022[order(SCBI_a_2022$DOY), ]
SCBI_a_2023 <- SCBI_a_2023[order(SCBI_a_2023$DOY), ]


##ODIN MODELS
full_odin <- odin({
  
  
  ##back transform parameters
  dln_min_BLAN = exp(log_dln_BLAN)
  dln_min_HARV = exp(log_dln_HARV)
  dln_min_SCBI = exp(log_dln_SCBI)
  log_dln_BLAN <- parameter()
  log_dln_HARV <- parameter()
  log_dln_SCBI <- parameter()
  
  a_n_BLAN <- parameter(30.5)
  a_n_HARV <- parameter(30.5)
  a_n_SCBI <- parameter(30.5)
  
  a_a_BLAN = exp(log_aa_BLAN)
  a_a_HARV = exp(log_aa_HARV)
  a_a_SCBI = exp(log_aa_SCBI)
  
  log_aa_BLAN <- parameter()
  log_aa_HARV <- parameter()
  log_aa_SCBI <- parameter()
  
  a_l_BLAN = exp(log_al_BLAN)
  a_l_HARV = exp(log_al_HARV)
  a_l_SCBI = exp(log_al_SCBI)
  
  log_al_BLAN <- parameter()
  log_al_HARV <- parameter()
  log_al_SCBI <- parameter()
  
  ##parameters
  y <- parameter(1)
  z <- parameter(21)
  r <- parameter(3)
  u <- parameter(5)
  w <- parameter(10)
  
  
  theta_i_norm <- parameter()
  theta_n_norm <- parameter()
  theta_a_norm <- parameter()
  
  ###fixed parameters across validation sets
  
  mu_e <- parameter(0.002)
  mu_hl <- parameter(0.006)
  mu_ql <- parameter(0.006)
  mu_qn <- parameter(0.006)
  mu_qa <- parameter(0.006)
  mu_el <- parameter(0.003)
  mu_en <- parameter(0.002)
  mu_ea <- parameter(0.0001)
  R <- parameter(2000)
  D <- parameter(20)
  p <- parameter(3000)
  b_l <- parameter(0.4)
  b_n <- parameter(0.4)
  b_a <- parameter(0.4)
  
  eps_n <- parameter(0.15)
  
  lambda_l <- 0.0013 * R ^ 0.515
  lambda_n <- 0.0013 * R ^ 0.515
  lambda_a <- 0.086 * D ^ 0.515
  
  ##BLAN 
  humid_int_BLAN_2016 <- parameter()
  humid_amp_BLAN_2016 <- parameter()
  humid_ps_BLAN_2016 <- parameter()
  
  humid_int_BLAN_2017 <- parameter()
  humid_amp_BLAN_2017 <- parameter()
  humid_ps_BLAN_2017 <- parameter()
  
  humid_int_BLAN_2018 <- parameter()
  humid_amp_BLAN_2018 <- parameter()
  humid_ps_BLAN_2018 <- parameter()
  
  humid_int_BLAN_2019 <- parameter()
  humid_amp_BLAN_2019 <- parameter()
  humid_ps_BLAN_2019 <- parameter()
  
  
  
  ## derived parameters, validation specific
  humidity_BLAN_2016 <- humid_int_BLAN_2016  +  humid_amp_BLAN_2016 *sin((2*pi/365)*(time-humid_ps_BLAN_2016))
  humidity_BLAN_2017 <- humid_int_BLAN_2017  +  humid_amp_BLAN_2017 *sin((2*pi/365)*(time-humid_ps_BLAN_2017))
  humidity_BLAN_2018 <- humid_int_BLAN_2018  +  humid_amp_BLAN_2018 *sin((2*pi/365)*(time-humid_ps_BLAN_2018))
  humidity_BLAN_2019 <- humid_int_BLAN_2019  +  humid_amp_BLAN_2019 *sin((2*pi/365)*(time-humid_ps_BLAN_2019))
  
  #BLAN2016 model
  
  
  #initial conditions
  initial(ELA_BLAN_2016) <- 0
  initial(E_BLAN_2016) <- 5000
  initial(HL_BLAN_2016) <- 0
  initial(RFL_BLAN_2016) <- 0
  initial(FL_BLAN_2016) <- 0
  initial(EL_BLAN_2016) <- 0
  initial(RFN_BLAN_2016) <- 0
  initial(FNabove_BLAN_2016) <- 0
  initial(FNbelow_BLAN_2016) <- 0
  initial(EN_BLAN_2016) <- 0
  initial(RFA_BLAN_2016) <- 0
  initial(FA_BLAN_2016) <- 0
  initial(EA_BLAN_2016) <- 0
  
  #equations
  deriv(ELA_BLAN_2016) <- x_BLAN_2016 * EA_BLAN_2016 - 1/y * ELA_BLAN_2016
  deriv(E_BLAN_2016) <- ELA_BLAN_2016 * (1 - (0.01 + (0.04 * log(1.01 + FA_BLAN_2016)/D))) * p - (q_BLAN_2016 + mu_e) * E_BLAN_2016
  deriv(HL_BLAN_2016) <- q_BLAN_2016 * E_BLAN_2016 - (1/z + mu_hl) * HL_BLAN_2016
  deriv(RFL_BLAN_2016) <- 1/z * HL_BLAN_2016 - (lambda_l * theta_i_BLAN_2016 * 1 / (1 + exp(-(-(a_n_BLAN+(a_l_BLAN)) + b_l * humidity_BLAN_2016))) + mu_ql) * RFL_BLAN_2016
  deriv(FL_BLAN_2016) <- lambda_l * theta_i_BLAN_2016 * 1 / (1 + exp(-(-(a_n_BLAN+(a_l_BLAN)) + b_l * humidity_BLAN_2016))) * RFL_BLAN_2016 - (1/r + (0.65 + (0.049 * log(1.01 + FL_BLAN_2016)/R))) * FL_BLAN_2016
  deriv(EL_BLAN_2016) <- 1/r * FL_BLAN_2016 - (s_BLAN_2016 + mu_el) * EL_BLAN_2016
  deriv(RFN_BLAN_2016) <- s_BLAN_2016 * EL_BLAN_2016 -
    (lambda_n * theta_n_BLAN_2016 * 1 / (1 + exp(-(-a_n_BLAN + b_n * humidity_BLAN_2016))) + lambda_n * theta_n_BLAN_2016 * eps_n + mu_qn) * RFN_BLAN_2016
  deriv(FNabove_BLAN_2016) <- lambda_n * theta_n_BLAN_2016 * 1 / (1 + exp(-(-a_n_BLAN + b_n * humidity_BLAN_2016))) * RFN_BLAN_2016 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_BLAN_2016 + FNbelow_BLAN_2016))/R))) * FNabove_BLAN_2016
  deriv(FNbelow_BLAN_2016) <- lambda_n * theta_n_BLAN_2016 * eps_n * RFN_BLAN_2016 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_BLAN_2016 + FNbelow_BLAN_2016))/R))) * FNbelow_BLAN_2016
  deriv(EN_BLAN_2016) <- 1/u * (FNabove_BLAN_2016 + FNbelow_BLAN_2016) - (v_BLAN_2016 + mu_en) * EN_BLAN_2016
  deriv(RFA_BLAN_2016) <- v_BLAN_2016 * EN_BLAN_2016 - (lambda_a * theta_a_BLAN_2016 *  1 / (1 + exp(-(-(a_n_BLAN-(a_a_BLAN)) + b_a * humidity_BLAN_2016))) + mu_qa) * RFA_BLAN_2016
  deriv(FA_BLAN_2016) <- 0.5 * lambda_a * theta_a_BLAN_2016 *  1 / (1 + exp(-(-(a_n_BLAN-(a_a_BLAN)) + b_a * humidity_BLAN_2016))) * RFA_BLAN_2016 - (1/w + (0.5 + (0.049 * log(1.01 + FA_BLAN_2016)/D))) * FA_BLAN_2016
  deriv(EA_BLAN_2016) <- 1/w * FA_BLAN_2016 - (x_BLAN_2016 + mu_ea) * EA_BLAN_2016
  
  
  theta_i_raw_BLAN_2016 <- if (temp_BLAN_2016 < 26.6) (-0.2912015 + 0.001641181 * temp_BLAN_2016^2) else (3.111262 - 0.08418802 * temp_BLAN_2016)
  theta_i_unscaled_BLAN_2016 <- if (theta_i_raw_BLAN_2016 > 0) theta_i_raw_BLAN_2016 else 0
  theta_i_BLAN_2016 <- if (temp_BLAN_2016 <= 0) 0 else theta_i_unscaled_BLAN_2016 / theta_i_norm
  #
  theta_n_raw_BLAN_2016 <- (-0.06003227 + 0.1140132 * temp_BLAN_2016 - 0.003410052 * temp_BLAN_2016^2)
  theta_n_BLAN_2016 <- if (theta_n_raw_BLAN_2016 > 0) theta_n_raw_BLAN_2016 / theta_n_norm else 0
  #
  theta_a_raw_BLAN_2016 <- (-0.6370415 + 0.3792334 * temp_BLAN_2016 - 0.02169659 * temp_BLAN_2016^2)
  theta_a_BLAN_2016 <- if (theta_a_raw_BLAN_2016 > 0) theta_a_raw_BLAN_2016 / theta_a_norm else 0
  
  q_BLAN_2016 <- if( temp_BLAN_2016 > 0) temp_BLAN_2016^2.27/34234 else 0
  #
  s_BLAN_2016 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2016 > 0) max((dln_min_BLAN/2000), temp_BLAN_2016^2.55/101181) else (dln_min_BLAN/2000)
  else 0
  #
  #
  v_BLAN_2016 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2016 > 0) temp_BLAN_2016^1.21/1596 else 0
  else 0
  #
  x_BLAN_2016 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2016 > 0) temp_BLAN_2016^1.42/1300 else 0
  else 0
  
  
  ##interpolated functions
  temp_BLAN_2016 <- interpolate(DOY_BLAN_2016, temperature_BLAN_2016, "spline")
  DOY_BLAN_2016 <- parameter(constant = TRUE)
  temperature_BLAN_2016 <- parameter(constant = TRUE)
  
  
  #dimesions
  dim(DOY_BLAN_2016, temperature_BLAN_2016) <- parameter(rank = 1)
  
  #BLAN2017
  
  initial(ELA_BLAN_2017) <- 0
  initial(E_BLAN_2017) <- 5000
  initial(HL_BLAN_2017) <- 0
  initial(RFL_BLAN_2017) <- 0
  initial(FL_BLAN_2017) <- 0
  initial(EL_BLAN_2017) <- 0
  initial(RFN_BLAN_2017) <- 0
  initial(FNabove_BLAN_2017) <- 0
  initial(FNbelow_BLAN_2017) <- 0
  initial(EN_BLAN_2017) <- 0
  initial(RFA_BLAN_2017) <- 0
  initial(FA_BLAN_2017) <- 0
  initial(EA_BLAN_2017) <- 0
  
  #equations
  deriv(ELA_BLAN_2017) <- x_BLAN_2017 * EA_BLAN_2017 - 1/y * ELA_BLAN_2017
  deriv(E_BLAN_2017) <- ELA_BLAN_2017 * (1 - (0.01 + (0.04 * log(1.01 + FA_BLAN_2017)/D))) * p - (q_BLAN_2017 + mu_e) * E_BLAN_2017
  deriv(HL_BLAN_2017) <- q_BLAN_2017 * E_BLAN_2017 - (1/z + mu_hl) * HL_BLAN_2017
  deriv(RFL_BLAN_2017) <- 1/z * HL_BLAN_2017 - (lambda_l * theta_i_BLAN_2017 * 1 / (1 + exp(-(-(a_n_BLAN+(a_l_BLAN)) + b_l * humidity_BLAN_2017))) + mu_ql) * RFL_BLAN_2017
  deriv(FL_BLAN_2017) <- lambda_l * theta_i_BLAN_2017 * 1 / (1 + exp(-(-(a_n_BLAN+(a_l_BLAN)) + b_l * humidity_BLAN_2017))) * RFL_BLAN_2017 - (1/r + (0.65 + (0.049 * log(1.01 + FL_BLAN_2017)/R))) * FL_BLAN_2017
  deriv(EL_BLAN_2017) <- 1/r * FL_BLAN_2017 - (s_BLAN_2017 + mu_el) * EL_BLAN_2017
  deriv(RFN_BLAN_2017) <- s_BLAN_2017 * EL_BLAN_2017 -
    (lambda_n * theta_n_BLAN_2017 * 1 / (1 + exp(-(-a_n_BLAN + b_n * humidity_BLAN_2017))) + lambda_n * theta_n_BLAN_2017 * eps_n + mu_qn) * RFN_BLAN_2017
  deriv(FNabove_BLAN_2017) <- lambda_n * theta_n_BLAN_2017 * 1 / (1 + exp(-(-a_n_BLAN + b_n * humidity_BLAN_2017))) * RFN_BLAN_2017 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_BLAN_2017 + FNbelow_BLAN_2017))/R))) * FNabove_BLAN_2017
  deriv(FNbelow_BLAN_2017) <- lambda_n * theta_n_BLAN_2017 * eps_n * RFN_BLAN_2017 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_BLAN_2017 + FNbelow_BLAN_2017))/R))) * FNbelow_BLAN_2017
  deriv(EN_BLAN_2017) <- 1/u * (FNabove_BLAN_2017 + FNbelow_BLAN_2017) - (v_BLAN_2017 + mu_en) * EN_BLAN_2017
  deriv(RFA_BLAN_2017) <- v_BLAN_2017 * EN_BLAN_2017 - (lambda_a * theta_a_BLAN_2017 *  1 / (1 + exp(-(-(a_n_BLAN-(a_a_BLAN)) + b_a * humidity_BLAN_2017))) + mu_qa) * RFA_BLAN_2017
  deriv(FA_BLAN_2017) <- 0.5 * lambda_a * theta_a_BLAN_2017 *  1 / (1 + exp(-(-(a_n_BLAN-(a_a_BLAN)) + b_a * humidity_BLAN_2017))) * RFA_BLAN_2017 - (1/w + (0.5 + (0.049 * log(1.01 + FA_BLAN_2017)/D))) * FA_BLAN_2017
  deriv(EA_BLAN_2017) <- 1/w * FA_BLAN_2017 - (x_BLAN_2017 + mu_ea) * EA_BLAN_2017
  
  theta_i_raw_BLAN_2017 <- if (temp_BLAN_2017 < 26.6) (-0.2912015 + 0.001641181 * temp_BLAN_2017^2) else (3.111262 - 0.08418802 * temp_BLAN_2017)
  theta_i_unscaled_BLAN_2017 <- if (theta_i_raw_BLAN_2017 > 0) theta_i_raw_BLAN_2017 else 0
  theta_i_BLAN_2017 <- if (temp_BLAN_2017 <= 0) 0 else theta_i_unscaled_BLAN_2017 / theta_i_norm
  #
  theta_n_raw_BLAN_2017 <- (-0.06003227 + 0.1140132 * temp_BLAN_2017 - 0.003410052 * temp_BLAN_2017^2)
  theta_n_BLAN_2017 <- if (theta_n_raw_BLAN_2017 > 0) theta_n_raw_BLAN_2017 / theta_n_norm else 0
  #
  theta_a_raw_BLAN_2017 <- (-0.6370415 + 0.3792334 * temp_BLAN_2017 - 0.02169659 * temp_BLAN_2017^2)
  theta_a_BLAN_2017 <- if (theta_a_raw_BLAN_2017 > 0) theta_a_raw_BLAN_2017 / theta_a_norm else 0
  
  q_BLAN_2017 <- if( temp_BLAN_2017 > 0) temp_BLAN_2017^2.27/34234 else 0
  #
  s_BLAN_2017 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2017 > 0) max((dln_min_BLAN/2000), temp_BLAN_2017^2.55/101181) else (dln_min_BLAN/2000)
  else 0
  #
  #
  v_BLAN_2017 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2017 > 0) temp_BLAN_2017^1.21/1596 else 0
  else 0
  #
  x_BLAN_2017 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2017 > 0) temp_BLAN_2017^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_BLAN_2017 <- interpolate(DOY_BLAN_2017, temperature_BLAN_2017, "spline")
  DOY_BLAN_2017 <- parameter(constant = TRUE)
  temperature_BLAN_2017 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_BLAN_2017, temperature_BLAN_2017) <- parameter(rank = 1)
  
  #BLAN2018
  
  initial(ELA_BLAN_2018) <- 0
  initial(E_BLAN_2018) <- 5000
  initial(HL_BLAN_2018) <- 0
  initial(RFL_BLAN_2018) <- 0
  initial(FL_BLAN_2018) <- 0
  initial(EL_BLAN_2018) <- 0
  initial(RFN_BLAN_2018) <- 0
  initial(FNabove_BLAN_2018) <- 0
  initial(FNbelow_BLAN_2018) <- 0
  initial(EN_BLAN_2018) <- 0
  initial(RFA_BLAN_2018) <- 0
  initial(FA_BLAN_2018) <- 0
  initial(EA_BLAN_2018) <- 0
  
  #equations
  deriv(ELA_BLAN_2018) <- x_BLAN_2018 * EA_BLAN_2018 - 1/y * ELA_BLAN_2018
  deriv(E_BLAN_2018) <- ELA_BLAN_2018 * (1 - (0.01 + (0.04 * log(1.01 + FA_BLAN_2018)/D))) * p - (q_BLAN_2018 + mu_e) * E_BLAN_2018
  deriv(HL_BLAN_2018) <- q_BLAN_2018 * E_BLAN_2018 - (1/z + mu_hl) * HL_BLAN_2018
  deriv(RFL_BLAN_2018) <- 1/z * HL_BLAN_2018 - (lambda_l * theta_i_BLAN_2018 * 1 / (1 + exp(-(-(a_n_BLAN+(a_l_BLAN)) + b_l * humidity_BLAN_2018))) + mu_ql) * RFL_BLAN_2018
  deriv(FL_BLAN_2018) <- lambda_l * theta_i_BLAN_2018 * 1 / (1 + exp(-(-(a_n_BLAN+(a_l_BLAN)) + b_l * humidity_BLAN_2018))) * RFL_BLAN_2018 - (1/r + (0.65 + (0.049 * log(1.01 + FL_BLAN_2018)/R))) * FL_BLAN_2018
  deriv(EL_BLAN_2018) <- 1/r * FL_BLAN_2018 - (s_BLAN_2018 + mu_el) * EL_BLAN_2018
  deriv(RFN_BLAN_2018) <- s_BLAN_2018 * EL_BLAN_2018 -
    (lambda_n * theta_n_BLAN_2018 * 1 / (1 + exp(-(-a_n_BLAN + b_n * humidity_BLAN_2018))) + lambda_n * theta_n_BLAN_2018 * eps_n + mu_qn) * RFN_BLAN_2018
  deriv(FNabove_BLAN_2018) <- lambda_n * theta_n_BLAN_2018 * 1 / (1 + exp(-(-a_n_BLAN + b_n * humidity_BLAN_2018))) * RFN_BLAN_2018 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_BLAN_2018 + FNbelow_BLAN_2018))/R))) * FNabove_BLAN_2018
  deriv(FNbelow_BLAN_2018) <- lambda_n * theta_n_BLAN_2018 * eps_n * RFN_BLAN_2018 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_BLAN_2018 + FNbelow_BLAN_2018))/R))) * FNbelow_BLAN_2018
  deriv(EN_BLAN_2018) <- 1/u * (FNabove_BLAN_2018 + FNbelow_BLAN_2018) - (v_BLAN_2018 + mu_en) * EN_BLAN_2018
  deriv(RFA_BLAN_2018) <- v_BLAN_2018 * EN_BLAN_2018 - (lambda_a * theta_a_BLAN_2018 *  1 / (1 + exp(-(-(a_n_BLAN-(a_a_BLAN)) + b_a * humidity_BLAN_2018))) + mu_qa) * RFA_BLAN_2018
  deriv(FA_BLAN_2018) <- 0.5 * lambda_a * theta_a_BLAN_2018 *  1 / (1 + exp(-(-(a_n_BLAN-(a_a_BLAN)) + b_a * humidity_BLAN_2018))) * RFA_BLAN_2018 - (1/w + (0.5 + (0.049 * log(1.01 + FA_BLAN_2018)/D))) * FA_BLAN_2018
  deriv(EA_BLAN_2018) <- 1/w * FA_BLAN_2018 - (x_BLAN_2018 + mu_ea) * EA_BLAN_2018
  
  theta_i_raw_BLAN_2018 <- if (temp_BLAN_2018 < 26.6) (-0.2912015 + 0.001641181 * temp_BLAN_2018^2) else (3.111262 - 0.08418802 * temp_BLAN_2018)
  theta_i_unscaled_BLAN_2018 <- if (theta_i_raw_BLAN_2018 > 0) theta_i_raw_BLAN_2018 else 0
  theta_i_BLAN_2018 <- if (temp_BLAN_2018 <= 0) 0 else theta_i_unscaled_BLAN_2018 / theta_i_norm
  #
  theta_n_raw_BLAN_2018 <- (-0.06003227 + 0.1140132 * temp_BLAN_2018 - 0.003410052 * temp_BLAN_2018^2)
  theta_n_BLAN_2018 <- if (theta_n_raw_BLAN_2018 > 0) theta_n_raw_BLAN_2018 / theta_n_norm else 0
  #
  theta_a_raw_BLAN_2018 <- (-0.6370415 + 0.3792334 * temp_BLAN_2018 - 0.02169659 * temp_BLAN_2018^2)
  theta_a_BLAN_2018 <- if (theta_a_raw_BLAN_2018 > 0) theta_a_raw_BLAN_2018 / theta_a_norm else 0
  
  q_BLAN_2018 <- if( temp_BLAN_2018 > 0) temp_BLAN_2018^2.27/34234 else 0
  #
  s_BLAN_2018 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2018 > 0) max((dln_min_BLAN/2000), temp_BLAN_2018^2.55/101181) else (dln_min_BLAN/2000)
  else 0
  #
  #
  v_BLAN_2018 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2018 > 0) temp_BLAN_2018^1.21/1596 else 0
  else 0
  #
  x_BLAN_2018 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2018 > 0) temp_BLAN_2018^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_BLAN_2018 <- interpolate(DOY_BLAN_2018, temperature_BLAN_2018, "spline")
  DOY_BLAN_2018 <- parameter(constant = TRUE)
  temperature_BLAN_2018 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_BLAN_2018, temperature_BLAN_2018) <- parameter(rank = 1)
  
  #BLAN2018
  
  initial(ELA_BLAN_2019) <- 0
  initial(E_BLAN_2019) <- 5000
  initial(HL_BLAN_2019) <- 0
  initial(RFL_BLAN_2019) <- 0
  initial(FL_BLAN_2019) <- 0
  initial(EL_BLAN_2019) <- 0
  initial(RFN_BLAN_2019) <- 0
  initial(FNabove_BLAN_2019) <- 0
  initial(FNbelow_BLAN_2019) <- 0
  initial(EN_BLAN_2019) <- 0
  initial(RFA_BLAN_2019) <- 0
  initial(FA_BLAN_2019) <- 0
  initial(EA_BLAN_2019) <- 0
  
  #equations
  deriv(ELA_BLAN_2019) <- x_BLAN_2019 * EA_BLAN_2019 - 1/y * ELA_BLAN_2019
  deriv(E_BLAN_2019) <- ELA_BLAN_2019 * (1 - (0.01 + (0.04 * log(1.01 + FA_BLAN_2019)/D))) * p - (q_BLAN_2019 + mu_e) * E_BLAN_2019
  deriv(HL_BLAN_2019) <- q_BLAN_2019 * E_BLAN_2019 - (1/z + mu_hl) * HL_BLAN_2019
  deriv(RFL_BLAN_2019) <- 1/z * HL_BLAN_2019 - (lambda_l * theta_i_BLAN_2019 * 1 / (1 + exp(-(-(a_n_BLAN+(a_l_BLAN)) + b_l * humidity_BLAN_2019))) + mu_ql) * RFL_BLAN_2019
  deriv(FL_BLAN_2019) <- lambda_l * theta_i_BLAN_2019 * 1 / (1 + exp(-(-(a_n_BLAN+(a_l_BLAN)) + b_l * humidity_BLAN_2019))) * RFL_BLAN_2019 - (1/r + (0.65 + (0.049 * log(1.01 + FL_BLAN_2019)/R))) * FL_BLAN_2019
  deriv(EL_BLAN_2019) <- 1/r * FL_BLAN_2019 - (s_BLAN_2019 + mu_el) * EL_BLAN_2019
  deriv(RFN_BLAN_2019) <- s_BLAN_2019 * EL_BLAN_2019 -
    (lambda_n * theta_n_BLAN_2019 * 1 / (1 + exp(-(-a_n_BLAN + b_n * humidity_BLAN_2019))) + lambda_n * theta_n_BLAN_2019 * eps_n + mu_qn) * RFN_BLAN_2019
  deriv(FNabove_BLAN_2019) <- lambda_n * theta_n_BLAN_2019 * 1 / (1 + exp(-(-a_n_BLAN + b_n * humidity_BLAN_2019))) * RFN_BLAN_2019 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_BLAN_2019 + FNbelow_BLAN_2019))/R))) * FNabove_BLAN_2019
  deriv(FNbelow_BLAN_2019) <- lambda_n * theta_n_BLAN_2019 * eps_n * RFN_BLAN_2019 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_BLAN_2019 + FNbelow_BLAN_2019))/R))) * FNbelow_BLAN_2019
  deriv(EN_BLAN_2019) <- 1/u * (FNabove_BLAN_2019 + FNbelow_BLAN_2019) - (v_BLAN_2019 + mu_en) * EN_BLAN_2019
  deriv(RFA_BLAN_2019) <- v_BLAN_2019 * EN_BLAN_2019 - (lambda_a * theta_a_BLAN_2019 *  1 / (1 + exp(-(-(a_n_BLAN-(a_a_BLAN)) + b_a * humidity_BLAN_2019))) + mu_qa) * RFA_BLAN_2019
  deriv(FA_BLAN_2019) <- 0.5 * lambda_a * theta_a_BLAN_2019 *  1 / (1 + exp(-(-(a_n_BLAN-(a_a_BLAN)) + b_a * humidity_BLAN_2019))) * RFA_BLAN_2019 - (1/w + (0.5 + (0.049 * log(1.01 + FA_BLAN_2019)/D))) * FA_BLAN_2019
  deriv(EA_BLAN_2019) <- 1/w * FA_BLAN_2019 - (x_BLAN_2019 + mu_ea) * EA_BLAN_2019
  
  theta_i_raw_BLAN_2019 <- if (temp_BLAN_2019 < 26.6) (-0.2912015 + 0.001641181 * temp_BLAN_2019^2) else (3.111262 - 0.08418802 * temp_BLAN_2019)
  theta_i_unscaled_BLAN_2019 <- if (theta_i_raw_BLAN_2019 > 0) theta_i_raw_BLAN_2019 else 0
  theta_i_BLAN_2019 <- if (temp_BLAN_2019 <= 0) 0 else theta_i_unscaled_BLAN_2019 / theta_i_norm
  #
  theta_n_raw_BLAN_2019 <- (-0.06003227 + 0.1140132 * temp_BLAN_2019 - 0.003410052 * temp_BLAN_2019^2)
  theta_n_BLAN_2019 <- if (theta_n_raw_BLAN_2019 > 0) theta_n_raw_BLAN_2019 / theta_n_norm else 0
  #
  theta_a_raw_BLAN_2019 <- (-0.6370415 + 0.3792334 * temp_BLAN_2019 - 0.02169659 * temp_BLAN_2019^2)
  theta_a_BLAN_2019 <- if (theta_a_raw_BLAN_2019 > 0) theta_a_raw_BLAN_2019 / theta_a_norm else 0
  
  q_BLAN_2019 <- if( temp_BLAN_2019 > 0) temp_BLAN_2019^2.27/34234 else 0
  #
  s_BLAN_2019 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2019 > 0) max((dln_min_BLAN/2000), temp_BLAN_2019^2.55/101181) else (dln_min_BLAN/2000)
  else 0
  #
  #
  v_BLAN_2019 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2019 > 0) temp_BLAN_2019^1.21/1596 else 0
  else 0
  #
  x_BLAN_2019 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_BLAN_2019 > 0) temp_BLAN_2019^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_BLAN_2019 <- interpolate(DOY_BLAN_2019, temperature_BLAN_2019, "spline")
  DOY_BLAN_2019 <- parameter(constant = TRUE)
  temperature_BLAN_2019 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_BLAN_2019, temperature_BLAN_2019) <- parameter(rank = 1)
  
  
  
  
  ##HARV 
  
  
  humid_int_HARV_2017 <- parameter()
  humid_amp_HARV_2017 <- parameter()
  humid_ps_HARV_2017 <- parameter()
  
  humid_int_HARV_2018 <- parameter()
  humid_amp_HARV_2018 <- parameter()
  humid_ps_HARV_2018 <- parameter()
  
  humid_int_HARV_2021 <- parameter()
  humid_amp_HARV_2021 <- parameter()
  humid_ps_HARV_2021 <- parameter()
  
  humid_int_HARV_2023 <- parameter()
  humid_amp_HARV_2023 <- parameter()
  humid_ps_HARV_2023 <- parameter()
  
  ## derived parameters, validation specific
  humidity_HARV_2017 <- humid_int_HARV_2017  +  humid_amp_HARV_2017 *sin((2*pi/365)*(time-humid_ps_HARV_2017))
  humidity_HARV_2018 <- humid_int_HARV_2018  +  humid_amp_HARV_2018 *sin((2*pi/365)*(time-humid_ps_HARV_2018))
  humidity_HARV_2021 <- humid_int_HARV_2021  +  humid_amp_HARV_2021 *sin((2*pi/365)*(time-humid_ps_HARV_2021))
  humidity_HARV_2023 <- humid_int_HARV_2023  +  humid_amp_HARV_2023 *sin((2*pi/365)*(time-humid_ps_HARV_2023))
  
  
  #HARV2017
  
  initial(ELA_HARV_2017) <- 0
  initial(E_HARV_2017) <- 5000
  initial(HL_HARV_2017) <- 0
  initial(RFL_HARV_2017) <- 0
  initial(FL_HARV_2017) <- 0
  initial(EL_HARV_2017) <- 0
  initial(RFN_HARV_2017) <- 0
  initial(FNabove_HARV_2017) <- 0
  initial(FNbelow_HARV_2017) <- 0
  initial(EN_HARV_2017) <- 0
  initial(RFA_HARV_2017) <- 0
  initial(FA_HARV_2017) <- 0
  initial(EA_HARV_2017) <- 0
  
  #equations
  deriv(ELA_HARV_2017) <- x_HARV_2017 * EA_HARV_2017 - 1/y * ELA_HARV_2017
  deriv(E_HARV_2017) <- ELA_HARV_2017 * (1 - (0.01 + (0.04 * log(1.01 + FA_HARV_2017)/D))) * p - (q_HARV_2017 + mu_e) * E_HARV_2017
  deriv(HL_HARV_2017) <- q_HARV_2017 * E_HARV_2017 - (1/z + mu_hl) * HL_HARV_2017
  deriv(RFL_HARV_2017) <- 1/z * HL_HARV_2017 - (lambda_l * theta_i_HARV_2017 * 1 / (1 + exp(-(-(a_n_HARV+(a_l_HARV)) + b_l * humidity_HARV_2017))) + mu_ql) * RFL_HARV_2017
  deriv(FL_HARV_2017) <- lambda_l * theta_i_HARV_2017 * 1 / (1 + exp(-(-(a_n_HARV+(a_l_HARV)) + b_l * humidity_HARV_2017))) * RFL_HARV_2017 - (1/r + (0.65 + (0.049 * log(1.01 + FL_HARV_2017)/R))) * FL_HARV_2017
  deriv(EL_HARV_2017) <- 1/r * FL_HARV_2017 - (s_HARV_2017 + mu_el) * EL_HARV_2017
  deriv(RFN_HARV_2017) <- s_HARV_2017 * EL_HARV_2017 -
    (lambda_n * theta_n_HARV_2017 * 1 / (1 + exp(-(-a_n_HARV + b_n * humidity_HARV_2017))) + lambda_n * theta_n_HARV_2017 * eps_n + mu_qn) * RFN_HARV_2017
  deriv(FNabove_HARV_2017) <- lambda_n * theta_n_HARV_2017 * 1 / (1 + exp(-(-a_n_HARV + b_n * humidity_HARV_2017))) * RFN_HARV_2017 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_HARV_2017 + FNbelow_HARV_2017))/R))) * FNabove_HARV_2017
  deriv(FNbelow_HARV_2017) <- lambda_n * theta_n_HARV_2017 * eps_n * RFN_HARV_2017 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_HARV_2017 + FNbelow_HARV_2017))/R))) * FNbelow_HARV_2017
  deriv(EN_HARV_2017) <- 1/u * (FNabove_HARV_2017 + FNbelow_HARV_2017) - (v_HARV_2017 + mu_en) * EN_HARV_2017
  deriv(RFA_HARV_2017) <- v_HARV_2017 * EN_HARV_2017 - (lambda_a * theta_a_HARV_2017 *  1 / (1 + exp(-(-(a_n_HARV-(a_a_HARV)) + b_a * humidity_HARV_2017))) + mu_qa) * RFA_HARV_2017
  deriv(FA_HARV_2017) <- 0.5 * lambda_a * theta_a_HARV_2017 *  1 / (1 + exp(-(-(a_n_HARV-(a_a_HARV)) + b_a * humidity_HARV_2017))) * RFA_HARV_2017 - (1/w + (0.5 + (0.049 * log(1.01 + FA_HARV_2017)/D))) * FA_HARV_2017
  deriv(EA_HARV_2017) <- 1/w * FA_HARV_2017 - (x_HARV_2017 + mu_ea) * EA_HARV_2017
  
  theta_i_raw_HARV_2017 <- if (temp_HARV_2017 < 26.6) (-0.2912015 + 0.001641181 * temp_HARV_2017^2) else (3.111262 - 0.08418802 * temp_HARV_2017)
  theta_i_unscaled_HARV_2017 <- if (theta_i_raw_HARV_2017 > 0) theta_i_raw_HARV_2017 else 0
  theta_i_HARV_2017 <- if (temp_HARV_2017 <= 0) 0 else theta_i_unscaled_HARV_2017 / theta_i_norm
  #
  theta_n_raw_HARV_2017 <- (-0.06003227 + 0.1140132 * temp_HARV_2017 - 0.003410052 * temp_HARV_2017^2)
  theta_n_HARV_2017 <- if (theta_n_raw_HARV_2017 > 0) theta_n_raw_HARV_2017 / theta_n_norm else 0
  #
  theta_a_raw_HARV_2017 <- (-0.6370415 + 0.3792334 * temp_HARV_2017 - 0.02169659 * temp_HARV_2017^2)
  theta_a_HARV_2017 <- if (theta_a_raw_HARV_2017 > 0) theta_a_raw_HARV_2017 / theta_a_norm else 0
  
  q_HARV_2017 <- if( temp_HARV_2017 > 0) temp_HARV_2017^2.27/34234 else 0
  #
  s_HARV_2017 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2017 > 0) max((dln_min_HARV/2000), temp_HARV_2017^2.55/101181) else (dln_min_HARV/2000)
  else 0
  #
  #
  v_HARV_2017 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2017 > 0) temp_HARV_2017^1.21/1596 else 0
  else 0
  #
  x_HARV_2017 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2017 > 0) temp_HARV_2017^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_HARV_2017 <- interpolate(DOY_HARV_2017, temperature_HARV_2017, "spline")
  DOY_HARV_2017 <- parameter(constant = TRUE)
  temperature_HARV_2017 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_HARV_2017, temperature_HARV_2017) <- parameter(rank = 1)
  
  #HARV2018
  
  initial(ELA_HARV_2018) <- 0
  initial(E_HARV_2018) <- 5000
  initial(HL_HARV_2018) <- 0
  initial(RFL_HARV_2018) <- 0
  initial(FL_HARV_2018) <- 0
  initial(EL_HARV_2018) <- 0
  initial(RFN_HARV_2018) <- 0
  initial(FNabove_HARV_2018) <- 0
  initial(FNbelow_HARV_2018) <- 0
  initial(EN_HARV_2018) <- 0
  initial(RFA_HARV_2018) <- 0
  initial(FA_HARV_2018) <- 0
  initial(EA_HARV_2018) <- 0
  
  #equations
  deriv(ELA_HARV_2018) <- x_HARV_2018 * EA_HARV_2018 - 1/y * ELA_HARV_2018
  deriv(E_HARV_2018) <- ELA_HARV_2018 * (1 - (0.01 + (0.04 * log(1.01 + FA_HARV_2018)/D))) * p - (q_HARV_2018 + mu_e) * E_HARV_2018
  deriv(HL_HARV_2018) <- q_HARV_2018 * E_HARV_2018 - (1/z + mu_hl) * HL_HARV_2018
  deriv(RFL_HARV_2018) <- 1/z * HL_HARV_2018 - (lambda_l * theta_i_HARV_2018 * 1 / (1 + exp(-(-(a_n_HARV+(a_l_HARV)) + b_l * humidity_HARV_2018))) + mu_ql) * RFL_HARV_2018
  deriv(FL_HARV_2018) <- lambda_l * theta_i_HARV_2018 * 1 / (1 + exp(-(-(a_n_HARV+(a_l_HARV)) + b_l * humidity_HARV_2018))) * RFL_HARV_2018 - (1/r + (0.65 + (0.049 * log(1.01 + FL_HARV_2018)/R))) * FL_HARV_2018
  deriv(EL_HARV_2018) <- 1/r * FL_HARV_2018 - (s_HARV_2018 + mu_el) * EL_HARV_2018
  deriv(RFN_HARV_2018) <- s_HARV_2018 * EL_HARV_2018 -
    (lambda_n * theta_n_HARV_2018 * 1 / (1 + exp(-(-a_n_HARV + b_n * humidity_HARV_2018))) + lambda_n * theta_n_HARV_2018 * eps_n + mu_qn) * RFN_HARV_2018
  deriv(FNabove_HARV_2018) <- lambda_n * theta_n_HARV_2018 * 1 / (1 + exp(-(-a_n_HARV + b_n * humidity_HARV_2018))) * RFN_HARV_2018 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_HARV_2018 + FNbelow_HARV_2018))/R))) * FNabove_HARV_2018
  deriv(FNbelow_HARV_2018) <- lambda_n * theta_n_HARV_2018 * eps_n * RFN_HARV_2018 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_HARV_2018 + FNbelow_HARV_2018))/R))) * FNbelow_HARV_2018
  deriv(EN_HARV_2018) <- 1/u * (FNabove_HARV_2018 + FNbelow_HARV_2018) - (v_HARV_2018 + mu_en) * EN_HARV_2018
  deriv(RFA_HARV_2018) <- v_HARV_2018 * EN_HARV_2018 - (lambda_a * theta_a_HARV_2018 *  1 / (1 + exp(-(-(a_n_HARV-(a_a_HARV)) + b_a * humidity_HARV_2018))) + mu_qa) * RFA_HARV_2018
  deriv(FA_HARV_2018) <- 0.5 * lambda_a * theta_a_HARV_2018 *  1 / (1 + exp(-(-(a_n_HARV-(a_a_HARV)) + b_a * humidity_HARV_2018))) * RFA_HARV_2018 - (1/w + (0.5 + (0.049 * log(1.01 + FA_HARV_2018)/D))) * FA_HARV_2018
  deriv(EA_HARV_2018) <- 1/w * FA_HARV_2018 - (x_HARV_2018 + mu_ea) * EA_HARV_2018
  
  theta_i_raw_HARV_2018 <- if (temp_HARV_2018 < 26.6) (-0.2912015 + 0.001641181 * temp_HARV_2018^2) else (3.111262 - 0.08418802 * temp_HARV_2018)
  theta_i_unscaled_HARV_2018 <- if (theta_i_raw_HARV_2018 > 0) theta_i_raw_HARV_2018 else 0
  theta_i_HARV_2018 <- if (temp_HARV_2018 <= 0) 0 else theta_i_unscaled_HARV_2018 / theta_i_norm
  #
  theta_n_raw_HARV_2018 <- (-0.06003227 + 0.1140132 * temp_HARV_2018 - 0.003410052 * temp_HARV_2018^2)
  theta_n_HARV_2018 <- if (theta_n_raw_HARV_2018 > 0) theta_n_raw_HARV_2018 / theta_n_norm else 0
  #
  theta_a_raw_HARV_2018 <- (-0.6370415 + 0.3792334 * temp_HARV_2018 - 0.02169659 * temp_HARV_2018^2)
  theta_a_HARV_2018 <- if (theta_a_raw_HARV_2018 > 0) theta_a_raw_HARV_2018 / theta_a_norm else 0
  
  q_HARV_2018 <- if( temp_HARV_2018 > 0) temp_HARV_2018^2.27/34234 else 0
  #
  s_HARV_2018 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2018 > 0) max((dln_min_HARV/2000), temp_HARV_2018^2.55/101181) else (dln_min_HARV/2000)
  else 0
  #
  #
  v_HARV_2018 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2018 > 0) temp_HARV_2018^1.21/1596 else 0
  else 0
  #
  x_HARV_2018 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2018 > 0) temp_HARV_2018^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_HARV_2018 <- interpolate(DOY_HARV_2018, temperature_HARV_2018, "spline")
  DOY_HARV_2018 <- parameter(constant = TRUE)
  temperature_HARV_2018 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_HARV_2018, temperature_HARV_2018) <- parameter(rank = 1)
  
  
  #HARV2021
  
  initial(ELA_HARV_2021) <- 0
  initial(E_HARV_2021) <- 5000
  initial(HL_HARV_2021) <- 0
  initial(RFL_HARV_2021) <- 0
  initial(FL_HARV_2021) <- 0
  initial(EL_HARV_2021) <- 0
  initial(RFN_HARV_2021) <- 0
  initial(FNabove_HARV_2021) <- 0
  initial(FNbelow_HARV_2021) <- 0
  initial(EN_HARV_2021) <- 0
  initial(RFA_HARV_2021) <- 0
  initial(FA_HARV_2021) <- 0
  initial(EA_HARV_2021) <- 0
  
  #equations
  deriv(ELA_HARV_2021) <- x_HARV_2021 * EA_HARV_2021 - 1/y * ELA_HARV_2021
  deriv(E_HARV_2021) <- ELA_HARV_2021 * (1 - (0.01 + (0.04 * log(1.01 + FA_HARV_2021)/D))) * p - (q_HARV_2021 + mu_e) * E_HARV_2021
  deriv(HL_HARV_2021) <- q_HARV_2021 * E_HARV_2021 - (1/z + mu_hl) * HL_HARV_2021
  deriv(RFL_HARV_2021) <- 1/z * HL_HARV_2021 - (lambda_l * theta_i_HARV_2021 * 1 / (1 + exp(-(-(a_n_HARV+(a_l_HARV)) + b_l * humidity_HARV_2021))) + mu_ql) * RFL_HARV_2021
  deriv(FL_HARV_2021) <- lambda_l * theta_i_HARV_2021 * 1 / (1 + exp(-(-(a_n_HARV+(a_l_HARV)) + b_l * humidity_HARV_2021))) * RFL_HARV_2021 - (1/r + (0.65 + (0.049 * log(1.01 + FL_HARV_2021)/R))) * FL_HARV_2021
  deriv(EL_HARV_2021) <- 1/r * FL_HARV_2021 - (s_HARV_2021 + mu_el) * EL_HARV_2021
  deriv(RFN_HARV_2021) <- s_HARV_2021 * EL_HARV_2021 -
    (lambda_n * theta_n_HARV_2021 * 1 / (1 + exp(-(-a_n_HARV + b_n * humidity_HARV_2021))) + lambda_n * theta_n_HARV_2021 * eps_n + mu_qn) * RFN_HARV_2021
  deriv(FNabove_HARV_2021) <- lambda_n * theta_n_HARV_2021 * 1 / (1 + exp(-(-a_n_HARV + b_n * humidity_HARV_2021))) * RFN_HARV_2021 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_HARV_2021 + FNbelow_HARV_2021))/R))) * FNabove_HARV_2021
  deriv(FNbelow_HARV_2021) <- lambda_n * theta_n_HARV_2021 * eps_n * RFN_HARV_2021 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_HARV_2021 + FNbelow_HARV_2021))/R))) * FNbelow_HARV_2021
  deriv(EN_HARV_2021) <- 1/u * (FNabove_HARV_2021 + FNbelow_HARV_2021) - (v_HARV_2021 + mu_en) * EN_HARV_2021
  deriv(RFA_HARV_2021) <- v_HARV_2021 * EN_HARV_2021 - (lambda_a * theta_a_HARV_2021 *  1 / (1 + exp(-(-(a_n_HARV-(a_a_HARV)) + b_a * humidity_HARV_2021))) + mu_qa) * RFA_HARV_2021
  deriv(FA_HARV_2021) <- 0.5 * lambda_a * theta_a_HARV_2021 *  1 / (1 + exp(-(-(a_n_HARV-(a_a_HARV)) + b_a * humidity_HARV_2021))) * RFA_HARV_2021 - (1/w + (0.5 + (0.049 * log(1.01 + FA_HARV_2021)/D))) * FA_HARV_2021
  deriv(EA_HARV_2021) <- 1/w * FA_HARV_2021 - (x_HARV_2021 + mu_ea) * EA_HARV_2021
  
  theta_i_raw_HARV_2021 <- if (temp_HARV_2021 < 26.6) (-0.2912015 + 0.001641181 * temp_HARV_2021^2) else (3.111262 - 0.08418802 * temp_HARV_2021)
  theta_i_unscaled_HARV_2021 <- if (theta_i_raw_HARV_2021 > 0) theta_i_raw_HARV_2021 else 0
  theta_i_HARV_2021 <- if (temp_HARV_2021 <= 0) 0 else theta_i_unscaled_HARV_2021 / theta_i_norm
  #
  theta_n_raw_HARV_2021 <- (-0.06003227 + 0.1140132 * temp_HARV_2021 - 0.003410052 * temp_HARV_2021^2)
  theta_n_HARV_2021 <- if (theta_n_raw_HARV_2021 > 0) theta_n_raw_HARV_2021 / theta_n_norm else 0
  #
  theta_a_raw_HARV_2021 <- (-0.6370415 + 0.3792334 * temp_HARV_2021 - 0.02169659 * temp_HARV_2021^2)
  theta_a_HARV_2021 <- if (theta_a_raw_HARV_2021 > 0) theta_a_raw_HARV_2021 / theta_a_norm else 0
  
  q_HARV_2021 <- if( temp_HARV_2021 > 0) temp_HARV_2021^2.27/34234 else 0
  #
  s_HARV_2021 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2021 > 0) max((dln_min_HARV/2000), temp_HARV_2021^2.55/101181) else (dln_min_HARV/2000)
  else 0
  #
  #
  v_HARV_2021 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2021 > 0) temp_HARV_2021^1.21/1596 else 0
  else 0
  #
  x_HARV_2021 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2021 > 0) temp_HARV_2021^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_HARV_2021 <- interpolate(DOY_HARV_2021, temperature_HARV_2021, "spline")
  DOY_HARV_2021 <- parameter(constant = TRUE)
  temperature_HARV_2021 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_HARV_2021, temperature_HARV_2021) <- parameter(rank = 1)
  
  
  
  
  #HARV2023
  
  initial(ELA_HARV_2023) <- 0
  initial(E_HARV_2023) <- 5000
  initial(HL_HARV_2023) <- 0
  initial(RFL_HARV_2023) <- 0
  initial(FL_HARV_2023) <- 0
  initial(EL_HARV_2023) <- 0
  initial(RFN_HARV_2023) <- 0
  initial(FNabove_HARV_2023) <- 0
  initial(FNbelow_HARV_2023) <- 0
  initial(EN_HARV_2023) <- 0
  initial(RFA_HARV_2023) <- 0
  initial(FA_HARV_2023) <- 0
  initial(EA_HARV_2023) <- 0
  
  #equations
  deriv(ELA_HARV_2023) <- x_HARV_2023 * EA_HARV_2023 - 1/y * ELA_HARV_2023
  deriv(E_HARV_2023) <- ELA_HARV_2023 * (1 - (0.01 + (0.04 * log(1.01 + FA_HARV_2023)/D))) * p - (q_HARV_2023 + mu_e) * E_HARV_2023
  deriv(HL_HARV_2023) <- q_HARV_2023 * E_HARV_2023 - (1/z + mu_hl) * HL_HARV_2023
  deriv(RFL_HARV_2023) <- 1/z * HL_HARV_2023 - (lambda_l * theta_i_HARV_2023 * 1 / (1 + exp(-(-(a_n_HARV+(a_l_HARV)) + b_l * humidity_HARV_2023))) + mu_ql) * RFL_HARV_2023
  deriv(FL_HARV_2023) <- lambda_l * theta_i_HARV_2023 * 1 / (1 + exp(-(-(a_n_HARV+(a_l_HARV)) + b_l * humidity_HARV_2023))) * RFL_HARV_2023 - (1/r + (0.65 + (0.049 * log(1.01 + FL_HARV_2023)/R))) * FL_HARV_2023
  deriv(EL_HARV_2023) <- 1/r * FL_HARV_2023 - (s_HARV_2023 + mu_el) * EL_HARV_2023
  deriv(RFN_HARV_2023) <- s_HARV_2023 * EL_HARV_2023 -
    (lambda_n * theta_n_HARV_2023 * 1 / (1 + exp(-(-a_n_HARV + b_n * humidity_HARV_2023))) + lambda_n * theta_n_HARV_2023 * eps_n + mu_qn) * RFN_HARV_2023
  deriv(FNabove_HARV_2023) <- lambda_n * theta_n_HARV_2023 * 1 / (1 + exp(-(-a_n_HARV + b_n * humidity_HARV_2023))) * RFN_HARV_2023 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_HARV_2023 + FNbelow_HARV_2023))/R))) * FNabove_HARV_2023
  deriv(FNbelow_HARV_2023) <- lambda_n * theta_n_HARV_2023 * eps_n * RFN_HARV_2023 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_HARV_2023 + FNbelow_HARV_2023))/R))) * FNbelow_HARV_2023
  deriv(EN_HARV_2023) <- 1/u * (FNabove_HARV_2023 + FNbelow_HARV_2023) - (v_HARV_2023 + mu_en) * EN_HARV_2023
  deriv(RFA_HARV_2023) <- v_HARV_2023 * EN_HARV_2023 - (lambda_a * theta_a_HARV_2023 *  1 / (1 + exp(-(-(a_n_HARV-(a_a_HARV)) + b_a * humidity_HARV_2023))) + mu_qa) * RFA_HARV_2023
  deriv(FA_HARV_2023) <- 0.5 * lambda_a * theta_a_HARV_2023 *  1 / (1 + exp(-(-(a_n_HARV-(a_a_HARV)) + b_a * humidity_HARV_2023))) * RFA_HARV_2023 - (1/w + (0.5 + (0.049 * log(1.01 + FA_HARV_2023)/D))) * FA_HARV_2023
  deriv(EA_HARV_2023) <- 1/w * FA_HARV_2023 - (x_HARV_2023 + mu_ea) * EA_HARV_2023
  
  theta_i_raw_HARV_2023 <- if (temp_HARV_2023 < 26.6) (-0.2912015 + 0.001641181 * temp_HARV_2023^2) else (3.111262 - 0.08418802 * temp_HARV_2023)
  theta_i_unscaled_HARV_2023 <- if (theta_i_raw_HARV_2023 > 0) theta_i_raw_HARV_2023 else 0
  theta_i_HARV_2023 <- if (temp_HARV_2023 <= 0) 0 else theta_i_unscaled_HARV_2023 / theta_i_norm
  #
  theta_n_raw_HARV_2023 <- (-0.06003227 + 0.1140132 * temp_HARV_2023 - 0.003410052 * temp_HARV_2023^2)
  theta_n_HARV_2023 <- if (theta_n_raw_HARV_2023 > 0) theta_n_raw_HARV_2023 / theta_n_norm else 0
  #
  theta_a_raw_HARV_2023 <- (-0.6370415 + 0.3792334 * temp_HARV_2023 - 0.02169659 * temp_HARV_2023^2)
  theta_a_HARV_2023 <- if (theta_a_raw_HARV_2023 > 0) theta_a_raw_HARV_2023 / theta_a_norm else 0
  
  q_HARV_2023 <- if( temp_HARV_2023 > 0) temp_HARV_2023^2.27/34234 else 0
  #
  s_HARV_2023 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2023 > 0) max((dln_min_HARV/2000), temp_HARV_2023^2.55/101181) else (dln_min_HARV/2000)
  else 0
  #
  #
  v_HARV_2023 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2023 > 0) temp_HARV_2023^1.21/1596 else 0
  else 0
  #
  x_HARV_2023 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_HARV_2023 > 0) temp_HARV_2023^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_HARV_2023 <- interpolate(DOY_HARV_2023, temperature_HARV_2023, "spline")
  DOY_HARV_2023 <- parameter(constant = TRUE)
  temperature_HARV_2023 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_HARV_2023, temperature_HARV_2023) <- parameter(rank = 1) 
  
  ##SCBI 
  humid_int_SCBI_2016 <- parameter()
  humid_amp_SCBI_2016 <- parameter()
  humid_ps_SCBI_2016 <- parameter()
  
  humid_int_SCBI_2017 <- parameter()
  humid_amp_SCBI_2017 <- parameter()
  humid_ps_SCBI_2017 <- parameter()
  
  humid_int_SCBI_2018 <- parameter()
  humid_amp_SCBI_2018 <- parameter()
  humid_ps_SCBI_2018 <- parameter()
  
  humid_int_SCBI_2019 <- parameter()
  humid_amp_SCBI_2019 <- parameter()
  humid_ps_SCBI_2019 <- parameter()
  
  humid_int_SCBI_2023 <- parameter()
  humid_amp_SCBI_2023 <- parameter()
  humid_ps_SCBI_2023 <- parameter()
  
  ## derived parameters, validation specific
  humidity_SCBI_2016 <- humid_int_SCBI_2016  +  humid_amp_SCBI_2016 *sin((2*pi/365)*(time-humid_ps_SCBI_2016))
  humidity_SCBI_2017 <- humid_int_SCBI_2017  +  humid_amp_SCBI_2017 *sin((2*pi/365)*(time-humid_ps_SCBI_2017))
  humidity_SCBI_2018 <- humid_int_SCBI_2018  +  humid_amp_SCBI_2018 *sin((2*pi/365)*(time-humid_ps_SCBI_2018))
  humidity_SCBI_2019 <- humid_int_SCBI_2019  +  humid_amp_SCBI_2019 *sin((2*pi/365)*(time-humid_ps_SCBI_2019))
  humidity_SCBI_2023 <- humid_int_SCBI_2023  +  humid_amp_SCBI_2023 *sin((2*pi/365)*(time-humid_ps_SCBI_2023))
  
  
  #SCBI2016 model
  
  
  #initial conditions
  initial(ELA_SCBI_2016) <- 0
  initial(E_SCBI_2016) <- 5000
  initial(HL_SCBI_2016) <- 0
  initial(RFL_SCBI_2016) <- 0
  initial(FL_SCBI_2016) <- 0
  initial(EL_SCBI_2016) <- 0
  initial(RFN_SCBI_2016) <- 0
  initial(FNabove_SCBI_2016) <- 0
  initial(FNbelow_SCBI_2016) <- 0
  initial(EN_SCBI_2016) <- 0
  initial(RFA_SCBI_2016) <- 0
  initial(FA_SCBI_2016) <- 0
  initial(EA_SCBI_2016) <- 0
  
  #equations
  deriv(ELA_SCBI_2016) <- x_SCBI_2016 * EA_SCBI_2016 - 1/y * ELA_SCBI_2016
  deriv(E_SCBI_2016) <- ELA_SCBI_2016 * (1 - (0.01 + (0.04 * log(1.01 + FA_SCBI_2016)/D))) * p - (q_SCBI_2016 + mu_e) * E_SCBI_2016
  deriv(HL_SCBI_2016) <- q_SCBI_2016 * E_SCBI_2016 - (1/z + mu_hl) * HL_SCBI_2016
  deriv(RFL_SCBI_2016) <- 1/z * HL_SCBI_2016 - (lambda_l * theta_i_SCBI_2016 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2016))) + mu_ql) * RFL_SCBI_2016
  deriv(FL_SCBI_2016) <- lambda_l * theta_i_SCBI_2016 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2016))) * RFL_SCBI_2016 - (1/r + (0.65 + (0.049 * log(1.01 + FL_SCBI_2016)/R))) * FL_SCBI_2016
  deriv(EL_SCBI_2016) <- 1/r * FL_SCBI_2016 - (s_SCBI_2016 + mu_el) * EL_SCBI_2016
  deriv(RFN_SCBI_2016) <- s_SCBI_2016 * EL_SCBI_2016 -
    (lambda_n * theta_n_SCBI_2016 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2016))) + lambda_n * theta_n_SCBI_2016 * eps_n + mu_qn) * RFN_SCBI_2016
  deriv(FNabove_SCBI_2016) <- lambda_n * theta_n_SCBI_2016 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2016))) * RFN_SCBI_2016 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2016 + FNbelow_SCBI_2016))/R))) * FNabove_SCBI_2016
  deriv(FNbelow_SCBI_2016) <- lambda_n * theta_n_SCBI_2016 * eps_n * RFN_SCBI_2016 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2016 + FNbelow_SCBI_2016))/R))) * FNbelow_SCBI_2016
  deriv(EN_SCBI_2016) <- 1/u * (FNabove_SCBI_2016 + FNbelow_SCBI_2016) - (v_SCBI_2016 + mu_en) * EN_SCBI_2016
  deriv(RFA_SCBI_2016) <- v_SCBI_2016 * EN_SCBI_2016 - (lambda_a * theta_a_SCBI_2016 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2016))) + mu_qa) * RFA_SCBI_2016
  deriv(FA_SCBI_2016) <- 0.5 * lambda_a * theta_a_SCBI_2016 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2016))) * RFA_SCBI_2016 - (1/w + (0.5 + (0.049 * log(1.01 + FA_SCBI_2016)/D))) * FA_SCBI_2016
  deriv(EA_SCBI_2016) <- 1/w * FA_SCBI_2016 - (x_SCBI_2016 + mu_ea) * EA_SCBI_2016
  
  
  theta_i_raw_SCBI_2016 <- if (temp_SCBI_2016 < 26.6) (-0.2912015 + 0.001641181 * temp_SCBI_2016^2) else (3.111262 - 0.08418802 * temp_SCBI_2016)
  theta_i_unscaled_SCBI_2016 <- if (theta_i_raw_SCBI_2016 > 0) theta_i_raw_SCBI_2016 else 0
  theta_i_SCBI_2016 <- if (temp_SCBI_2016 <= 0) 0 else theta_i_unscaled_SCBI_2016 / theta_i_norm
  #
  theta_n_raw_SCBI_2016 <- (-0.06003227 + 0.1140132 * temp_SCBI_2016 - 0.003410052 * temp_SCBI_2016^2)
  theta_n_SCBI_2016 <- if (theta_n_raw_SCBI_2016 > 0) theta_n_raw_SCBI_2016 / theta_n_norm else 0
  #
  theta_a_raw_SCBI_2016 <- (-0.6370415 + 0.3792334 * temp_SCBI_2016 - 0.02169659 * temp_SCBI_2016^2)
  theta_a_SCBI_2016 <- if (theta_a_raw_SCBI_2016 > 0) theta_a_raw_SCBI_2016 / theta_a_norm else 0
  
  q_SCBI_2016 <- if( temp_SCBI_2016 > 0) temp_SCBI_2016^2.27/34234 else 0
  #
  s_SCBI_2016 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2016 > 0) max((dln_min_SCBI/2000), temp_SCBI_2016^2.55/101181) else (dln_min_SCBI/2000)
  else 0
  #
  #
  v_SCBI_2016 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2016 > 0) temp_SCBI_2016^1.21/1596 else 0
  else 0
  #
  x_SCBI_2016 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2016 > 0) temp_SCBI_2016^1.42/1300 else 0
  else 0
  
  
  ##interpolated functions
  temp_SCBI_2016 <- interpolate(DOY_SCBI_2016, temperature_SCBI_2016, "spline")
  DOY_SCBI_2016 <- parameter(constant = TRUE)
  temperature_SCBI_2016 <- parameter(constant = TRUE)
  
  
  #dimesions
  dim(DOY_SCBI_2016, temperature_SCBI_2016) <- parameter(rank = 1)
  
  #SCBI2017
  
  initial(ELA_SCBI_2017) <- 0
  initial(E_SCBI_2017) <- 5000
  initial(HL_SCBI_2017) <- 0
  initial(RFL_SCBI_2017) <- 0
  initial(FL_SCBI_2017) <- 0
  initial(EL_SCBI_2017) <- 0
  initial(RFN_SCBI_2017) <- 0
  initial(FNabove_SCBI_2017) <- 0
  initial(FNbelow_SCBI_2017) <- 0
  initial(EN_SCBI_2017) <- 0
  initial(RFA_SCBI_2017) <- 0
  initial(FA_SCBI_2017) <- 0
  initial(EA_SCBI_2017) <- 0
  
  #equations
  deriv(ELA_SCBI_2017) <- x_SCBI_2017 * EA_SCBI_2017 - 1/y * ELA_SCBI_2017
  deriv(E_SCBI_2017) <- ELA_SCBI_2017 * (1 - (0.01 + (0.04 * log(1.01 + FA_SCBI_2017)/D))) * p - (q_SCBI_2017 + mu_e) * E_SCBI_2017
  deriv(HL_SCBI_2017) <- q_SCBI_2017 * E_SCBI_2017 - (1/z + mu_hl) * HL_SCBI_2017
  deriv(RFL_SCBI_2017) <- 1/z * HL_SCBI_2017 - (lambda_l * theta_i_SCBI_2017 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2017))) + mu_ql) * RFL_SCBI_2017
  deriv(FL_SCBI_2017) <- lambda_l * theta_i_SCBI_2017 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2017))) * RFL_SCBI_2017 - (1/r + (0.65 + (0.049 * log(1.01 + FL_SCBI_2017)/R))) * FL_SCBI_2017
  deriv(EL_SCBI_2017) <- 1/r * FL_SCBI_2017 - (s_SCBI_2017 + mu_el) * EL_SCBI_2017
  deriv(RFN_SCBI_2017) <- s_SCBI_2017 * EL_SCBI_2017 -
    (lambda_n * theta_n_SCBI_2017 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2017))) + lambda_n * theta_n_SCBI_2017 * eps_n + mu_qn) * RFN_SCBI_2017
  deriv(FNabove_SCBI_2017) <- lambda_n * theta_n_SCBI_2017 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2017))) * RFN_SCBI_2017 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2017 + FNbelow_SCBI_2017))/R))) * FNabove_SCBI_2017
  deriv(FNbelow_SCBI_2017) <- lambda_n * theta_n_SCBI_2017 * eps_n * RFN_SCBI_2017 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2017 + FNbelow_SCBI_2017))/R))) * FNbelow_SCBI_2017
  deriv(EN_SCBI_2017) <- 1/u * (FNabove_SCBI_2017 + FNbelow_SCBI_2017) - (v_SCBI_2017 + mu_en) * EN_SCBI_2017
  deriv(RFA_SCBI_2017) <- v_SCBI_2017 * EN_SCBI_2017 - (lambda_a * theta_a_SCBI_2017 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2017))) + mu_qa) * RFA_SCBI_2017
  deriv(FA_SCBI_2017) <- 0.5 * lambda_a * theta_a_SCBI_2017 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2017))) * RFA_SCBI_2017 - (1/w + (0.5 + (0.049 * log(1.01 + FA_SCBI_2017)/D))) * FA_SCBI_2017
  deriv(EA_SCBI_2017) <- 1/w * FA_SCBI_2017 - (x_SCBI_2017 + mu_ea) * EA_SCBI_2017
  
  theta_i_raw_SCBI_2017 <- if (temp_SCBI_2017 < 26.6) (-0.2912015 + 0.001641181 * temp_SCBI_2017^2) else (3.111262 - 0.08418802 * temp_SCBI_2017)
  theta_i_unscaled_SCBI_2017 <- if (theta_i_raw_SCBI_2017 > 0) theta_i_raw_SCBI_2017 else 0
  theta_i_SCBI_2017 <- if (temp_SCBI_2017 <= 0) 0 else theta_i_unscaled_SCBI_2017 / theta_i_norm
  #
  theta_n_raw_SCBI_2017 <- (-0.06003227 + 0.1140132 * temp_SCBI_2017 - 0.003410052 * temp_SCBI_2017^2)
  theta_n_SCBI_2017 <- if (theta_n_raw_SCBI_2017 > 0) theta_n_raw_SCBI_2017 / theta_n_norm else 0
  #
  theta_a_raw_SCBI_2017 <- (-0.6370415 + 0.3792334 * temp_SCBI_2017 - 0.02169659 * temp_SCBI_2017^2)
  theta_a_SCBI_2017 <- if (theta_a_raw_SCBI_2017 > 0) theta_a_raw_SCBI_2017 / theta_a_norm else 0
  
  q_SCBI_2017 <- if( temp_SCBI_2017 > 0) temp_SCBI_2017^2.27/34234 else 0
  #
  s_SCBI_2017 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2017 > 0) max((dln_min_SCBI/2000), temp_SCBI_2017^2.55/101181) else (dln_min_SCBI/2000)
  else 0
  #
  #
  v_SCBI_2017 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2017 > 0) temp_SCBI_2017^1.21/1596 else 0
  else 0
  #
  x_SCBI_2017 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2017 > 0) temp_SCBI_2017^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_SCBI_2017 <- interpolate(DOY_SCBI_2017, temperature_SCBI_2017, "spline")
  DOY_SCBI_2017 <- parameter(constant = TRUE)
  temperature_SCBI_2017 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_SCBI_2017, temperature_SCBI_2017) <- parameter(rank = 1)
  
  #SCBI2018
  
  initial(ELA_SCBI_2018) <- 0
  initial(E_SCBI_2018) <- 5000
  initial(HL_SCBI_2018) <- 0
  initial(RFL_SCBI_2018) <- 0
  initial(FL_SCBI_2018) <- 0
  initial(EL_SCBI_2018) <- 0
  initial(RFN_SCBI_2018) <- 0
  initial(FNabove_SCBI_2018) <- 0
  initial(FNbelow_SCBI_2018) <- 0
  initial(EN_SCBI_2018) <- 0
  initial(RFA_SCBI_2018) <- 0
  initial(FA_SCBI_2018) <- 0
  initial(EA_SCBI_2018) <- 0
  
  #equations
  deriv(ELA_SCBI_2018) <- x_SCBI_2018 * EA_SCBI_2018 - 1/y * ELA_SCBI_2018
  deriv(E_SCBI_2018) <- ELA_SCBI_2018 * (1 - (0.01 + (0.04 * log(1.01 + FA_SCBI_2018)/D))) * p - (q_SCBI_2018 + mu_e) * E_SCBI_2018
  deriv(HL_SCBI_2018) <- q_SCBI_2018 * E_SCBI_2018 - (1/z + mu_hl) * HL_SCBI_2018
  deriv(RFL_SCBI_2018) <- 1/z * HL_SCBI_2018 - (lambda_l * theta_i_SCBI_2018 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2018))) + mu_ql) * RFL_SCBI_2018
  deriv(FL_SCBI_2018) <- lambda_l * theta_i_SCBI_2018 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2018))) * RFL_SCBI_2018 - (1/r + (0.65 + (0.049 * log(1.01 + FL_SCBI_2018)/R))) * FL_SCBI_2018
  deriv(EL_SCBI_2018) <- 1/r * FL_SCBI_2018 - (s_SCBI_2018 + mu_el) * EL_SCBI_2018
  deriv(RFN_SCBI_2018) <- s_SCBI_2018 * EL_SCBI_2018 -
    (lambda_n * theta_n_SCBI_2018 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2018))) + lambda_n * theta_n_SCBI_2018 * eps_n + mu_qn) * RFN_SCBI_2018
  deriv(FNabove_SCBI_2018) <- lambda_n * theta_n_SCBI_2018 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2018))) * RFN_SCBI_2018 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2018 + FNbelow_SCBI_2018))/R))) * FNabove_SCBI_2018
  deriv(FNbelow_SCBI_2018) <- lambda_n * theta_n_SCBI_2018 * eps_n * RFN_SCBI_2018 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2018 + FNbelow_SCBI_2018))/R))) * FNbelow_SCBI_2018
  deriv(EN_SCBI_2018) <- 1/u * (FNabove_SCBI_2018 + FNbelow_SCBI_2018) - (v_SCBI_2018 + mu_en) * EN_SCBI_2018
  deriv(RFA_SCBI_2018) <- v_SCBI_2018 * EN_SCBI_2018 - (lambda_a * theta_a_SCBI_2018 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2018))) + mu_qa) * RFA_SCBI_2018
  deriv(FA_SCBI_2018) <- 0.5 * lambda_a * theta_a_SCBI_2018 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2018))) * RFA_SCBI_2018 - (1/w + (0.5 + (0.049 * log(1.01 + FA_SCBI_2018)/D))) * FA_SCBI_2018
  deriv(EA_SCBI_2018) <- 1/w * FA_SCBI_2018 - (x_SCBI_2018 + mu_ea) * EA_SCBI_2018
  
  theta_i_raw_SCBI_2018 <- if (temp_SCBI_2018 < 26.6) (-0.2912015 + 0.001641181 * temp_SCBI_2018^2) else (3.111262 - 0.08418802 * temp_SCBI_2018)
  theta_i_unscaled_SCBI_2018 <- if (theta_i_raw_SCBI_2018 > 0) theta_i_raw_SCBI_2018 else 0
  theta_i_SCBI_2018 <- if (temp_SCBI_2018 <= 0) 0 else theta_i_unscaled_SCBI_2018 / theta_i_norm
  #
  theta_n_raw_SCBI_2018 <- (-0.06003227 + 0.1140132 * temp_SCBI_2018 - 0.003410052 * temp_SCBI_2018^2)
  theta_n_SCBI_2018 <- if (theta_n_raw_SCBI_2018 > 0) theta_n_raw_SCBI_2018 / theta_n_norm else 0
  #
  theta_a_raw_SCBI_2018 <- (-0.6370415 + 0.3792334 * temp_SCBI_2018 - 0.02169659 * temp_SCBI_2018^2)
  theta_a_SCBI_2018 <- if (theta_a_raw_SCBI_2018 > 0) theta_a_raw_SCBI_2018 / theta_a_norm else 0
  
  q_SCBI_2018 <- if( temp_SCBI_2018 > 0) temp_SCBI_2018^2.27/34234 else 0
  #
  s_SCBI_2018 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2018 > 0) max((dln_min_SCBI/2000), temp_SCBI_2018^2.55/101181) else (dln_min_SCBI/2000)
  else 0
  #
  #
  v_SCBI_2018 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2018 > 0) temp_SCBI_2018^1.21/1596 else 0
  else 0
  #
  x_SCBI_2018 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2018 > 0) temp_SCBI_2018^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_SCBI_2018 <- interpolate(DOY_SCBI_2018, temperature_SCBI_2018, "spline")
  DOY_SCBI_2018 <- parameter(constant = TRUE)
  temperature_SCBI_2018 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_SCBI_2018, temperature_SCBI_2018) <- parameter(rank = 1)
  
  #SCBI2018
  
  initial(ELA_SCBI_2019) <- 0
  initial(E_SCBI_2019) <- 5000
  initial(HL_SCBI_2019) <- 0
  initial(RFL_SCBI_2019) <- 0
  initial(FL_SCBI_2019) <- 0
  initial(EL_SCBI_2019) <- 0
  initial(RFN_SCBI_2019) <- 0
  initial(FNabove_SCBI_2019) <- 0
  initial(FNbelow_SCBI_2019) <- 0
  initial(EN_SCBI_2019) <- 0
  initial(RFA_SCBI_2019) <- 0
  initial(FA_SCBI_2019) <- 0
  initial(EA_SCBI_2019) <- 0
  
  #equations
  deriv(ELA_SCBI_2019) <- x_SCBI_2019 * EA_SCBI_2019 - 1/y * ELA_SCBI_2019
  deriv(E_SCBI_2019) <- ELA_SCBI_2019 * (1 - (0.01 + (0.04 * log(1.01 + FA_SCBI_2019)/D))) * p - (q_SCBI_2019 + mu_e) * E_SCBI_2019
  deriv(HL_SCBI_2019) <- q_SCBI_2019 * E_SCBI_2019 - (1/z + mu_hl) * HL_SCBI_2019
  deriv(RFL_SCBI_2019) <- 1/z * HL_SCBI_2019 - (lambda_l * theta_i_SCBI_2019 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2019))) + mu_ql) * RFL_SCBI_2019
  deriv(FL_SCBI_2019) <- lambda_l * theta_i_SCBI_2019 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2019))) * RFL_SCBI_2019 - (1/r + (0.65 + (0.049 * log(1.01 + FL_SCBI_2019)/R))) * FL_SCBI_2019
  deriv(EL_SCBI_2019) <- 1/r * FL_SCBI_2019 - (s_SCBI_2019 + mu_el) * EL_SCBI_2019
  deriv(RFN_SCBI_2019) <- s_SCBI_2019 * EL_SCBI_2019 -
    (lambda_n * theta_n_SCBI_2019 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2019))) + lambda_n * theta_n_SCBI_2019 * eps_n + mu_qn) * RFN_SCBI_2019
  deriv(FNabove_SCBI_2019) <- lambda_n * theta_n_SCBI_2019 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2019))) * RFN_SCBI_2019 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2019 + FNbelow_SCBI_2019))/R))) * FNabove_SCBI_2019
  deriv(FNbelow_SCBI_2019) <- lambda_n * theta_n_SCBI_2019 * eps_n * RFN_SCBI_2019 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2019 + FNbelow_SCBI_2019))/R))) * FNbelow_SCBI_2019
  deriv(EN_SCBI_2019) <- 1/u * (FNabove_SCBI_2019 + FNbelow_SCBI_2019) - (v_SCBI_2019 + mu_en) * EN_SCBI_2019
  deriv(RFA_SCBI_2019) <- v_SCBI_2019 * EN_SCBI_2019 - (lambda_a * theta_a_SCBI_2019 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2019))) + mu_qa) * RFA_SCBI_2019
  deriv(FA_SCBI_2019) <- 0.5 * lambda_a * theta_a_SCBI_2019 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2019))) * RFA_SCBI_2019 - (1/w + (0.5 + (0.049 * log(1.01 + FA_SCBI_2019)/D))) * FA_SCBI_2019
  deriv(EA_SCBI_2019) <- 1/w * FA_SCBI_2019 - (x_SCBI_2019 + mu_ea) * EA_SCBI_2019
  
  theta_i_raw_SCBI_2019 <- if (temp_SCBI_2019 < 26.6) (-0.2912015 + 0.001641181 * temp_SCBI_2019^2) else (3.111262 - 0.08418802 * temp_SCBI_2019)
  theta_i_unscaled_SCBI_2019 <- if (theta_i_raw_SCBI_2019 > 0) theta_i_raw_SCBI_2019 else 0
  theta_i_SCBI_2019 <- if (temp_SCBI_2019 <= 0) 0 else theta_i_unscaled_SCBI_2019 / theta_i_norm
  #
  theta_n_raw_SCBI_2019 <- (-0.06003227 + 0.1140132 * temp_SCBI_2019 - 0.003410052 * temp_SCBI_2019^2)
  theta_n_SCBI_2019 <- if (theta_n_raw_SCBI_2019 > 0) theta_n_raw_SCBI_2019 / theta_n_norm else 0
  #
  theta_a_raw_SCBI_2019 <- (-0.6370415 + 0.3792334 * temp_SCBI_2019 - 0.02169659 * temp_SCBI_2019^2)
  theta_a_SCBI_2019 <- if (theta_a_raw_SCBI_2019 > 0) theta_a_raw_SCBI_2019 / theta_a_norm else 0
  
  q_SCBI_2019 <- if( temp_SCBI_2019 > 0) temp_SCBI_2019^2.27/34234 else 0
  #
  s_SCBI_2019 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2019 > 0) max((dln_min_SCBI/2000), temp_SCBI_2019^2.55/101181) else (dln_min_SCBI/2000)
  else 0
  #
  #
  v_SCBI_2019 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2019 > 0) temp_SCBI_2019^1.21/1596 else 0
  else 0
  #
  x_SCBI_2019 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2019 > 0) temp_SCBI_2019^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_SCBI_2019 <- interpolate(DOY_SCBI_2019, temperature_SCBI_2019, "spline")
  DOY_SCBI_2019 <- parameter(constant = TRUE)
  temperature_SCBI_2019 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_SCBI_2019, temperature_SCBI_2019) <- parameter(rank = 1)
  
  
  #SCBI2023
  
  initial(ELA_SCBI_2023) <- 0
  initial(E_SCBI_2023) <- 5000
  initial(HL_SCBI_2023) <- 0
  initial(RFL_SCBI_2023) <- 0
  initial(FL_SCBI_2023) <- 0
  initial(EL_SCBI_2023) <- 0
  initial(RFN_SCBI_2023) <- 0
  initial(FNabove_SCBI_2023) <- 0
  initial(FNbelow_SCBI_2023) <- 0
  initial(EN_SCBI_2023) <- 0
  initial(RFA_SCBI_2023) <- 0
  initial(FA_SCBI_2023) <- 0
  initial(EA_SCBI_2023) <- 0
  
  #equations
  deriv(ELA_SCBI_2023) <- x_SCBI_2023 * EA_SCBI_2023 - 1/y * ELA_SCBI_2023
  deriv(E_SCBI_2023) <- ELA_SCBI_2023 * (1 - (0.01 + (0.04 * log(1.01 + FA_SCBI_2023)/D))) * p - (q_SCBI_2023 + mu_e) * E_SCBI_2023
  deriv(HL_SCBI_2023) <- q_SCBI_2023 * E_SCBI_2023 - (1/z + mu_hl) * HL_SCBI_2023
  deriv(RFL_SCBI_2023) <- 1/z * HL_SCBI_2023 - (lambda_l * theta_i_SCBI_2023 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2023))) + mu_ql) * RFL_SCBI_2023
  deriv(FL_SCBI_2023) <- lambda_l * theta_i_SCBI_2023 * 1 / (1 + exp(-(-(a_n_SCBI+(a_l_SCBI)) + b_l * humidity_SCBI_2023))) * RFL_SCBI_2023 - (1/r + (0.65 + (0.049 * log(1.01 + FL_SCBI_2023)/R))) * FL_SCBI_2023
  deriv(EL_SCBI_2023) <- 1/r * FL_SCBI_2023 - (s_SCBI_2023 + mu_el) * EL_SCBI_2023
  deriv(RFN_SCBI_2023) <- s_SCBI_2023 * EL_SCBI_2023 -
    (lambda_n * theta_n_SCBI_2023 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2023))) + lambda_n * theta_n_SCBI_2023 * eps_n + mu_qn) * RFN_SCBI_2023
  deriv(FNabove_SCBI_2023) <- lambda_n * theta_n_SCBI_2023 * 1 / (1 + exp(-(-a_n_SCBI + b_n * humidity_SCBI_2023))) * RFN_SCBI_2023 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2023 + FNbelow_SCBI_2023))/R))) * FNabove_SCBI_2023
  deriv(FNbelow_SCBI_2023) <- lambda_n * theta_n_SCBI_2023 * eps_n * RFN_SCBI_2023 -
    (1/u + (0.55 + (0.049 * log(1.01 + (FNabove_SCBI_2023 + FNbelow_SCBI_2023))/R))) * FNbelow_SCBI_2023
  deriv(EN_SCBI_2023) <- 1/u * (FNabove_SCBI_2023 + FNbelow_SCBI_2023) - (v_SCBI_2023 + mu_en) * EN_SCBI_2023
  deriv(RFA_SCBI_2023) <- v_SCBI_2023 * EN_SCBI_2023 - (lambda_a * theta_a_SCBI_2023 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2023))) + mu_qa) * RFA_SCBI_2023
  deriv(FA_SCBI_2023) <- 0.5 * lambda_a * theta_a_SCBI_2023 *  1 / (1 + exp(-(-(a_n_SCBI-(a_a_SCBI)) + b_a * humidity_SCBI_2023))) * RFA_SCBI_2023 - (1/w + (0.5 + (0.049 * log(1.01 + FA_SCBI_2023)/D))) * FA_SCBI_2023
  deriv(EA_SCBI_2023) <- 1/w * FA_SCBI_2023 - (x_SCBI_2023 + mu_ea) * EA_SCBI_2023
  
  theta_i_raw_SCBI_2023 <- if (temp_SCBI_2023 < 26.6) (-0.2912015 + 0.001641181 * temp_SCBI_2023^2) else (3.111262 - 0.08418802 * temp_SCBI_2023)
  theta_i_unscaled_SCBI_2023 <- if (theta_i_raw_SCBI_2023 > 0) theta_i_raw_SCBI_2023 else 0
  theta_i_SCBI_2023 <- if (temp_SCBI_2023 <= 0) 0 else theta_i_unscaled_SCBI_2023 / theta_i_norm
  #
  theta_n_raw_SCBI_2023 <- (-0.06003227 + 0.1140132 * temp_SCBI_2023 - 0.003410052 * temp_SCBI_2023^2)
  theta_n_SCBI_2023 <- if (theta_n_raw_SCBI_2023 > 0) theta_n_raw_SCBI_2023 / theta_n_norm else 0
  #
  theta_a_raw_SCBI_2023 <- (-0.6370415 + 0.3792334 * temp_SCBI_2023 - 0.02169659 * temp_SCBI_2023^2)
  theta_a_SCBI_2023 <- if (theta_a_raw_SCBI_2023 > 0) theta_a_raw_SCBI_2023 / theta_a_norm else 0
  
  q_SCBI_2023 <- if( temp_SCBI_2023 > 0) temp_SCBI_2023^2.27/34234 else 0
  #
  s_SCBI_2023 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2023 > 0) max((dln_min_SCBI/2000), temp_SCBI_2023^2.55/101181) else (dln_min_SCBI/2000)
  else 0
  #
  #
  v_SCBI_2023 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2023 > 0) temp_SCBI_2023^1.21/1596 else 0
  else 0
  #
  x_SCBI_2023 <- if (time - floor(time / 365) * 365 < 181)
    if(temp_SCBI_2023 > 0) temp_SCBI_2023^1.42/1300 else 0
  else 0
  
  ##interpolated functions
  temp_SCBI_2023 <- interpolate(DOY_SCBI_2023, temperature_SCBI_2023, "spline")
  DOY_SCBI_2023 <- parameter(constant = TRUE)
  temperature_SCBI_2023 <- parameter(constant = TRUE)
  
  #dimesions
  dim(DOY_SCBI_2023, temperature_SCBI_2023) <- parameter(rank = 1)
  
})



## create file to save results 
B_L_16_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_B_L16_preds.csv",sep="")
B_L_17_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_B_L17_preds.csv",sep="")
B_L_18_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_B_L18_preds.csv",sep="")
B_L_19_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_B_L19_preds.csv",sep="")


B_N_16_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_B_N16_preds.csv",sep="")
B_N_17_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_B_N17_preds.csv",sep="")
B_N_18_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_B_N18_preds.csv",sep="")
B_N_19_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_B_N19_preds.csv",sep="")



H_L_17_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_H_L17_preds.csv",sep="")
H_L_18_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_H_L18_preds.csv",sep="")
H_L_21_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_H_L21_preds.csv",sep="")
H_L_23_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_H_L23_preds.csv",sep="")



H_N_17_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_H_N17_preds.csv",sep="")
H_N_18_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_H_N18_preds.csv",sep="")
H_N_21_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_H_N21_preds.csv",sep="")
H_N_23_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_H_N23_preds.csv",sep="")


S_L_16_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_L16_preds.csv",sep="")
S_L_17_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_L17_preds.csv",sep="")
S_L_18_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_L18_preds.csv",sep="")
S_L_19_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_L19_preds.csv",sep="")
S_L_23_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_L23_preds.csv",sep="")


S_N_16_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_N16_preds.csv",sep="")
S_N_17_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_N17_preds.csv",sep="")
S_N_18_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_N18_preds.csv",sep="")
S_N_19_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_N19_preds.csv",sep="")
S_N_23_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_S_N23_preds.csv",sep="")


all_overlaps_file<-paste("~/Desktop/CH2_FEB26/posterioreffect_unc_overlap.csv",sep="")

B_L_16_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
B_L_17_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
B_L_18_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
B_L_19_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)

B_N_16_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
B_N_17_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
B_N_18_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
B_N_19_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)



H_L_17_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
H_L_18_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
H_L_21_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
H_L_23_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)




H_N_17_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
H_N_18_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
H_N_21_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
H_N_23_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)


S_L_16_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
S_L_17_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
S_L_18_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
S_L_19_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
S_L_23_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)



S_N_16_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
S_N_17_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
S_N_18_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
S_N_19_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)
S_N_23_mat<-matrix(0, nrow=365,ncol=n_posterior_samples)


all_overlaps<-matrix(0, nrow=13,ncol=n_posterior_samples)


for (samp_idx in 1: n_posterior_samples)
{
  print(samp_idx)
  full_pars <- list( theta_i_norm = theta_i_norm,
                     theta_n_norm = theta_n_norm,
                     theta_a_norm = theta_a_norm,
                     humid_int_BLAN_2016 = BLAN_humid_int_2016,
                     humid_amp_BLAN_2016 = BLAN_humid_amp_2016,
                     humid_ps_BLAN_2016 = BLAN_phaseshift_2016,
                     humid_int_BLAN_2017 = BLAN_humid_int_2017,
                     humid_amp_BLAN_2017 = BLAN_humid_amp_2017,
                     humid_ps_BLAN_2017 = BLAN_phaseshift_2017,
                     humid_int_BLAN_2018 = BLAN_humid_int_2018,
                     humid_amp_BLAN_2018 = BLAN_humid_amp_2018,
                     humid_ps_BLAN_2018 = BLAN_phaseshift_2018,
                     
                     humid_int_BLAN_2019 = BLAN_humid_int_2019,
                     humid_amp_BLAN_2019 = BLAN_humid_amp_2019,
                     humid_ps_BLAN_2019 = BLAN_phaseshift_2019,
                     
                     
                     humid_int_HARV_2017 = HARV_humid_int_2017,
                     humid_amp_HARV_2017 = HARV_humid_amp_2017,
                     humid_ps_HARV_2017 = HARV_phaseshift_2017,
                     humid_int_HARV_2018 = HARV_humid_int_2018,
                     humid_amp_HARV_2018 = HARV_humid_amp_2018,
                     humid_ps_HARV_2018 = HARV_phaseshift_2018,
                     humid_int_HARV_2021 = HARV_humid_int_2021,
                     humid_amp_HARV_2021 = HARV_humid_amp_2021,
                     humid_ps_HARV_2021 = HARV_phaseshift_2021,
                     humid_int_HARV_2023 = HARV_humid_int_2023,
                     humid_amp_HARV_2023 = HARV_humid_amp_2023,
                     humid_ps_HARV_2023 = HARV_phaseshift_2023,
                     
                     
                     humid_int_SCBI_2016 = SCBI_humid_int_2016,
                     humid_amp_SCBI_2016 = SCBI_humid_amp_2016,
                     humid_ps_SCBI_2016 = SCBI_phaseshift_2016,
                     humid_int_SCBI_2017 = SCBI_humid_int_2017,
                     humid_amp_SCBI_2017 = SCBI_humid_amp_2017,
                     humid_ps_SCBI_2017 = SCBI_phaseshift_2017,
                     humid_int_SCBI_2018 = SCBI_humid_int_2018,
                     humid_amp_SCBI_2018 = SCBI_humid_amp_2018,
                     humid_ps_SCBI_2018 = SCBI_phaseshift_2018,
                     humid_int_SCBI_2019 = SCBI_humid_int_2019,
                     humid_amp_SCBI_2019 = SCBI_humid_amp_2019,
                     humid_ps_SCBI_2019 = SCBI_phaseshift_2019,
                     humid_int_SCBI_2023 = SCBI_humid_int_2023,
                     humid_amp_SCBI_2023 = SCBI_humid_amp_2023,
                     humid_ps_SCBI_2023 = SCBI_phaseshift_2023,
                     
                     
                     DOY_BLAN_2016 = DOY_BLAN_2016,
                     temperature_BLAN_2016 = temperature_BLAN_2016,
                     DOY_BLAN_2017 = DOY_BLAN_2017,
                     temperature_BLAN_2017 = temperature_BLAN_2017,
                     DOY_BLAN_2018 = DOY_BLAN_2018,
                     temperature_BLAN_2018 = temperature_BLAN_2018,
                     DOY_BLAN_2019 = DOY_BLAN_2019,
                     temperature_BLAN_2019 = temperature_BLAN_2019,
                     
                     
                     
                     DOY_HARV_2017 = DOY_HARV_2017,
                     temperature_HARV_2017 = temperature_HARV_2017,
                     DOY_HARV_2018 = DOY_HARV_2018,
                     temperature_HARV_2018 = temperature_HARV_2018,
                     
                     DOY_HARV_2021 = DOY_HARV_2021,
                     temperature_HARV_2021 = temperature_HARV_2021,
                     
                     DOY_HARV_2023 = DOY_HARV_2023,
                     temperature_HARV_2023 = temperature_HARV_2023,
                     
                     DOY_SCBI_2016 = DOY_SCBI_2016,
                     temperature_SCBI_2016 = temperature_SCBI_2016,
                     DOY_SCBI_2017 = DOY_SCBI_2017,
                     temperature_SCBI_2017 = temperature_SCBI_2017,
                     DOY_SCBI_2018 = DOY_SCBI_2018,
                     temperature_SCBI_2018 = temperature_SCBI_2018,
                     
                     DOY_SCBI_2019 = DOY_SCBI_2019,
                     temperature_SCBI_2019 = temperature_SCBI_2019,
                     DOY_SCBI_2021 = DOY_SCBI_2021,
                     
                     DOY_SCBI_2023 = DOY_SCBI_2023,
                     temperature_SCBI_2023 = temperature_SCBI_2023,
                     
                     
                     
                     log_dln_BLAN = posterior_pars[samp_idx,"log_dln_BLAN"],
                     log_dln_SCBI = posterior_pars[samp_idx,"log_dln_SCBI"],
                     log_dln_HARV = posterior_pars[samp_idx,"log_dln_HARV"],
                     log_aa_BLAN = posterior_pars[samp_idx,"log_aa_BLAN"],
                     log_aa_HARV = posterior_pars[samp_idx,"log_aa_HARV"],
                     log_aa_SCBI = posterior_pars[samp_idx,"log_aa_SCBI"],
                     log_al_BLAN = posterior_pars[samp_idx,"log_al_BLAN"],
                     log_al_HARV = posterior_pars[samp_idx,"log_al_HARV"],
                     log_al_SCBI = posterior_pars[samp_idx,"log_al_SCBI"],

                     a_n_BLAN = posterior_pars[samp_idx,"a_n_BLAN"],
                     a_n_HARV = posterior_pars[samp_idx,"a_n_HARV"],
                     a_n_SCBI = posterior_pars[samp_idx,"a_n_SCBI"]
                     
                     # log_dln_BLAN = median(posterior_pars[,"log_dln_BLAN"]),
                     # log_dln_SCBI = median(posterior_pars[,"log_dln_SCBI"]),
                     # log_dln_HARV = median(posterior_pars[,"log_dln_HARV"]),
                     # log_aa_BLAN = median(posterior_pars[,"log_aa_BLAN"]),
                     # log_aa_HARV = median(posterior_pars[,"log_aa_HARV"]),
                     # log_aa_SCBI = median(posterior_pars[,"log_aa_SCBI"]),
                     # log_al_BLAN = median(posterior_pars[,"log_al_BLAN"]),
                     # log_al_HARV = median(posterior_pars[,"log_al_HARV"]),
                     # log_al_SCBI = median(posterior_pars[,"log_al_SCBI"]),
                     # 
                     # a_n_BLAN = median(posterior_pars[,"a_n_BLAN"]),
                     # a_n_HARV = median(posterior_pars[,"a_n_HARV"]),
                     # a_n_SCBI = median(posterior_pars[,"a_n_SCBI"])
                     
                     
  )
  
  
  ode_ctrl <- dust_ode_control(
    max_steps = 10000000,
    step_size_min = .000000001
  )
  
  full_sys <- dust_system_create(full_odin(), 
                                 pars = full_pars,
                                 n_particles = 1,
                                 ode_control = ode_ctrl)
  
  
  
  
  ##initializes system
  dust_system_set_state_initial(full_sys)
  
  
  ##retrieves the current state vector of the model system
  full_sys_s<-dust_system_state(full_sys)
  
  
  ##simulate model over t time 
  full_sys_y <- dust_system_simulate(full_sys, t)
  
  ##save outputs
  
  
  L_BLAN_2016_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_BLAN_2016[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_BLAN_2016[i:(i+364)])) * sample_window_BLAN_2016
  N_BLAN_2016_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_BLAN_2016[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_BLAN_2016[i:(i+364)])) * sample_window_BLAN_2016
  
  L_BLAN_2017_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_BLAN_2017[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_BLAN_2017[i:(i+364)])) * sample_window_BLAN_2017
  N_BLAN_2017_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_BLAN_2017[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_BLAN_2017[i:(i+364)])) * sample_window_BLAN_2017
  
  L_BLAN_2018_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_BLAN_2018[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_BLAN_2018[i:(i+364)])) * sample_window_BLAN_2018
  N_BLAN_2018_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_BLAN_2018[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_BLAN_2018[i:(i+364)])) * sample_window_BLAN_2018
  
  L_BLAN_2019_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_BLAN_2019[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_BLAN_2019[i:(i+364)])) * sample_window_BLAN_2019
  N_BLAN_2019_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_BLAN_2019[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_BLAN_2019[i:(i+364)])) * sample_window_BLAN_2019
  
  
  
  L_HARV_2017_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_HARV_2017[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_HARV_2017[i:(i+364)])) * sample_window_HARV_2017
  N_HARV_2017_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_HARV_2017[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_HARV_2017[i:(i+364)])) * sample_window_HARV_2017
  
  L_HARV_2018_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_HARV_2018[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_HARV_2018[i:(i+364)])) * sample_window_HARV_2018
  N_HARV_2018_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_HARV_2018[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_HARV_2018[i:(i+364)])) * sample_window_HARV_2018
  
  L_HARV_2021_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_HARV_2021[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_HARV_2021[i:(i+364)])) * sample_window_HARV_2021
  N_HARV_2021_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_HARV_2021[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_HARV_2021[i:(i+364)])) * sample_window_HARV_2021
  
  L_HARV_2023_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_HARV_2023[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_HARV_2023[i:(i+364)])) * sample_window_HARV_2023
  N_HARV_2023_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_HARV_2023[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_HARV_2023[i:(i+364)])) * sample_window_HARV_2023
  
  
  L_SCBI_2016_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2016[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2016[i:(i+364)])) * sample_window_SCBI_2016
  N_SCBI_2016_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2016[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2016[i:(i+364)])) * sample_window_SCBI_2016
  
  L_SCBI_2017_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2017[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2017[i:(i+364)])) * sample_window_SCBI_2017
  N_SCBI_2017_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2017[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2017[i:(i+364)])) * sample_window_SCBI_2017
  
  L_SCBI_2018_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2018[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2018[i:(i+364)])) * sample_window_SCBI_2018
  N_SCBI_2018_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2018[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2018[i:(i+364)])) * sample_window_SCBI_2018
  
  L_SCBI_2019_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2019[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2019[i:(i+364)])) * sample_window_SCBI_2019
  N_SCBI_2019_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2019[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2019[i:(i+364)])) * sample_window_SCBI_2019
  
  L_SCBI_2023_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2023[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FL_SCBI_2023[i:(i+364)])) * sample_window_SCBI_2023
  N_SCBI_2023_out<-as.vector(dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2023[i:(i+364)]/AUC(t[i:(i+364)],dust_unpack_state(full_sys, full_sys_y)$FNabove_SCBI_2023[i:(i+364)])) * sample_window_SCBI_2023
  
  
  
  ##normalize based on sampling window
  L_BLAN_2016_out<-L_BLAN_2016_out/sum(L_BLAN_2016_out)
  N_BLAN_2016_out<-N_BLAN_2016_out/sum(N_BLAN_2016_out)
  
  L_BLAN_2017_out<-L_BLAN_2017_out/sum(L_BLAN_2017_out)
  N_BLAN_2017_out<-N_BLAN_2017_out/sum(N_BLAN_2017_out)
  
  L_BLAN_2018_out<-L_BLAN_2018_out/sum(L_BLAN_2018_out)
  N_BLAN_2018_out<-N_BLAN_2018_out/sum(N_BLAN_2018_out)
  
  L_BLAN_2019_out<-L_BLAN_2019_out/sum(L_BLAN_2019_out)
  N_BLAN_2019_out<-N_BLAN_2019_out/sum(N_BLAN_2019_out)
  
  
  
  
  L_HARV_2017_out<-L_HARV_2017_out/sum(L_HARV_2017_out)
  N_HARV_2017_out<-N_HARV_2017_out/sum(N_HARV_2017_out)
  
  L_HARV_2018_out<-L_HARV_2018_out/sum(L_HARV_2018_out)
  N_HARV_2018_out<-N_HARV_2018_out/sum(N_HARV_2018_out)
  
  L_HARV_2021_out<-L_HARV_2021_out/sum(L_HARV_2021_out)
  N_HARV_2021_out<-N_HARV_2021_out/sum(N_HARV_2021_out)
  
  L_HARV_2023_out<-L_HARV_2023_out/sum(L_HARV_2023_out)
  N_HARV_2023_out<-N_HARV_2023_out/sum(N_HARV_2023_out)
  
  
  
  L_SCBI_2016_out<-L_SCBI_2016_out/sum(L_SCBI_2016_out)
  N_SCBI_2016_out<-N_SCBI_2016_out/sum(N_SCBI_2016_out)
  
  L_SCBI_2017_out<-L_SCBI_2017_out/sum(L_SCBI_2017_out)
  N_SCBI_2017_out<-N_SCBI_2017_out/sum(N_SCBI_2017_out)
  
  L_SCBI_2018_out<-L_SCBI_2018_out/sum(L_SCBI_2018_out)
  N_SCBI_2018_out<-N_SCBI_2018_out/sum(N_SCBI_2018_out)
  
  L_SCBI_2019_out<-L_SCBI_2019_out/sum(L_SCBI_2019_out)
  N_SCBI_2019_out<-N_SCBI_2019_out/sum(N_SCBI_2019_out)
  
  L_SCBI_2023_out<-L_SCBI_2023_out/sum(L_SCBI_2023_out)
  N_SCBI_2023_out<-N_SCBI_2023_out/sum(N_SCBI_2023_out)
  
  
  
  
  
  
  
  ##calculate overlap based on GAMS
  time <- c(i:(i+364))
  new_data_BLAN_2016 <- data.frame(time = c(i:(i+364)),keep=sample_window_BLAN_2016)
  new_data_BLAN_2016<-subset(new_data_BLAN_2016,new_data_BLAN_2016$keep==1)
  new_data_BLAN_2016$keep <- NULL
  
  new_data_BLAN_2017 <- data.frame(time = c(i:(i+364)),keep=sample_window_BLAN_2017)
  new_data_BLAN_2017<-subset(new_data_BLAN_2017,new_data_BLAN_2017$keep==1)
  new_data_BLAN_2017$keep <- NULL
  
  new_data_BLAN_2018 <- data.frame(time = c(i:(i+364)),keep=sample_window_BLAN_2018)
  new_data_BLAN_2018<-subset(new_data_BLAN_2018,new_data_BLAN_2018$keep==1)
  new_data_BLAN_2018$keep <- NULL
  
  new_data_BLAN_2019 <- data.frame(time = c(i:(i+364)),keep=sample_window_BLAN_2019)
  new_data_BLAN_2019<-subset(new_data_BLAN_2019,new_data_BLAN_2019$keep==1)
  new_data_BLAN_2019$keep <- NULL
  
  
  
  
  new_data_HARV_2017 <- data.frame(time = c(i:(i+364)),keep=sample_window_HARV_2017)
  new_data_HARV_2017<-subset(new_data_HARV_2017,new_data_HARV_2017$keep==1)
  new_data_HARV_2017$keep <- NULL
  
  new_data_HARV_2018 <- data.frame(time = c(i:(i+364)),keep=sample_window_HARV_2018)
  new_data_HARV_2018<-subset(new_data_HARV_2018,new_data_HARV_2018$keep==1)
  new_data_HARV_2018$keep <- NULL
  
  
  new_data_HARV_2021 <- data.frame(time = c(i:(i+364)),keep=sample_window_HARV_2021)
  new_data_HARV_2021<-subset(new_data_HARV_2021,new_data_HARV_2021$keep==1)
  new_data_HARV_2021$keep <- NULL
  
  
  new_data_HARV_2023 <- data.frame(time = c(i:(i+364)),keep=sample_window_HARV_2023)
  new_data_HARV_2023<-subset(new_data_HARV_2023,new_data_HARV_2023$keep==1)
  new_data_HARV_2023$keep <- NULL
  
  new_data_SCBI_2016 <- data.frame(time = c(i:(i+364)),keep=sample_window_SCBI_2016)
  new_data_SCBI_2016<-subset(new_data_SCBI_2016,new_data_SCBI_2016$keep==1)
  new_data_SCBI_2016$keep <- NULL
  
  new_data_SCBI_2017 <- data.frame(time = c(i:(i+364)),keep=sample_window_SCBI_2017)
  new_data_SCBI_2017<-subset(new_data_SCBI_2017,new_data_SCBI_2017$keep==1)
  new_data_SCBI_2017$keep <- NULL
  
  new_data_SCBI_2018 <- data.frame(time = c(i:(i+364)),keep=sample_window_SCBI_2018)
  new_data_SCBI_2018<-subset(new_data_SCBI_2018,new_data_SCBI_2018$keep==1)
  new_data_SCBI_2018$keep <- NULL
  
  new_data_SCBI_2019 <- data.frame(time = c(i:(i+364)),keep=sample_window_SCBI_2019)
  new_data_SCBI_2019<-subset(new_data_SCBI_2019,new_data_SCBI_2019$keep==1)
  new_data_SCBI_2019$keep <- NULL
  
  
  new_data_SCBI_2023 <- data.frame(time = c(i:(i+364)),keep=sample_window_SCBI_2023)
  new_data_SCBI_2023<-subset(new_data_SCBI_2023,new_data_SCBI_2023$keep==1)
  new_data_SCBI_2023$keep <- NULL
  
  
  
  
  
  
  
  gam_rate = 0.5
  
  yl_BLAN_2016 <- log(L_BLAN_2016_out + .0001)
  yn_BLAN_2016 <- log(N_BLAN_2016_out + .0001)
  
  gam_model_L_BLAN_2016 <- gam(yl_BLAN_2016 ~ s(time,k=29), gamma = gam_rate)
  gam_model_N_BLAN_2016 <- gam(yn_BLAN_2016  ~ s(time,k=29), gamma = gam_rate)
  
  
  preds_L_BLAN_2016 <- exp(predict(gam_model_L_BLAN_2016,new_data_BLAN_2016 ))
  preds_N_BLAN_2016 <- exp(predict(gam_model_N_BLAN_2016,new_data_BLAN_2016 ))
  
  preds_L_BLAN_2016 <- exp(predict(gam_model_L_BLAN_2016,new_data_BLAN_2016 ))/sum(preds_L_BLAN_2016)
  preds_N_BLAN_2016 <- exp(predict(gam_model_N_BLAN_2016,new_data_BLAN_2016 ))/sum(preds_N_BLAN_2016)
  
  f_BLAN_2016 <- approxfun(new_data_BLAN_2016$time, pmin(preds_L_BLAN_2016,preds_N_BLAN_2016))
  synch_BLAN_2016<-integrate( f_BLAN_2016, min(new_data_BLAN_2016$time), max(new_data_BLAN_2016$time),subdivisions = 1000)
  overlap_BLAN_2016<-synch_BLAN_2016$value/(AUC(new_data_BLAN_2016$time,preds_L_BLAN_2016)+AUC(new_data_BLAN_2016$time,preds_N_BLAN_2016)-synch_BLAN_2016$value)
  
  yl_BLAN_2017 <- log(L_BLAN_2017_out + .0001)
  yn_BLAN_2017 <- log(N_BLAN_2017_out + .0001)
  
  gam_model_L_BLAN_2017 <- gam(yl_BLAN_2017 ~ s(time,k=29), gamma = 0.55)
  gam_model_N_BLAN_2017 <- gam(yn_BLAN_2017  ~ s(time,k=29), gamma = 0.55)
  
  preds_L_BLAN_2017 <- exp(predict(gam_model_L_BLAN_2017,new_data_BLAN_2017 ))
  preds_N_BLAN_2017 <- exp(predict(gam_model_N_BLAN_2017,new_data_BLAN_2017 ))
  
  preds_L_BLAN_2017 <- exp(predict(gam_model_L_BLAN_2017,new_data_BLAN_2017 ))/sum(preds_L_BLAN_2017)
  preds_N_BLAN_2017 <- exp(predict(gam_model_N_BLAN_2017,new_data_BLAN_2017 ))/sum(preds_N_BLAN_2017)
  
  f_BLAN_2017 <- approxfun(new_data_BLAN_2017$time, pmin(preds_L_BLAN_2017,preds_N_BLAN_2017))
  synch_BLAN_2017<-integrate( f_BLAN_2017, min(new_data_BLAN_2017$time), max(new_data_BLAN_2017$time),subdivisions = 1000)
  overlap_BLAN_2017<-synch_BLAN_2017$value/(AUC(new_data_BLAN_2017$time,preds_L_BLAN_2017)+AUC(new_data_BLAN_2017$time,preds_N_BLAN_2017)-synch_BLAN_2017$value)
  
  yl_BLAN_2018 <- log(L_BLAN_2018_out + .0001)
  yn_BLAN_2018 <- log(N_BLAN_2018_out + .0001)
  
  gam_model_L_BLAN_2018 <- gam(yl_BLAN_2018 ~ s(time,k=29), gamma = gam_rate)
  gam_model_N_BLAN_2018 <- gam(yn_BLAN_2018  ~ s(time,k=29), gamma = gam_rate)
  
  preds_L_BLAN_2018 <- exp(predict(gam_model_L_BLAN_2018,new_data_BLAN_2018 ))
  preds_N_BLAN_2018 <- exp(predict(gam_model_N_BLAN_2018,new_data_BLAN_2018 ))
  
  
  preds_L_BLAN_2018 <- exp(predict(gam_model_L_BLAN_2018,new_data_BLAN_2018 ))/sum(preds_L_BLAN_2018)
  preds_N_BLAN_2018 <- exp(predict(gam_model_N_BLAN_2018,new_data_BLAN_2018 ))/sum(preds_N_BLAN_2018)
  
  f_BLAN_2018 <- approxfun(new_data_BLAN_2018$time, pmin(preds_L_BLAN_2018,preds_N_BLAN_2018))
  synch_BLAN_2018<-integrate( f_BLAN_2018, min(new_data_BLAN_2018$time), max(new_data_BLAN_2018$time),subdivisions = 1000)
  overlap_BLAN_2018<-synch_BLAN_2018$value/(AUC(new_data_BLAN_2018$time,preds_L_BLAN_2018)+AUC(new_data_BLAN_2018$time,preds_N_BLAN_2018)-synch_BLAN_2018$value)
  
  yl_BLAN_2019 <- log(L_BLAN_2019_out + .0001)
  yn_BLAN_2019 <- log(N_BLAN_2019_out + .0001)
  
  gam_model_L_BLAN_2019 <- gam(yl_BLAN_2019 ~ s(time,k=29), gamma = gam_rate)
  gam_model_N_BLAN_2019 <- gam(yn_BLAN_2019  ~ s(time,k=29), gamma = gam_rate)
  
  preds_L_BLAN_2019 <- exp(predict(gam_model_L_BLAN_2019,new_data_BLAN_2019 ))
  preds_N_BLAN_2019 <- exp(predict(gam_model_N_BLAN_2019,new_data_BLAN_2019 ))
  
  
  preds_L_BLAN_2019 <- exp(predict(gam_model_L_BLAN_2019,new_data_BLAN_2019 ))/sum(preds_L_BLAN_2019)
  preds_N_BLAN_2019 <- exp(predict(gam_model_N_BLAN_2019,new_data_BLAN_2019 ))/sum(preds_N_BLAN_2019)
  
  f_BLAN_2019 <- approxfun(new_data_BLAN_2019$time, pmin(preds_L_BLAN_2019,preds_N_BLAN_2019))
  synch_BLAN_2019<-integrate( f_BLAN_2019, min(new_data_BLAN_2019$time), max(new_data_BLAN_2019$time),subdivisions = 1000)
  overlap_BLAN_2019<-synch_BLAN_2019$value/(AUC(new_data_BLAN_2019$time,preds_L_BLAN_2019)+AUC(new_data_BLAN_2019$time,preds_N_BLAN_2019)-synch_BLAN_2019$value)
  
  
  
  yl_HARV_2017 <- log(L_HARV_2017_out + .0001)
  yn_HARV_2017 <- log(N_HARV_2017_out + .0001)
  
  gam_model_L_HARV_2017 <- gam(yl_HARV_2017 ~ s(time,k=29), gamma = 0.55)
  gam_model_N_HARV_2017 <- gam(yn_HARV_2017  ~ s(time,k=29), gamma = 0.55)
  
  preds_L_HARV_2017 <- exp(predict(gam_model_L_HARV_2017,new_data_HARV_2017 ))
  preds_N_HARV_2017 <- exp(predict(gam_model_N_HARV_2017,new_data_HARV_2017 ))
  
  preds_L_HARV_2017 <- exp(predict(gam_model_L_HARV_2017,new_data_HARV_2017 ))/sum(preds_L_HARV_2017)
  preds_N_HARV_2017 <- exp(predict(gam_model_N_HARV_2017,new_data_HARV_2017 ))/sum(preds_N_HARV_2017)
  
  f_HARV_2017 <- approxfun(new_data_HARV_2017$time, pmin(preds_L_HARV_2017,preds_N_HARV_2017))
  synch_HARV_2017<-integrate( f_HARV_2017, min(new_data_HARV_2017$time), max(new_data_HARV_2017$time),subdivisions = 1000)
  overlap_HARV_2017<-synch_HARV_2017$value/(AUC(new_data_HARV_2017$time,preds_L_HARV_2017)+AUC(new_data_HARV_2017$time,preds_N_HARV_2017)-synch_HARV_2017$value)
  
  yl_HARV_2018 <- log(L_HARV_2018_out + .0001)
  yn_HARV_2018 <- log(N_HARV_2018_out + .0001)
  
  gam_model_L_HARV_2018 <- gam(yl_HARV_2018 ~ s(time,k=29), gamma = gam_rate)
  gam_model_N_HARV_2018 <- gam(yn_HARV_2018  ~ s(time,k=29), gamma = gam_rate)
  
  preds_L_HARV_2018 <- exp(predict(gam_model_L_HARV_2018,new_data_HARV_2018 ))
  preds_N_HARV_2018 <- exp(predict(gam_model_N_HARV_2018,new_data_HARV_2018 ))
  
  
  preds_L_HARV_2018 <- exp(predict(gam_model_L_HARV_2018,new_data_HARV_2018 ))/sum(preds_L_HARV_2018)
  preds_N_HARV_2018 <- exp(predict(gam_model_N_HARV_2018,new_data_HARV_2018 ))/sum(preds_N_HARV_2018)
  
  f_HARV_2018 <- approxfun(new_data_HARV_2018$time, pmin(preds_L_HARV_2018,preds_N_HARV_2018))
  synch_HARV_2018<-integrate( f_HARV_2018, min(new_data_HARV_2018$time), max(new_data_HARV_2018$time),subdivisions = 1000)
  overlap_HARV_2018<-synch_HARV_2018$value/(AUC(new_data_HARV_2018$time,preds_L_HARV_2018)+AUC(new_data_HARV_2018$time,preds_N_HARV_2018)-synch_HARV_2018$value)
  
  
  yl_HARV_2021 <- log(L_HARV_2021_out + .0001)
  yn_HARV_2021 <- log(N_HARV_2021_out + .0001)
  
  gam_model_L_HARV_2021 <- gam(yl_HARV_2021 ~ s(time,k=29), gamma = gam_rate)
  gam_model_N_HARV_2021 <- gam(yn_HARV_2021  ~ s(time,k=29), gamma = gam_rate)
  
  preds_L_HARV_2021 <- exp(predict(gam_model_L_HARV_2021,new_data_HARV_2021 ))
  preds_N_HARV_2021 <- exp(predict(gam_model_N_HARV_2021,new_data_HARV_2021 ))
  
  
  preds_L_HARV_2021 <- exp(predict(gam_model_L_HARV_2021,new_data_HARV_2021 ))/sum(preds_L_HARV_2021)
  preds_N_HARV_2021 <- exp(predict(gam_model_N_HARV_2021,new_data_HARV_2021 ))/sum(preds_N_HARV_2021)
  
  f_HARV_2021 <- approxfun(new_data_HARV_2021$time, pmin(preds_L_HARV_2021,preds_N_HARV_2021))
  synch_HARV_2021<-integrate( f_HARV_2021, min(new_data_HARV_2021$time), max(new_data_HARV_2021$time),subdivisions = 1000)
  overlap_HARV_2021<-synch_HARV_2021$value/(AUC(new_data_HARV_2021$time,preds_L_HARV_2021)+AUC(new_data_HARV_2021$time,preds_N_HARV_2021)-synch_HARV_2021$value)
  
  
  yl_HARV_2023 <- log(L_HARV_2023_out + .0001)
  yn_HARV_2023 <- log(N_HARV_2023_out + .0001)
  
  gam_model_L_HARV_2023 <- gam(yl_HARV_2023 ~ s(time,k=29), gamma = gam_rate)
  gam_model_N_HARV_2023 <- gam(yn_HARV_2023  ~ s(time,k=29), gamma = gam_rate)
  
  preds_L_HARV_2023 <- exp(predict(gam_model_L_HARV_2023,new_data_HARV_2023 ))
  preds_N_HARV_2023 <- exp(predict(gam_model_N_HARV_2023,new_data_HARV_2023 ))
  
  
  preds_L_HARV_2023 <- exp(predict(gam_model_L_HARV_2023,new_data_HARV_2023 ))/sum(preds_L_HARV_2023)
  preds_N_HARV_2023 <- exp(predict(gam_model_N_HARV_2023,new_data_HARV_2023 ))/sum(preds_N_HARV_2023)
  
  f_HARV_2023 <- approxfun(new_data_HARV_2023$time, pmin(preds_L_HARV_2023,preds_N_HARV_2023))
  synch_HARV_2023<-integrate( f_HARV_2023, min(new_data_HARV_2023$time), max(new_data_HARV_2023$time),subdivisions = 1000)
  overlap_HARV_2023<-synch_HARV_2023$value/(AUC(new_data_HARV_2023$time,preds_L_HARV_2023)+AUC(new_data_HARV_2023$time,preds_N_HARV_2023)-synch_HARV_2023$value)
  
  yl_SCBI_2016 <- log(L_SCBI_2016_out + .0001)
  yn_SCBI_2016 <- log(N_SCBI_2016_out + .0001)
  
  gam_model_L_SCBI_2016 <- gam(yl_SCBI_2016 ~ s(time,k=29), gamma = gam_rate)
  gam_model_N_SCBI_2016 <- gam(yn_SCBI_2016  ~ s(time,k=29), gamma = gam_rate)
  
  
  preds_L_SCBI_2016 <- exp(predict(gam_model_L_SCBI_2016,new_data_SCBI_2016 ))
  preds_N_SCBI_2016 <- exp(predict(gam_model_N_SCBI_2016,new_data_SCBI_2016 ))
  
  preds_L_SCBI_2016 <- exp(predict(gam_model_L_SCBI_2016,new_data_SCBI_2016 ))/sum(preds_L_SCBI_2016)
  preds_N_SCBI_2016 <- exp(predict(gam_model_N_SCBI_2016,new_data_SCBI_2016 ))/sum(preds_N_SCBI_2016)
  
  f_SCBI_2016 <- approxfun(new_data_SCBI_2016$time, pmin(preds_L_SCBI_2016,preds_N_SCBI_2016))
  synch_SCBI_2016<-integrate( f_SCBI_2016, min(new_data_SCBI_2016$time), max(new_data_SCBI_2016$time),subdivisions = 1000)
  overlap_SCBI_2016<-synch_SCBI_2016$value/(AUC(new_data_SCBI_2016$time,preds_L_SCBI_2016)+AUC(new_data_SCBI_2016$time,preds_N_SCBI_2016)-synch_SCBI_2016$value)
  
  yl_SCBI_2017 <- log(L_SCBI_2017_out + .0001)
  yn_SCBI_2017 <- log(N_SCBI_2017_out + .0001)
  
  gam_model_L_SCBI_2017 <- gam(yl_SCBI_2017 ~ s(time,k=29), gamma = 0.55)
  gam_model_N_SCBI_2017 <- gam(yn_SCBI_2017  ~ s(time,k=29), gamma = 0.55)
  
  preds_L_SCBI_2017 <- exp(predict(gam_model_L_SCBI_2017,new_data_SCBI_2017 ))
  preds_N_SCBI_2017 <- exp(predict(gam_model_N_SCBI_2017,new_data_SCBI_2017 ))
  
  preds_L_SCBI_2017 <- exp(predict(gam_model_L_SCBI_2017,new_data_SCBI_2017 ))/sum(preds_L_SCBI_2017)
  preds_N_SCBI_2017 <- exp(predict(gam_model_N_SCBI_2017,new_data_SCBI_2017 ))/sum(preds_N_SCBI_2017)
  
  f_SCBI_2017 <- approxfun(new_data_SCBI_2017$time, pmin(preds_L_SCBI_2017,preds_N_SCBI_2017))
  synch_SCBI_2017<-integrate( f_SCBI_2017, min(new_data_SCBI_2017$time), max(new_data_SCBI_2017$time),subdivisions = 1000)
  overlap_SCBI_2017<-synch_SCBI_2017$value/(AUC(new_data_SCBI_2017$time,preds_L_SCBI_2017)+AUC(new_data_SCBI_2017$time,preds_N_SCBI_2017)-synch_SCBI_2017$value)
  
  yl_SCBI_2018 <- log(L_SCBI_2018_out + .0001)
  yn_SCBI_2018 <- log(N_SCBI_2018_out + .0001)
  
  gam_model_L_SCBI_2018 <- gam(yl_SCBI_2018 ~ s(time,k=29), gamma = gam_rate)
  gam_model_N_SCBI_2018 <- gam(yn_SCBI_2018  ~ s(time,k=29), gamma = gam_rate)
  
  preds_L_SCBI_2018 <- exp(predict(gam_model_L_SCBI_2018,new_data_SCBI_2018 ))
  preds_N_SCBI_2018 <- exp(predict(gam_model_N_SCBI_2018,new_data_SCBI_2018 ))
  
  
  preds_L_SCBI_2018 <- exp(predict(gam_model_L_SCBI_2018,new_data_SCBI_2018 ))/sum(preds_L_SCBI_2018)
  preds_N_SCBI_2018 <- exp(predict(gam_model_N_SCBI_2018,new_data_SCBI_2018 ))/sum(preds_N_SCBI_2018)
  
  f_SCBI_2018 <- approxfun(new_data_SCBI_2018$time, pmin(preds_L_SCBI_2018,preds_N_SCBI_2018))
  synch_SCBI_2018<-integrate( f_SCBI_2018, min(new_data_SCBI_2018$time), max(new_data_SCBI_2018$time),subdivisions = 1000)
  overlap_SCBI_2018<-synch_SCBI_2018$value/(AUC(new_data_SCBI_2018$time,preds_L_SCBI_2018)+AUC(new_data_SCBI_2018$time,preds_N_SCBI_2018)-synch_SCBI_2018$value)
  
  yl_SCBI_2019 <- log(L_SCBI_2019_out + .0001)
  yn_SCBI_2019 <- log(N_SCBI_2019_out + .0001)
  
  gam_model_L_SCBI_2019 <- gam(yl_SCBI_2019 ~ s(time,k=27), gamma = gam_rate)
  gam_model_N_SCBI_2019 <- gam(yn_SCBI_2019  ~ s(time,k=27), gamma = gam_rate)
  
  preds_L_SCBI_2019 <- exp(predict(gam_model_L_SCBI_2019,new_data_SCBI_2019 ))
  preds_N_SCBI_2019 <- exp(predict(gam_model_N_SCBI_2019,new_data_SCBI_2019 ))
  
  
  preds_L_SCBI_2019 <- exp(predict(gam_model_L_SCBI_2019,new_data_SCBI_2019 ))/sum(preds_L_SCBI_2019)
  preds_N_SCBI_2019 <- exp(predict(gam_model_N_SCBI_2019,new_data_SCBI_2019 ))/sum(preds_N_SCBI_2019)
  
  f_SCBI_2019 <- approxfun(new_data_SCBI_2019$time, pmin(preds_L_SCBI_2019,preds_N_SCBI_2019))
  synch_SCBI_2019<-integrate( f_SCBI_2019, min(new_data_SCBI_2019$time), max(new_data_SCBI_2019$time),subdivisions = 1000)
  overlap_SCBI_2019<-synch_SCBI_2019$value/(AUC(new_data_SCBI_2019$time,preds_L_SCBI_2019)+AUC(new_data_SCBI_2019$time,preds_N_SCBI_2019)-synch_SCBI_2019$value)
  
  
  
  yl_SCBI_2023 <- log(L_SCBI_2023_out + .0001)
  yn_SCBI_2023 <- log(N_SCBI_2023_out + .0001)
  
  gam_model_L_SCBI_2023 <- gam(yl_SCBI_2023 ~ s(time,k=29), gamma = gam_rate)
  gam_model_N_SCBI_2023 <- gam(yn_SCBI_2023  ~ s(time,k=29), gamma = gam_rate)
  
  preds_L_SCBI_2023 <- exp(predict(gam_model_L_SCBI_2023,new_data_SCBI_2023 ))
  preds_N_SCBI_2023 <- exp(predict(gam_model_N_SCBI_2023,new_data_SCBI_2023 ))
  
  
  preds_L_SCBI_2023 <- exp(predict(gam_model_L_SCBI_2023,new_data_SCBI_2023 ))/sum(preds_L_SCBI_2023)
  preds_N_SCBI_2023 <- exp(predict(gam_model_N_SCBI_2023,new_data_SCBI_2023 ))/sum(preds_N_SCBI_2023)
  
  f_SCBI_2023 <- approxfun(new_data_SCBI_2023$time, pmin(preds_L_SCBI_2023,preds_N_SCBI_2023))
  synch_SCBI_2023<-integrate( f_SCBI_2023, min(new_data_SCBI_2023$time), max(new_data_SCBI_2023$time),subdivisions = 1000)
  overlap_SCBI_2023<-synch_SCBI_2023$value/(AUC(new_data_SCBI_2023$time,preds_L_SCBI_2023)+AUC(new_data_SCBI_2023$time,preds_N_SCBI_2023)-synch_SCBI_2023$value)
  
  
  
  
  
  ##combine all the overlaps
  all_overlaps_new<-c(overlap_BLAN_2016,
                      overlap_BLAN_2017,
                      overlap_BLAN_2018,
                      overlap_BLAN_2019,
                      
                      overlap_HARV_2017,
                      overlap_HARV_2018,
                      
                      overlap_HARV_2021,
                      
                      overlap_HARV_2023,
                      overlap_SCBI_2016,
                      overlap_SCBI_2017,
                      overlap_SCBI_2018,
                      overlap_SCBI_2019,
                      
                      overlap_SCBI_2023)
  
  
  
  
  
  B_L_16_mat[,samp_idx]<-L_BLAN_2016_out
  B_L_17_mat[,samp_idx]<-L_BLAN_2017_out
  B_L_18_mat[,samp_idx]<-L_BLAN_2018_out
  B_L_19_mat[,samp_idx]<-L_BLAN_2019_out
  
  
  B_N_16_mat[,samp_idx]<-N_BLAN_2016_out
  B_N_17_mat[,samp_idx]<-N_BLAN_2017_out
  B_N_18_mat[,samp_idx]<-N_BLAN_2018_out
  B_N_19_mat[,samp_idx]<-N_BLAN_2019_out
 

  H_L_17_mat[,samp_idx]<-L_HARV_2017_out
  H_L_18_mat[,samp_idx]<-L_HARV_2018_out
  
  H_L_21_mat[,samp_idx]<-L_HARV_2021_out
  
  H_L_23_mat[,samp_idx]<-L_HARV_2023_out
  
  
  H_N_17_mat[,samp_idx]<-N_HARV_2017_out
  H_N_18_mat[,samp_idx]<-N_HARV_2018_out

  H_N_21_mat[,samp_idx]<-N_HARV_2021_out
  
  H_N_23_mat[,samp_idx]<-N_HARV_2023_out
  
  S_L_16_mat[,samp_idx]<-L_SCBI_2016_out
  S_L_17_mat[,samp_idx]<-L_SCBI_2017_out
  S_L_18_mat[,samp_idx]<-L_SCBI_2018_out
  S_L_19_mat[,samp_idx]<-L_SCBI_2019_out
  
  S_L_23_mat[,samp_idx]<-L_SCBI_2023_out
  
  S_N_16_mat[,samp_idx]<-N_SCBI_2016_out
  S_N_17_mat[,samp_idx]<-N_SCBI_2017_out
  S_N_18_mat[,samp_idx]<-N_SCBI_2018_out
  S_N_19_mat[,samp_idx]<-N_SCBI_2019_out
  
  S_N_23_mat[,samp_idx]<-N_SCBI_2023_out
  
  all_overlaps[,samp_idx]<-all_overlaps_new
  
}

all_overlaps <- all_overlaps[, colSums(all_overlaps != 0) > 0]
B_L_16_mat <- B_L_16_mat[, colSums(B_L_16_mat != 0) > 0]
B_N_16_mat <- B_N_16_mat[, colSums(B_N_16_mat != 0) > 0]
B_L_17_mat <- B_L_17_mat[, colSums(B_L_17_mat != 0) > 0]
B_N_17_mat <- B_N_17_mat[, colSums(B_N_17_mat != 0) > 0]
B_L_18_mat <- B_L_18_mat[, colSums(B_L_18_mat != 0) > 0]
B_N_18_mat <- B_N_18_mat[, colSums(B_N_18_mat != 0) > 0]
B_L_19_mat <- B_L_19_mat[, colSums(B_L_19_mat != 0) > 0]
B_N_19_mat <- B_N_19_mat[, colSums(B_N_19_mat != 0) > 0]



H_L_17_mat <- H_L_17_mat[, colSums(H_L_17_mat != 0) > 0]
H_N_17_mat <- H_N_17_mat[, colSums(H_N_17_mat != 0) > 0]
H_L_18_mat <- H_L_18_mat[, colSums(H_L_18_mat != 0) > 0]
H_N_18_mat <- H_N_18_mat[, colSums(H_N_18_mat != 0) > 0]


H_L_21_mat <- H_L_21_mat[, colSums(H_L_21_mat != 0) > 0]
H_N_21_mat <- H_N_21_mat[, colSums(H_N_21_mat != 0) > 0]
H_L_23_mat <- H_L_23_mat[, colSums(H_L_23_mat != 0) > 0]
H_N_23_mat <- H_N_23_mat[, colSums(H_N_23_mat != 0) > 0]

S_L_16_mat <- S_L_16_mat[, colSums(S_L_16_mat != 0) > 0]
S_N_16_mat <- S_N_16_mat[, colSums(S_N_16_mat != 0) > 0]
S_L_17_mat <- S_L_17_mat[, colSums(S_L_17_mat != 0) > 0]
S_N_17_mat <- S_N_17_mat[, colSums(S_N_17_mat != 0) > 0]
S_L_18_mat <- S_L_18_mat[, colSums(S_L_18_mat != 0) > 0]
S_N_18_mat <- S_N_18_mat[, colSums(S_N_18_mat != 0) > 0]

S_L_19_mat <- S_L_19_mat[, colSums(S_L_19_mat != 0) > 0]
S_N_19_mat <- S_N_19_mat[, colSums(S_N_19_mat != 0) > 0]
S_L_23_mat <- S_L_23_mat[, colSums(S_L_23_mat != 0) > 0]
S_N_23_mat <- S_N_23_mat[, colSums(S_N_23_mat != 0) > 0]


write.csv(all_overlaps[,1:100], paste0("~/Desktop/CH2_FEB26/all_overlaps_posterior_unc.csv"))


write.csv(B_L_16_mat[,1:100], paste0("~/Desktop/CH2_FEB26/B_L_16_mat_posterior_unc.csv"))
write.csv(B_L_17_mat[,1:100], paste0("~/Desktop/CH2_FEB26/B_L_17_mat_posterior_unc.csv"))
write.csv(B_L_18_mat[,1:100], paste0("~/Desktop/CH2_FEB26/B_L_18_mat_posterior_unc.csv"))
write.csv(B_L_19_mat[,1:100], paste0("~/Desktop/CH2_FEB26/B_L_19_mat_posterior_unc.csv"))



write.csv(H_L_17_mat[,1:100], paste0("~/Desktop/CH2_FEB26/H_L_17_mat_posterior_unc.csv"))
write.csv(H_L_18_mat[,1:100], paste0("~/Desktop/CH2_FEB26/H_L_18_mat_posterior_unc.csv"))
write.csv(H_L_21_mat[,1:100], paste0("~/Desktop/CH2_FEB26/H_L_21_mat_posterior_unc.csv"))
write.csv(H_L_23_mat[,1:100], paste0("~/Desktop/CH2_FEB26/H_L_23_mat_posterior_unc.csv"))

write.csv(S_L_16_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_L_16_mat_posterior_unc.csv"))
write.csv(S_L_17_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_L_17_mat_posterior_unc.csv"))
write.csv(S_L_18_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_L_18_mat_posterior_unc.csv"))
write.csv(S_L_19_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_L_19_mat_posterior_unc.csv"))
write.csv(S_L_23_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_L_23_mat_posterior_unc.csv"))


# 
# 
write.csv(B_N_16_mat[,1:100], paste0("~/Desktop/CH2_FEB26/B_N_16_mat_posterior_unc.csv"))
write.csv(B_N_17_mat[,1:100], paste0("~/Desktop/CH2_FEB26/B_N_17_mat_posterior_unc.csv"))
write.csv(B_N_18_mat[,1:100], paste0("~/Desktop/CH2_FEB26/B_N_18_mat_posterior_unc.csv"))
write.csv(B_N_19_mat[,1:100], paste0("~/Desktop/CH2_FEB26/B_N_19_mat_posterior_unc.csv"))



write.csv(H_N_17_mat[,1:100], paste0("~/Desktop/CH2_FEB26/H_N_17_mat_posterior_unc.csv"))
write.csv(H_N_18_mat[,1:100], paste0("~/Desktop/CH2_FEB26/H_N_18_mat_posterior_unc.csv"))
write.csv(H_N_21_mat[,1:100], paste0("~/Desktop/CH2_FEB26/H_N_21_mat_posterior_unc.csv"))
write.csv(H_N_23_mat[,1:100], paste0("~/Desktop/CH2_FEB26/H_N_23_mat_posterior_unc.csv"))



write.csv(S_N_16_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_N_16_mat_posterior_unc.csv"))
write.csv(S_N_17_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_N_17_mat_posterior_unc.csv"))
write.csv(S_N_18_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_N_18_mat_posterior_unc.csv"))
write.csv(S_N_19_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_N_19_mat_posterior_unc.csv"))
write.csv(S_N_23_mat[,1:100], paste0("~/Desktop/CH2_FEB26/S_N_23_mat_posterior_unc.csv"))




##make bar graph of model compared to data 

BLAN_16L_data_from_posterior<-matrix(0,nrow=nrow(BLAN_a_2016), ncol = 100)
BLAN_16N_data_from_posterior<-matrix(0,nrow=nrow(BLAN_a_2016), ncol = 100)
BLAN_17L_data_from_posterior<-matrix(0,nrow=nrow(BLAN_a_2017), ncol = 100)
BLAN_17N_data_from_posterior<-matrix(0,nrow=nrow(BLAN_a_2017), ncol = 100)
BLAN_18L_data_from_posterior<-matrix(0,nrow=nrow(BLAN_a_2018), ncol = 100)
BLAN_18N_data_from_posterior<-matrix(0,nrow=nrow(BLAN_a_2018), ncol = 100)
BLAN_19L_data_from_posterior<-matrix(0,nrow=nrow(BLAN_a_2019), ncol = 100)
BLAN_19N_data_from_posterior<-matrix(0,nrow=nrow(BLAN_a_2019), ncol = 100)




HARV_17L_data_from_posterior<-matrix(0,nrow=nrow(HARV_a_2017), ncol = 100)
HARV_17N_data_from_posterior<-matrix(0,nrow=nrow(HARV_a_2017), ncol = 100)
HARV_18L_data_from_posterior<-matrix(0,nrow=nrow(HARV_a_2018), ncol = 100)
HARV_18N_data_from_posterior<-matrix(0,nrow=nrow(HARV_a_2018), ncol = 100)
HARV_21L_data_from_posterior<-matrix(0,nrow=nrow(HARV_a_2021), ncol = 100)
HARV_21N_data_from_posterior<-matrix(0,nrow=nrow(HARV_a_2021), ncol = 100)
HARV_23L_data_from_posterior<-matrix(0,nrow=nrow(HARV_a_2023), ncol = 100)
HARV_23N_data_from_posterior<-matrix(0,nrow=nrow(HARV_a_2023), ncol = 100)


SCBI_16L_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2016), ncol = 100)
SCBI_16N_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2016), ncol = 100)
SCBI_17L_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2017), ncol = 100)
SCBI_17N_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2017), ncol = 100)
SCBI_18L_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2018), ncol = 100)
SCBI_18N_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2018), ncol = 100)
SCBI_19L_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2019), ncol = 100)
SCBI_19N_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2019), ncol = 100)
SCBI_23L_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2023), ncol = 100)
SCBI_23N_data_from_posterior<-matrix(0,nrow=nrow(SCBI_a_2023), ncol = 100)

##Pull out sampling days

for (pred_idx in 1:100)
{
  
 
  se_B_L_16<-B_L_16_mat[unique(BLAN_a_2016$DOY), pred_idx] * (BLAN_a_2016$area[seq(1, length(BLAN_a_2016$area), by = 2)] /max(BLAN_a_2016$area))
  se_B_L_17<-B_L_17_mat[unique(BLAN_a_2017$DOY), pred_idx] * (BLAN_a_2017$area[seq(1, length(BLAN_a_2017$area), by = 2)] /max(BLAN_a_2017$area))
  se_B_L_18<-B_L_18_mat[unique(BLAN_a_2018$DOY), pred_idx] * (BLAN_a_2018$area[seq(1, length(BLAN_a_2018$area), by = 2)] /max(BLAN_a_2018$area))
  se_B_L_19<-B_L_19_mat[unique(BLAN_a_2019$DOY), pred_idx] * (BLAN_a_2019$area[seq(1, length(BLAN_a_2019$area), by = 2)] /max(BLAN_a_2019$area))
  
  se_B_N_16<-B_N_16_mat[unique(BLAN_a_2016$DOY), pred_idx] * (BLAN_a_2016$area[seq(1, length(BLAN_a_2016$area), by = 2)] /max(BLAN_a_2016$area))
  se_B_N_17<-B_N_17_mat[unique(BLAN_a_2017$DOY), pred_idx] * (BLAN_a_2017$area[seq(1, length(BLAN_a_2017$area), by = 2)] /max(BLAN_a_2017$area))
  se_B_N_18<-B_N_18_mat[unique(BLAN_a_2018$DOY), pred_idx] * (BLAN_a_2018$area[seq(1, length(BLAN_a_2018$area), by = 2)] /max(BLAN_a_2018$area))
  se_B_N_19<-B_N_19_mat[unique(BLAN_a_2019$DOY), pred_idx] * (BLAN_a_2019$area[seq(1, length(BLAN_a_2019$area), by = 2)] /max(BLAN_a_2019$area))
  
  se_H_L_17<-H_L_17_mat[unique(HARV_a_2017$DOY), pred_idx] * (HARV_a_2017$area[seq(1, length(HARV_a_2017$area), by = 2)] /max(HARV_a_2017$area))
  se_H_L_18<-H_L_18_mat[unique(HARV_a_2018$DOY), pred_idx] * (HARV_a_2018$area[seq(1, length(HARV_a_2018$area), by = 2)] /max(HARV_a_2018$area))
  se_H_L_21<-H_L_21_mat[unique(HARV_a_2021$DOY), pred_idx] * (HARV_a_2021$area[seq(1, length(HARV_a_2021$area), by = 2)] /max(HARV_a_2021$area))
  se_H_L_23<-H_L_23_mat[unique(HARV_a_2023$DOY), pred_idx] * (HARV_a_2023$area[seq(1, length(HARV_a_2023$area), by = 2)] /max(HARV_a_2023$area))
  
  se_H_N_17<-H_N_17_mat[unique(HARV_a_2017$DOY), pred_idx] * (HARV_a_2017$area[seq(1, length(HARV_a_2017$area), by = 2)] /max(HARV_a_2017$area))
  se_H_N_18<-H_N_18_mat[unique(HARV_a_2018$DOY), pred_idx] * (HARV_a_2018$area[seq(1, length(HARV_a_2018$area), by = 2)] /max(HARV_a_2018$area))
  se_H_N_21<-H_N_21_mat[unique(HARV_a_2021$DOY), pred_idx] * (HARV_a_2021$area[seq(1, length(HARV_a_2021$area), by = 2)] /max(HARV_a_2021$area))
  se_H_N_23<-H_N_23_mat[unique(HARV_a_2023$DOY), pred_idx] * (HARV_a_2023$area[seq(1, length(HARV_a_2023$area), by = 2)] /max(HARV_a_2023$area))
  
  se_S_L_16<-S_L_16_mat[unique(SCBI_a_2016$DOY), pred_idx] * (SCBI_a_2016$area[seq(1, length(SCBI_a_2016$area), by = 2)] /max(SCBI_a_2016$area))
  se_S_L_17<-S_L_17_mat[unique(SCBI_a_2017$DOY), pred_idx] * (SCBI_a_2017$area[seq(1, length(SCBI_a_2017$area), by = 2)] /max(SCBI_a_2017$area))
  se_S_L_18<-S_L_18_mat[unique(SCBI_a_2018$DOY), pred_idx] * (SCBI_a_2018$area[seq(1, length(SCBI_a_2018$area), by = 2)] /max(SCBI_a_2018$area))
  se_S_L_19<-S_L_19_mat[unique(SCBI_a_2019$DOY), pred_idx] * (SCBI_a_2019$area[seq(1, length(SCBI_a_2019$area), by = 2)] /max(SCBI_a_2019$area))
  se_S_L_23<-S_L_23_mat[unique(SCBI_a_2023$DOY), pred_idx] * (SCBI_a_2023$area[seq(1, length(SCBI_a_2023$area), by = 2)] /max(SCBI_a_2023$area))
  
  
  se_S_N_16<-S_N_16_mat[unique(SCBI_a_2016$DOY), pred_idx] * (SCBI_a_2016$area[seq(1, length(SCBI_a_2016$area), by = 2)] /max(SCBI_a_2016$area))
  se_S_N_17<-S_N_17_mat[unique(SCBI_a_2017$DOY), pred_idx] * (SCBI_a_2017$area[seq(1, length(SCBI_a_2017$area), by = 2)] /max(SCBI_a_2017$area))
  se_S_N_18<-S_N_18_mat[unique(SCBI_a_2018$DOY), pred_idx] * (SCBI_a_2018$area[seq(1, length(SCBI_a_2018$area), by = 2)] /max(SCBI_a_2018$area))
  se_S_N_19<-S_N_19_mat[unique(SCBI_a_2019$DOY), pred_idx] * (SCBI_a_2019$area[seq(1, length(SCBI_a_2019$area), by = 2)] /max(SCBI_a_2019$area))
  se_S_N_23<-S_N_23_mat[unique(SCBI_a_2023$DOY), pred_idx] * (SCBI_a_2023$area[seq(1, length(SCBI_a_2023$area), by = 2)] /max(SCBI_a_2023$area))
  
 
  
  

  

  
  
  p_B_L_16<-se_B_L_16/sum(se_B_L_16)
  p_B_L_17<-se_B_L_17/sum(se_B_L_17)
  p_B_L_18<-se_B_L_18/sum(se_B_L_18)
  p_B_L_19<-se_B_L_19/sum(se_B_L_19)
  
  
  p_B_N_16<-se_B_N_16/sum(se_B_N_16)
  p_B_N_17<-se_B_N_17/sum(se_B_N_17)
  p_B_N_18<-se_B_N_18/sum(se_B_N_18)
  p_B_N_19<-se_B_N_19/sum(se_B_N_19)
  
  
  
  p_H_L_17<-se_H_L_17/sum(se_H_L_17)
  p_H_L_18<-se_H_L_18/sum(se_H_L_18)
  
  p_H_L_21<-se_H_L_21/sum(se_H_L_21)
  
  p_H_L_23<-se_H_L_23/sum(se_H_L_23)
  
  
  p_H_N_17<-se_H_N_17/sum(se_H_N_17)
  p_H_N_18<-se_H_N_18/sum(se_H_N_18)
  
  p_H_N_21<-se_H_N_21/sum(se_H_N_21)
  
  p_H_N_23<-se_H_N_23/sum(se_H_N_23)
  
  p_S_L_16<-se_S_L_16/sum(se_S_L_16)
  p_S_L_17<-se_S_L_17/sum(se_S_L_17)
  p_S_L_18<-se_S_L_18/sum(se_S_L_18)
  p_S_L_19<-se_S_L_19/sum(se_S_L_19)
  
  p_S_L_23<-se_S_L_23/sum(se_S_L_23)
  
  p_S_N_16<-se_S_N_16/sum(se_S_N_16)
  p_S_N_17<-se_S_N_17/sum(se_S_N_17)
  p_S_N_18<-se_S_N_18/sum(se_S_N_18)
  p_S_N_19<-se_S_N_19/sum(se_S_N_19)
  
  p_S_N_23<-se_S_N_23/sum(se_S_N_23)
  
  
  
  
  current_B_L_2016<-rdirmnom(1,sum(subset(BLAN_a_2016, BLAN_a_2016$stage == "L")$ticks),p_B_L_16*exp(log_alpha0_B_L_2016[pred_idx]))
  current_B_L_2017<-rdirmnom(1,sum(subset(BLAN_a_2017, BLAN_a_2017$stage == "L")$ticks),p_B_L_17*exp(log_alpha0_B_L_2017[pred_idx]))
  current_B_L_2018<-rdirmnom(1,sum(subset(BLAN_a_2018, BLAN_a_2018$stage == "L")$ticks),p_B_L_18*exp(log_alpha0_B_L_2018[pred_idx]))
  current_B_L_2019<-rdirmnom(1,sum(subset(BLAN_a_2019, BLAN_a_2019$stage == "L")$ticks),p_B_L_19*exp(log_alpha0_B_L_2019[pred_idx]))
  
  current_B_N_2016<-rdirmnom(1,sum(subset(BLAN_a_2016, BLAN_a_2016$stage == "N")$ticks),p_B_N_16*exp(log_alpha0_B_N_2016[pred_idx]))
  current_B_N_2017<-rdirmnom(1,sum(subset(BLAN_a_2017, BLAN_a_2017$stage == "N")$ticks),p_B_N_17*exp(log_alpha0_B_N_2017[pred_idx]))
  current_B_N_2018<-rdirmnom(1,sum(subset(BLAN_a_2018, BLAN_a_2018$stage == "N")$ticks),p_B_N_18*exp(log_alpha0_B_N_2018[pred_idx]))
  current_B_N_2019<-rdirmnom(1,sum(subset(BLAN_a_2019, BLAN_a_2019$stage == "N")$ticks),p_B_N_19*exp(log_alpha0_B_N_2019[pred_idx]))
  
  
  current_H_L_2017<-rdirmnom(1,sum(subset(HARV_a_2017, HARV_a_2017$stage == "L")$ticks),p_H_L_17*exp(log_alpha0_H_L_2017[pred_idx]))
  current_H_L_2018<-rdirmnom(1,sum(subset(HARV_a_2018, HARV_a_2018$stage == "L")$ticks),p_H_L_18*exp(log_alpha0_H_L_2018[pred_idx]))
  current_H_L_2021<-rdirmnom(1,sum(subset(HARV_a_2021, HARV_a_2021$stage == "L")$ticks),p_H_L_21*exp(log_alpha0_H_L_2021[pred_idx]))
  current_H_L_2023<-rdirmnom(1,sum(subset(HARV_a_2023, HARV_a_2023$stage == "L")$ticks),p_H_L_23*exp(log_alpha0_H_L_2023[pred_idx]))
  
  
  current_H_N_2017<-rdirmnom(1,sum(subset(HARV_a_2017, HARV_a_2017$stage == "N")$ticks),p_H_N_17*exp(log_alpha0_H_N_2017[pred_idx]))
  current_H_N_2018<-rdirmnom(1,sum(subset(HARV_a_2018, HARV_a_2018$stage == "N")$ticks),p_H_N_18*exp(log_alpha0_H_N_2018[pred_idx]))
  current_H_N_2021<-rdirmnom(1,sum(subset(HARV_a_2021, HARV_a_2021$stage == "N")$ticks),p_H_N_21*exp(log_alpha0_H_N_2021[pred_idx]))
  current_H_N_2023<-rdirmnom(1,sum(subset(HARV_a_2023, HARV_a_2023$stage == "N")$ticks),p_H_N_23*exp(log_alpha0_H_N_2023[pred_idx]))
  
  
  
  current_S_L_2016<-rdirmnom(1,sum(subset(SCBI_a_2016, SCBI_a_2016$stage == "L")$ticks),p_S_L_16*exp(log_alpha0_S_L_2016[pred_idx]))
  current_S_L_2017<-rdirmnom(1,sum(subset(SCBI_a_2017, SCBI_a_2017$stage == "L")$ticks),p_S_L_17*exp(log_alpha0_S_L_2017[pred_idx]))
  current_S_L_2018<-rdirmnom(1,sum(subset(SCBI_a_2018, SCBI_a_2018$stage == "L")$ticks),p_S_L_18*exp(log_alpha0_S_L_2018[pred_idx]))
  current_S_L_2019<-rdirmnom(1,sum(subset(SCBI_a_2019, SCBI_a_2019$stage == "L")$ticks),p_S_L_19*exp(log_alpha0_S_L_2019[pred_idx]))
  current_S_L_2023<-rdirmnom(1,sum(subset(SCBI_a_2023, SCBI_a_2023$stage == "L")$ticks),p_S_L_23*exp(log_alpha0_S_L_2023[pred_idx]))
  
  current_S_N_2016<-rdirmnom(1,sum(subset(SCBI_a_2016, SCBI_a_2016$stage == "N")$ticks),p_S_N_16*exp(log_alpha0_S_N_2016[pred_idx]))
  current_S_N_2017<-rdirmnom(1,sum(subset(SCBI_a_2017, SCBI_a_2017$stage == "N")$ticks),p_S_N_17*exp(log_alpha0_S_N_2017[pred_idx]))
  current_S_N_2018<-rdirmnom(1,sum(subset(SCBI_a_2018, SCBI_a_2018$stage == "N")$ticks),p_S_N_18*exp(log_alpha0_S_N_2018[pred_idx]))
  current_S_N_2019<-rdirmnom(1,sum(subset(SCBI_a_2019, SCBI_a_2019$stage == "N")$ticks),p_S_N_19*exp(log_alpha0_S_N_2019[pred_idx]))
  current_S_N_2023<-rdirmnom(1,sum(subset(SCBI_a_2023, SCBI_a_2023$stage == "N")$ticks),p_S_N_23*exp(log_alpha0_S_N_2023[pred_idx]))
  
  
  
  ##bind preds to DOY
  sort_current_B_L_2016<-cbind(BLAN_a_2016$DOY,as.vector(current_B_L_2016))
  sort_current_B_L_2017<-cbind(BLAN_a_2017$DOY,as.vector(current_B_L_2017))
  sort_current_B_L_2018<-cbind(BLAN_a_2018$DOY,as.vector(current_B_L_2018))
  sort_current_B_L_2019<-cbind(BLAN_a_2019$DOY,as.vector(current_B_L_2019))
  
  sort_current_B_N_2016<-cbind(BLAN_a_2016$DOY,as.vector(current_B_N_2016))
  sort_current_B_N_2017<-cbind(BLAN_a_2017$DOY,as.vector(current_B_N_2017))
  sort_current_B_N_2018<-cbind(BLAN_a_2018$DOY,as.vector(current_B_N_2018))
  sort_current_B_N_2019<-cbind(BLAN_a_2019$DOY,as.vector(current_B_N_2019))
  
  
  sort_current_H_L_2017<-cbind(HARV_a_2017$DOY,as.vector(current_H_L_2017))
  sort_current_H_L_2018<-cbind(HARV_a_2018$DOY,as.vector(current_H_L_2018))
  
  sort_current_H_L_2021<-cbind(HARV_a_2021$DOY,as.vector(current_H_L_2021))
  
  sort_current_H_L_2023<-cbind(HARV_a_2023$DOY,as.vector(current_H_L_2023))
  
  
  sort_current_H_N_2017<-cbind(HARV_a_2017$DOY,as.vector(current_H_N_2017))
  sort_current_H_N_2018<-cbind(HARV_a_2018$DOY,as.vector(current_H_N_2018))
  
  sort_current_H_N_2021<-cbind(HARV_a_2021$DOY,as.vector(current_H_N_2021))
  
  sort_current_H_N_2023<-cbind(HARV_a_2023$DOY,as.vector(current_H_N_2023))
  
  sort_current_S_L_2016<-cbind(SCBI_a_2016$DOY,as.vector(current_S_L_2016))
  sort_current_S_L_2017<-cbind(SCBI_a_2017$DOY,as.vector(current_S_L_2017))
  sort_current_S_L_2018<-cbind(SCBI_a_2018$DOY,as.vector(current_S_L_2018))
  sort_current_S_L_2019<-cbind(SCBI_a_2019$DOY,as.vector(current_S_L_2019))
  
  sort_current_S_L_2023<-cbind(SCBI_a_2023$DOY,as.vector(current_S_L_2023))
  
  sort_current_S_N_2016<-cbind(SCBI_a_2016$DOY,as.vector(current_S_N_2016))
  sort_current_S_N_2017<-cbind(SCBI_a_2017$DOY,as.vector(current_S_N_2017))
  sort_current_S_N_2018<-cbind(SCBI_a_2018$DOY,as.vector(current_S_N_2018))
  sort_current_S_N_2019<-cbind(SCBI_a_2019$DOY,as.vector(current_S_N_2019))
  
  sort_current_S_N_2023<-cbind(SCBI_a_2023$DOY,as.vector(current_S_N_2023))
  
  save_sort_current_B_L_2016<-sort_current_B_L_2016[order(sort_current_B_L_2016[,1]),][,2]
  save_sort_current_B_L_2017<-sort_current_B_L_2017[order(sort_current_B_L_2017[,1]),][,2]
  save_sort_current_B_L_2018<-sort_current_B_L_2018[order(sort_current_B_L_2018[,1]),][,2]
  save_sort_current_B_L_2019<-sort_current_B_L_2019[order(sort_current_B_L_2019[,1]),][,2]
  
  save_sort_current_B_N_2016<-sort_current_B_N_2016[order(sort_current_B_N_2016[,1]),][,2]
  save_sort_current_B_N_2017<-sort_current_B_N_2017[order(sort_current_B_N_2017[,1]),][,2]
  save_sort_current_B_N_2018<-sort_current_B_N_2018[order(sort_current_B_N_2018[,1]),][,2]
  save_sort_current_B_N_2019<-sort_current_B_N_2019[order(sort_current_B_N_2019[,1]),][,2]
  
  
  save_sort_current_H_L_2017<-sort_current_H_L_2017[order(sort_current_H_L_2017[,1]),][,2]
  save_sort_current_H_L_2018<-sort_current_H_L_2018[order(sort_current_H_L_2018[,1]),][,2]
  save_sort_current_H_L_2021<-sort_current_H_L_2021[order(sort_current_H_L_2021[,1]),][,2]
  save_sort_current_H_L_2023<-sort_current_H_L_2023[order(sort_current_H_L_2023[,1]),][,2]
  
  
  save_sort_current_H_N_2017<-sort_current_H_N_2017[order(sort_current_H_N_2017[,1]),][,2]
  save_sort_current_H_N_2018<-sort_current_H_N_2018[order(sort_current_H_N_2018[,1]),][,2]
  save_sort_current_H_N_2021<-sort_current_H_N_2021[order(sort_current_H_N_2021[,1]),][,2]
  save_sort_current_H_N_2023<-sort_current_H_N_2023[order(sort_current_H_N_2023[,1]),][,2]
  
  save_sort_current_S_L_2016<-sort_current_S_L_2016[order(sort_current_S_L_2016[,1]),][,2]
  save_sort_current_S_L_2017<-sort_current_S_L_2017[order(sort_current_S_L_2017[,1]),][,2]
  save_sort_current_S_L_2018<-sort_current_S_L_2018[order(sort_current_S_L_2018[,1]),][,2]
  save_sort_current_S_L_2019<-sort_current_S_L_2019[order(sort_current_S_L_2019[,1]),][,2]
  save_sort_current_S_L_2023<-sort_current_S_L_2023[order(sort_current_S_L_2023[,1]),][,2]
  
  save_sort_current_S_N_2016<-sort_current_S_N_2016[order(sort_current_S_N_2016[,1]),][,2]
  save_sort_current_S_N_2017<-sort_current_S_N_2017[order(sort_current_S_N_2017[,1]),][,2]
  save_sort_current_S_N_2018<-sort_current_S_N_2018[order(sort_current_S_N_2018[,1]),][,2]
  save_sort_current_S_N_2019<-sort_current_S_N_2019[order(sort_current_S_N_2019[,1]),][,2]
  save_sort_current_S_N_2023<-sort_current_S_N_2023[order(sort_current_S_N_2023[,1]),][,2]
  
  
  
  BLAN_16L_data_from_posterior[,pred_idx]<-save_sort_current_B_L_2016
  BLAN_17L_data_from_posterior[,pred_idx]<-save_sort_current_B_L_2017
  BLAN_18L_data_from_posterior[,pred_idx]<-save_sort_current_B_L_2018
  BLAN_19L_data_from_posterior[,pred_idx]<-save_sort_current_B_L_2019

  
  BLAN_16N_data_from_posterior[,pred_idx]<-save_sort_current_B_N_2016
  BLAN_17N_data_from_posterior[,pred_idx]<-save_sort_current_B_N_2017
  BLAN_18N_data_from_posterior[,pred_idx]<-save_sort_current_B_N_2018
  BLAN_19N_data_from_posterior[,pred_idx]<-save_sort_current_B_N_2019

  
  HARV_17L_data_from_posterior[,pred_idx]<-save_sort_current_H_L_2017
  HARV_18L_data_from_posterior[,pred_idx]<-save_sort_current_H_L_2018
  HARV_21L_data_from_posterior[,pred_idx]<-save_sort_current_H_L_2021
  HARV_23L_data_from_posterior[,pred_idx]<-save_sort_current_H_L_2023
  
  
  
  HARV_17N_data_from_posterior[,pred_idx]<-save_sort_current_H_N_2017
  HARV_18N_data_from_posterior[,pred_idx]<-save_sort_current_H_N_2018
  HARV_21N_data_from_posterior[,pred_idx]<-save_sort_current_H_N_2021
  HARV_23N_data_from_posterior[,pred_idx]<-save_sort_current_H_N_2023
  
  
  SCBI_16L_data_from_posterior[,pred_idx]<-save_sort_current_S_L_2016
  SCBI_17L_data_from_posterior[,pred_idx]<-save_sort_current_S_L_2017
  SCBI_18L_data_from_posterior[,pred_idx]<-save_sort_current_S_L_2018
  SCBI_19L_data_from_posterior[,pred_idx]<-save_sort_current_S_L_2019
  SCBI_23L_data_from_posterior[,pred_idx]<-save_sort_current_S_L_2023
  
  
  SCBI_16N_data_from_posterior[,pred_idx]<-save_sort_current_S_N_2016
  SCBI_17N_data_from_posterior[,pred_idx]<-save_sort_current_S_N_2017
  SCBI_18N_data_from_posterior[,pred_idx]<-save_sort_current_S_N_2018
  SCBI_19N_data_from_posterior[,pred_idx]<-save_sort_current_S_N_2019
  SCBI_23N_data_from_posterior[,pred_idx]<-save_sort_current_S_N_2023
  
  
}


##write out estimates
write.csv(BLAN_16L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_BLAN_L_16.csv")
write.csv(BLAN_17L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_BLAN_L_17.csv")
write.csv(BLAN_18L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_BLAN_L_18.csv")
write.csv(BLAN_19L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_BLAN_L_19.csv")

write.csv(BLAN_16N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_BLAN_N_16.csv")
write.csv(BLAN_17N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_BLAN_N_17.csv")
write.csv(BLAN_18N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_BLAN_N_18.csv")
write.csv(BLAN_19N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_BLAN_N_19.csv")


write.csv(HARV_17L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_HARV_L_17.csv")
write.csv(HARV_18L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_HARV_L_18.csv")
write.csv(HARV_21L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_HARV_L_21.csv")
write.csv(HARV_23L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_HARV_L_23.csv")


write.csv(HARV_17N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_HARV_N_17.csv")
write.csv(HARV_18N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_HARV_N_18.csv")
write.csv(HARV_21N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_HARV_N_21.csv")
write.csv(HARV_23N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_HARV_N_23.csv")

write.csv(SCBI_16L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_L_16.csv")
write.csv(SCBI_17L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_L_17.csv")
write.csv(SCBI_18L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_L_18.csv")
write.csv(SCBI_19L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_L_19.csv")
write.csv(SCBI_23L_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_L_23.csv")



write.csv(SCBI_16N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_N_16.csv")
write.csv(SCBI_17N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_N_17.csv")
write.csv(SCBI_18N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_N_18.csv")
write.csv(SCBI_19N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_N_19.csv")
write.csv(SCBI_23N_data_from_posterior,"~/Desktop/CH2_FEB26/annual_posterior_preds_SCBI_N_23.csv")
