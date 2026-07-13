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

##load in parameter estimates
someyrs_10k_out<-readRDS("model_estimates.rds")
someyrs_10k_pars <- getSample(someyrs_10k_out, parametersOnly = TRUE)


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




 

##LOAD TICK ABUNDANCE DATA

stackByTable("NEON_count-ticks.zip")

tckabun <- readTableNEON(
  dataFile="NEON_count-ticks/stackedFiles/tck_fielddata.csv",
  varFile="NEON_count-ticks/stackedFiles/variables_10093.csv")


ticktaxon<-read.csv("NEON_count-ticks/stackedFiles/tck_taxonomyProcessed.csv")


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
#For each day calculate: total nymphs, total larvae, total adults, total area flagged

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




BLAN_a_2016$time<- ((365*7) + BLAN_a_2016$DOY - 1)

BLAN_a_2017$time<- ((365*7) + BLAN_a_2017$DOY - 1)

BLAN_a_2019$time<- ((365*7) + BLAN_a_2019$DOY - 1)

BLAN_a_2019$time<- ((365*7) + BLAN_a_2019$DOY - 1)

BLAN_a_2020$time<- ((365*7) + BLAN_a_2020$DOY - 1)

BLAN_a_2021$time<- ((365*7) + BLAN_a_2021$DOY - 1)

BLAN_a_2022$time<- ((365*7) + BLAN_a_2022$DOY - 1)

BLAN_a_2023$time<- ((365*7) + BLAN_a_2023$DOY - 1)

BLAN_a_2024$time<- ((365*7) + BLAN_a_2024$DOY - 1)

# ##HARV

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
BLAN<-read.csv("BLAN.csv")



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
HARV<-read.csv("~/HARV.csv")

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
SCBI<-read.csv("~/SCBI.csv")


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
BLAN_humidity<-read.csv("~/BLAN_RH.csv")


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
HARV_humidity<-read.csv("~/HARV_RH.csv")


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


##calculate normalizing parameters for questing

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




##ODIN MODELS
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


##creates infrastructure to interat with model in dust2
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
                   
                   
                   
                   # ##no effect
                   # log_dln_BLAN = log(0.00001),
                   # log_dln_SCBI =log(0.00001),
                   # log_dln_HARV =log(0.00001),
                   # log_aa_BLAN = log(0.00001),
                   # log_aa_HARV = log(0.00001),
                   # log_aa_SCBI = log(0.00001),
                   # log_al_BLAN = log(0.00001),
                   # log_al_HARV = log(0.00001),
                   # log_al_SCBI = log(0.00001),
                   # # 
                   # a_n_BLAN = 0,
                   # a_n_HARV = 0,
                   # a_n_SCBI =0
                   # 
                   
                   ##prior effect
                   #  log_dln_BLAN = log(0.00001),
                   #  log_dln_SCBI =log(0.00001),
                   #  log_dln_HARV =log(0.00001),
                   #  log_aa_BLAN = log(1),
                   #  log_aa_HARV = log(1),
                   #  log_aa_SCBI = log(1),
                   #  log_al_BLAN = log(1),
                   #  log_al_HARV = log(1),
                   #  log_al_SCBI = log(1),
                   # #
                   #  a_n_BLAN = 30.5,
                   #  a_n_HARV = 30.5,
                   #  a_n_SCBI =30.5
                   
                   
                   
                   log_dln_BLAN = median(someyrs_10k_pars[,"log_dln_BLAN"]),
                   log_dln_SCBI = median(someyrs_10k_pars[,"log_dln_SCBI"]),
                   log_dln_HARV = median(someyrs_10k_pars[,"log_dln_HARV"]),
                   log_aa_BLAN = median(someyrs_10k_pars[,"log_aa_BLAN"]),
                   log_aa_HARV = median(someyrs_10k_pars[,"log_aa_HARV"]),
                   log_aa_SCBI = median(someyrs_10k_pars[,"log_aa_SCBI"]),
                   log_al_BLAN = median(someyrs_10k_pars[,"log_al_BLAN"]),
                   log_al_HARV = median(someyrs_10k_pars[,"log_al_HARV"]),
                   log_al_SCBI = median(someyrs_10k_pars[,"log_al_SCBI"]),
                   
                   a_n_BLAN = median(someyrs_10k_pars[,"a_n_BLAN"]),
                   a_n_HARV = median(someyrs_10k_pars[,"a_n_HARV"]),
                   a_n_SCBI = median(someyrs_10k_pars[,"a_n_SCBI"])
                   
                   
                   
                   
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

# ##plot model outputs
years<-list(
  BLAN= c(2016,2017,2018,2019),
  HARV = c(2017, 2018, 2021, 2023),
  SCBI = c(2016, 2017, 2018, 2019, 2023))



idx_BLAN_2016 <- BLAN_a_2016$DOY
idx_BLAN_2017 <- BLAN_a_2017$DOY
idx_BLAN_2018 <- BLAN_a_2018$DOY
idx_BLAN_2019 <- BLAN_a_2019$DOY


yL_BLAN_2016  <- subset(BLAN_a_2016, BLAN_a_2016$stage == "L")$ticks
yN_BLAN_2016  <- subset(BLAN_a_2016, BLAN_a_2016$stage == "N")$ticks
yL_BLAN_2017  <- subset(BLAN_a_2017, BLAN_a_2017$stage == "L")$ticks
yN_BLAN_2017  <- subset(BLAN_a_2017, BLAN_a_2017$stage == "N")$ticks
yL_BLAN_2018  <- subset(BLAN_a_2018, BLAN_a_2018$stage == "L")$ticks
yN_BLAN_2018  <- subset(BLAN_a_2018, BLAN_a_2018$stage == "N")$ticks
yL_BLAN_2019  <- subset(BLAN_a_2019, BLAN_a_2019$stage == "L")$ticks
yN_BLAN_2019  <- subset(BLAN_a_2019, BLAN_a_2019$stage == "N")$ticks

ef_BLAN_2016  <- BLAN_a_2016$area
ef_BLAN_2017  <- BLAN_a_2017$area
ef_BLAN_2018  <- BLAN_a_2018$area
ef_BLAN_2019  <- BLAN_a_2019$area


idx_HARV_2017 <- HARV_a_2017$DOY
idx_HARV_2018 <- HARV_a_2018$DOY
idx_HARV_2021 <- HARV_a_2021$DOY
idx_HARV_2023 <- HARV_a_2023$DOY


yL_HARV_2017  <- subset(HARV_a_2017, HARV_a_2017$stage == "L")$ticks
yN_HARV_2017  <- subset(HARV_a_2017, HARV_a_2017$stage == "N")$ticks
yL_HARV_2018  <- subset(HARV_a_2018, HARV_a_2018$stage == "L")$ticks
yN_HARV_2018  <- subset(HARV_a_2018, HARV_a_2018$stage == "N")$ticks
yL_HARV_2021  <- subset(HARV_a_2021, HARV_a_2021$stage == "L")$ticks
yN_HARV_2021  <- subset(HARV_a_2021, HARV_a_2021$stage == "N")$ticks
yL_HARV_2023  <- subset(HARV_a_2023, HARV_a_2023$stage == "L")$ticks
yN_HARV_2023  <- subset(HARV_a_2023, HARV_a_2023$stage == "N")$ticks


ef_HARV_2017  <- HARV_a_2017$area
ef_HARV_2018  <- HARV_a_2018$area
ef_HARV_2021  <- HARV_a_2021$area
ef_HARV_2023  <- HARV_a_2023$area

idx_SCBI_2016 <- SCBI_a_2016$DOY
idx_SCBI_2017 <- SCBI_a_2017$DOY
idx_SCBI_2018 <- SCBI_a_2018$DOY
idx_SCBI_2019 <- SCBI_a_2019$DOY
idx_SCBI_2023 <- SCBI_a_2023$DOY

yL_SCBI_2016  <- subset(SCBI_a_2016, SCBI_a_2016$stage == "L")$ticks
yN_SCBI_2016  <- subset(SCBI_a_2016, SCBI_a_2016$stage == "N")$ticks
yL_SCBI_2017  <- subset(SCBI_a_2017, SCBI_a_2017$stage == "L")$ticks
yN_SCBI_2017  <- subset(SCBI_a_2017, SCBI_a_2017$stage == "N")$ticks
yL_SCBI_2018  <- subset(SCBI_a_2018, SCBI_a_2018$stage == "L")$ticks
yN_SCBI_2018  <- subset(SCBI_a_2018, SCBI_a_2018$stage == "N")$ticks
yL_SCBI_2019  <- subset(SCBI_a_2019, SCBI_a_2019$stage == "L")$ticks
yN_SCBI_2019  <- subset(SCBI_a_2019, SCBI_a_2019$stage == "N")$ticks
yL_SCBI_2023  <- subset(SCBI_a_2023, SCBI_a_2023$stage == "L")$ticks
yN_SCBI_2023  <- subset(SCBI_a_2023, SCBI_a_2023$stage == "N")$ticks

ef_SCBI_2016  <- SCBI_a_2016$area
ef_SCBI_2017  <- SCBI_a_2017$area
ef_SCBI_2018  <- SCBI_a_2018$area
ef_SCBI_2019  <- SCBI_a_2019$area
ef_SCBI_2023  <- SCBI_a_2023$area


##likelihood 

#parameter vector
param_names <- c(
  "log_dln_bar","log_dln_sigma",
  "a_n_bar","a_n_sigma",
  "log_al_bar","log_al_sigma",
  "log_aa_bar","log_aa_sigma",
  # BLAN
  "log_dln_BLAN","a_n_BLAN","log_al_BLAN","log_aa_BLAN",
  "log_alpha0_L_BLAN_2016","log_alpha0_N_BLAN_2016",
  "log_alpha0_L_BLAN_2017","log_alpha0_N_BLAN_2017",
  "log_alpha0_L_BLAN_2018","log_alpha0_N_BLAN_2018",
  "log_alpha0_L_BLAN_2019","log_alpha0_N_BLAN_2019",
  
  # HARV
  "log_dln_HARV","a_n_HARV","log_al_HARV","log_aa_HARV",
  "log_alpha0_L_HARV_2017","log_alpha0_N_HARV_2017",
  "log_alpha0_L_HARV_2018","log_alpha0_N_HARV_2018",
  "log_alpha0_L_HARV_2021","log_alpha0_N_HARV_2021",
  "log_alpha0_L_HARV_2023","log_alpha0_N_HARV_2023",
  # SCBI
  "log_dln_SCBI","a_n_SCBI","log_al_SCBI","log_aa_SCBI",
  "log_alpha0_L_SCBI_2016","log_alpha0_N_SCBI_2016",
  "log_alpha0_L_SCBI_2017","log_alpha0_N_SCBI_2017",
  "log_alpha0_L_SCBI_2018","log_alpha0_N_SCBI_2018",
  "log_alpha0_L_SCBI_2019","log_alpha0_N_SCBI_2019",
  "log_alpha0_L_SCBI_2023","log_alpha0_N_SCBI_2023"
)

## ----- Hyperprior constants (from your DSL) -----
m_log_dln_bar <- 0.8
s_log_dln_bar <- 0.7
rate_log_dln_sigma <- 1/0.7

m_a_n_bar <- 30.5
s_a_n_bar <- 2.0
rate_a_n_sigma <- 1/2

m_log_al_bar <- -0.11
s_log_al_bar <- 0.48
rate_log_al_sigma <- 1/0.5

m_log_aa_bar <- -0.11
s_log_aa_bar <- 0.48
rate_log_aa_sigma <- 1/0.5

m_log_alpha0 <- 2;
s_log_alpha0 <- 1   # all log-alpha0 priors



fn <- function(par) {
  names(par) <- param_names
  
  # pull parameters
  
  
  log_dln_BLAN <- par["log_dln_BLAN"]
  a_n_BLAN     <- par["a_n_BLAN"]
  log_al_BLAN  <- par["log_al_BLAN"]
  log_aa_BLAN  <- par["log_aa_BLAN"]
  
  log_alpha0_L_BLAN_2016 <- par["log_alpha0_L_BLAN_2016"]
  log_alpha0_N_BLAN_2016 <- par["log_alpha0_N_BLAN_2016"]
  log_alpha0_L_BLAN_2017 <- par["log_alpha0_L_BLAN_2017"]
  log_alpha0_N_BLAN_2017 <- par["log_alpha0_N_BLAN_2017"]
  log_alpha0_L_BLAN_2018 <- par["log_alpha0_L_BLAN_2018"]
  log_alpha0_N_BLAN_2018 <- par["log_alpha0_N_BLAN_2018"]
  log_alpha0_L_BLAN_2019 <- par["log_alpha0_L_BLAN_2019"]
  log_alpha0_N_BLAN_2019 <- par["log_alpha0_N_BLAN_2019"]
  
  
  log_dln_HARV <- par["log_dln_HARV"]
  a_n_HARV     <- par["a_n_HARV"]
  log_al_HARV  <- par["log_al_HARV"]
  log_aa_HARV  <- par["log_aa_HARV"]
  
  
  log_alpha0_L_HARV_2017 <- par["log_alpha0_L_HARV_2017"]
  log_alpha0_N_HARV_2017 <- par["log_alpha0_N_HARV_2017"]
  log_alpha0_L_HARV_2018 <- par["log_alpha0_L_HARV_2018"]
  log_alpha0_N_HARV_2018 <- par["log_alpha0_N_HARV_2018"]
  log_alpha0_L_HARV_2021 <- par["log_alpha0_L_HARV_2021"]
  log_alpha0_N_HARV_2021 <- par["log_alpha0_N_HARV_2021"]
  log_alpha0_L_HARV_2023 <- par["log_alpha0_L_HARV_2023"]
  log_alpha0_N_HARV_2023 <- par["log_alpha0_N_HARV_2023"]
  
  
  log_dln_SCBI <- par["log_dln_SCBI"]
  a_n_SCBI     <- par["a_n_SCBI"]
  log_al_SCBI  <- par["log_al_SCBI"]
  log_aa_SCBI  <- par["log_aa_SCBI"]
  
  log_alpha0_L_SCBI_2016 <- par["log_alpha0_L_SCBI_2016"]
  log_alpha0_N_SCBI_2016 <- par["log_alpha0_N_SCBI_2016"]
  log_alpha0_L_SCBI_2017 <- par["log_alpha0_L_SCBI_2017"]
  log_alpha0_N_SCBI_2017 <- par["log_alpha0_N_SCBI_2017"]
  log_alpha0_L_SCBI_2018 <- par["log_alpha0_L_SCBI_2018"]
  log_alpha0_N_SCBI_2018 <- par["log_alpha0_N_SCBI_2018"]
  log_alpha0_L_SCBI_2019 <- par["log_alpha0_L_SCBI_2019"]
  log_alpha0_N_SCBI_2019 <- par["log_alpha0_N_SCBI_2019"]
  log_alpha0_L_SCBI_2023 <- par["log_alpha0_L_SCBI_2023"]
  log_alpha0_N_SCBI_2023 <- par["log_alpha0_N_SCBI_2023"]
  
  
  full_pars <- list(
    theta_i_norm = theta_i_norm, theta_n_norm = theta_n_norm, theta_a_norm = theta_a_norm,
    
    
    
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
    
    
    
    DOY_BLAN_2016 = DOY_BLAN_2016, temperature_BLAN_2016 = temperature_BLAN_2016,
    DOY_BLAN_2017 = DOY_BLAN_2017, temperature_BLAN_2017 = temperature_BLAN_2017,
    DOY_BLAN_2018 = DOY_BLAN_2018, temperature_BLAN_2018 = temperature_BLAN_2018,
    DOY_BLAN_2019 = DOY_BLAN_2019, temperature_BLAN_2019 = temperature_BLAN_2019,
    
    log_dln_BLAN = log_dln_BLAN, 
    a_n_BLAN = a_n_BLAN,
    log_aa_BLAN = log_aa_BLAN,
    log_al_BLAN = log_al_BLAN,
    
    
    
    
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
    
    
    DOY_HARV_2017 = DOY_HARV_2017, temperature_HARV_2017 = temperature_HARV_2017,
    DOY_HARV_2018 = DOY_HARV_2018, temperature_HARV_2018 = temperature_HARV_2018,
    DOY_HARV_2021 = DOY_HARV_2021, temperature_HARV_2021 = temperature_HARV_2021,
    DOY_HARV_2023 = DOY_HARV_2023, temperature_HARV_2023 = temperature_HARV_2023,
    
    log_dln_HARV = log_dln_HARV, a_n_HARV = a_n_HARV,
    log_aa_HARV = log_aa_HARV, log_al_HARV = log_al_HARV,
    
    
    
    
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
    
    
    DOY_SCBI_2016 = DOY_SCBI_2016, temperature_SCBI_2016 = temperature_SCBI_2016,
    DOY_SCBI_2017 = DOY_SCBI_2017, temperature_SCBI_2017 = temperature_SCBI_2017,
    DOY_SCBI_2018 = DOY_SCBI_2018, temperature_SCBI_2018 = temperature_SCBI_2018,
    DOY_SCBI_2019 = DOY_SCBI_2019, temperature_SCBI_2019 = temperature_SCBI_2019,
    DOY_SCBI_2023 = DOY_SCBI_2023, temperature_SCBI_2023 = temperature_SCBI_2023,
    
    log_dln_SCBI = log_dln_SCBI, a_n_SCBI = a_n_SCBI,
    log_aa_SCBI = log_aa_SCBI, log_al_SCBI = log_al_SCBI
  )
  
  # simulate odin model
  total <- 365*10 - 1
  yr <- 7
  i <- yr*365
  t <- 0:(total-1)
  ode_ctrl <- dust_ode_control(max_steps = 1e7, step_size_min = 1e-9)
  sys <- dust_system_create(full_odin(), pars = full_pars, n_particles = 1, ode_control = ode_ctrl)
  dust_system_set_state_initial(sys)
  y <- dust_system_simulate(sys, t)
  
  # normalized model predictions
  
  BLAN_N16 <- as.vector(dust_unpack_state(sys, y)$FNabove_BLAN_2016[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_BLAN_2016[i:(i+364)]))[idx_BLAN_2016]
  BLAN_L16 <- as.vector(dust_unpack_state(sys, y)$FL_BLAN_2016[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_BLAN_2016[i:(i+364)]))[idx_BLAN_2016]
  BLAN_N17 <- as.vector(dust_unpack_state(sys, y)$FNabove_BLAN_2017[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_BLAN_2017[i:(i+364)]))[idx_BLAN_2017]
  BLAN_L17 <- as.vector(dust_unpack_state(sys, y)$FL_BLAN_2017[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_BLAN_2017[i:(i+364)]))[idx_BLAN_2017]
  BLAN_N18 <- as.vector(dust_unpack_state(sys, y)$FNabove_BLAN_2018[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_BLAN_2018[i:(i+364)]))[idx_BLAN_2018]
  BLAN_L18 <- as.vector(dust_unpack_state(sys, y)$FL_BLAN_2018[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_BLAN_2018[i:(i+364)]))[idx_BLAN_2018]
  BLAN_N19 <- as.vector(dust_unpack_state(sys, y)$FNabove_BLAN_2019[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_BLAN_2019[i:(i+364)]))[idx_BLAN_2019]
  BLAN_L19 <- as.vector(dust_unpack_state(sys, y)$FL_BLAN_2019[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_BLAN_2019[i:(i+364)]))[idx_BLAN_2019]
  
  
  
  HARV_N17 <- as.vector(dust_unpack_state(sys, y)$FNabove_HARV_2017[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_HARV_2017[i:(i+364)]))[idx_HARV_2017]
  HARV_L17 <- as.vector(dust_unpack_state(sys, y)$FL_HARV_2017[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_HARV_2017[i:(i+364)]))[idx_HARV_2017]
  HARV_N18 <- as.vector(dust_unpack_state(sys, y)$FNabove_HARV_2018[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_HARV_2018[i:(i+364)]))[idx_HARV_2018]
  HARV_L18 <- as.vector(dust_unpack_state(sys, y)$FL_HARV_2018[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_HARV_2018[i:(i+364)]))[idx_HARV_2018]
  HARV_N21 <- as.vector(dust_unpack_state(sys, y)$FNabove_HARV_2021[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_HARV_2021[i:(i+364)]))[idx_HARV_2021]
  HARV_L21 <- as.vector(dust_unpack_state(sys, y)$FL_HARV_2021[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_HARV_2021[i:(i+364)]))[idx_HARV_2021]
  HARV_N23 <- as.vector(dust_unpack_state(sys, y)$FNabove_HARV_2023[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_HARV_2023[i:(i+364)]))[idx_HARV_2023]
  HARV_L23 <- as.vector(dust_unpack_state(sys, y)$FL_HARV_2023[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_HARV_2023[i:(i+364)]))[idx_HARV_2023]
  
  SCBI_N16 <- as.vector(dust_unpack_state(sys, y)$FNabove_SCBI_2016[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_SCBI_2016[i:(i+364)]))[idx_SCBI_2016]
  SCBI_L16 <- as.vector(dust_unpack_state(sys, y)$FL_SCBI_2016[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_SCBI_2016[i:(i+364)]))[idx_SCBI_2016]
  SCBI_N17 <- as.vector(dust_unpack_state(sys, y)$FNabove_SCBI_2017[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_SCBI_2017[i:(i+364)]))[idx_SCBI_2017]
  SCBI_L17 <- as.vector(dust_unpack_state(sys, y)$FL_SCBI_2017[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_SCBI_2017[i:(i+364)]))[idx_SCBI_2017]
  SCBI_N18 <- as.vector(dust_unpack_state(sys, y)$FNabove_SCBI_2018[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_SCBI_2018[i:(i+364)]))[idx_SCBI_2018]
  SCBI_L18 <- as.vector(dust_unpack_state(sys, y)$FL_SCBI_2018[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_SCBI_2018[i:(i+364)]))[idx_SCBI_2018]
  SCBI_N19 <- as.vector(dust_unpack_state(sys, y)$FNabove_SCBI_2019[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_SCBI_2019[i:(i+364)]))[idx_SCBI_2019]
  SCBI_L19 <- as.vector(dust_unpack_state(sys, y)$FL_SCBI_2019[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_SCBI_2019[i:(i+364)]))[idx_SCBI_2019]
  SCBI_N23 <- as.vector(dust_unpack_state(sys, y)$FNabove_SCBI_2023[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FNabove_SCBI_2023[i:(i+364)]))[idx_SCBI_2023]
  SCBI_L23 <- as.vector(dust_unpack_state(sys, y)$FL_SCBI_2023[i:(i+364)] / AUC(t[i:(i+364)], dust_unpack_state(sys, y)$FL_SCBI_2023[i:(i+364)]))[idx_SCBI_2023]
  
  # normalize based on effort 
  
  BLAN_pL16 <- BLAN_L16 * ef_BLAN_2016 / max(ef_BLAN_2016)
  BLAN_pN16 <- BLAN_N16 * ef_BLAN_2016 / max(ef_BLAN_2016)
  BLAN_pL17 <- BLAN_L17 * ef_BLAN_2017 / max(ef_BLAN_2017)
  BLAN_pN17 <- BLAN_N17 * ef_BLAN_2017 / max(ef_BLAN_2017)
  BLAN_pL18 <- BLAN_L18 * ef_BLAN_2018 / max(ef_BLAN_2018)
  BLAN_pN18 <- BLAN_N18 * ef_BLAN_2018 / max(ef_BLAN_2018)
  BLAN_pL19 <- BLAN_L19 * ef_BLAN_2019 / max(ef_BLAN_2019)
  BLAN_pN19 <- BLAN_N19 * ef_BLAN_2019 / max(ef_BLAN_2019)
  
  HARV_pL17 <- HARV_L17 * ef_HARV_2017 / max(ef_HARV_2017)
  HARV_pN17 <- HARV_N17 * ef_HARV_2017 / max(ef_HARV_2017)
  HARV_pL18 <- HARV_L18 * ef_HARV_2018 / max(ef_HARV_2018)
  HARV_pN18 <- HARV_N18 * ef_HARV_2018 / max(ef_HARV_2018)
  HARV_pL21 <- HARV_L21 * ef_HARV_2021 / max(ef_HARV_2021)
  HARV_pN21 <- HARV_N21 * ef_HARV_2021 / max(ef_HARV_2021)
  HARV_pL23 <- HARV_L23 * ef_HARV_2023 / max(ef_HARV_2023)
  HARV_pN23 <- HARV_N23 * ef_HARV_2023 / max(ef_HARV_2023)
  
  SCBI_pL16 <- SCBI_L16 * ef_SCBI_2016 / max(ef_SCBI_2016)
  SCBI_pN16 <- SCBI_N16 * ef_SCBI_2016 / max(ef_SCBI_2016)
  SCBI_pL17 <- SCBI_L17 * ef_SCBI_2017 / max(ef_SCBI_2017)
  SCBI_pN17 <- SCBI_N17 * ef_SCBI_2017 / max(ef_SCBI_2017)
  SCBI_pL18 <- SCBI_L18 * ef_SCBI_2018 / max(ef_SCBI_2018)
  SCBI_pN18 <- SCBI_N18 * ef_SCBI_2018 / max(ef_SCBI_2018)
  SCBI_pL19 <- SCBI_L19 * ef_SCBI_2019 / max(ef_SCBI_2019)
  SCBI_pN19 <- SCBI_N19 * ef_SCBI_2019 / max(ef_SCBI_2019)
  SCBI_pL23 <- SCBI_L23 * ef_SCBI_2023 / max(ef_SCBI_2023)
  SCBI_pN23 <- SCBI_N23 * ef_SCBI_2023 / max(ef_SCBI_2023)
  
  #add in alpha for DM
  
  BLAN_aL16 <- (BLAN_pL16 / sum(BLAN_pL16)) * exp(log_alpha0_L_BLAN_2016)
  BLAN_aN16 <- (BLAN_pN16 / sum(BLAN_pN16)) * exp(log_alpha0_N_BLAN_2016)
  BLAN_aL17 <- (BLAN_pL17 / sum(BLAN_pL17)) * exp(log_alpha0_L_BLAN_2017)
  BLAN_aN17 <- (BLAN_pN17 / sum(BLAN_pN17)) * exp(log_alpha0_N_BLAN_2017)
  BLAN_aL18 <- (BLAN_pL18 / sum(BLAN_pL18)) * exp(log_alpha0_L_BLAN_2018)
  BLAN_aN18 <- (BLAN_pN18 / sum(BLAN_pN18)) * exp(log_alpha0_N_BLAN_2018)
  BLAN_aL19 <- (BLAN_pL19 / sum(BLAN_pL19)) * exp(log_alpha0_L_BLAN_2019)
  BLAN_aN19 <- (BLAN_pN19 / sum(BLAN_pN19)) * exp(log_alpha0_N_BLAN_2019)
  
  
  HARV_aL17 <- (HARV_pL17 / sum(HARV_pL17)) * exp(log_alpha0_L_HARV_2017)
  HARV_aN17 <- (HARV_pN17 / sum(HARV_pN17)) * exp(log_alpha0_N_HARV_2017)
  HARV_aL18 <- (HARV_pL18 / sum(HARV_pL18)) * exp(log_alpha0_L_HARV_2018)
  HARV_aN18 <- (HARV_pN18 / sum(HARV_pN18)) * exp(log_alpha0_N_HARV_2018)
  HARV_aL21 <- (HARV_pL21 / sum(HARV_pL21)) * exp(log_alpha0_L_HARV_2021)
  HARV_aN21 <- (HARV_pN21 / sum(HARV_pN21)) * exp(log_alpha0_N_HARV_2021)
  HARV_aL23 <- (HARV_pL23 / sum(HARV_pL23)) * exp(log_alpha0_L_HARV_2023)
  HARV_aN23 <- (HARV_pN23 / sum(HARV_pN23)) * exp(log_alpha0_N_HARV_2023)
  
  SCBI_aL16 <- (SCBI_pL16 / sum(SCBI_pL16)) * exp(log_alpha0_L_SCBI_2016)
  SCBI_aN16 <- (SCBI_pN16 / sum(SCBI_pN16)) * exp(log_alpha0_N_SCBI_2016)
  SCBI_aL17 <- (SCBI_pL17 / sum(SCBI_pL17)) * exp(log_alpha0_L_SCBI_2017)
  SCBI_aN17 <- (SCBI_pN17 / sum(SCBI_pN17)) * exp(log_alpha0_N_SCBI_2017)
  SCBI_aL18 <- (SCBI_pL18 / sum(SCBI_pL18)) * exp(log_alpha0_L_SCBI_2018)
  SCBI_aN18 <- (SCBI_pN18 / sum(SCBI_pN18)) * exp(log_alpha0_N_SCBI_2018)
  SCBI_aL19 <- (SCBI_pL19 / sum(SCBI_pL19)) * exp(log_alpha0_L_SCBI_2019)
  SCBI_aN19 <- (SCBI_pN19 / sum(SCBI_pN19)) * exp(log_alpha0_N_SCBI_2019)
  SCBI_aL23 <- (SCBI_pL23 / sum(SCBI_pL23)) * exp(log_alpha0_L_SCBI_2023)
  SCBI_aN23 <- (SCBI_pN23 / sum(SCBI_pN23)) * exp(log_alpha0_N_SCBI_2023)
  
  
  
  # log-likelihood (Dirichlet-multinomial)
  ll <- 0
  
  ll <- ll + extraDistr::ddirmnom(yL_BLAN_2016, sum(yL_BLAN_2016), BLAN_aL16, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_BLAN_2016, sum(yN_BLAN_2016), BLAN_aN16, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_BLAN_2017, sum(yL_BLAN_2017), BLAN_aL17, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_BLAN_2017, sum(yN_BLAN_2017), BLAN_aN17, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_BLAN_2018, sum(yL_BLAN_2018), BLAN_aL18, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_BLAN_2018, sum(yN_BLAN_2018), BLAN_aN18, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_BLAN_2019, sum(yL_BLAN_2019), BLAN_aL19, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_BLAN_2019, sum(yN_BLAN_2019), BLAN_aN19, log = TRUE)
  
  
  
  ll <- ll + extraDistr::ddirmnom(yL_HARV_2017, sum(yL_HARV_2017), HARV_aL17, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_HARV_2017, sum(yN_HARV_2017), HARV_aN17, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_HARV_2018, sum(yL_HARV_2018), HARV_aL18, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_HARV_2018, sum(yN_HARV_2018), HARV_aN18, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_HARV_2021, sum(yL_HARV_2021), HARV_aL21, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_HARV_2021, sum(yN_HARV_2021), HARV_aN21, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_HARV_2023, sum(yL_HARV_2023), HARV_aL23, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_HARV_2023, sum(yN_HARV_2023), HARV_aN23, log = TRUE)
  
  ll <- ll + extraDistr::ddirmnom(yL_SCBI_2016, sum(yL_SCBI_2016), SCBI_aL16, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_SCBI_2016, sum(yN_SCBI_2016), SCBI_aN16, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_SCBI_2017, sum(yL_SCBI_2017), SCBI_aL17, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_SCBI_2017, sum(yN_SCBI_2017), SCBI_aN17, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_SCBI_2018, sum(yL_SCBI_2018), SCBI_aL18, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_SCBI_2018, sum(yN_SCBI_2018), SCBI_aN18, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_SCBI_2019, sum(yL_SCBI_2019), SCBI_aL19, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_SCBI_2019, sum(yN_SCBI_2019), SCBI_aN19, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yL_SCBI_2023, sum(yL_SCBI_2023), SCBI_aL23, log = TRUE)
  ll <- ll + extraDistr::ddirmnom(yN_SCBI_2023, sum(yN_SCBI_2023), SCBI_aN23, log = TRUE)
  
  
  
  
  return(ll)
  
}



dprior <- function(par, log = TRUE){
  names(par) <- param_names
  
  # extract hypers
  log_dln_bar   <- par["log_dln_bar"]
  log_dln_sigma <- par["log_dln_sigma"]
  a_n_bar       <- par["a_n_bar"]
  a_n_sigma     <- par["a_n_sigma"]
  log_al_bar    <- par["log_al_bar"]
  log_al_sigma  <- par["log_al_sigma"]
  log_aa_bar    <- par["log_aa_bar"]
  log_aa_sigma  <- par["log_aa_sigma"]
  
  
  ll <- 0
  ## hyperpriors
  ll <- ll + dnorm(log_dln_bar, m_log_dln_bar, s_log_dln_bar, log=TRUE)
  ll <- ll + dexp( log_dln_sigma, rate = rate_log_dln_sigma, log=TRUE)
  
  ll <- ll + dnorm(a_n_bar, m_a_n_bar, s_a_n_bar, log=TRUE)
  ll <- ll + dexp( a_n_sigma, rate = rate_a_n_sigma, log=TRUE)
  
  ll <- ll + dnorm(log_al_bar, m_log_al_bar, s_log_al_bar, log=TRUE)
  ll <- ll + dexp( log_al_sigma, rate = rate_log_al_sigma, log=TRUE)
  
  ll <- ll + dnorm(log_aa_bar, m_log_aa_bar, s_log_aa_bar, log=TRUE)
  ll <- ll + dexp( log_aa_sigma, rate = rate_log_aa_sigma, log=TRUE)
  
  ## BLAN conditionals
  ll <- ll + dnorm(par["log_dln_BLAN"], mean=log_dln_bar, sd=log_dln_sigma, log=TRUE)
  ll <- ll + dnorm(par["a_n_BLAN"],     mean=a_n_bar,     sd=a_n_sigma,     log=TRUE)
  ll <- ll + dnorm(par["log_al_BLAN"],  mean=log_al_bar,  sd=log_al_sigma,  log=TRUE)
  ll <- ll + dnorm(par["log_aa_BLAN"],  mean=log_aa_bar,  sd=log_aa_sigma,  log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_BLAN_2016"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_BLAN_2016"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_BLAN_2017"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_BLAN_2017"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_BLAN_2018"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_BLAN_2018"], m_log_alpha0, s_log_alpha0, log=TRUE)
  
  ll <- ll + dnorm(par["log_alpha0_L_BLAN_2019"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_BLAN_2019"], m_log_alpha0, s_log_alpha0, log=TRUE)
  
  
  
  
  ## HARV conditionals
  ll <- ll + dnorm(par["log_dln_HARV"], mean=log_dln_bar, sd=log_dln_sigma, log=TRUE)
  ll <- ll + dnorm(par["a_n_HARV"],     mean=a_n_bar,     sd=a_n_sigma,     log=TRUE)
  ll <- ll + dnorm(par["log_al_HARV"],  mean=log_al_bar,  sd=log_al_sigma,  log=TRUE)
  ll <- ll + dnorm(par["log_aa_HARV"],  mean=log_aa_bar,  sd=log_aa_sigma,  log=TRUE)
  
  ll <- ll + dnorm(par["log_alpha0_L_HARV_2017"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_HARV_2017"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_HARV_2018"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_HARV_2018"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_HARV_2021"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_HARV_2021"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_HARV_2023"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_HARV_2023"], m_log_alpha0, s_log_alpha0, log=TRUE)
  
  
  
  ## SCBI conditionals
  ll <- ll + dnorm(par["log_dln_SCBI"], mean=log_dln_bar, sd=log_dln_sigma, log=TRUE)
  ll <- ll + dnorm(par["a_n_SCBI"],     mean=a_n_bar,     sd=a_n_sigma,     log=TRUE)
  ll <- ll + dnorm(par["log_al_SCBI"],  mean=log_al_bar,  sd=log_al_sigma,  log=TRUE)
  ll <- ll + dnorm(par["log_aa_SCBI"],  mean=log_aa_bar,  sd=log_aa_sigma,  log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_SCBI_2016"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_SCBI_2016"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_SCBI_2017"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_SCBI_2017"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_SCBI_2018"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_SCBI_2018"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_SCBI_2019"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_SCBI_2019"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_L_SCBI_2023"], m_log_alpha0, s_log_alpha0, log=TRUE)
  ll <- ll + dnorm(par["log_alpha0_N_SCBI_2023"], m_log_alpha0, s_log_alpha0, log=TRUE)
  
  
  if (log) ll else exp(ll)
}

## ----- Prior sampler (top-down) -----
rprior <- function(n = 1){
  L <- length(param_names)
  out <- matrix(NA_real_, nrow = n, ncol = L, dimnames = list(NULL, param_names))
  
  for (k in 1:n){
    # hypers
    log_dln_bar   <- rnorm(1, m_log_dln_bar, s_log_dln_bar)
    log_dln_sigma <- rexp(1, rate_log_dln_sigma)
    a_n_bar       <- rnorm(1, m_a_n_bar, s_a_n_bar)
    a_n_sigma     <- rexp(1, rate_a_n_sigma)
    log_al_bar    <- rnorm(1, m_log_al_bar, s_log_al_bar)
    log_al_sigma  <- rexp(1, rate_log_al_sigma)
    log_aa_bar    <- rnorm(1, m_log_aa_bar, s_log_aa_bar)
    log_aa_sigma  <- rexp(1, rate_log_aa_sigma)
    
    draw <- numeric(L); names(draw) <- param_names
    draw[c("log_dln_bar","log_dln_sigma","a_n_bar","a_n_sigma",
           "log_al_bar","log_al_sigma","log_aa_bar","log_aa_sigma")] <-
      c(log_dln_bar, log_dln_sigma, a_n_bar, a_n_sigma,
        log_al_bar, log_al_sigma, log_aa_bar, log_aa_sigma)
    
    # BLAN
    draw["log_dln_BLAN"] <- rnorm(1, log_dln_bar, log_dln_sigma)
    draw["a_n_BLAN"]     <- rnorm(1, a_n_bar,     a_n_sigma)
    draw["log_al_BLAN"]  <- rnorm(1, log_al_bar,  log_al_sigma)
    draw["log_aa_BLAN"]  <- rnorm(1, log_aa_bar,  log_aa_sigma)
    draw["log_alpha0_L_BLAN_2016"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_BLAN_2016"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_BLAN_2017"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_BLAN_2017"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_BLAN_2018"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_BLAN_2018"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_BLAN_2019"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_BLAN_2019"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    
    
    # HARV
    draw["log_dln_HARV"] <- rnorm(1, log_dln_bar, log_dln_sigma)
    draw["a_n_HARV"]     <- rnorm(1, a_n_bar,     a_n_sigma)
    draw["log_al_HARV"]  <- rnorm(1, log_al_bar,  log_al_sigma)
    draw["log_aa_HARV"]  <- rnorm(1, log_aa_bar,  log_aa_sigma)
    draw["log_alpha0_L_HARV_2017"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_HARV_2017"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_HARV_2018"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_HARV_2018"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_HARV_2021"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_HARV_2021"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_HARV_2023"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_HARV_2023"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    
    
    # SCBI
    draw["log_dln_SCBI"] <- rnorm(1, log_dln_bar, log_dln_sigma)
    draw["a_n_SCBI"]     <- rnorm(1, a_n_bar,     a_n_sigma)
    draw["log_al_SCBI"]  <- rnorm(1, log_al_bar,  log_al_sigma)
    draw["log_aa_SCBI"]  <- rnorm(1, log_aa_bar,  log_aa_sigma)
    draw["log_alpha0_L_SCBI_2016"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_SCBI_2016"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_SCBI_2017"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_SCBI_2017"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_SCBI_2018"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_SCBI_2018"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_SCBI_2019"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_SCBI_2019"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_L_SCBI_2023"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    draw["log_alpha0_N_SCBI_2023"] <- rnorm(1, m_log_alpha0, s_log_alpha0)
    
    
    out[k, ] <- draw
  }
  out
}

##set domains 
lower <- rep(-Inf, length(param_names))
upper <- rep( Inf, length(param_names))
names(lower) <- names(upper) <- param_names

# hyper-sd’s must be > 0
lower["log_dln_sigma"] <- 0
lower["a_n_sigma"]     <- 0
lower["log_al_sigma"]  <- 0
lower["log_aa_sigma"]  <- 0

# your site-level bounds (on the log / natural scales as specified)
upper["log_dln_BLAN"] <- log(18)
upper["log_dln_HARV"] <- log(18)
upper["log_dln_SCBI"] <- log(18)

lower["a_n_BLAN"] <- 22; upper["a_n_BLAN"] <- 42
lower["a_n_HARV"] <- 22; upper["a_n_HARV"] <- 42
lower["a_n_SCBI"] <- 22; upper["a_n_SCBI"] <- 42

upper["log_al_BLAN"] <- log(22)
upper["log_al_HARV"] <- log(22)
upper["log_al_SCBI"] <- log(22)

upper["log_aa_BLAN"] <- log(22)
upper["log_aa_HARV"] <- log(22)
upper["log_aa_SCBI"] <- log(22)




prior <- BayesianTools::createPrior(density = dprior, sampler = rprior,
                                    lower = lower, upper = upper, best = NULL)

bs  <- createBayesianSetup(likelihood = fn, prior = prior, names = param_names)

#run sampler
set <- list(iterations = 10000, thin = 1, burnin = 2000, nrChains = 3)   



##load in parameter estimates
someyrs_10k_out<-readRDS("model_estimates.rds")

setup<-someyrs_10k_out[[1]]$setup

# 1. Get the raw likelihood function
lik_func <- someyrs_10k_out[[1]]$setup$likelihood$density

# 2. Get your samples and calculate the mean of the parameters
samps <- getSample(someyrs_10k_out)
param_means <- colMeans(samps)

# 3. Calculate D_bar (Mean of the deviances)
# We take the log-likelihoods already stored in the object
logLiks <- getSample(someyrs_10k_out, parametersOnly = FALSE)[, "Llikelihood"]
D_bar <- mean(-2 * logLiks)

# 4. Calculate D_hat (Deviance at the mean)
# This is the step that was likely failing
logLik_at_mean <- lik_func(param_means)
D_hat <- -2 * logLik_at_mean

# 5. Final DIC assembly
pD <- D_bar - D_hat
manual_DIC <- D_bar + pD

# Print results
cat("D_bar:", D_bar, "\nD_hat:", D_hat, "\npD:", pD, "\nDIC:", manual_DIC)






# 3. Get your likelihood function from the model setup
# Note: You need the 'bayesianSetup' object you used to run the model
# Let's assume it's named 'setup'
likelihood_func <- setup$likelihood$density

# 4. Calculate log-likelihood at the parameter means
log_lik_at_mean <- likelihood_func(param_means)

# 5. Calculate D_hat (Deviance at the mean)
D_hat <- -2 * log_lik_at_mean

