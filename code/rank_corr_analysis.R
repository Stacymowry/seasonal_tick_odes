##rank corerlation with unc 

##MAKE PLOT OF DATA VS PRIOR VS POSTERIOR VS NO EFFECT ESTIMATES

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(lubridate)
library(neonUtilities)
library(neonOS)
library(stringr)
library(purrr)



##load in NEON data and sort
options(stringsAsFactors=F)

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

#calculate the area sampled on each day
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
BLAN_a<-cbind(BLAN_DOY, BLAN_year, BLAN_area[,2],BLAN_larvae_pc,BLAN_nymphs_pc,BLAN_adults_pc,BLAN_larvae[,2],BLAN_nymphs[,2],BLAN_adults[,2])
colnames(BLAN_a)<-c("DOY","year","area","larvae_pc","nymphs_pc","adults_pc", "larvae", "nymphs", "adults")
BLAN_a<-as.data.frame(BLAN_a)

#Cchange format into long form
BLAN_a_larv<-BLAN_a[,c(1,2,3,4,7)]
BLAN_a_larv$stage<-rep("L",nrow(BLAN_a_larv))
colnames(BLAN_a_larv)<-c("DOY","year","area","ticks_pc","ticks","stage")

BLAN_a_nymphs<-BLAN_a[,c(1,2,3,5,8)]
BLAN_a_nymphs$stage<-rep("N",nrow(BLAN_a_nymphs))
colnames(BLAN_a_nymphs)<-c("DOY","year","area","ticks_pc","ticks","stage")

BLAN_a_adult<-BLAN_a[,c(1,2,3,6,9)]
BLAN_a_adult$stage<-rep("A",nrow(BLAN_a_adult))
colnames(BLAN_a_adult)<-c("DOY","year","area","ticks_pc","ticks","stage")

#get rid of adult rows 
BLAN_a<-as.data.frame(rbind(BLAN_a_larv,BLAN_a_nymphs,BLAN_a_adult))
BLAN_a<-subset(BLAN_a,BLAN_a$stage != "A")


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
HARV_a<-cbind(HARV_DOY, HARV_year, HARV_area[,2],HARV_larvae_pc,HARV_nymphs_pc,HARV_adults_pc,HARV_larvae[,2],HARV_nymphs[,2],HARV_adults[,2])
colnames(HARV_a)<-c("DOY","year","area","larvae_pc","nymphs_pc","adults_pc","larvae","nymphs","adults")
HARV_a<-as.data.frame(HARV_a)


HARV_a_larv<-HARV_a[,c(1,2,3,4,7)]
HARV_a_larv$stage<-rep("L",nrow(HARV_a_larv))
colnames(HARV_a_larv)<-c("DOY","year","area","ticks_pc","ticks","stage")
HARV_a_nymphs<-HARV_a[,c(1,2,3,5,8)]
HARV_a_nymphs$stage<-rep("N",nrow(HARV_a_nymphs))
colnames(HARV_a_nymphs)<-c("DOY","year","area","ticks_pc","ticks","stage")
HARV_a_adult<-HARV_a[,c(1,2,3,6,9)]
HARV_a_adult$stage<-rep("A",nrow(HARV_a_adult))
colnames(HARV_a_adult)<-c("DOY","year","area","ticks_pc","ticks","stage")
HARV_a<-as.data.frame(rbind(HARV_a_larv,HARV_a_nymphs,HARV_a_adult))
HARV_a<-subset(HARV_a,HARV_a$stage != "A")


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
SCBI_a<-cbind(SCBI_DOY, SCBI_year, SCBI_area[,2],SCBI_larvae_pc,SCBI_nymphs_pc,SCBI_adults_pc, SCBI_larvae[,2],SCBI_nymphs[,2],SCBI_adults[,2])
colnames(SCBI_a)<-c("DOY","year","area","larvae_pc","nymphs_pc","adults_pc", "larvae","nymphs","adults")
SCBI_a<-as.data.frame(SCBI_a)
SCBI_a_larv<-SCBI_a[,c(1,2,3,4,7)]
SCBI_a_larv$stage<-rep("L",nrow(SCBI_a_larv))
colnames(SCBI_a_larv)<-c("DOY","year","area","ticks_pc","ticks","stage")
SCBI_a_nymphs<-SCBI_a[,c(1,2,3,5,8)]
SCBI_a_nymphs$stage<-rep("N",nrow(SCBI_a_nymphs))
colnames(SCBI_a_nymphs)<-c("DOY","year","area","ticks_pc","ticks","stage")
SCBI_a_adult<-SCBI_a[,c(1,2,3,6,9)]
SCBI_a_adult$stage<-rep("A",nrow(SCBI_a_adult))
colnames(SCBI_a_adult)<-c("DOY","year","area","ticks_pc","ticks","stage")
SCBI_a<-as.data.frame(rbind(SCBI_a_larv,SCBI_a_nymphs,SCBI_a_adult))
SCBI_a<-subset(SCBI_a,SCBI_a$stage != "A")

##Subset for each year

BLAN_a_2016<-subset(BLAN_a,BLAN_a$year == 2016)
BLAN_a_2017<-subset(BLAN_a,BLAN_a$year == 2017)
BLAN_a_2018<-subset(BLAN_a,BLAN_a$year == 2018)
BLAN_a_2019<-subset(BLAN_a,BLAN_a$year == 2019)


HARV_a_2017<-subset(HARV_a,HARV_a$year == 2017)
HARV_a_2018<-subset(HARV_a,HARV_a$year == 2018)
HARV_a_2021<-subset(HARV_a,HARV_a$year == 2021)
HARV_a_2023<-subset(HARV_a,HARV_a$year == 2023)

SCBI_a_2016<-subset(SCBI_a,SCBI_a$year == 2016)
SCBI_a_2017<-subset(SCBI_a,SCBI_a$year == 2017)
SCBI_a_2018<-subset(SCBI_a,SCBI_a$year == 2018)
SCBI_a_2019<-subset(SCBI_a,SCBI_a$year == 2019)
SCBI_a_2023<-subset(SCBI_a,SCBI_a$year == 2023)

##format data for fitting
BLAN_a_2016 <- BLAN_a_2016 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )



BLAN_a_2016$time<- ((365*7) + BLAN_a_2016$DOY - 1)


BLAN_a_2017 <- BLAN_a_2017 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )


BLAN_a_2017$time<- ((365*7) + BLAN_a_2017$DOY - 1)

BLAN_a_2018 <- BLAN_a_2018 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )


BLAN_a_2018$time<- ((365*7) + BLAN_a_2018$DOY - 1)


BLAN_a_2019 <- BLAN_a_2019 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )


BLAN_a_2019$time<- ((365*7) + BLAN_a_2019$DOY - 1)


##HARV


HARV_a_2017 <- HARV_a_2017 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )


HARV_a_2017$time<- ((365*7) + HARV_a_2017$DOY - 1)

HARV_a_2018 <- HARV_a_2018 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )



HARV_a_2018$time<- ((365*7) + HARV_a_2018$DOY - 1)


HARV_a_2021 <- HARV_a_2021 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )



HARV_a_2021$time<- ((365*7) + HARV_a_2021$DOY - 1)


HARV_a_2023 <- HARV_a_2023 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )



HARV_a_2023$time<- ((365*7) + HARV_a_2023$DOY - 1)



###SCBI
SCBI_a_2016 <- SCBI_a_2016 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )


SCBI_a_2016$time<- ((365*7) + SCBI_a_2016$DOY - 1)

SCBI_a_2017 <- SCBI_a_2017 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )


SCBI_a_2017$time<- ((365*7) + SCBI_a_2017$DOY - 1)

SCBI_a_2018 <- SCBI_a_2018 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )


SCBI_a_2018$time<- ((365*7) + SCBI_a_2018$DOY - 1)

SCBI_a_2019 <- SCBI_a_2019 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )


SCBI_a_2019$time<- ((365*7) + SCBI_a_2019$DOY - 1)

SCBI_a_2023 <- SCBI_a_2023 %>%
  pivot_wider(
    id_cols = c(DOY, year, area),
    names_from = stage,
    values_from = ticks,
    names_prefix = "ticks_"
  )


SCBI_a_2023$time<- ((365*7) + SCBI_a_2023$DOY - 1)


##sort data by DOY

BLAN_a_2016 <- BLAN_a_2016[order(BLAN_a_2016$DOY), ]
BLAN_a_2017 <- BLAN_a_2017[order(BLAN_a_2017$DOY), ]
BLAN_a_2018 <- BLAN_a_2018[order(BLAN_a_2018$DOY), ]
BLAN_a_2019 <- BLAN_a_2019[order(BLAN_a_2019$DOY), ]



HARV_a_2017 <- HARV_a_2017[order(HARV_a_2017$DOY), ]
HARV_a_2018 <- HARV_a_2018[order(HARV_a_2018$DOY), ]
HARV_a_2021 <- HARV_a_2021[order(HARV_a_2021$DOY), ]
HARV_a_2023 <- HARV_a_2023[order(HARV_a_2023$DOY), ]

SCBI_a_2016 <- SCBI_a_2016[order(SCBI_a_2016$DOY), ]
SCBI_a_2017 <- SCBI_a_2017[order(SCBI_a_2017$DOY), ]
SCBI_a_2018 <- SCBI_a_2018[order(SCBI_a_2018$DOY), ]
SCBI_a_2019 <- SCBI_a_2019[order(SCBI_a_2019$DOY), ]
SCBI_a_2023 <- SCBI_a_2018[order(SCBI_a_2023$DOY), ]


##read in prior estimates to add to plot

# Helper function to read and trim duplicated rows
read_half <- function(path) {
  data <- read.csv(path)
  # Keep rows 1 through (Total Rows / 2)
  return(data[1:(nrow(data)/2), ])
}


BLAN_16L_data_from_prior<-read_half("/annual_prior_preds_BLAN_L_16.csv")
BLAN_17L_data_from_prior<-read_half("/annual_prior_preds_BLAN_L_17.csv")
BLAN_18L_data_from_prior<-read_half("annual_prior_preds_BLAN_L_18.csv")
BLAN_19L_data_from_prior<-read_half("/annual_prior_preds_BLAN_L_19.csv")

BLAN_16N_data_from_prior<-read_half("/annual_prior_preds_BLAN_N_16.csv")
BLAN_17N_data_from_prior<-read_half("/annual_prior_preds_BLAN_N_17.csv")
BLAN_18N_data_from_prior<-read_half("/annual_prior_preds_BLAN_N_18.csv")
BLAN_19N_data_from_prior<-read_half("/annual_prior_preds_BLAN_N_19.csv")


HARV_17L_data_from_prior<-read_half("/annual_prior_preds_HARV_L_17.csv")
HARV_18L_data_from_prior<-read_half("/annual_prior_preds_HARV_L_18.csv")
HARV_21L_data_from_prior<-read_half("/annual_prior_preds_HARV_L_21.csv")
HARV_23L_data_from_prior<-read_half("/annual_prior_preds_HARV_L_23.csv")

HARV_17N_data_from_prior<-read_half("/annual_prior_preds_HARV_N_17.csv")
HARV_18N_data_from_prior<-read_half("/annual_prior_preds_HARV_N_18.csv")
HARV_21N_data_from_prior<-read_half("/annual_prior_preds_HARV_N_21.csv")
HARV_23N_data_from_prior<-read_half("/annual_prior_preds_HARV_N_23.csv")

SCBI_16L_data_from_prior<-read_half("/annual_prior_preds_SCBI_L_16.csv")
SCBI_17L_data_from_prior<-read_half("/annual_prior_preds_SCBI_L_17.csv")
SCBI_18L_data_from_prior<-read_half("/annual_prior_preds_SCBI_L_18.csv")
SCBI_19L_data_from_prior<-read_half("/annual_prior_preds_SCBI_L_19.csv")
SCBI_23L_data_from_prior<-read_half("/annual_prior_preds_SCBI_L_23.csv")

SCBI_16N_data_from_prior<-read_half("/annual_prior_preds_SCBI_N_16.csv")
SCBI_17N_data_from_prior<-read_half("/annual_prior_preds_SCBI_N_17.csv")
SCBI_18N_data_from_prior<-read_half("/annual_prior_preds_SCBI_N_18.csv")
SCBI_19N_data_from_prior<-read_half("/annual_prior_preds_SCBI_N_19.csv")
SCBI_23N_data_from_prior<-read_half("/annual_prior_preds_SCBI_N_23.csv")


##read in posterior predictions
BLAN_16L_data_from_posterior<-read_half("/annual_posterior_preds_BLAN_L_16.csv")
BLAN_17L_data_from_posterior<-read_half("/annual_posterior_preds_BLAN_L_17.csv")
BLAN_18L_data_from_posterior<-read_half("/annual_posterior_preds_BLAN_L_18.csv")
BLAN_19L_data_from_posterior<-read_half("/annual_posterior_preds_BLAN_L_19.csv")

BLAN_16N_data_from_posterior<-read_half("/annual_posterior_preds_BLAN_N_16.csv")
BLAN_17N_data_from_posterior<-read_half("/annual_posterior_preds_BLAN_N_17.csv")
BLAN_18N_data_from_posterior<-read_half("/annual_posterior_preds_BLAN_N_18.csv")
BLAN_19N_data_from_posterior<-read_half("/annual_posterior_preds_BLAN_N_19.csv")


HARV_17L_data_from_posterior<-read_half("/annual_posterior_preds_HARV_L_17.csv")
HARV_18L_data_from_posterior<-read_half("/annual_posterior_preds_HARV_L_18.csv")
HARV_21L_data_from_posterior<-read_half("/annual_posterior_preds_HARV_L_21.csv")
HARV_23L_data_from_posterior<-read_half("/annual_posterior_preds_HARV_L_23.csv")

HARV_17N_data_from_posterior<-read_half("/annual_posterior_preds_HARV_N_17.csv")
HARV_18N_data_from_posterior<-read_half("/annual_posterior_preds_HARV_N_18.csv")
HARV_21N_data_from_posterior<-read_half("/annual_posterior_preds_HARV_N_21.csv")
HARV_23N_data_from_posterior<-read_half("/annual_posterior_preds_HARV_N_23.csv")

SCBI_16L_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_L_16.csv")
SCBI_17L_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_L_17.csv")
SCBI_18L_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_L_18.csv")
SCBI_19L_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_L_19.csv")
SCBI_23L_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_L_23.csv")

SCBI_16N_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_N_16.csv")
SCBI_17N_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_N_17.csv")
SCBI_18N_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_N_18.csv")
SCBI_19N_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_N_19.csv")
SCBI_23N_data_from_posterior<-read_half("/annual_posterior_preds_SCBI_N_23.csv")


##read in no effect predictions
BLAN_16L_data_from_NE<-read_half("/annual_NE_preds_BLAN_L_16.csv")
BLAN_17L_data_from_NE<-read_half("/annual_NE_preds_BLAN_L_17.csv")
BLAN_18L_data_from_NE<-read_half("/annual_NE_preds_BLAN_L_18.csv")
BLAN_19L_data_from_NE<-read_half("/annual_NE_preds_BLAN_L_19.csv")

BLAN_16N_data_from_NE<-read_half("/annual_NE_preds_BLAN_N_16.csv")
BLAN_17N_data_from_NE<-read_half("/annual_NE_preds_BLAN_N_17.csv")
BLAN_18N_data_from_NE<-read_half("/annual_NE_preds_BLAN_N_18.csv")
BLAN_19N_data_from_NE<-read_half("/annual_NE_preds_BLAN_N_19.csv")


HARV_17L_data_from_NE<-read_half("/annual_NE_preds_HARV_L_17.csv")
HARV_18L_data_from_NE<-read_half("/annual_NE_preds_HARV_L_18.csv")
HARV_21L_data_from_NE<-read_half("/annual_NE_preds_HARV_L_21.csv")
HARV_23L_data_from_NE<-read_half("/annual_NE_preds_HARV_L_23.csv")

HARV_17N_data_from_NE<-read_half("/annual_NE_preds_HARV_N_17.csv")
HARV_18N_data_from_NE<-read_half("/annual_NE_preds_HARV_N_18.csv")
HARV_21N_data_from_NE<-read_half("/annual_NE_preds_HARV_N_21.csv")
HARV_23N_data_from_NE<-read_half("/annual_NE_preds_HARV_N_23.csv")

SCBI_16L_data_from_NE<-read_half("/annual_NE_preds_SCBI_L_16.csv")
SCBI_17L_data_from_NE<-read_half("/annual_NE_preds_SCBI_L_17.csv")
SCBI_18L_data_from_NE<-read_half("/annual_NE_preds_SCBI_L_18.csv")
SCBI_19L_data_from_NE<-read_half("/annual_NE_preds_SCBI_L_19.csv")
SCBI_23L_data_from_NE<-read_half("/annual_NE_preds_SCBI_L_23.csv")

SCBI_16N_data_from_NE<-read_half("/annual_NE_preds_SCBI_N_16.csv")
SCBI_17N_data_from_NE<-read_half("/annual_NE_preds_SCBI_N_17.csv")
SCBI_18N_data_from_NE<-read_half("/annual_NE_preds_SCBI_N_18.csv")
SCBI_19N_data_from_NE<-read_half("/annual_NE_preds_SCBI_N_19.csv")
SCBI_23N_data_from_NE<-read_half("/annual_NE_preds_SCBI_N_23.csv")


# -------------------------
# 1. Identify prediction objects
# -------------------------
objs <- ls(pattern = "^[A-Z]+_[0-9]{2}[LN]_data_from_(posterior|prior|NE)$")

# Parse names into metadata table
meta <- tibble(name = objs) %>%
  mutate(
    SITE = str_extract(name, "^[A-Z]+"),
    YR   = str_extract(name, "(?<=_)[0-9]{2}"),
    L    = str_extract(name, "(?<=[0-9]{2})[LN]"),
    model = str_extract(name, "(posterior|prior|NE)$")
  )


# -------------------------
# 2. Function to compute correlations for one object
# -------------------------
get_corrs <- function(name, SITE, YR, L, model){
  
  preds <- get(name)
  obsname <- paste0(SITE, "_a_20", YR)
  obsdf <- get(obsname)
  
  obs <- if (L == "L") obsdf$ticks_L else obsdf$ticks_N
  
  cors <- apply(preds, 2, function(col)
    cor(col, obs, method="spearman", use="complete.obs")
  )
  
  tibble(
    SITE = SITE,
    YR = YR,
    L = L,
    model = model,
    cors = list(cors)
  )
}

allcors <- pmap_dfr(meta, get_corrs)


