# seasonal_tick_odes

##DATA
#NEON tick abundance data provided in original downloaded format ("NEON_count-ticks.zip") and a cleaned csv format ("NEON_tick_abundance.csv). 

##Zip files for temperature (i.e., HARV.zip) and relative humidity (i.e., HARV_RH.zip) provided but need to be unzipped to run in code. 

##CODE
#fit_HDQ_model estimates parameters for full model 
#noeffect_predictions, prior_preictions, and posterior_predictions generate predictions for each model 
#DIC_HDQ_model ad DIC_nohumidit_model calculate DIC for the full humidity dependent and temperature only model, respectively. 
##likelihood_HDQ_no_priors calculates marginal parameter likelihood profiles independent of the prior distributions. 
##Rank_corr_analysis compares model predictions against data for all 3 model 
##Sync_analysis calculates Weibull distributions for empirical phenology data (and uncertainty) and calculates best estimates of empirical synchrony, then compares these to estimates of sychrony in data. 

