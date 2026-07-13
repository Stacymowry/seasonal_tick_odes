library(MASS) # For mvrnorm

# 1. LOAD DATA
# ==============================================================================
# Ensure this path is correct for your machine
x <- read.csv('/NEON_tick_abundance.csv', header=TRUE)

# ==============================================================================
# UPDATED PLOTTING CODE: WITH SCALED RAW DATA POINTS
# ==============================================================================
fit_tick_model <- function(doy, counts, areas) {
  
  if(sum(counts > 0) < 3) return(NULL)
  
  nll_poisson <- function(par) {
    C <- par[1]; k <- par[2]; w_scale <- par[3]
    mu <- areas * C * dweibull(doy, shape = k, scale = w_scale)
    mu[mu < 1e-9] <- 1e-9
    if(any(is.infinite(mu)) || any(is.na(mu))) return(1e12)
    -sum(dpois(counts, lambda = mu, log = TRUE))
  }
  
  guess_k <- 4.0; guess_w_scale <- 210
  start_C <- sum(counts) / sum(areas * dweibull(doy, guess_k, guess_w_scale))
  if(start_C <= 0) start_C <- 100
  
  tryCatch({
    # Fit Once
    opt <- optim(par = c(start_C, guess_k, guess_w_scale), fn = nll_poisson, 
                 method = "L-BFGS-B", lower = c(0.1, 1.1, 50), upper = c(1e8, 30, 400),
                 hessian = TRUE)
    
    fisher_info <- solve(opt$hessian)
    best_pars   <- opt$par            
    
    # SAMPLE 100 TIMES
    n_samples <- 100
    param_samples <- mvrnorm(n = n_samples, mu = best_pars, Sigma = fisher_info)
    
    # Generate 100 Curves (1:365)
    pred_doy  <- 1:365
    posterior_curves <- matrix(NA, nrow = n_samples, ncol = 365)
    
    for(i in 1:n_samples) {
      p <- param_samples[i, ]
      raw_curve <- p[1] * dweibull(pred_doy, shape = p[2], scale = p[3])
      
      # Normalize to Sum = 1
      if(sum(raw_curve) > 0) {
        posterior_curves[i, ] <- raw_curve / sum(raw_curve)
      } else {
        posterior_curves[i, ] <- rep(0, 365)
      }
    }
    
    # --- RESTORED CODE START ---
    # Calculate the mean of the 100 curves for plotting
    mean_curve <- colMeans(posterior_curves, na.rm = TRUE)
    # --- RESTORED CODE END ---
    
    return(list(
      success = TRUE, 
      posterior_curves = posterior_curves,
      mean = mean_curve  # <--- Now available for your plotting code
    ))
    
  }, error = function(e) { return(NULL) })
}
# B. Efficient Synchrony Function
# Calculates Jaccard-style overlap: Intersection / Union
calc_sync_vec <- function(curve_L, curve_N) {
  # Intersection = Area where both curves overlap
  intersection <- sum(pmin(curve_L, curve_N))
  
  # Union = AreaL + AreaN - Intersection
  # Since curves sum to 1, Union is approx (1 + 1 - intersection)
  area_L <- sum(curve_L)
  area_N <- sum(curve_N)
  union  <- area_L + area_N - intersection
  
  # Handle division by zero edge case
  if(union == 0) return(0)
  
  return(intersection / union)
}


# 1. HELPER FUNCTION
calc_sync_vec <- function(curve_L, curve_N) {
  intersection <- sum(pmin(curve_L, curve_N))
  area_L <- sum(curve_L)
  area_N <- sum(curve_N)
  union  <- area_L + area_N - intersection
  if(union == 0) return(0)
  return(intersection / union)
}

# 2. SETUP GRID
target_sites <- c("BLAN", "HARV", "SCBI")
target_years <- c(2016, 2017, 2018, 2019, 2021, 2023)

excluded_pairs <- c(
  "HARV_2016", "HARV_2019", "HARV_2022", "HARV_2023",
  "BLAN_2021", "BLAN_2022", "BLAN_2023",
  "SCBI_2021", "SCBI_2022"
)

par(mfrow = c(length(target_years), length(target_sites)), 
    mar = c(0.5, 0.5, 0.5, 0.5), 
    oma = c(4, 5, 3, 1)) # Bottom, Left, Top, Right

# ==============================================================================
# ADDED: Initialize a list to capture synchrony values across iterations
# ==============================================================================
synchrony_storage <- list()

# 3. LOOP THROUGH GRID
for(y in target_years) {
  for(s in target_sites) {
    
    current_id <- paste(s, y, sep = "_")
    is_excluded <- current_id %in% excluded_pairs
    
    fit_L <- NULL
    fit_N <- NULL
    sub_data <- subset(x, site == s & year == y & area > 0 & !is.na(DOY))
    
    # Fit Models (if not excluded)
    if(!is_excluded && nrow(sub_data) > 0) {
      if(nrow(subset(sub_data, !is.na(larvae))) > 0) {
        fit_L <- fit_tick_model(sub_data$DOY[!is.na(sub_data$larvae)], 
                                sub_data$larvae[!is.na(sub_data$larvae)], 
                                sub_data$area[!is.na(sub_data$larvae)])
      }
      if(nrow(subset(sub_data, !is.na(nymphs))) > 0) {
        fit_N <- fit_tick_model(sub_data$DOY[!is.na(sub_data$nymphs)], 
                                sub_data$nymphs[!is.na(sub_data$nymphs)], 
                                sub_data$area[!is.na(sub_data$nymphs)])
      }
    }
    
    # 4. DRAW PLOT
    if(!is_excluded && !is.null(fit_L) && !is.null(fit_N)) {
      
      # Calculate Synchrony
      sync_values <- numeric(100)
      for(k in 1:100) {
        sync_values[k] <- calc_sync_vec(fit_L$posterior_curves[k,], 
                                        fit_N$posterior_curves[k,])
      }
      mean_sync <- mean(sync_values, na.rm=TRUE)
      
      # ==============================================================================
      # ADDED: Save the 100 uncertainty values for this specific site and year
      # ==============================================================================
      synchrony_storage[[current_id]] <- data.frame(
        site = s,
        year = y,
        iteration = 1:100,
        synchrony = sync_values
      )
      
      # Determine Max Y (Curves)
      max_val_L <- max(fit_L$mean, na.rm=TRUE)
      max_val_N <- max(fit_N$mean, na.rm=TRUE)
      y_max <- max(c(max_val_L, max_val_N)) * 1.2
      if(y_max == 0) y_max <- 0.02
      
      # Plot Setup
      plot(1, type="n", xlim=c(1, 365), ylim=c(0, y_max), 
           xaxt="n", yaxt="n", xlab="", ylab="")
      
      # Grid & Labels
      abline(v = seq(50, 350, 50), col = "lightgray", lty = 3)
      if(y == max(target_years)) axis(1, cex.axis=0.8)        
      if(s == target_sites[1])   axis(2, las=1, cex.axis=0.8) 
      if(y == target_years[1]) mtext(s, side=3, line=0.2, cex=0.9, font=2)
      if(s == target_sites[1]) mtext(y, side=2, line=3.5, cex=0.9, font=2, las=0)
      
      # --- PLOT CURVES ---
      # Nymphs (Purple)
      for(k in 1:100) lines(1:365, fit_N$posterior_curves[k, ], col=rgb(0.13, 0.55, 0.13, 0.1), lwd=0.25)
      lines(1:365, fit_N$mean, col="forestgreen", lwd=1) 
      
      # Larvae (Orange)
      for(k in 1:100) lines(1:365, fit_L$posterior_curves[k, ], col=rgb(0, 0, 0, 0.1), lwd=0.25)
      lines(1:365, fit_L$mean, col="black", lwd=1.5, lty=1) 
      
      # --- PLOT SCALED DATA POINTS ---
      
      # Nymphs Points
      n_data <- subset(sub_data, !is.na(nymphs) & nymphs > 0)
      if(nrow(n_data) > 0) {
        # Scaling Factor: Max Curve / Max Data
        scale_N <- max_val_N / max(n_data$nymphs, na.rm=TRUE)
        points(n_data$DOY, n_data$nymphs * scale_N, 
               pch=21, bg="forestgreen", col="black", cex=1)
      }
      
      # Larvae Points
      l_data <- subset(sub_data, !is.na(larvae) & larvae > 0)
      if(nrow(l_data) > 0) {
        # Scaling Factor: Max Curve / Max Data
        scale_L <- max_val_L / max(l_data$larvae, na.rm=TRUE)
        points(l_data$DOY, l_data$larvae * scale_L, 
               pch=24, bg="black", col="black", cex=1) # pch=24 is a triangle
      }
      
      # Shade Overlap
      overlap_curve <- pmin(fit_L$mean, fit_N$mean)
      polygon(c(1:365, 365:1), c(overlap_curve, rep(0, 365)), 
              col = rgb(0.2, 0.2, 0.2, 0.3), border = NA)
      
      # Synchrony Legend
      # legend("topleft", legend=sprintf("%.2f", mean_sync), 
      #    bty="n", cex=1, inset=0.02)
      
    } else {
      # Blank Panel
      plot(1, type="n", axes=FALSE, xlab="", ylab="")
      if(y == target_years[1]) mtext(s, side=3, line=0.2, cex=0.9, font=2)
      if(s == target_sites[1]) mtext(y, side=2, line=3.5, cex=0.9, font=2, las=0)
    }
  }
}

# Global Labels
mtext("Day of Year", side=1, outer=TRUE, line=2.5, cex=1)

# ==============================================================================
# ADDED: Compile all stored lists into a master table and save to disk
# ==============================================================================
if (length(synchrony_storage) > 0) {
  master_synchrony_df <- do.call(rbind, synchrony_storage)
  write.csv(master_synchrony_df, "~/Desktop/weibull_synchrony_uncertainty.csv", row.names = FALSE)
  message("Success: Synchrony uncertainty data frames compiled and saved to '~/Desktop/weibull_synchrony_uncertainty.csv'!")
}


##load in prior, posterior and NE mats
no<-read.csv("/noeffect_all_mat_nosamplewindow.csv")


##load in all estimates
prior_B_L_16<-read.csv("/B_L_16_mat_prior_unc_nosamplewindow.csv")
prior_B_N_16<-read.csv("/B_N_16_mat_prior_unc_nosamplewindow.csv")
prior_B_L_17<-read.csv("/B_L_17_mat_prior_unc_nosamplewindow.csv")
prior_B_N_17<-read.csv("/B_N_17_mat_prior_unc_nosamplewindow.csv")
prior_B_L_18<-read.csv("/B_L_18_mat_prior_unc_nosamplewindow.csv")
prior_B_N_18<-read.csv("/B_N_18_mat_prior_unc_nosamplewindow.csv")
prior_B_L_19<-read.csv("/B_L_19_mat_prior_unc_nosamplewindow.csv")
prior_B_N_19<-read.csv("/B_N_19_mat_prior_unc_nosamplewindow.csv")


prior_H_L_17<-read.csv("/H_L_17_mat_prior_unc_nosamplewindow.csv")
prior_H_N_17<-read.csv("/H_N_17_mat_prior_unc_nosamplewindow.csv")
prior_H_L_18<-read.csv("/H_L_18_mat_prior_unc_nosamplewindow.csv")
prior_H_N_18<-read.csv("/H_N_18_mat_prior_unc_nosamplewindow.csv")
prior_H_L_21<-read.csv("/H_L_21_mat_prior_unc_nosamplewindow.csv")
prior_H_N_21<-read.csv("/H_N_21_mat_prior_unc_nosamplewindow.csv")
prior_H_L_23<-read.csv("/H_L_23_mat_prior_unc_nosamplewindow.csv")
prior_H_N_23<-read.csv("/H_N_23_mat_prior_unc_nosamplewindow.csv")

prior_S_L_16<-read.csv("/S_L_16_mat_prior_unc_nosamplewindow.csv")
prior_S_N_16<-read.csv("/S_N_16_mat_prior_unc_nosamplewindow.csv")
prior_S_L_17<-read.csv("/S_L_17_mat_prior_unc_nosamplewindow.csv")
prior_S_N_17<-read.csv("/S_N_17_mat_prior_unc_nosamplewindow.csv")
prior_S_L_18<-read.csv("/S_L_18_mat_prior_unc_nosamplewindow.csv")
prior_S_N_18<-read.csv("/S_N_18_mat_prior_unc_nosamplewindow.csv")
prior_S_L_19<-read.csv("/S_L_19_mat_prior_unc_nosamplewindow.csv")
prior_S_N_19<-read.csv("/S_N_19_mat_prior_unc_nosamplewindow.csv")
prior_S_L_23<-read.csv("/S_L_23_mat_prior_unc_nosamplewindow.csv")
prior_S_N_23<-read.csv("/S_N_23_mat_prior_unc_nosamplewindow.csv")


posterior_B_L_16<-read.csv("/B_L_16_mat_posterior_unc_nosamplewindow.csv")
posterior_B_N_16<-read.csv("/B_N_16_mat_posterior_unc_nosamplewindow.csv")
posterior_B_L_17<-read.csv("/B_L_17_mat_posterior_unc_nosamplewindow.csv")
posterior_B_N_17<-read.csv("/B_N_17_mat_posterior_unc_nosamplewindow.csv")
posterior_B_L_18<-read.csv("/B_L_18_mat_posterior_unc_nosamplewindow.csv")
posterior_B_N_18<-read.csv("/B_N_18_mat_posterior_unc_nosamplewindow.csv")
posterior_B_L_19<-read.csv("/B_L_19_mat_posterior_unc_nosamplewindow.csv")
posterior_B_N_19<-read.csv("/B_N_19_mat_posterior_unc_nosamplewindow.csv")


posterior_H_L_17<-read.csv("/H_L_17_mat_posterior_unc_nosamplewindow.csv")
posterior_H_N_17<-read.csv("/H_N_17_mat_posterior_unc_nosamplewindow.csv")
posterior_H_L_18<-read.csv("/H_L_18_mat_posterior_unc_nosamplewindow.csv")
posterior_H_N_18<-read.csv("/H_N_18_mat_posterior_unc_nosamplewindow.csv")
posterior_H_L_21<-read.csv("/H_L_21_mat_posterior_unc_nosamplewindow.csv")
posterior_H_N_21<-read.csv("/H_N_21_mat_posterior_unc_nosamplewindow.csv")
posterior_H_L_23<-read.csv("/H_L_23_mat_posterior_unc_nosamplewindow.csv")
posterior_H_N_23<-read.csv("/H_N_23_mat_posterior_unc_nosamplewindow.csv")

posterior_S_L_16<-read.csv("/S_L_16_mat_posterior_unc_nosamplewindow.csv")
posterior_S_N_16<-read.csv("/S_N_16_mat_posterior_unc_nosamplewindow.csv")
posterior_S_L_17<-read.csv("/S_L_17_mat_posterior_unc_nosamplewindow.csv")
posterior_S_N_17<-read.csv("/S_N_17_mat_posterior_unc_nosamplewindow.csv")
posterior_S_L_18<-read.csv("/S_L_18_mat_posterior_unc_nosamplewindow.csv")
posterior_S_N_18<-read.csv("/S_N_18_mat_posterior_unc_nosamplewindow.csv")
posterior_S_L_19<-read.csv("/S_L_19_mat_posterior_unc_nosamplewindow.csv")
posterior_S_N_19<-read.csv("/S_N_19_mat_posterior_unc_nosamplewindow.csv")
posterior_S_L_23<-read.csv("/S_L_23_mat_posterior_unc_nosamplewindow.csv")
posterior_S_N_23<-read.csv("/S_N_23_mat_posterior_unc_nosamplewindow.csv")



# 1. Grab all object names
model_names <- ls(pattern = "^(prior|posterior)_")

# 1. Combine and filter HARV 2023 out of model_df
model_df <- map_dfr(model_names, function(name) {
  obj <- get(name)
  parts <- str_split(name, "_")[[1]]
  
  obj %>%
    as.data.frame() %>%
    rename(DOY_col = 1) %>%
    pivot_longer(cols = -DOY_col, names_to = "iteration", values_to = "value") %>%
    mutate(
      DOY = DOY_col,
      model = parts[1],
      SITE  = case_when(parts[2] == "B" ~ "BLAN", parts[2] == "S" ~ "SCBI", parts[2] == "H" ~ "HARV"),
      L     = parts[3],
      YR    = paste0("20", parts[4])
    ) %>%
    dplyr::select(-DOY_col)
}) %>% 
  # DIRECT FIX: Remove HARV 2023
  filter(!(SITE == "HARV" & YR == "2023"))

# 2. Prepare and filter HARV 2023 out of x_long
x_long <- x %>%
  dplyr::select(SITE = site, YR = year, DOY, larvae_perarea, nymphs_perarea) %>%
  mutate(YR = as.character(YR)) %>%
  pivot_longer(cols = ends_with("perarea"), names_to = "L", values_to = "obs_density") %>%
  mutate(L = if_else(str_detect(L, "larvae"), "L", "N")) %>%
  # DIRECT FIX: Remove HARV 2023
  filter(!(SITE == "HARV" & YR == "2023"))






# ... [Your previous data loading and model_df / x_long creation code] ...

# 1. CALCULATE THE MEDIAN CURVE
# Group by day, site, year, life stage, and model to get the median across all iterations
median_df <- model_df %>%
  group_by(SITE, YR, L, model, DOY) %>%
  summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")

# 2. UPDATE SCALING FACTORS
# Find the maximum value of the *median curve* instead of the absolute max of all curves
scaling_factors <- median_df %>%
  group_by(SITE, YR, L, model) %>%
  summarise(max_median = max(median_value, na.rm = TRUE), .groups = "drop")

# 3. SCALE RAW OBSERVATIONS
# Scale the data points so the max observation equals the max of the median curve
x_scaled <- x_long %>%
  inner_join(scaling_factors, by = c("SITE", "YR", "L")) %>%
  group_by(SITE, YR, L, model) %>%
  mutate(scaled_obs = (obs_density / max(obs_density, na.rm = TRUE)) * max_median) %>%
  ungroup()

# 4. UPDATED PLOTTING FUNCTION
plot_site_eval <- function(target_site) {
  site_models <- model_df %>% filter(SITE == target_site)
  site_median <- median_df %>% filter(SITE == target_site) # Grab median data
  site_obs    <- x_scaled %>% filter(SITE == target_site)
  
  ggplot() +
    # Faint lines for all posterior/prior iterations
    geom_line(data = site_models, 
              aes(x = DOY, y = value, group = interaction(iteration, L), color = L), 
              alpha = 0.1) +
    # BOLD line for the median curve
    geom_line(data = site_median,
              aes(x = DOY, y = median_value, group = L, color = L),
              linewidth = 1.2) + # use linewidth (or size if on an older ggplot2 version)
    # Scaled observation points
    geom_point(data = site_obs, 
               aes(x = DOY, y = scaled_obs, color = L), 
               size = 1.5) +
    facet_grid(model ~ YR, scales = "free_y") +
    scale_color_manual(values = c("L" = "forestgreen", "N" = "black"),
                       labels = c("L" = "Larvae", "N" = "Nymphs")) +
    theme_bw() +
    labs(title = target_site, y = "Tick Activity", x = "Day of Year") +
    theme(legend.position = "bottom")
}

# Run the plot
plot_site_eval("HARV")












# --- Proceed with scaling (no changes needed to these steps) ---

scaling_factors <- model_df %>%
  group_by(SITE, YR, L, model) %>%
  summarise(max_pred = max(value, na.rm = TRUE), .groups = "drop")

x_scaled <- x_long %>%
  inner_join(scaling_factors, by = c("SITE", "YR", "L")) %>%
  group_by(SITE, YR, L, model) %>%
  mutate(scaled_obs = (obs_density / max(obs_density, na.rm = TRUE)) * max_pred) %>%
  ungroup()

# --- Proceed with plotting function (no changes needed) ---

plot_site_eval <- function(target_site) {
  site_models <- model_df %>% filter(SITE == target_site)
  site_obs    <- x_scaled %>% filter(SITE == target_site)
  
  ggplot() +
    geom_line(data = site_models, 
              aes(x = DOY, y = value, group = interaction(iteration, L), color = L), 
              alpha = 0.1) +
    geom_point(data = site_obs, 
               aes(x = DOY, y = scaled_obs, color = L), 
               size = 1.5) +
    facet_grid(model ~ YR, scales = "free_y") +
    scale_color_manual(values = c("L" = "forestgreen", "N" = "black"),
                       labels = c("L" = "Larvae", "N" = "Nymphs")) +
    theme_bw() +
    labs(title = target_site, y = "Tick Activity", x = "Day of Year") +
    theme(legend.position = "bottom")
}

plot_site_eval("HARV")