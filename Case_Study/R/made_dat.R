####### Make data ########
# 1. remove all counties with less that 40K people
# 2. select covs (while removing other indexes of mobility) and rescale weeks
# 3. create dummy for state
# 4. Lag data by two weeks
# 5. remove late time points with a lot of missing data
# 6. Impute data using prior history with single imputation
# 7. scale data
# 8. Select random dups
# 9. Switch to wide format
# 10. Save data
##########################

rm(list= ls())
# load packages #
library(here)
library(dplyr)
library(mice)
library(fastDummies)
library(tidyr)
library(tibble)
# load data #
load(here("Case_Study/data/2021-01-25_covid_wk.RData"))


# 1. select covs and rescale weeks
dat <- dfwk %>%
  filter(population >= 40000, !is.na(dex_a)) %>% # remove all counties with less that 40K people
  select(
    # id, week, exposure, outcome #
    week, fips, dex_a, cum_cases_percap,
    new_deaths_1wks_ahead_percap, new_deaths_2wks_ahead_percap, new_deaths_3wks_ahead_percap, new_deaths_4wks_ahead_percap,
    new_cases_1wks_behind_percap,
    
    # covariates geography #
    state_abb, log_pop_dens, log_pop, pct_rural_2010, air_pollution,
    
    # pop covariates #
    ppl_below_poverty, hh_below_poverty, pct_black, pct_hispanic, pct_asian,
    pct_under_18, pct_65plus, some_college, pct_female, inc_ineq, pct_dem_2016,
    meatpack_cat4, meatpack_cat5, meatpack_4or5, drive_alone_to_work, prison,
    pct_crowded_housing, svi_ses, avg_hh_size, unemployment,
    
    # health covariates #
    uninsured_adults, adult_smoke_rate, log_med_hh_inc, obesity_pct, diabetes_pct, mask_level
  ) %>%
  mutate(week = week - 20,
         dex_a = dex_a,
         cum_cases_percap = cum_cases_percap
  )




# 2. create dummy for state
dat_dum <- dat %>%
  dummy_cols(select_columns = "state_abb", remove_most_frequent_dummy = F,
             remove_first_dummy = TRUE,
             remove_selected_columns = TRUE)



# 3. don't lag data
lag_dat_2 <- dat_dum

# Impute data by sequential time point #
weeks <- unique(lag_dat_2$week)

# Initialize an empty list to store imputed datasets for each week
imputed_list <- list()
temp_to_imp <- lag_dat_2
# Loop over weeks, imputing each sequentially using data up to that week
for(i in weeks) {
  print(i)
  # Subset the data up to the current week
  subset_dat <- temp_to_imp %>% filter(week <= i)
  
  # Perform mice imputation on the subset
  temp_imputed <- mice(subset_dat, m = 1, method = 'pmm', maxit = 1, seed = 500)
  
  # Store the completed data
  imputed_list[[i]] <- complete(temp_imputed, 1)
  temp_to_imp <- rbind(imputed_list[[i]], temp_to_imp %>% filter(week > i))
}
dat_comp <- temp_to_imp


# covariates baseline vs time varying #
baseline_vars <- time_vars <- c()

for(i in unique(dat_comp$fips)){
  tmp <- dat_comp %>%
    filter(fips == i) %>%
    summarise(across(.cols = -week, .fns = ~ n_distinct(.)))
  
  baseline_vars <- c(baseline_vars, colnames(tmp)[which(tmp == 1)])
  time_vars <- c(time_vars, colnames(tmp)[which(tmp > 1)])
}

time_vars <- unique(time_vars)
baseline_vars <- unique(baseline_vars)
baseline_vars <- baseline_vars[!baseline_vars %in% time_vars]

confounder_vars <- list("time_vars"=time_vars, "baseline_vars"=baseline_vars)
#save(confounder_vars, file = here("data","confounder_vars.RData"))


# 6. Select random dups
# keep only fips ids that are in all weeks #
keep_ops <- unique(dat_comp$fips)
for(i in unique(dat_comp$week)){
  keep_ops <- intersect(keep_ops,unique(dat_comp$fips[dat_comp$week == i]))
}
dat_comp <- dat_comp %>% filter(fips %in% keep_ops)

# keep one of each week at random... #
df_random_sample <- dat_comp %>%
  group_by(fips, week) %>%
  sample_n(size = 1) %>%
  ungroup()



# 8. Switch to wide format
# Switch to wide format #
# baseline_vars <- c("meatpack_cat5","meatpack_cat4", "ppl_below_poverty", "prison", "log_pop", "log_pop_dens",
#                    colnames(df_random_sample)[grepl("state_abb", colnames(df_random_sample))])
dat_wide <- df_random_sample %>%
  pivot_wider(
    id_cols = fips,
    names_from = week,
    names_sep = "_",
    values_from = colnames(df_random_sample)[-c(1,2)]
    # values_from = colnames(df_random_sample)[-c(c(1,2),match(baseline_vars, colnames(df_random_sample)))], unused_fn = match(baseline_vars, colnames(df_random_sample))
  ) %>% ungroup()

# 9. Save data
data <- dat_wide
save(data, file = here("Case_Study/data","paper_data_full.RData"))
