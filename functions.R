#' Run a nimue model using online fits
create_vacc_fit <- function(country_example = NULL,
                            continue_rt = TRUE,
                            forecast = 180,
                            prop_immune = 0,
                            increase_rt = FALSE){
  
  iso3c <- countrycode(country_example, origin = 'country.name', destination = 'iso3c')
  
  # get parameters based on desired coverage
  coverage <- seq(from=0, to = 0.8, by=0.1)
  
  # get country fit
  # grab the json from the data exports
  country_fit <- get_fit_output(iso3c)
  
  # get inputs for forward sims
  model_inputs <- get_parameters_byCov(country_fit, country = country_example, cov_vector = coverage, continue_rt = continue_rt, prop_immune = prop_immune, increase_rt = increase_rt)
  
  ## get country specific params
  contact_matrix = squire::get_mixing_matrix(country_example)
  population = squire::get_population(country_example)
  
  # get inputs from json
  tt_R0 <- as.numeric(country_fit$tt_R0)
  dates <- country_fit$dates
  deaths <- as.numeric(country_fit$deaths)
  Rts <- as.numeric(country_fit$rt)
  max_vaccine <- as.numeric(country_fit$max_vaccine)
  betas <- as.numeric(country_fit$betas)
  
  # future Rts and vaccines based on user inputs
  Rt_vec <- c(Rts, model_inputs$rt)
  tt_s <- c(as.numeric(country_fit$tt_R0), model_inputs$tt_R0)
  
  # future betas based on user inputs
  durs <- nimue:::default_durations()
  probs <- nimue:::default_probs()
  future_beta <- nimue::beta_est_infectiousness(dur_IMild = durs$dur_IMild,
                                                  dur_ICase = durs$dur_ICase,
                                                  prob_hosp = probs$prob_hosp,
                                                  rel_infectiousness = rep(1, 17),
                                                  mixing_matrix = squire:::process_contact_matrix_scaled_age(squire:::get_mixing_matrix(iso3c = iso3c), population$n),
                                                  R0 = model_inputs$rt)
    
  new_betas <- c(betas, future_beta)

  vacc_days_vec <- seq(tail(as.numeric(country_fit$tt_R0),1)+1, tail(model_inputs$tt_R0,1))
  vacc_days <- length(vacc_days_vec)
  new_vaccines <- c(as.numeric(country_fit$max_vaccine), rep(model_inputs$max_vaccine[1], vacc_days))
  tt_vacc <- c(as.numeric(country_fit$tt_R0), vacc_days_vec)

  # efficacy
  vaccine_efficacy_infection <- c(country_fit$eff_infection, model_inputs$eff_infection)
  vaccine_efficacy_disease <- c(country_fit$eff_disease, model_inputs$eff_disease)
  # Vaccine strategy
  vacc_json <- paste0("https://github.com/mrc-ide/nimue_js/releases/download/v1.0.10/", iso3c, ".json")
  vacc_strat_json <- jsonlite::read_json(vacc_json)
  vaccine_uptake <- 0.8
  vaccine_available <- 1
  strategy <- "HCW, Elderly and High-Risk"
  cov_mat <- matrix(unlist(vacc_strat_json$etagePriority), ncol = 17)  * vaccine_uptake
  
  # allow children to be vaccinated
  cov_mat_sub <- tail(cov_mat,1)
  cov_mat_sub2 <- matrix(rbind(cov_mat_sub, cov_mat_sub, cov_mat_sub), nrow = 3)
  cov_mat_sub2[3,1] <- vaccine_uptake
  cov_mat_sub2[2,2] <- vaccine_uptake
  cov_mat_sub2[1,3] <- vaccine_uptake
  cov_mat <- matrix(rbind(cov_mat, cov_mat_sub2), nrow = nrow(cov_mat) + 3)
  cov_mat <- scale_cov_mat(cov_mat, vaccine_available, population$n)
  
  # NIMUE RUN
  det_out_vac <- nimue::run(
    country = country_example,
    beta_set = new_betas,
    dur_R = 365,
    use_dde = TRUE,
    seeding_cases = 5,
    seeding_age_order = 6:10,
    tt_R0 = tt_s,
    R0 = new_betas,
    max_vaccine = new_vaccines,
    tt_vaccine = tt_vacc,
    time_period = length(tt_R0)+forecast,
    dur_V = Inf,
    vaccine_efficacy_infection = as.list(as.numeric(vaccine_efficacy_infection)),
    tt_vaccine_efficacy_infection = seq_along(as.numeric(vaccine_efficacy_infection)),
    vaccine_efficacy_disease = as.list(as.numeric(vaccine_efficacy_disease)),
    tt_vaccine_efficacy_disease = seq_along(as.numeric(vaccine_efficacy_infection)),
    rel_infectiousness_vaccinated = 0.5,
    vaccine_coverage_mat = cov_mat)
  
  # get results
  index <- squire:::odin_index(det_out_vac$model)
  D_index <- index$D
  inf_cumu_index <- index$infections_cumu
  hosp_demand_index <- index$hospital_demand
  icu_demand_index <- index$ICU_demand
  vacc_cumu_index <- index$vaccines_cumu
  recovered_index_1 <- index$R1[,]
  recovered_index_2 <- index$R2[,]
  
  df_rec <- data.frame(rec = rowSums(det_out_vac$output[,recovered_index_1, 1]) + rowSums(det_out_vac$output[,recovered_index_2, 1]))
  pr_i <- df_rec/sum(population$n)*100
  
  # build data frame for main plot outputs
  df <- data.frame(date = as.Date(dates), real = deaths)
  df2 <- data.frame(deaths = diff(rowSums(det_out_vac$output[,D_index,1])),
                    infections = diff(rowSums(det_out_vac$output[,inf_cumu_index,1])),
                    hospitalisations = rowSums(det_out_vac$output[-1,hosp_demand_index,1]),
                    critical = rowSums(det_out_vac$output[-1,icu_demand_index,1]),
                    vaccines = (rowSums(det_out_vac$output[-1,vacc_cumu_index,1])))
  df2$date <- seq.Date(as.Date(dates[2]), as.Date(dates[1]) + length(df2$deaths), 1)
  df <- dplyr::left_join(df2, df, by = "date")
  df$country <- country_example
  df$iso3c <- iso3c
  df$continue_rt <- continue_rt
  
  # join with Rt
  df$timestep <- 1:nrow(df)
  Rdf <- data.frame(timestep = tt_s, Rt = Rt_vec)
  df <- left_join(df, Rdf)
  
  # join with vaccinated people
  Vdf <- data.frame(timestep = tt_vacc, new_vaccines = new_vaccines)
  df <- left_join(df, Vdf)

  return(df)
}

# scale vaccine coverage for availability function
scale_cov_mat <- function(cov_mat, vaccine_available, pop) {
  
  # total vaccs available
  tot_vaccines <- sum(pop*vaccine_available)
  
  # step 1, find when max allocation exceeds capacity
  step <- 1
  step_found <- FALSE
  tot_vaccs_steps <- 0
  cov_mat_dup_ex <- rbind(0, cov_mat)
  
  while(!step_found && step <= nrow(cov_mat)) {
    
    if(nrow(cov_mat) == 1) {
      step_found <- TRUE
    }
    
    vaccs_in_step <- sum((cov_mat_dup_ex[step+1, ] - cov_mat_dup_ex[step, ]) * pop)
    tot_vaccs_steps <- tot_vaccs_steps + vaccs_in_step
    if(tot_vaccs_steps > tot_vaccines) {
      step_found <- TRUE
    } else {
      step <- step+1
    }
  }
  
  # if we have enough vaccine return now
  if(step > nrow(cov_mat)) {
    return(cov_mat)
  }
  
  # set steps after max available reached to 0
  if(step < nrow(cov_mat)) {
    cov_mat[(step+1):nrow(cov_mat),] <- 0
  }
  
  # now set this step to be correct for available
  tots_given <- sum(cov_mat[step-1,] %*% pop)
  tots_tried <- sum(cov_mat[step,] %*% pop)
  remaining <- tot_vaccines - tots_given
  
  # next_group
  next_group <- cov_mat[step,]-cov_mat[step-1,]
  poss_to_vacc <- (next_group[which(next_group > 0)] * pop[which(next_group > 0)])
  new_cov <- (remaining/sum(poss_to_vacc)) * cov_mat[step, which(next_group > 0)]
  cov_mat[step, which(next_group > 0)] <- new_cov
  return(cov_mat)
}
