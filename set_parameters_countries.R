library(squire)
library(dplyr)

# Get country fits 
get_fit_output <- function (country_id){
  
  today <- Sys.Date()
  
  # Pull latest country fit
  iso3c <- country_id
  country <- squire::population$country[squire::population$iso3c==iso3c][1]
  json_path <- file.path("https://raw.githubusercontent.com/mrc-ide/global-lmic-reports/master/",iso3c,"input_params.json")
  json <- jsonlite::read_json(json_path)
  
  date_deaths <- unlist(lapply(json, function(x){
    if("deaths" %in% names(x)) {
      x["date"]
    } else {
      NULL
    }
  }))
  dates<- c(unlist(lapply(json, "[[", "date")))
  tt_R0 <- c(unlist(lapply(json, "[[", "tt_beta")) + 1)
  rt <- c(unlist(lapply(json, "[[", "Rt")))
  betas <- as.numeric(unlist(lapply(json, "[[", "beta_set")))
  max_vaccine <- c(unlist(lapply(json, "[[", "max_vaccine")))
  eff_infection <-  c(unlist(lapply(json, "[[", "vaccine_efficacy_infection")))
  eff_disease <-  c(unlist(lapply(json, "[[", "vaccine_efficacy_disease")))
  deaths <- c(unlist(lapply(json, "[[", "deaths")))
  deaths_df <- data.frame(dates = date_deaths, deaths = deaths)
  
  rt <- rt[which(dates <= max(date_deaths))]
  tt_R0 <- tt_R0[which(dates <= max(date_deaths))]
  betas <- betas[which(dates <= max(date_deaths))]
  max_vaccine <- max_vaccine[which(dates <= max(date_deaths))]
  eff_infection <- eff_infection[which(dates <= max(date_deaths))]
  eff_disease <- eff_disease[which(dates <= max(date_deaths))]
  dates <- dates[which(dates <= max(date_deaths))]
  
  out <- data.frame(cbind(dates,tt_R0,rt,betas,max_vaccine,eff_infection, eff_disease)) %>% filter(dates <= today ) %>% left_join(deaths_df, by = "dates")
  
}

# Get times for lifting intervention based on Rt vector

time_estimation <- function (Rt,efficacy,daily_vaccine,target_pop, current_cov){
  
  coverage_needed <- (1-1/Rt)*(1/efficacy)   # Coverage needed to achieve that Rt 
  
  coverage_reamining <- coverage_needed - current_cov 
  
  times <- round((coverage_reamining*target_pop)/daily_vaccine) # Number of days to achieve that coverage 
  
  return <- list(times,coverage_needed)
}

# Get parameters to run the model based on fitting results and vaccine efficacy 
get_parameters_byRt <- function(fit_output,country,last_rt,continue_rt = TRUE, steps=10){
  
  # Formatting data 
  fit_output$dates <- as.Date(fit_output$dates)
  fit_output[,2:6] <- as.numeric(unlist(fit_output[,2:6]))
  
  
  # Vaccination rate parameters
  pop <- squire::get_population(country = country)$n
  target_pop <- sum(pop)
  
  daily_vaccine <- round(mean(tail(fit_output$max_vaccine, n=30)))   # Daily vaccine assumed to continue as an average of last 30 days reported 
  efficacy_disease <- mean(tail(as.numeric(fit_output$eff_disease), n=30))      #  Vaccine efficacy assumed to remain the same as an average of the last 30 days estimated 
  efficacy_infection <- mean(tail(as.numeric(fit_output$eff_infection), n=30))
  
  # Day 0 of simulation
  last_t <- last(fit_output$tt_R0)
  
  
  # current coverage:
  current_cov<- sum(fit_output$max_vaccine)/sum(squire::get_population(country)$n)
  current_ideal_rt <- 1/(1- current_cov*efficacy_disease)
  
  # current Rt:
  current_Rt <-last(fit_output$rt)
  
  # Rt vector
  if(!continue_rt){
    r_vacc <- seq(from=current_ideal_rt, to=last_rt,length.out = steps)}  
  else{
    r_vacc <- seq(from=current_Rt, to=last_rt,length.out =steps) 
  }
  
  
  
  # Time vector that represents the days Rt is increased based on coverage 
  t_estimation <-  time_estimation(Rt= r_vacc, efficacy = efficacy_disease,
                                   daily_vaccine =daily_vaccine ,target_pop = target_pop,
                                   current_cov = current_cov)
  
  t_needed <- t_estimation[[1]]
  t_needed <- t_needed + last_t +1
  t_dates <- fit_output$dates[1] + t_needed
  
  # Binding new parameters
  out <- data.frame(cbind (t_needed,r_vacc))
  colnames(out) <- c("tt_R0","rt")
  out$max_vaccine <- daily_vaccine
  out$dates <- as.Date(fit_output$dates[1] + t_needed)
  out$eff_infection <- efficacy_infection
  out$eff_disease <- efficacy_disease
  out$coverage <- t_estimation [[2]]
  return (out)
  
}

# Get times for lifting intervention based on coverage vector: #### 

time_estimation_cov <- function (cov_vector, efficacy, daily_vaccine, target_pop, current_cov){
  
  coverage_reamining <- cov_vector - current_cov 
  
  times <- round((coverage_reamining*target_pop)/daily_vaccine) # Number of days to achieve that coverage 
  
  return <- times
  
}


# Get parameters to run the model based on fitting results and vaccine efficacy 
get_parameters_byCov <- function(fit_output,country,cov_vector,continue_rt = TRUE, steps=10, prop_immune = 0, increase_rt = FALSE){
  
  # Formatting data 
  fit_output$dates <- as.Date(fit_output$dates)
  fit_output[,2:6] <- as.numeric(unlist(fit_output[,2:6]))
  
  
  # Vaccination rate parameters
  pop <- squire::get_population(country = country)$n
  target_pop <- sum(pop)
  
  daily_vaccine <- round(mean(tail(fit_output$max_vaccine, n=30)))   # Daily vaccine assumed to continue as an average of last 30 days reported 
  efficacy_disease <- mean(tail(as.numeric(fit_output$eff_disease), n=30))      #  Vaccine efficacy assumed to remain the same as an average of the last 30 days estimated 
  efficacy_infection <- mean(tail(as.numeric(fit_output$eff_infection), n=30))
  
  # Day 0 of simulation
  last_t <- last(fit_output$tt_R0)
  
  
  # current coverage:
  current_cov<- sum(fit_output$max_vaccine)/sum(squire::get_population(country)$n)
  current_ideal_rt <- 1/(1- current_cov*efficacy_infection)
  
  if(continue_rt == FALSE){
  cov_togo <- which(current_cov < cov_vector)
  cov_steps <- c(current_cov, cov_vector[cov_togo])
  } else {
  # current Rt:
  current_Rt <-last(fit_output$rt)
  current_ideal_cov <- (1-1/current_Rt)*(1/efficacy_infection) - prop_immune
  cov_togo <- which(current_ideal_cov < cov_vector)
  cov_steps <- c(current_ideal_cov, cov_vector[cov_togo])
  }
  
  r_vacc<- 1/(1-cov_steps*efficacy_infection)
  
  
  # Time vector that represents the days Rt is increased based on coverage 
  t_needed <- time_estimation_cov(cov_vector= cov_steps, efficacy = efficacy_infection,
                                  daily_vaccine =daily_vaccine,
                                  target_pop = target_pop,
                                  current_cov = current_cov)
  t_needed <- t_needed + last_t +1
  t_dates <- fit_output$dates[1] + t_needed
  
  # Binding new parameters
  out <- data.frame(cbind (t_needed,r_vacc))
  colnames(out) <- c("tt_R0","rt")
  out$max_vaccine <- daily_vaccine
  out$dates <- as.Date(fit_output$dates[1] + t_needed)
  out$eff_infection <- efficacy_infection
  out$eff_disease <- efficacy_disease
  out$coverage <- cov_steps
  
  if (increase_rt == TRUE) {
    out_final <- data.frame(tt_R0 = tail(out$tt_R0, 1) + 7,
                            rt = 5,
                            max_vaccine = 0,
                            dates = tail(out$dates,1)+7,
                            eff_infection = tail(out$eff_infection,1),
                            eff_disease = tail(out$eff_disease,1),
                            coverage = tail(out$coverage,1))
    out <- rbind(out, out_final)
  }
  return (out)
  
}




