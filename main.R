library(squire)
library(nimue)
library(purrr)
library(furrr)
library(tidyverse)
library(zoo)
library(countrycode)
library(patchwork)

col1 <- "#5da0b5"
col2 <- "#c59e96"
col3 <- "#747473"
col4 <- "#5c8e72"
col5 <- "#2a73bb" # Reff color

source("set_parameters_countries.R")
source("functions.R")

country_example <- "Singapore"
continue_rt <- TRUE
increase_rt <- TRUE

population = sum(squire::get_population(country_example)$n)
prop_R_dat <- read.csv("prop_R.csv")
prop_R = prop_R_dat[which(prop_R_dat$country == country_example),]$prop_R

out <- create_vacc_fit(country_example, continue_rt, forecast = 365, prop_immune = prop_R, increase_rt = increase_rt)

Rt_new <- na.locf(out$Rt, fromLast = TRUE)
Rt_new <- c(Rt_new, rep(tail(Rt_new,1), nrow(out) - length(Rt_new)))
out <- out %>%
  mutate(Rt = Rt_new)

p3 <- ggplot(data = out, aes(x = date, y = deaths)) +
  geom_bar(aes(x = as.Date(date), y = real, fill = "Reported"),
           stat = "identity",
           fill = "#c59e96") +
  geom_line() +
  geom_line(aes(date, zoo::rollmean(real, 7, na.pad = TRUE), color = "7-day Weekly Mean"), lwd = 1) +
  geom_line(aes(date, deaths, color = "Modelled deaths"), lwd = 1) +
  ylab("Daily deaths") +
  scale_fill_manual(values = "#c59e96") +
  scale_color_manual(values = c("black", "darkgrey")) +
  xlab("") +
  scale_y_continuous(expand = c(0,0)) +
  ggpubr::theme_pubclean() +
  theme(axis.line = element_line(), legend.title = element_blank(), legend.key = element_blank())
  
p4 <- ggplot(data = out) +
  geom_line(aes(x = date, y = vaccines/population*100, color = "Vaccine coverage", linetype = "Vaccine coverage"), size = 1) +
  geom_line(aes(x = date, y = Rt*12, color = "Rt", linetype = "Rt"), size = 0.8) +
  scale_linetype_manual(values = c("solid","dashed"), labels = c("Rt", "Vaccine coverage"), name = "Measure") +
  geom_hline(yintercept = 5.*12, linetype = "dotted", size = 1) +
  geom_hline(yintercept = 1.*12, linetype = "dotted", size = 1) +
  scale_color_manual(values = c(col5, col4), labels = c("Rt", "Vaccine coverage"), name = "Measure") +
  ggpubr::theme_pubclean() +
  theme(axis.line = element_line(),
        axis.line.y.right = element_line(color = col5), 
        axis.ticks.y.right = element_line(color = col5),
        axis.text.y.right = element_text(color = col5),
        legend.title = element_blank(),
        legend.key = element_blank()) +
  #ggtitle(country_example) +
  labs(x = "Date") +
  scale_y_continuous(    # Features of the first axis
    name = "Vaccine coverage (%)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./12, name="Rt"))


plots_out <- (p2 + p4) / (p1 + p3) + plot_layout(nrow = 2, byrow = TRUE) + plot_annotation(tag_levels = "A") +
  plot_annotation(title = country_example) +
  plot_layout(guides = "collect")& theme(legend.position = 'bottom')


ggsave(paste0("plots_out_", country_example, ".png"), plots_out, height = 8, width = 10)

