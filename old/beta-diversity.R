library(tidyverse)
library(lubridate)
library(readxl)
library(betapart)
setwd('/Users/raphael/Library/CloudStorage/Box-Box/BMM-US')

d <- readxl::read_xlsx('data/eBird/data_diversity.xlsx')

ds <- d %>% 
  filter(name=="Tompkins") %>% 
  filter(year(obs_dt_day)==2020) %>% 
  filter(month(obs_dt_day)<8) %>% 
  mutate(pres = ifelse(is.na(fun1_obs_count),1,0)) %>% 
  pivot_wider(names_from="species_code",values_from="pres", values_fill = 0) %>% 
  select(-c(obs_dt_day,GroupCount,name, fun1_obs_count))


core <- betapart.core(ds)


a=beta.pair(core)
