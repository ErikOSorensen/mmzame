source(here::here("packages.R"))
source(here::here("functions.R"))
NP_S = 5
SELFISH_REQUIREMENT = 0.05

plan <- drake_plan(
  mmzame_decisions = readRDS(file_in(here("data/mmzame_decisions.rds"))),
  background = readRDS( file_in(here("data/background.rds"))),
  
  prop3list = prepare_decisions(mmzame_decisions, c("moral", "risk")),
  prop4list = prepare_decisions(mmzame_decisions, c("moral", "dictator")),
  p3_00 = foreach(d = prop3list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                          np=NP_S, p_Bronars=0, rFOSD=TRUE),
  p3_05 = foreach(d = prop3list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                            np=NP_S, p_Bronars=0.05, rFOSD=TRUE),
  p3_10 = foreach(d = prop3list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                            np=NP_S, p_Bronars=0.10, rFOSD=TRUE),
  p3_15 = foreach(d = prop3list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                            np=NP_S, p_Bronars=0.15, rFOSD=TRUE),
  p3_20 = foreach(d = prop3list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                            np=NP_S, p_Bronars=0.20, rFOSD=TRUE),
  p4_00 = foreach(d = prop4list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                            np=NP_S, p_Bronars=0, rFOSD=TRUE),
  p4_05 = foreach(d = prop4list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                            np=NP_S, p_Bronars=0.05, rFOSD=TRUE),
  p4_10 = foreach(d = prop4list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                            np=NP_S, p_Bronars=0.10, rFOSD=TRUE),
  p4_15 = foreach(d = prop4list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                            np=NP_S, p_Bronars=0.15, rFOSD=TRUE),
  p4_20 = foreach(d = prop4list) %do% p_permutations_bronars(d, c("moral","risk"),
                                                                   np=NP_S, p_Bronars=0.20, rFOSD=TRUE),
  prop3 = list(p3_00, p3_05, p3_10, p3_15, p3_20),
  prop4 = list(p4_00, p4_05, p4_10, p4_15, p4_20),
  symmetric = mmzame_decisions %>% 
    filter(treatment=="dictator") %>% 
    mutate(yshare=y/(x+y)) %>% 
    group_by(id) %>%
    summarize(mean_yshare = mean(yshare)) %>%
    filter( abs((mean_yshare - 0.50)) < SELFISH_REQUIREMENT) ,
  selfish = mmzame_decisions %>%
    filter(treatment=="dictator") %>% 
    mutate(yshare=y/(x+y)) %>% 
    group_by(id) %>%
    summarize(mean_yshare = mean(yshare)) %>%
    filter( mean_yshare > (1.0 - SELFISH_REQUIREMENT) ),
  hypotheses_data = prepare_hypothesis_data(prop3, selfish, prop4, symmetric)
)