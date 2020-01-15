source("packages.R")
source("functions.R")

NP_S = 100
N_STANDARD_BRONARS = 100
SELFISH_REQUIREMENT = 0.05

mmzame_plan <- drake_plan(
  mmzame_decisions = readRDS(file_in(here("data/mmzame_decisions.rds"))),
  background = readRDS( file_in(here("data/background.rds"))),
  
  prop3list = prepare_decisions(mmzame_decisions, c("moral", "risk")),
  prop4list = prepare_decisions(mmzame_decisions, c("moral", "dictator")),
  p3_00 = purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.00, rFOSD=TRUE),
  p3_05 = purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.05, rFOSD=TRUE),
  p3_10 = purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.10, rFOSD=TRUE),
  p3_15 = purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.15, rFOSD=TRUE),
  p3_20 = purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.20, rFOSD=TRUE),
  p4_00 = purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.00, rFOSD=TRUE),
  p4_05 = purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.05, rFOSD=TRUE),
  p4_10 = purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.10, rFOSD=TRUE),
  p4_15 = purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.15, rFOSD=TRUE),
  p4_20 = purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.20, rFOSD=TRUE),
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
  hypotheses_data = prepare_hypothesis_data(prop3, selfish, prop4, symmetric),
  bronars_budgets = bronars_datasets(mmzame_decisions, 50, N_STANDARD_BRONARS),
  pure_bronars = purrr::map_dbl(bronars_budgets, ccei_on_bronars_budgets_df),
  Aaggregate_behavior = rmarkdown::render(
    knitr_in('vignettes/aggregate_behavior.Rmd'),
    output_file = file_out('aggregate_behavior.html'),
    quiet = TRUE),
  Bindividual_behavior = rmarkdown::render(
    knitr_in("vignettes/individual_behavior.Rmd"),
    output_file = file_out("individual_behavior.html"),
    quiet = TRUE),
  Ctesting_rationality = rmarkdown::render(
    knitr_in("vignettes/testing_rationality.Rmd"),
    output_file = file_out("testing_rationality.html"),
    quiet = TRUE),
  Dtesting_theory = rmarkdown::render(
    knitr_in("vignettes/testing_theory.Rmd"),
    output_file = file_out("testing_theory.html"),
    quiet = TRUE)
)