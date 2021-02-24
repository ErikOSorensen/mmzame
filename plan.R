source("packages.R")
source("functions.R")

NP_S = 10000
N_STANDARD_BRONARS = 10000
SELFISH_REQUIREMENT = 0.05

mmzame_plan <- drake_plan(
  mmzame_decisions = readRDS(file_in("data/mmzame_decisions.rds")),
  background = readRDS( file_in("data/background.rds")),
  prop3list = prepare_decisions(mmzame_decisions, c("moral", "risk")),
  prop4list = prepare_decisions(mmzame_decisions, c("moral", "dictator")),
  prop2list = prepare_decisions(mmzame_decisions, c("risk", "dictator")),
  all_domains = df_3way(mmzame_decisions),
  pall_domains = calculate_3waytest(all_domains, np=NP_S),
  p2_00 = purrr::map(prop2list, p_permutations_bronars, np=NP_S, p_Bronars=0.00, rFOSD=TRUE),
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
  sym_moral = prepare_decisions(mmzame_decisions, c("moral")),
  sym_dict  = prepare_decisions(mmzame_decisions, c("dictator")),
  sym_risk  = prepare_decisions(mmzame_decisions, c("risk")),
  symmetricp_dict = purrr::map(sym_dict, p_symmetric, np=NP_S),
  symmetricp_moral = purrr::map(sym_moral, p_symmetric, np=NP_S),
  symmetricp_risk = purrr::map(sym_risk, p_symmetric, np=NP_S),
  CES_dfs = make_CES_dfs(mmzame_decisions), 
  CES_estimates = make_CES_estimates(CES_dfs),
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
    output_file = file_out("testing_rationality.html" ),
    quiet = TRUE),
  Dtesting_theory = rmarkdown::render(
    knitr_in("vignettes/testing_theory.Rmd"),
    output_file = file_out("testing_theory.html"),
    quiet = TRUE), 
  Revision = rmarkdown::render(
    knitr_in("vignettes/revision.Rmd"),
    output_file = file_out("revision.html"),
    quiet = TRUE)
)