library(targets)
library(tarchetypes)
library(future)
library(visNetwork)
library(future.callr)
future::plan(callr)
source("functions.R")

# Set target-specific options such as packages.
tar_option_set(
  packages = c(
    "dplyr",
    "tidyr",
    "here",
    "revealedPrefs",
    "patchwork",
    "ggplot2",
    "showtext",
    "extrafont",
    "numDeriv",
    "forcats",
    "purrr"
  )
)

NP_S = 10000
N_STANDARD_BRONARS = 100000
SELFISH_REQUIREMENT = 0.05
EPSILON = 0.0001
# Now for the id and server to download data from.
DATA_FILE_ID = 5446287
BACKGROUND_FILE_ID = 5446286
DATA_SERVER = "dataverse.harvard.edu"

# End this file with a list of target objects.
list(
  tar_target(mmzame_decisions_alltreatments, dataverse::get_dataframe_by_id(DATA_FILE_ID, 
                                                                            .f = readr::read_tsv, 
                                                                            server = DATA_SERVER) |>
               filter(id<200|id>300)),
  tar_target(mmzame_background, dataverse::get_dataframe_by_id(BACKGROUND_FILE_ID,
                                                               .f = readr::read_tsv,
                                                               server = DATA_SERVER) |>
               filter(id<200|id>300)),
  tar_target(mmzame_decisions, mmzame_decisions_alltreatments %>%
               filter(treatment %in% c("moral", "risk", "dictator"))),
  tar_target(prop3list, prepare_decisions(mmzame_decisions, c("moral","risk"))),
  tar_target(prop4list, prepare_decisions(mmzame_decisions, c("moral","dictator"))),
  tar_target(prop2list, prepare_decisions(mmzame_decisions, c("risk","dictator"))),
  tar_target(all_domains, df_3way(mmzame_decisions)),
  tar_target(pall_domains, calculate_3waytest(all_domains, np=NP_S)),
  tar_target(p2_00, purrr::map(prop2list, p_permutations_bronars, np=NP_S, p_Bronars=0.00, rFOSD=TRUE)),
  tar_target(p3_00, purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.00, rFOSD=TRUE)),
  tar_target(p3_05, purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.05, rFOSD=TRUE)),
  tar_target(p3_10, purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.10, rFOSD=TRUE)),
  tar_target(p3_15, purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.15, rFOSD=TRUE)),
  tar_target(p3_20, purrr::map(prop3list, p_permutations_bronars, np=NP_S, p_Bronars=0.20, rFOSD=TRUE)),
  tar_target(p4_00, purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.00, rFOSD=TRUE)),
  tar_target(p4_05, purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.05, rFOSD=TRUE)),
  tar_target(p4_10, purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.10, rFOSD=TRUE)),
  tar_target(p4_15, purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.15, rFOSD=TRUE)),
  tar_target(p4_20, purrr::map(prop4list, p_permutations_bronars, np=NP_S, p_Bronars=0.20, rFOSD=TRUE)),
  tar_target(prop3, list(p3_00, p3_05, p3_10, p3_15, p3_20)),
  tar_target(prop4, list(p4_00, p4_05, p4_10, p4_15, p4_20)),
  tar_target(symmetric, mmzame_decisions %>% filter(treatment=="dictator") %>%
               mutate(yshare = y/(x+y)) %>%
               group_by(id) %>%
               summarize(mean_yshare = mean(yshare)) %>%
               filter( abs(mean_yshare - 0.5) < SELFISH_REQUIREMENT)),
  tar_target(selfish, mmzame_decisions %>% filter(treatment=="dictator") %>%
               mutate(yshare=y/(x+y)) %>%
               group_by(id) %>%
               summarize(mean_yshare = mean(yshare)) %>%
               filter( mean_yshare > (1.0 - SELFISH_REQUIREMENT))),
  tar_target(hypotheses_data, prepare_hypothesis_data(prop3, selfish, prop4, symmetric)),
  tar_target(bronars_budgets, bronars_datasets(mmzame_decisions, 50, N_STANDARD_BRONARS)),
  tar_target(pure_bronars, purrr::map_dbl(bronars_budgets, ccei_on_bronars_budgets_df)),
  tar_target(sym_moral, prepare_decisions(mmzame_decisions, c("moral"))),
  tar_target(sym_dict, prepare_decisions(mmzame_decisions, c("dictator"))),
  tar_target(sym_risk, prepare_decisions(mmzame_decisions, c("risk"))),
  tar_target(symmetricp_dict, purrr::map(sym_dict, p_symmetric, np=NP_S)),
  tar_target(symmetricp_moral, purrr::map(sym_moral, p_symmetric, np=NP_S)),
  tar_target(symmetricp_risk, purrr::map(sym_risk, p_symmetric, np=NP_S)),
  tar_target(CES_dfs, make_CES_dfs(mmzame_decisions)),
  tar_target(CES_estimates, make_CES_estimates(CES_dfs)),
  tar_render(Aaggregate_behavior, here::here("vignettes","aggregate_behavior.Rmd")),
  tar_render(Bindividual_behavior, here::here("vignettes","individual_behavior.Rmd")),
  tar_render(Ctesting_rationality, here::here("vignettes","testing_rationality.Rmd")),
  tar_render(Dtesting_theory, here::here("vignettes", "testing_theory.Rmd")),
  tar_render(Revision, here::here("vignettes","revision.Rmd")),
  tar_render(TableA1, here::here("vignettes","background_table.Rmd")),
  tar_render(long_table, here::here("vignettes","long_table_individuals.Rmd"))
)

