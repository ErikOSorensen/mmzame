source(here::here("plan.R"))
future::plan(future.callr::callr, workers = 10L)
make(plan)
