source(here::here("plan.R"))
future::plan(multicore)
getwd()
make(mmzame_plan, parallelism = "future", jobs = 30)
