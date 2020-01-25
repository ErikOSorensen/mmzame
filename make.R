source(here::here("plan.R"))
future::plan(multiprocess)
getwd()
make(mmzame_plan, parallelism = "future", jobs = 10)
