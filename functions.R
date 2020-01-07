
# Julian Boelaert has created the R package [revealedPrefs](http://cran.r-project.org/web/packages/revealedPrefs/index.html) 
# which, among other things,
# calculate a test for whether a data set is consistent with GARP. All the
# heavy lifting is done by C++ code, the only R part is the interfacing
# commands, so the code should be reasonably fast. The library is also minimally
# ambitious: The GARP-test only checks whether there exists a transitive closure
# of the preferences, and stops as soon as it has found a cycle - it does not
# identify where all these cycles are, and so are not slowed down so much by 
# very inconsistent individuals.
# 
# My idea is that it is fairly easy to add a binary search algorithm to 
# calculate the CCEI on top of that. A default tolerance of 0.0001 shold be 
# accurate enough for our purposes. 
ccei <- function(x, p, tol = 1e-4, print.progress = FALSE) {
  if ( !checkGarp(x,p)$violation ) 
    return(1.0)
  l <- 0.0
  u <- 1.0
  f <- function(l, u) { # This function does recursive binary search for CCEI on interval (l,u)
    h <- (u+l)/2.0
    if ( print.progress )
      print( c(l,h,u) )
    if ( checkGarp(x, p, afriat.par = h)$violation ) {
      u <- h
    } else {
      l <- h
    }
    if ( abs(u-l) < 2 * tol) {
      return( (u+l)/2 )
    } else {
      return( f(l, u) )
    }
  }
  return( f(l, u) )
}



# I want the simulation code code to allow for adding a certain fraction
# of Bronars choices, i.e. uniform choices over the budget sets.
# First, I define a function that takes the internal data from the 
# permutation file and draws a Bronars dataset. This is based
# on an actual dataset (one data frame for choices, one for prices),
# and it returns a dataset of simulated Bronars data of choices,
# with choices being made uniformly on the budget line.
# 
# 
# I include a parameter for whether the simulated choices should 
# be on the non-dominated (FOSD) par of the budget line.
simulated_bronars <- function(x, p, respect_FOSD=FALSE) {
  # The x and the p arguments are datasets of equal size. 
  # I assume that all x are on the budget line.
  ns<-nrow(x)
  M <- p$px*x$x + p$py*x$y
  maxx <- M/p$px
  if (respect_FOSD) {
    xe <- M/(p$px + p$py)
    sxl <- runif(ns, 0, xe)
    sxr <- runif(ns, xe, maxx)
    left_random <- p$px > p$py
    right_random <- !left_random
    sx <- sxl * left_random + sxr * right_random
  }
  else {
    sx <- runif(ns, 0, maxx)  # Uniform selection of x on (0,maxx)
  }
  sy <- (M - p$px*sx)/p$py  # Residual determination of y, spending the budget.
  data.frame(x=sx, y=sy)    # The prices do not change.
}

p_permutations_bronars <- function(decisions_individual, ns, np = 100000, p_Bronars=0, rFOSD=FALSE) {
  # Arranged with n budget sets (sorted same way) for each of the two treatments being tested.
  id <- decisions_individual$id[1]
  n <- nrow(decisions_individual)/2
  x <- decisions_individual %>% 
    select(x,y)
  p <- decisions_individual %>%
    select(px,py)
  ccei1 <- ccei(head(x,n), head(p,n))
  ccei2 <- ccei(tail(x,n), tail(p,n))
  cceis <- numeric(np)
  for (i in 1:np) { 
    s <- (rbinom(n,1,0.5)==1)
    s <- c(s,!s) # This is the simulated mix of actual data.
    xactual <- x[s,]
    pactual <- p[s,] # Really only to cut down length.
    xbronars <- simulated_bronars(xactual,pactual, respect_FOSD = rFOSD)
    # Now mixing the actual and the bronars.
    bronars_selection <- rbinom(n, 1, p_Bronars)
    xtest <- (1-bronars_selection)*xactual + bronars_selection*xbronars
    cceis[i] <- ccei(xtest, pactual)
  }
  ps <- ccei_p_value(ccei1, ccei2, cceis)
  list(id=id, treatments=ns, ccei1=ccei1, ccei2=ccei2, cceis_permuted=cceis, 
       p_min=ps[1], p_max=ps[2], p_com=min(min(ps)*2, 1),
       p_Bronars=p_Bronars)
}

ccei_p_value <- function(ccei1, ccei2, cceis_permuted) {
  tmin <- min(ccei1, ccei2)
  tmax <- max(ccei1, ccei2)
  f <- ecdf(cceis_permuted)
  p_min <- (1-f(tmin-0.00001))^2
  p_max <- 1 - f(tmax-0.00001)^2
  c(p_min, p_max)
}

prepare_decisions <- function(decisions, treatments) {
  df <- decisions %>% 
    filter(treatment %in% treatments) %>%
    group_by(id, bset) %>% 
    filter(n()==2) %>%
    mutate( t = as.factor(treatment), px = maxy/maxx, py = 1.0) %>%
    select(id, t, bset, px, py, x, y) %>%
    ungroup() %>%
    arrange(id, t, bset)
  split(df, df$id)
}

prepare_hypothesis_data <- function(prop3, selfish_ids, prop4, symmetric_ids) {
  hyp <- list()
  j <- 1
  for (s in prop4) {
    for (o in s) {
      if (o$id %in% symmetric_ids$id) {
        hyp[[j]] <- tibble(id=o$id, comparison="Social Risk = Social", 
                           p_com=o$p_com, p_Bronars=o$p_Bronars)
        j <- j+1
      }
    }
  }
  for (s in prop3) {
    for (o in s) {
      if (o$id %in% selfish_ids$id) {
        hyp[[j]] <- tibble(id=o$id, comparison="Social Risk = Risk", 
                           p_com=o$p_com, p_Bronars=o$p_Bronars)
        j <- j+1
      }
    }
  }
  do.call(rbind,hyp)
}

bronars_datasets <- function(decisions_df, n_per_datasets, n_datasets, smallest_x=0) {
  n_total <- n_per_datasets * n_datasets
  source <- decisions_df %>% select(maxx, maxy)
  all_budgets <- source %>% sample_n(n_total, replace=TRUE)
  
  all_budgets$x <- runif(n_total, min=smallest_x, max=1-smallest_x) * all_budgets$maxx
  all_budgets$y <- all_budgets$maxy - all_budgets$x * (all_budgets$maxy/all_budgets$maxx)
  all_budgets$px <- 1
  all_budgets$py <- all_budgets$maxx / all_budgets$maxy
  all_budgets <- all_budgets %>% select(px,py,x,y)
  ab_list <- split(all_budgets, rep(1:n_datasets, each=n_per_datasets))
  ab_list
}

ccei_on_bronars_budgets_df <- function(df) {
  df_x <- df %>% select(x,y)
  df_p <- df %>% select(px,py)
  ccei(df_x, df_p)
}