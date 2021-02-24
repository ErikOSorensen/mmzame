
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

p_permutations_bronars <- function(decisions_individual, np = 10000, p_Bronars=0, rFOSD=FALSE) {
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
  list(id=id, ccei1=ccei1, ccei2=ccei2, cceis_permuted=cceis, 
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
    filter(n()==length(treatments)) %>%
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

p_symmetric <- function(decisions_df, np = 99) {
  id <- decisions_df$id[1]
  n <- nrow(decisions_df)
  x <- decisions_df %>% 
    select(x,y)
  p <- decisions_df %>%
    select(px,py)
  ccei_actual <- ccei(x,p)
  cceis <- numeric(np) 
  xm <- decisions_df %>%
    mutate(xm = y, ym = x) %>% select(xm, ym) %>% rename(x = xm, y = ym)
  pm <- decisions_df %>% # 
    mutate(pxm = py, pym = px) %>% select(pxm, pym) %>% rename(px=pxm, py=pym) 
  for (i in 1:np) {
    s <- (rbinom(n,1,0.5)==1)
    x1 <- x[s,]
    p1 <- p[s,]
    x2 <- xm[!s,]
    p2 <- pm[!s,]
    xp <- rbind(x1,x2)
    pp <- rbind(p1,p2)
    cceis[i] <- ccei(xp,pp) 
  }
  f <- ecdf(cceis)
  p_upper <- 1 - f(ccei_actual - 0.000001)
  list(id=id, ccei=ccei_actual, cceis_permuted=cceis, 
       p = p_upper)
}

cesinv <- function(x) {
  ta <- x[['ta']]
  tr <- x[['tr']]
  alpha <- 1 / (1 + exp(ta))
  rho <- 1 - exp(tr)
  c("rho"=rho, "alpha"=alpha)
}
delta.method <- function(func, vcov, x) {
  params <- func(x)
  stopifnot(all(names(x) == colnames(vcov)))
  stopifnot(names(x)==c("ta","tr"))
  d <- numDeriv::jacobian(func, x)
  new.vcov <- d %*% vcov %*% t(d)
  colnames(new.vcov) <- names(params)
  rownames(new.vcov) <- names(params)
  
  list(parameters = params, vcov = new.vcov)
}


pred_bud_share <- function(price,talpha,trho) {
  alpha <- max(min( 1 / (1 + exp(talpha)), 1-1E-14), 1E-14)
  rho <- min(1 - exp(trho), 1-1E-14)
  r <- -rho/(1-rho)
  g <- (alpha/(1-alpha))^(1/(1-rho))
  
  g/(price^r + g)
}

alpha_inv <- function(alpha) {
  talpha <- log(1/alpha - 1)
  talpha
}

rho_inv <- function(rho) {
  trho <- log( 1 - rho)
  trho
}

# Possible starting values:
starting_values <- list( 
  c(alpha = 0.5, rho = -1),
  c(alpha = 0.5, rho = 0),
  c(alpha = 0.5, rho = 0.8),
  c(alpha = 0.7, rho = -1),
  c(alpha = 0.7, rho = 0),
  c(alpha = 0.7, rho = 0.8),
  c(alpha = 0.9, rho = -1),
  c(alpha = 0.9, rho = 0),
  c(alpha = 0.9, rho = 0.8),
  c(alpha = 0.5, rho = -10),
  c(alpha = 0.5, rho = -100)
)

estimate_ces <- function(df, sv = starting_values) {
  lowdev <- 1e10
  est <- NULL
  id <- min(df$id)
  for (v in sv) {
    estc <- try(nls( bshare_self ~ pred_bud_share(p, ta, tr), 
                     data=df, 
                     start=list( ta= alpha_inv(v[['alpha']]), 
                                 tr= rho_inv(v[['rho']])), 
                     trace=FALSE),
                silent = TRUE)
    if (class(estc)!="try-error") {
      if (deviance(estc)<lowdev) {
        est <- estc
        lowdev <- deviance(est)
      }
    }
  }
  if(!is.null(est)) { 
    orgest <- delta.method(cesinv, vcov(est), coef(est))
    tibble( id = id,
            alpha = orgest$parameters[['alpha']],
            rho = orgest$parameters[['rho']],
            alpha_se = sqrt(orgest$vcov['alpha','alpha']),
            rho_se = sqrt(orgest$vcov['rho','rho']),
            deviance = lowdev,
            stopMessage = est$convInfo$stopMessage)
  } else {
    tibble(id = id,
           alpha = NA, rho=NA, alpha_se=NA, rho_se=NA, deviance = NA, stopMessage="Failed!")
  }
}
make_CES_dfs <- function(mmzame_decisions) {
  mmzame_decisions %>% 
    filter(treatment=="dictator") %>%
    mutate( px = 1/maxx, 
            py = 1/maxy,
            m = x*px + y*py,
            bshare_self = py*y/m,
            p = px/py) %>%
    group_by(id) %>%
    mutate(mean_bshare_self = mean(bshare_self)) %>%
    filter(mean_bshare_self <0.95) %>%
    ungroup() %>%
    select(id, bshare_self, p)
}
make_CES_estimates <- function(dfs) {  
  dfs %>% 
    group_by(id) %>%
    group_split() %>% 
    map(estimate_ces) %>%
    bind_rows()
}

p_permutationsK <- function(id, p_df, x_list, np = 99) {
  # Inputs should be 
  # id: id of individual
  # p_df: NxG dataframe of prices, (N choices, G goods). 
  # x_list: K-list of dataframes, each of size NxG.
  # The order of budget sets in each of the x_list datasets must match p_df !!!
  K <- length(x_list)
  nobs <- nrow(p_df)
  ngoods <- ncol(p_df)
  for (i in 1:K) {
    stopifnot( nobs == nrow(x_list[[i]]) & (ngoods == ncol(x_list[[i]])))
  }
  actual_cceis <- numeric(K)
  permuted_cceis <- numeric(np)
  for (i in 1:K) {
    actual_cceis[i] <- ccei(x_list[[i]], p_df)
  }
  xperm <- tibble(x_list[[1]]) # Allocates new memory to the copy
  for (i in 1:np) {
    xchoice <- sample(seq_along(x_list), nobs, TRUE)
    for (j in 1:nobs) {
      xperm[j,] <- x_list[[xchoice[j]]][j,]
    }
    permuted_cceis[i] <- ccei(xperm, p_df)
  }
  list(actual_cceis, permuted_cceis)
  tmin <- min(actual_cceis)
  tmax <- max(actual_cceis)
  f <- ecdf(permuted_cceis)
  p_min <- (1 - f(tmin-0.00001))^K
  p_max <- 1 - f(tmax-0.00001)^K
  c(p_min, p_max)
  list(id = id, ccei=actual_cceis, cceis_permuted = permuted_cceis,
       p_min =p_min, p_max=p_max, p_com = min( min(p_min,p_max)*2, 1))
}

ind_3way <- function(df) {
  # Set up data for estimation of a 3-way test
  id <- min(df$id)
  plist <- df %>% 
    group_by(t) %>%
    group_split() %>%
    map( function(df2) {
      df2 %>% arrange(bset) %>% select(px, py)
    })
  p_df <- plist[[1]]
  xlist <- df %>%
    group_by(t) %>%
    group_split() %>%
    map( function(df2) {
      df2 %>% arrange(bset) %>% select(x, y)
    })
  list(id=id, p_df=p_df, x_list=xlist)
}

df_3way <- function(df) {
  prepare_decisions(df, c("dictator","moral","risk")) %>%
    purrr::map(ind_3way)
}

calculate_3waytest <- function(data_3way_list, np=99) {
  out <- vector("list", length(data_3way_list))
  for(i in seq_along(data_3way_list)) {
    out[[i]] <- p_permutationsK(data_3way_list[[i]]$id,
                                data_3way_list[[i]]$p_df,
                                data_3way_list[[i]]$x_list, np=np)
  }
  out
}
