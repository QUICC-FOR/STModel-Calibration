support_limits = function (model, par, var, source_data, pdf, 
    delta = 100, slimit = 2, progress=TRUE) 
{

        par_lo <- par
        for (i in 1:length(par)) {
            par_lo[[i]] = rep(-.Machine$double.xmax, times = length(par[[i]]))
        }

        par_hi <- par
        for (i in 1:length(par)) {
            par_hi[[i]] = rep(.Machine$double.xmax, times = length(par[[i]]))
        }
    

    lower_limit <- par
    upper_limit <- par
    
    best_lh <- likeli(model, par, var, source_data, pdf)
    if (is.infinite(best_lh) || is.nan(best_lh) || is.na(best_lh)) {
        for (i in 1:length(lower_limit)) lower_limit[[i]] <- NA
        for (i in 1:length(upper_limit)) upper_limit[[i]] <- NA
        return(list(upper_limits = upper_limit, lower_limits = lower_limit))
    }
    
    bad_likeli <- best_lh - (slimit + 10)
    for (i in 1:length(par)) {
    if(progress) cat("varying paramter ", names(par)[i], "\n")
        for (j in 1:length(par[[i]])) {
            par_copy <- par[[i]][[j]]
            smallest_step <- par[[i]][[j]]/delta
            par[[i]][[j]] <- par_hi[[i]][[j]]
            lhood <- likeli(model, par, var, source_data, pdf)
            if (is.infinite(lhood) || is.nan(lhood) || is.na(lhood)) 
                lhood <- bad_likeli
            lhdiff <- best_lh - lhood
            if (!is.nan(lhdiff) && !is.na(lhdiff) && lhdiff < 
                slimit) {
                upper_limit[[i]][[j]] <- par_hi[[i]][[j]]
            }
            else {
                par[[i]][[j]] <- par_copy
                lhdiff <- 0
                biggest_step <- (par_hi[[i]][[j]] - par[[i]][[j]])/50
                if (is.infinite(biggest_step)) 
                  biggest_step <- 1e+306
                step <- biggest_step
                while (step > smallest_step) {
                  lhdiff <- 0
                  while (!is.nan(lhdiff) && !is.na(lhdiff) && 
                    (lhdiff < slimit) && (par[[i]][[j]] < par_hi[[i]][[j]])) {
                    par[[i]][[j]] <- par[[i]][[j]] + step
                    lhood <- likeli(model, par, var, source_data, 
                      pdf)
                    if (is.infinite(lhood) || is.nan(lhood) || 
                      is.na(lhood)) 
                      lhood <- bad_likeli
                    lhdiff <- best_lh - lhood
                  }
                  par[[i]][[j]] <- par[[i]][[j]] - step
                  step <- step/10
                  if (!((par[[i]][[j]] + step) > par[[i]][[j]])) {
                    break
                  }
                }
                step <- smallest_step
                if ((par[[i]][[j]] + step) > par[[i]][[j]]) {
                  lhdiff <- 0
                  while (!is.nan(lhdiff) && !is.na(lhdiff) && 
                    (lhdiff < slimit) && (par[[i]][[j]] < par_hi[[i]][[j]])) {
                    par[[i]][[j]] <- par[[i]][[j]] + step
                    lhood <- likeli(model, par, var, source_data, 
                      pdf)
                    if (is.infinite(lhood) || is.nan(lhood) || 
                      is.na(lhood)) 
                      lhood <- bad_likeli
                    lhdiff <- best_lh - lhood
                  }
                }
                upper_limit[[i]][[j]] <- par[[i]][[j]] - step
                par[[i]][[j]] <- par_copy
            }
            par[[i]][[j]] <- par_lo[[i]][[j]]
            lhood <- likeli(model, par, var, source_data, pdf)
            if (is.infinite(lhood) || is.nan(lhood) || is.na(lhood)) 
                lhood <- bad_likeli
            lhdiff <- best_lh - lhood
            if (!is.nan(lhdiff) && !is.na(lhdiff) && lhdiff < 
                slimit) {
                lower_limit[[i]][[j]] <- par_lo[[i]][[j]]
            }
            else {
                par[[i]][[j]] <- par_copy
                lhdiff <- 0
                biggest_step <- (par[[i]][[j]] - par_lo[[i]][[j]])/50
                step <- biggest_step
                if (is.infinite(biggest_step)) 
                  biggest_step <- 1e+306
                while (step > smallest_step) {
                  lhdiff <- 0
                  while (!is.nan(lhdiff) && !is.na(lhdiff) && 
                    (lhdiff < slimit) && (par_lo[[i]][[j]] < 
                    par[[i]][[j]])) {
                    par[[i]][[j]] <- par[[i]][[j]] - step
                    lhood <- likeli(model, par, var, source_data, 
                      pdf)
                    if (is.infinite(lhood) || is.nan(lhood) || 
                      is.na(lhood)) 
                      lhood <- bad_likeli
                    lhdiff <- best_lh - lhood
                  }
                  par[[i]][[j]] <- par[[i]][[j]] + step
                  step <- step/10
                  if (!((par[[i]][[j]] - step) < par[[i]][[j]])) {
                    break
                  }
                }
                step <- smallest_step
                if ((par[[i]][[j]] - step) < par[[i]][[j]]) {
                  lhdiff <- 0
                  while (!is.nan(lhdiff) && !is.na(lhdiff) && 
                    (lhdiff < slimit) && (par_lo[[i]][[j]] < 
                    par[[i]][[j]])) {
                    par[[i]][[j]] <- par[[i]][[j]] - step
                    lhood <- likeli(model, par, var, source_data, 
                      pdf)
                    if (is.infinite(lhood) || is.nan(lhood) || 
                      is.na(lhood)) 
                      lhood <- bad_likeli
                    lhdiff <- best_lh - lhood
                  }
                }
                lower_limit[[i]][[j]] <- par[[i]][[j]] + step
                par[[i]][[j]] <- par_copy
            }
        }
    }
    list(upper_limits = upper_limit, lower_limits = lower_limit)
}