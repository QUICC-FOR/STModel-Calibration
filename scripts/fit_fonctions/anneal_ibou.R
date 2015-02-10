## modified version by Isabelle BOULANGEAT
# 26-07-2014
#####################################################
# Performs simulated annealing to find maximum likelihood
# estimates for a set of parameters.
#
# Arguments:
# model = model function to parameterize. Arguments to this function will be
# provided from par, var, and source_data.
#
# par = list of parameters for which we are using simulated annealing to
# find maximum likelihood estimates.  Each par element name matches an
# argument in a function (either model or pdf). All values listed in par
# must be numeric vectors.  Vectors of length greater than one have each of
# their elements treated separately as individual parameters to estimate.
#
# var = list of other variables and data needed by the model and pdf
# functions, any type as needed.  These will be kept constant.
#
# source_data = data frame with dependent variable and associated
# independent variables.
#
# par_initStep = initial step for all the parameters. an element in par in both name and size
#
# par_testFct = function a given set of parameters that test if they are realistic (inside some boundaries themselves of in combination)
#
# pdf = probability density function. Make sure to use a function such as
# dnorm that can calculate log probability.
#
# dep_var = dependent variable label in source_data.
#
# initial_temp = temperature at which to start annealing.
#
# temp_red = interval by which to reduce temperature.
#
# max_iter = maximum iterations, where one iteration is one varying of each
# parameter to parameterize.  The annealing will run for this many
# iterations, unless alternate conditions have been specified (see min_change
# and min_drops).
#
# dep_var = label of model dependent variable
#
# ns = interval between changes in range
#
# nt = interval between drops in temperature
#
# min_change = an alternate way to specify quitting conditions. This is the
# minimum amount of change in likelihood in min_drop number of temperature
# drops.  If the change is less, execution stops.
#
# min_drops = the companion to min_change for alternate quitting conditions.
# This is the number of temperature drops over which the likelihood must
# have changed more than min_change for execution to continue.
#
# hessian = if TRUE, include the standard errors and the Hessian matrix
# in the output.  If FALSE, do not.
#
# delta = when calculating support limits, the number of pieces into which
# to divide each parameter.  This is the size of the "step" the function
# takes in trying to find the support limits.
#
# slimit = when calculating support limits, the number of likelihood
# units less than the optimum likelihood for which the support intervals
# will be calculated.  2 units is semi-standard.  1.92 units corresponds
# roughly to a 95% confidence interval.
#
# note = any note to self the user wants to make.  This will be
# included in the output.
#
# Author:  Lora Murphy, Institute of Ecosystem Studies
# murphyl@ecostudies.org
######################################################
# debug
#model = model
#par = pars
#var = var
#source_data = data
#dep_var = "st1"
#pdf = PDF
#par_initStep = par_initStep
#par_testFct = testBounds
#max_iter = 100
#initial_temp = 3
#note = ""
#progress = TRUE
#display=FALSE
#support = FALSE
#min_change = 1
#min_drops = 10
#temp_red = 0.95
#ns = 20
#nt = 100
#delta = 100 
#slimit = 2

######################################################

anneal<-function(model, 
par, 
var, 
source_data, 
par_initStep,
#par_testFct,  
pdf, 
dep_var, 
initial_temp = 3, 
temp_red = 0.95,
ns = 20, 
nt = 100, 
max_iter = 50000, 
min_change = 0, 
min_drops = 100,
#hessian = TRUE, 
delta = 100, 
slimit = 2, 
note = "", 
support = TRUE,
progress = TRUE,
display=TRUE) 
{

  ##
  ## Error checking
  ##

  # Check types
  if (!is.function(model)) {
    stop("anneal: model is not a function.\n")
  }
  if (!is.list(par)) {
    stop("anneal: par is not a list.\n")
  }
  if (!is.list(var)) {
    stop("anneal: var is not a list.\n")
  }
  if (!is.list(par_initStep)) {
    stop("anneal: par_initStep is not a list.\n")
  }
#  if (!is.function(par_testFct) ) {
#    stop("anneal: par_testFct is not a function.\n")
#  }
  if (!is.data.frame(source_data)) {
    stop("anneal: source_data is not a data frame.\n")
  }
  if (length(source_data) < 1) {
    stop("anneal: source_data contains no data.\n")
  }
  if (length(source_data[[1]]) < 1) {
    stop("anneal: source_data contains no data.\n")
  }
  
  if (!is.function(pdf)) {
    stop("anneal: pdf is not a function.\n")
  }

  # Make sure there is something in par
  parnames <- names(par)
  numpars <- length(par)
  if (numpars < 1) {
    stop("anneal: no parameters in par - nothing to anneal!\n")
  }

  # Make sure all values in par are numeric
  for (i in 1:numpars) {
    if (!is.numeric(par[[i]])) {
      stop("anneal: All values in par must be numeric.\n")
    }
    if (!is.vector(par[[i]])) {
      stop("anneal: All values in par must be vectors.\n")
    }
  }

  # Make sure the dependent variable is part of source_data
  if (!is.character(dep_var)) {
    stop("anneal: I cannot understand the dependent variable name.  I expect a string.\n")
  }
  if(!any(names(source_data)==dep_var)) {
    stop("anneal: the dependent variable must be part of source_data.\n")
  }

# Make sure that par and par_initStep have corresponding labels in par,
# and that the vectors are the same length

  for (i in 1:length(names(par_initStep))) {
    if(!any(parnames == names(par_initStep)[i])) {
      stop("anneal: all values in par_initStep must be present in par.  Can't find value ", names(par_initStep)[i], ".\n")
    }

   if (length(par_initStep[[i]]) != length( par[(parnames==names(par_initStep)[i])][[1]] )) {
     stop("anneal: par_initStep value ", names(par_initStep)[i], " must be of the same length as its corresponding element in par.\n")
     }
  }

  ##
  ## Prep the likelihood calculations
  ##

  # Here's where we'll put results of all pre-evaluations
  eval_results <- NULL

  # This is where to put the predicted
  predicted <- NULL

  # Combine par and var into one - we will put the varying values here
  par_template = par
  par <- c(par, var)

  # Set up the PDF call - set it up to use a copy of par
  datasets<-NULL
  datasets[[1]]<-list(value=source_data, varname="source_data")
  pdf_call<-analyze_function(pdf, list(value=par, varname="par_copy"), NULL, datasets)

  # Set up the model call - set it up to use a copy of par
  model_call<-analyze_function(model, list(value=par, varname="par_copy"), NULL, datasets)

  # Initialize support limits to zeroes before we "flatten" arrays
  upper_limit<-list()
  lower_limit<-list()
  for (i in 1:length(par_template)) {
    upper_limit[[names(par_template)[[i]]]]<-rep(0, length(par_template[[i]]))
    lower_limit[[names(par_template)[[i]]]]<-rep(0, length(par_template[[i]]))
  }

  # If there are any vectors of values in a varying parameter, "flatten"
  # into 1-d vectors; create an array of index vectors for
  # pinpointing just one parameter value in par, even if it's in
  # an array
  #
  # Plus make a flat array of the initial values for calculating step

# Make a placeholder for step
  par_step_copy <- par_template
#

#  initial_vals <- NULL
  nm <- NULL # new names
  par_index<-NULL
  for (i in 1:length(par_template)) {
    numpars <- numpars + (length(par_template[[i]]) - 1)
    # Create new names for display
    if (length(par_template[[i]]) == 1) nm<-c(nm, parnames[[i]])
      else nm <- c(nm, paste(parnames[[i]], c(1:length(par_template[[i]])), sep=""))


#    initial_vals<-c(initial_vals, par_template[[parnames[[i]]]])

    for (j in 1:length(names(par_template)))
      if (parnames[[i]] == names(par_template)[[j]]) ind <- j
    for (j in 1:length(par_template[[i]])) {
      par_index[[length(par_index)+1]] <- c(ind, j)
    }
  }

#  names(initial_vals) <- nm
  nm <- NULL

  ##
  ## Since I can't break the habits of C++ - a rundown of my variables...
  ##
  iterate <- TRUE             # flag for whether or not to keep cycling
  using <- 1                  # index for the parameter we're working with
  cycles <- 0                 # how many cycles of annealing we've done
  best_lh <- NULL             # best likelihood
  acc_lh <- NULL              # last accepted likelihood
  best_par <- par             # set of parameters providing the best likelihood
  best_cycle <- 1             # number when the current best likelihood was calculated
  best_predicted<-NULL        # Predicted values that go with the best parameters
  nacc <- 0                   # number of accepted likelihood values
  temp <- initial_temp        # annealing temperature
  lhist <- NULL               # array of lists with likelihood history
  #par_copy <- NULL           # copy of parameters that we will vary
  #parnames <- names(par_step) # names of parameters to vary - already defined above
  #numpars                    # number of parameters to vary - already defined above
  #par_index                  # index array for retrieving values from par - already defined above
  #upper_limit                # upper support limit for all varying parameters - already defined above
  #lower_limit                # lower support limit for all varying parameters - already defined above
  #par_step                   # search interval - defined below
  nacp <- rep(0, numpars)     # number of times we accepted changes to each parameter
  range_cycles <- numpars * ns # number of loops before we adjust range
  temp_cycles <- numpars * ns * nt # number of loops before we reduce temperature
  display_cycles<-numpars # regular update point
  numobs <- length(source_data[[dep_var]]) # number of observations
  aiccorr <- NULL             # value for AIC corr (not entirely sure what that is)

  # Calculate the initial search interval, which for each parameter is the
  # greater distance between the initial value and the upper and lower bound
    par_step = par_initStep

  ##
  ## Do the first likelihood calculation with initial values
  ##
  
  par_copy <- par
  # :  the model predicted values
  predicted <- do.call(model,model_call$call)
  if (any(is.infinite(predicted)) || any(is.nan(predicted)) || any(is.na(predicted))) {
    stop("anneal: initial conditions caused the model to produce math errors.\n")
  }
#  

  # : sum over the PDF
  best_lh <- sum (do.call(pdf,pdf_call$call)) # best likelihood
  best_predicted <- predicted
  if (is.infinite(best_lh) || is.nan(best_lh) || is.na(best_lh)) best_lh <- -1000000
  acc_lh <- best_lh
  lhistcycles <- cycles
  lhisthood <- best_lh
  lhisttemp <- temp
  lhistpar <- vector(mode="list")
  temppar <- list()
  for (q in 1:length(parnames)) temppar[[parnames[q]]] <- par[[parnames[q]]]
  lhistpar[[1]] <- temppar


  # Calculate AIC corr (corrected for small sample size)
  aiccorr <- (-2.0 * best_lh)+((2*numpars)*(numobs/(numobs - numpars - 1)))

  # Setup output display and show it for the first time
  # Create a window so likdisplay will default to showing the graphics
  if(display)
  {
  plot(x=1,y=1)
  layout(matrix(c(1,2), 2, 1, byrow = TRUE)) # this makes a two-part graph
  likdisplay(lhistcycles, lhisthood, slp, R2, aiccorr, temp, max_iter)
  }

   if(progress)	cat("Initial likelihood = ", best_lh, '\n')


	
  tryCatch ({
    ##
    ## DO ANNEALING
    ##

    # Keep cycling until we're done
    while (iterate) {

      cycles <- cycles + 1;

	if(progress) cat("Cycle = ", cycles / numpars, " Likelihood = ", best_lh, '\n')
 
      # Make a copy of parameters before we vary
      par_copy <- par

      # Choose a new value for the parameter being varied this time
     # Make sure the new value is in bounds
      
#      test=FALSE
#      
#      while(!test)
#      {
      par_copy[[par_index[[using]]]] <- par[[par_index[[using]]]] + ((runif(1)*2 - 1)* par_step[[using]]);

#      test = par_testFct(par_copy, source_data)
#      }
#      

      # :  the model predicted values
      predicted <- do.call(model,model_call$call)

      # Make sure that all values are defined
      if (!any(is.infinite(predicted)) && !any(is.nan(predicted)) && !any(is.na(predicted))) {

        # : sum over the PDF
        lh <- sum (do.call(pdf,pdf_call$call))
        print(lh)
        if (is.infinite(lh) || is.nan(lh) || is.na(lh)) lh <- -1000000

        # Accept the new values if likelihood increases or at least stays the same
        # since the last accepted likelihood
        if (lh >= acc_lh) {
          par <- par_copy
          acc_lh <- lh
          nacc <- nacc + 1
          nacp[[using]] <- nacp[[using]] + 1
          # If this is a new maximum, then update the maximum likelihood
          if (acc_lh > best_lh) {
            best_par <- par_copy
            best_lh <- acc_lh
            best_cycle <- cycles
            best_predicted <- predicted


            # Calculate AIC corrected for small sample size
            aiccorr <- (-2.0 * best_lh)+((2*numpars)*(numobs/(numobs - numpars - 1)))

            lhisttemp[length(lhisttemp)+1] <- temp
            lhistcycles[length(lhistcycles) + 1] <- cycles / numpars
            lhisthood[length(lhisthood) + 1] <- best_lh
            temppar <- list()
            for (q in 1:length(parnames)) temppar[[parnames[q]]] <- best_par[[parnames[q]]]
            lhistpar[[length(lhistpar) + 1]] <- temppar
            
            # Display
            if(display) likdisplay(lhistcycles, lhisthood, slp, R2, aiccorr, temp, max_iter)
          }
        }
        else {
          # Use Metropolis criteria to determine whether to accept a downhill move
          # lh < f, so the code below is a shortcut for exp(-1.0(abs(acc_lh-lh)/t)
          prob <- exp((lh-acc_lh)/temp);  # temp = current temperature
          if (runif(1) < prob) {
            par <- par_copy
            acc_lh <- lh
            nacp[[using]] <- nacp[[using]] + 1
          }
        }
        
        # After the user-specified interval, adjust step so that half of
        # evaluations are accepted
        if (cycles %% range_cycles == 0) {
          for (i in 1:numpars) {
            if (par_step[[i]] != 0) {
              ratio <- nacp[[i]] / ns
              # C controls the adjustment of range - references suggest
              # setting at 2.0.  I've hard-coded the 2 since Charlie doesn't
              # seem to allow setting it anywhere.
              if (ratio > 0.6) par_step[[i]] <- par_step[[i]] * (1.0 + 2*((ratio - 0.6)/0.4))
                else {if (ratio < 0.4) par_step[[i]] <- par_step[[i]] / (1.0+2*((0.4 - ratio)/0.4))}
              if (is.infinite(par_step[[i]])) par_step[[i]] = .Machine$double.xmax
#              if (par_step[[i]] > (par_hi[[i]] - par_lo[[i]])) par_step[[i]] <- par_hi[[i]] - par_lo[[i]]
            }
          }
          # Reset number of times parameters were accepted
          nacp <- rep(0, numpars)
        }
        
        # After numpars * ns * nt, reduce temperature
        if (cycles %% temp_cycles == 0) {
          temp <- temp_red * temp
          # NOTE:  Goffe et al. restart the search at the current best fit each
          # time the temperature drops, but I (CDC) generally don't do this
          # Store current maximum lhood in history list
          lhisttemp[length(lhisttemp)+1] <- temp
          lhistcycles[length(lhistcycles) + 1] <- cycles / numpars
          lhisthood[length(lhisthood) + 1] <- best_lh
          temppar <- list()
          for (q in 1:length(parnames)) temppar[[parnames[q]]] <- best_par[[parnames[q]]]
          lhistpar[[length(lhistpar) + 1]] <- temppar
          # Update display
          if(display) likdisplay(lhistcycles, lhisthood, slp, R2, aiccorr, temp, max_iter)

          # Check to see if likelihood has changed sufficiently to keep going,
          # if the user has specified alternate quitting conditions
          if (length(lhisthood) >= min_drops+1) {
            if (lhisthood[length(lhisthood)] - lhisthood[length(lhisthood) - min_drops] < min_change) {
              iterate <- FALSE
            }
          }
        }
        
        # Also update display every 100 cycles (but don't update if it just was
        # updated)
        if (cycles %% display_cycles == 0 && cycles %% temp_cycles != 0) {
          lhisttemp[length(lhisttemp)+1] <- temp
          lhistcycles[length(lhistcycles) + 1] <- cycles / numpars
          lhisthood[length(lhisthood) + 1] <- best_lh
          temppar <- list()
          for (q in 1:length(parnames)) temppar[[parnames[q]]] <- best_par[[parnames[q]]]
          lhistpar[[length(lhistpar) + 1]] <- temppar
          if(display) likdisplay(lhistcycles, lhisthood, slp, R2, aiccorr, temp, max_iter)
        }

        if (using == numpars) using <- 1 else using <- using + 1
          if ((cycles / numpars) >= max_iter) iterate <- FALSE

        # Here's where we skip to if there was a problem calculating likelihood
        # I don't want to update the cycle number yet -
        # that way we won't miss any scheduled temperature drops or likelihood updates
      }
      
    } #en iterate

#-----------------------------------------
if(support)
{
    # Calculate support limits
    par <- list()
    for (i in 1:length(parnames)) {
      par[[parnames[i]]] <- best_par[[parnames[i]]]
    }
    limits<-support_limits(model = model, par = par, var = var, source_data = source_data, pdf = pdf, delta = delta, slimit = slimit)
    upper_limit<-limits$upper_limits
    lower_limit<-limits$lower_limits
}

    # Update likelihood history one last time
    lhistcycles[length(lhistcycles)] <- cycles / numpars
    lhisthood[length(lhisthood)] <- best_lh
    temppar <- list()
    for (q in 1:length(parnames)) temppar[[parnames[q]]] <- best_par[[parnames[q]]]
    lhistpar[[length(lhistpar)]] <- temppar
    lhisttemp[length(lhisttemp)] <- temp
    if(display) likdisplay(lhistcycles, lhisthood, slp, R2, aiccorr, temp, max_iter)

  }, finally={

    if (exists("neighparam_cleanup")) {
      neighparam_cleanup()
    }

# END TRY CATCH------ 

    # Do output
    
 #---------------------
    
    # Update the par_step copy array with the values
    # from the "flattened" array
    for (i in 1:length(par_step_copy)) {
      if (length(par_step_copy[[i]]) == 1)
        par_step_copy[[i]] <- par_step[[parnames[[i]]]]
      else {
        for (j in 1:length(par_step_copy[[i]])) {
          par_step_copy[[i]][[j]] <- par_step[[paste(parnames[[i]], j, sep="")]]
        }
      }
    }

    # Make par once again just the values that were varied, with the best
    # values; we did this when we calculated support intervals, but do it
    # again in case the user aborted the run
    par <- list()
    for (i in 1:length(parnames)) {
      par[[parnames[i]]] <- best_par[[parnames[i]]]
    }

pdfname<-"User-defined function"

    # Find out what the user called the model
    model_name<-""
    base_all<-ls(.GlobalEnv)
    for (i in 1:length(base_all)) {
      if (identical(get(base_all[i], pos=.GlobalEnv),model)) {
        model_name<-base_all[i]
      }
    }


    # To reduce var to manageable size, remove any data frames
    for (i in 1:length(var)) {if (is.data.frame(var[[i]])) var[[i]]<-"[dataset removed]" }
       
    parhistory <- array(dim=c(length(lhistpar),length(par_template)))
    colnames(parhistory)<-names(par_template)
    for (i in 1:length(lhistpar)) {
      parvals<-vector()
      for (j in 1:length(lhistpar[[i]])) parvals <- c(parvals,lhistpar[[i]][[j]])
      parhistory[i,] <- parvals 
    }
  
      return(list(best_pars = par, var = var, iterations = cycles / numpars,
    source_data = data.frame(source_data, predicted=best_predicted),
    par_step = par_step_copy,
    support_interval_range = slimit,
    upper_limits = upper_limit, lower_limits = lower_limit,
    initial_temp = initial_temp, temp_red = temp_red, ns = ns, nt = nt,
    pdf = pdfname, note = note, model=model_name, 
    max_likeli = best_lh, aic_corr = aiccorr,
    aic = (-2.0*best_lh) + (2*numpars), 
    likeli_hist = data.frame(temp = lhisttemp, iter = lhistcycles, likeli = lhisthood, parhistory)))}
  )
}
