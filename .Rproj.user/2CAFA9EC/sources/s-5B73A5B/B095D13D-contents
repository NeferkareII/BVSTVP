#' @export
BVS_TVP <- function(formula,
                    data,
                    niter = 10000,
                    nburn = round(niter / 2),
                    nthin = 1,
                    hyperprior_param,
                    display_progress = TRUE,
                    sv = FALSE,
                    sv_param,
                    starting_vals){


  # Input checking ----------------------------------------------------------

  # default hyperparameter values
  default_hyper <- list(s0 = 1,
                        S0 = 1,
                        tau2 = 10,
                        c0 = 2.5,
                        g0 = 5,
                        G0 = 5 / (2.5 - 1))

  # default sv params
  default_hyper_sv <- list(Bsigma_sv = 1,
                           a0_sv = 5,
                           b0_sv = 1.5,
                           bmu = 0,
                           Bmu = 1)

  # Change hyperprior values if user overwrites them
  if (missing(hyperprior_param)){
    hyperprior_param <- default_hyper
  } else {
    hyperprior_param <- shrinkTVP:::list_merger(default_hyper, hyperprior_param)
  }


  # Same procedure for sv_param
  if (missing(sv_param) | sv == FALSE){
    sv_param <- default_hyper_sv
  } else {
    sv_param <- shrinkTVP:::list_merger(default_hyper_sv, sv_param)
  }

  # Check if all numeric inputs are correct
  to_test_num <- list()


  if (missing(hyperprior_param) == FALSE){
    to_test_num <- c(to_test_num, hyperprior_param)
  }

  if (missing(sv_param) == FALSE){
    to_test_num <- c(to_test_num, sv_param[names(sv_param) != "bmu"])
  }


  bad_inp <- sapply(to_test_num, shrinkTVP:::numeric_input_bad)


  if (any(bad_inp)){
    stand_names <- c(names(default_hyper), names(default_hyper_sv))
    bad_inp_names <- names(to_test_num)[bad_inp]
    bad_inp_names <- bad_inp_names[bad_inp_names %in% stand_names]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a real, positive number"))
  }

  # Check bmu seperately
  if (!is.numeric(sv_param$bmu) | !shrinkTVP:::is.scalar(sv_param$bmu)){
    stop("bmu has to be a real number")
  }


  # Check if all integer inputs are correct
  to_test_int <- c(niter = niter,
                   nburn = nburn,
                   nthin = nthin)

  bad_int_inp <- sapply(to_test_int, shrinkTVP:::int_input_bad)

  if (any(bad_int_inp)){
    bad_inp_names <- names(to_test_int)[bad_int_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single, positive integer"))

  }

  if ((niter - nburn) < 2){
    stop("niter has to be larger than or equal to nburn + 2")
  }

  if (nthin == 0){
    stop("nthin can not be 0")
  }

  if ((niter - nburn)/2 < nthin){
    stop("nthin can not be larger than (niter - nburn)/2")
  }

  # Check if all boolean inputs are correct
  to_test_bool <- c(display_progress = display_progress,
                    sv = sv)

  bad_bool_inp <- sapply(to_test_bool, shrinkTVP:::bool_input_bad)

  if (any(bad_bool_inp)){
    bad_inp_names <- names(to_test_bool)[bad_bool_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single logical value"))

  }

  # Check if formula is a formula
  if (inherits(formula, "formula") == FALSE){
    stop("formula is not of class formula")
  }




  # Formula interface -------------------------------------------------------


  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data"), table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  # Create Vector y
  y <- model.response(mf, "numeric")
  mt <- attr(x = mf, which = "terms")
  # Create Matrix X with dummies and transformations
  x <- model.matrix(object = mt, data = mf)

  # Check that there are no NAs in y and x
  if (any(is.na(y))) {
    stop("No NA values are allowed in response variable")
  }

  if (any(is.na(x))){
    stop("No NA values are allowed in covariates")
  }

  # Get the index
  if (missing(data)){
    index <- zoo::index(y)
  } else {
    index <- zoo::index(data)
  }

  p <- 0

  colnames(x)[colnames(x) == "(Intercept)"] <- "Intercept"
  d <- dim(x)[2]

  # Fuse user starting vals with standard ones
  default_starting_vals <- list(beta_mean_st = rep(0, d),
                                theta_st = rep(1, d),
                                sv_mu_st = -10,
                                sv_phi_st = 0.5,
                                sv_sigma2_st = 1,
                                C0_st = 1,
                                sigma2_st = 1,
                                h0_st = 0)

  if (sv == TRUE){
    default_starting_vals$sigma2_st <- rep(1, length(y))
  }


  # Change starting values of MCMC algorithm if user overwrites them
  if (missing(starting_vals)){
    starting_vals <- default_starting_vals
  } else {
    starting_vals <- shrinkTVP:::list_merger(default_starting_vals, starting_vals)
  }

  # Input check starting vals
  # Check length of vectors
  vec_valued <- c("beta_mean_st",
                  "theta_st")
  bad_length <- sapply(starting_vals[vec_valued], function(x) length(x) != d)

  if (any(bad_length)){
    bad_length_names <- vec_valued[bad_length]
    stop(paste0(paste(bad_length_names, collapse = ", "),
                ifelse(length(bad_length_names) == 1, " has", " have"),
                " to be of length ", d))

  }

  # check sigma2_st seperately
  if (sv == TRUE) {
    if (length(starting_vals$sigma2_st) != length(y)) {
      stop("sigma2_st has to be the same length as y if sv is TRUE")
    }

    num_input_bad <- sapply(starting_vals$sigma2_st, shrinkTVP:::numeric_input_bad)

    if (any(num_input_bad)) {
      stop("sigma2_st may only contain real, positive numbers")
    }
  } else if (shrinkTVP:::numeric_input_bad(starting_vals$sigma2_st)) {
    stop("sigma2_st has to be a real, positive number")
  }

  # Check content of vectors
  vec_valued_pos <- vec_valued[vec_valued != "beta_mean_st"]
  bad_content <- sapply(starting_vals[vec_valued_pos], function(x) any(sapply(x, shrinkTVP:::numeric_input_bad)))

  if (any(bad_content)) {
    bad_content_names <- vec_valued_pos[bad_content]
    stop(paste0(paste(bad_content_names, collapse = ", "), " may only contain real, positive numbers"))

  }

  # Check beta_mean_st seperately
  if (any(sapply(starting_vals$beta_mean_st, shrinkTVP:::numeric_input_bad_))) {
    stop("beta_mean_st may only contain real numbers")
  }

  # Check single values
  to_check_num_st <- names(starting_vals)[!names(starting_vals) %in% c(vec_valued, "h0_st", "sv_mu_st", "sigma2_st")]
  bad_num_st <- sapply(starting_vals[to_check_num_st], shrinkTVP:::numeric_input_bad)

  if (any(bad_num_st)) {
    bad_num_names <- to_check_num_st[bad_num_st]
    stop(paste0(paste(bad_num_names, collapse = ", "),
                ifelse(length(bad_num_names) == 1, " has", " have"),
                " to be a real, positive number"))

  }

  # Check h0_st and sv_mu_st seperately
  if (shrinkTVP:::numeric_input_bad_(starting_vals$h0_st)) {
    stop("h0_st has to be a real number")
  }

  if (shrinkTVP:::numeric_input_bad_(starting_vals$sv_mu_st)) {
    stop("sv_mu_st has to be a real number")
  }

  # Check that sv_phi_st falls between -1 and 1
  if (abs(starting_vals$sv_phi_st) >= 1) {
    stop("sv_phi_st has to be between -1 and 1")
  }


  # Run sampler -------------------------------------------------------------


  runtime <- system.time({
    suppressWarnings({
      res <- BVSTVP_cpp(y,
                        x,
                        niter,
                        nburn,
                        nthin,
                        hyperprior_param$c0,
                        hyperprior_param$g0,
                        hyperprior_param$G0,
                        hyperprior_param$tau2,
                        hyperprior_param$s0,
                        hyperprior_param$S0,
                        display_progress,
                        sv,
                        sv_param$Bsigma_sv,
                        sv_param$a0_sv,
                        sv_param$b0_sv,
                        sv_param$bmu,
                        sv_param$Bmu,
                        starting_vals)
    })
  })


  # Throw an error if the sampler failed
  if (res$internals$success_vals$success == FALSE){
    stop(paste0("The sampler failed at iteration ",
                res$internals$success_vals$fail_iter,
                " while trying to ",
                res$internals$success_vals$fail, ". ",
                "Try rerunning the model. ",
                "If the sampler fails again, try changing the prior to be more informative. ",
                "If the problem still persists, please contact the maintainer: ",
                maintainer("shrinkTVP")))
  } else {
    res$internals$success_vals <- NULL
  }

  # Post process sampler results --------------------------------------------


  if (display_progress == TRUE){
    cat("Timing (elapsed): ", file = stderr())
    cat(runtime["elapsed"], file = stderr())
    cat(" seconds.\n", file = stderr())
    cat(round( (niter + nburn) / runtime[3]), "iterations per second.\n\n", file = stderr())
    cat("Converting results to coda objects and summarizing draws... ", file = stderr())
  }


  # Hack for some methods to work
  res$theta_sr <- sqrt(res$theta)

  # Collapse sigma2 to single vector if sv=FALSE
  if (sv == FALSE){
    res$sigma2 <- matrix(res$sigma2[1, 1, ], ncol = 1)
  }



  # Remove empty storage elements
  res[sapply(res, function(x) 0 %in% dim(x))] <- NULL
  res$MH_diag[sapply(res$MH_diag, function(x) 0 %in% dim(x))] <- NULL


  # Create object to hold prior values
  res$priorvals <- c(hyperprior_param,
                     sv_param)

  # Add data to output
  res[["model"]] <- list()
  res$model$x <- x
  res$model$y <- y
  res$model$formula <- formula
  res$model$xlevels <- .getXlevels(mt, mf)
  res$model$terms <- mt

  res$summaries <- list()

  # add attributes to the individual objects if they are distributions or individual statistics
  nsave <- floor((niter - nburn)/nthin)
  for (i in names(res)){

    attr(res[[i]], "type") <- ifelse(nsave %in% dim(res[[i]]), "sample", "stat")

    # Name each individual sample for plotting frontend
    if (attr(res[[i]], "type") == "sample"){

      if (dim(res[[i]])[2] == d){

        colnames(res[[i]]) <- paste0(i, "_",  colnames(x))

      } else if (dim(res[[i]])[2] == 2 * d){

        colnames(res[[i]]) <- paste0(i, "_", rep( colnames(x), 2))

      } else {

        colnames(res[[i]]) <- i

      }
    }

    # Change objects to be coda compatible
    # Only apply to posterior samples
    if (attr(res[[i]], "type") == "sample"){

      # Differentiate between TVP and non TVP
      if (is.na(dim(res[[i]])[3]) == FALSE){

        # Create a sub list containing an mcmc object for each parameter in TVP case
        dat <- res[[i]]
        res[[i]] <- list()
        for (j in 1:dim(dat)[2]){
          res[[i]][[j]] <- as.mcmc(t(dat[, j, ]), start = niter - nburn, end = niter, thin = nthin)
          colnames(res[[i]][[j]]) <- paste0(i, "_", j, "_", 1:ncol(res[[i]][[j]]))

          # make it of class mcmc.tvp for custom plotting function
          class(res[[i]][[j]]) <- c("mcmc.tvp", "mcmc")

          attr(res[[i]][[j]], "type") <- "sample"

          # Imbue each mcmc.tvp object with index
          attr(res[[i]][[j]], "index") <- index
        }

        if (length(res[[i]]) == 1){
          res[[i]] <- res[[i]][[j]]
          attr(res[[i]][[j]], "index") <- index
        }

        # Make it of type 'sample' again
        attr(res[[i]], "type") <- "sample"

        # Rename
        if (dim(dat)[2] > 1){
          names(res[[i]]) <- colnames(dat)
        }


      } else {

        res[[i]] <- as.mcmc(res[[i]], start = niter - nburn, end = niter, thin = nthin)

      }
    }

    # Create summary of posterior
    if (is.list(res[[i]]) == FALSE & attr(res[[i]], "type") == "sample") {
      if (i != "theta_sr" & !(i == "sigma2" & sv == TRUE) & i != "beta") {
        res$summaries[[i]] <- t(apply(res[[i]], 2, function(x){

          obj <- as.mcmc(x, start = niter - nburn, end = niter, thin = nthin)
          ESS <- tryCatch(coda::effectiveSize(obj),
                          error = function(err) {
                            warning("Calculation of effective sample size failed for one or more variable. This can happen if the prior placed on the model induces extreme shrinkage.")
                            return(NA)
                          }, silent = TRUE)

          return(c("mean" = mean(obj),
                   "sd" = sd(obj),
                   "median" = median(obj),
                   "HPD" = HPDinterval(obj)[c(1, 2)],
                   "ESS" = round(ESS)))
        }))
      } else if (i == "theta_sr") {
        res$summaries[[i]] <- t(apply(res[[i]], 2, function(x){

          obj <- as.mcmc(abs(x), start = niter - nburn, end = niter, thin = nthin)
          ESS <- tryCatch(coda::effectiveSize(obj),
                          error = function(err) {
                            warning("Calculation of effective sample size failed for one or more variable. This can happen if the prior placed on the model induces extreme shrinkage.")
                            return(NA)
                          }, silent = TRUE)

          return(c("mean" = mean(obj),
                   "sd" = sd(obj),
                   "median" = median(obj),
                   "HPD" = HPDinterval(obj)[c(1, 2)],
                   "ESS" = round(ESS)))
        }))
      }
    }
  }


  if (display_progress == TRUE) {
    cat("Done!\n", file = stderr())
  }


  # add some attributes for the methods and plotting
  attr(res, "class") <- "shrinkTVP"
  attr(res, "niter") <- niter
  attr(res, "nburn") <- nburn
  attr(res, "nthin") <- nthin
  attr(res, "sv") <- sv
  attr(res, "colnames") <-  colnames(x)
  attr(res, "index") <- index
  attr(res, "p") <- p



  return(res)
}
