    # initial_par, the vector of the MLE of the remaining parameters, to 
    #           use as starting values for the optimization
    # fixed_par the value at the fixed parameter will be fixed
    # coord_fixed_par coordinate of the fixed parameter
    # loglik_fun the log likelihood function
    # max_loglik the maximum value of the log likelihood
    # returns the signed square root of the deviance and the profile mle of the remaining parameters

    #' @noRd 
    dev_fun_score <- function(initial_par, fixed_par, base_par, coord_fixed_par, loglik_fun, max_loglik, optim_method, optim_controls){
        n_par <- length(initial_par)

        # function of the non-fixed parameters 
        new_loglik <- function(nf_par){
                if(coord_fixed_par == 1){
                    new_par <- c(fixed_par, nf_par)
                } else if (coord_fixed_par == n_par+1){
                    new_par <- c(nf_par, fixed_par)
                } else{
                    new_par <- c(nf_par[1:(coord_fixed_par-1)], fixed_par, nf_par[coord_fixed_par:n_par])
                }             
               return(loglik_fun(new_par))
        }
        res <- optim(initial_par, 
                  new_loglik, method = optim_method,
                  control = optim_controls,
                  hessian = FALSE)             
        lik1 <- -res$val
        dev <- -2*(lik1 - max_loglik)
        dev_0 <- max(0, dev)
        return(list(score = sign(fixed_par - base_par) * sqrt(dev_0), par = res$par, deviance = dev, new_lik = lik1))
    }


# which_par can be a character vector with the names of the parameters, 
# if which par is numeric it must obey the order c(measurement_error, random effects parameters, fixed_effects parameters)

# if which_par contains ".fixed", this stands for all fixed parameters
# and ".random" contains stands for all the random effects.

# Inspired by lme4::profile.merMod

#' @title Profile method for `graph_lme` objects
#' @description Function providing a summary of several informations/characteristics of a metric graph object.
#' @aliases profile profile.graph_lme
#' @param fitted A fitted model using `graph_lme()`.
#' @param which_par A character vector indicating which parameters to profile. `NULL` means all parameters. `.fixed` will act as a vector of all fixed effects, `.random` will act as a vector of all random effects.
#' @param parameterization For Matern and Graph Laplacian models. Which parameterization to use? `matern` or `spde`? The default is `matern`.
#' @param alphamax	a number in `(0,1)` such that `1 - alphamax` is the maximum alpha value for likelihood ratio confidence regions; used to establish the range of values to be profiled.
#' @param maxpts maximum number of points (in each direction, for each parameter) to evaluate in attempting to construct the profile.
#' @param delta	stepping scale for deciding on next point to profile.
#' @param delta_cutoff stepping scale (see delta) expressed as a fraction of the target maximum value of the profile on the square-root-deviance scale.
#' @param minstep The minimum step to be taken.
#' @param verbose Should update messages be printed with updates?
#' @param optim_method Which optimization method to pass to `optim`.
#' @param optim_controls Additional controls to be passed to `optim`
#' @param parallel Should parallel optimization be used?
#' @param n_cores Number of cores to be used in case of parallel optimization.
#' @param maxmult maximum multiplier of the original step size allowed, defaults to 10.
#' @param ... Currently not used.
#' @return A `data.frame` containing the profiled parameters.
#' @method profile graph_lme
#' @noRd

profile.graph_lme <- function(fitted, which_par = NULL, parameterization = "matern", alphamax = 0.01, maxpts = 100,
            delta = NULL, delta_cutoff = 1/8, maxmult = 10, minstep = 1e-6, verbose = FALSE,
            optim_method = "L-BFGS-B", parallel = FALSE, n_cores = parallel::detectCores()-1,optim_controls = list(),...){

                max_loglik <- logLik(fitted)
                df_model <- attr(logLik(fitted),"df")
                par_names_fixed <- names(fitted$coeff$fixed_effects)
                loglikfun <- fitted$lik_fun

                change_par <- FALSE
                include_nu <- FALSE

                if(is.null(fitted$stationary)){
                    mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$coeff$random_effects), fitted$coeff$fixed_effects)
                    if(fitted$latent_model$type %in% c("graphLaplacian", "WhittleMatern")){
                        if(parameterization == "spde"){
                            mle_par[2] <- 1/mle_par[2]
                            par_names_random <- names(fitted$coeff$random_effects)
                        } else{
                            mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$matern_coeff$random_effects), fitted$coeff$fixed_effects)
                            change_par <- TRUE
                            par_names_random <- names(fitted$matern_coeff$random_effects)
                        }
                    }
                } else if(fitted$stationary){
                    if(fitted$estimate_nu){
                        if(parameterization == "spde"){
                            par_names_random <- names(fitted$coeff$random_effects)
                            mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$coeff$random_effects[1]-0.5), log(fitted$coeff$random_effects[2:length(fitted$coeff$random_effects)]), fitted$coeff$fixed_effects)
                        } else{
                            mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$matern_coeff$random_effects[1]-0.5), log(fitted$matern_coeff$random_effects[2:length(fitted$matern_coeff$random_effects)]), fitted$coeff$fixed_effects)
                            par_names_random <- names(fitted$matern_coeff$random_effects)                     
                            change_par <- TRUE
                            include_nu <- TRUE
                        }
                    } else{
                        if(parameterization == "spde"){
                            mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$coeff$random_effects), fitted$coeff$fixed_effects)
                            par_names_random <- names(fitted$coeff$random_effects)                            
                        } else{
                            mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$matern_coeff$random_effects), fitted$coeff$fixed_effects)
                            par_names_random <- names(fitted$matern_coeff$random_effects)           
                            change_par <- TRUE                                      
                        }
                    }                    
                } else{
                    par_names_random <- names(fitted$coeff$random_effects)                        
                    if(fitted$estimate_nu){
                        if(parameterization == "spde"){
                            mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$coeff$random_effects[1]-0.5), fitted$coeff$random_effects[2:length(fitted$coeff$random_effects)], fitted$coeff$fixed_effects)
                        } else{
                            mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$coeff$random_effects[1]), fitted$coeff$random_effects[2:length(fitted$coeff$random_effects)], fitted$coeff$fixed_effects)
                            par_names_random[1] <- "nu"
                        }
                    } else{
                        mle_par <- c(log(fitted$coeff$measurement_error), fitted$coeff$random_effects, fitted$coeff$fixed_effects)
                    }
                }

                if(change_par){
                    if(include_nu){
                        loglikfun <- function(theta){
                            nu <- exp(theta[2])
                            alpha <- nu + 0.5
                            theta[2] <- log(alpha)
                            sigma <- exp(theta[3])
                            range <- exp(theta[4])
                            kappa <- sqrt(8 * nu) / range
                            theta[4] <- log(kappa)
                            tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) * (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
                            theta[3] <- log(sigma)
                            return(loglikfun(theta))
                        }
                    } else{
                        loglikfun <- function(theta){
                            if(!is.null(fitted$nu)){
                                nu <- fitted$nu
                            } else{
                                alpha <- fitted$latent_model$alpha
                                nu <- alpha - 0.5
                            }
                            sigma <- exp(theta[2])
                            range <- exp(theta[3])
                            kappa <- sqrt(8 * nu) / range
                            theta[3] <- log(kappa)
                            tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) * (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
                            theta[2] <- log(sigma)
                            return(loglikfun(theta))                            
                        }
                    }
                }

                par_names <- c("std. dev", par_names_random, par_names_fixed)
                n_par_full <- length(par_names)

                cutoff <- sqrt(qchisq(1 - alphamax, df_model)) 

                if(is.null(delta)){
                    delta <- delta_cutoff * cutoff
                } 

                # Processing which_names

                if(is.null(which_par)){
                    which_par <- par_names
                } else {
                    if(is.numeric(which_par)){
                        which_par <- par_names[which_par]
                    }
                    if(".fixed" %in% which_par){
                        which_par <- setdiff(which_par, ".fixed")
                        which_par <- c(which_par, names(fitted$coeff$fixed_effects))
                    }
                    if(".random"%in% which_par){
                         which_par <- setdiff(which_par, ".random")
                        which_par <- c(which_par, par_names_random)
                            }

                    if(any(!(which_par %in% par_names))){
                        wrong_names <- which_par[!(which_par %in% par_names)]
                        which_par <- setdiff(which_par, wrong_names)
                        if(length(wrong_names)>1){
                            warning(paste(wrong_names, "are not valid parameters names and were, thus, excluded."))
                        } else{
                            warning(paste(wrong_names, "is not a valid parameter names and was, thus, excluded."))
                        }
                    }
                }
                
                # Creating initial matrix
                res_mat <- matrix(nrow = maxpts, ncol = length(par_names)+1)
                row_next_par <- 1
                col_par_names <- vector(mode = "character", length = maxpts)

                if(verbose){
                    message("Starting to compute the profiles...")
                }
                for(par_ in which_par){
                    if(verbose){
                        message(paste("Profiling parameter",par_))
                    }

                    
                    # We start by filling the first row
                    # the first column is the quantile, the last column is the parameter name
                    count_par <- 1
                    row_tmp <- c(0, mle_par)
                    if(row_next_par > nrow(res_mat)){
                        res_mat <- rbind(res_mat, row_tmp)
                        col_par_names <- c(col_par_names, par_)
                    } else{
                        res_mat[row_next_par,1:(length(par_names)+1)] <- row_tmp
                        col_par_names[row_next_par] <- par_
                    }
                    current_row <- row_next_par

                    base_par <- mle_par[[par_]]
                    col_num_par <- which(par_names == par_)
                    initial_par_bkp <- initial_par <- mle_par[-col_num_par]

                    # We start by adding the parameters in the positive direction
                    z <- 0
                    while((z < cutoff) && count_par < maxpts){
                        if(current_row == row_next_par){
                            current_par <- base_par
                            if(base_par == 0){
                                new_par <- 0.001
                            } else {
                                new_par <- base_par + sign(base_par) * 0.01 * base_par
                            }
                        } else{
                            numer <- res_mat[current_row, 1 + col_num_par] - res_mat[current_row - 1, 1 + col_num_par]
                            denom <- res_mat[current_row, 1] - res_mat[current_row - 1, 1]
                            current_par <- res_mat[current_row, 1 + col_num_par]
                            if(numer == 0){
                                numer <- minstep * sign(denom)
                            }
                            max_step <- abs(numer * maxmult)
                            if(denom == 0){
                                step <- max_step
                            } else { 
                                step <- delta * numer / denom
                            }
                             if(step < 0){
                                step <- minstep
                             }

                            if(abs(step) > max_step){
                                step <- max_step
                            }
                            new_par <- current_par + step
                        }

                        dev_list <- withCallingHandlers(tryCatch(dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls), error = function(e){return(NA)}), 
                                    warning = function(w){invokeRestart("muffleWarning")})

                        possible_methods <- c("Nelder-Mead", "L-BFGS-B", "BFGS", "CG")
                        possible_methods <- setdiff(possible_methods, optim_method)
                        new_method <- optim_method

                        while(all(is.na(dev_list)) && (length(possible_methods) >= 1)){
                            if(verbose){
                                message(paste(new_method, "failed for one profile. Trying", possible_methods[1]))
                            }
                            new_method <- possible_methods[1]
                            dev_list <- withCallingHandlers(tryCatch(dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls), error = function(e){return(NA)}), 
                                warning = function(w){invokeRestart("muffleWarning")})
                            possible_methods <- setdiff(possible_methods, new_method)
                        }

                        if(all(is.na(dev_list))){
                            if(verbose){
                                message("The optimizers failed, we will try again with a new step size.")
                            }
                                # Trying a very small step
                                new_par <- current_par + minstep

                                dev_list <- withCallingHandlers(tryCatch(dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls), error = function(e){return(NA)}), 
                                        warning = function(w){invokeRestart("muffleWarning")})

                                possible_methods <- c("Nelder-Mead", "L-BFGS-B", "BFGS", "CG")
                                possible_methods <- setdiff(possible_methods, optim_method)
                                new_method <- optim_method
                                while(all(is.na(dev_list)) && (length(possible_methods) >= 1)){
                                    if(verbose){
                                        message(paste(new_method, "failed for one profile. Trying", possible_methods[1]))
                                    }
                                    new_method <- possible_methods[1]
                                    dev_list <- withCallingHandlers(tryCatch(dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls), error = function(e){return(NA)}), 
                                    warning = function(w){invokeRestart("muffleWarning")})
                                    possible_methods <- setdiff(possible_methods, new_method)
                            }
                        }                        

                        if(all(is.na(dev_list))){
                            stop("The profiling failed. Try changing the controls parameters.")
                        } 

                        initial_par <- dev_list[["par"]]
                        z <- dev_list[["score"]]
                        if(col_num_par == 1){
                            row_par <- c(new_par, initial_par)
                        } else if (col_num_par == n_par_full){
                            row_par <- c(initial_par, new_par)
                        } else{
                            row_par <- c(initial_par[1:(col_num_par-1)], new_par, initial_par[col_num_par:(n_par_full-1)])
                        }
                        dev <- dev_list[["deviance"]]
                        if(dev < 0){
                            warning("Larger likelihood found. Updating. The larger likelihood will be saved as an attribute.")
                            mle_par <- row_par
                            max_loglik <- dev_list[["new_lik"]]
                        }
                        row_tmp <- c(z, row_par)
                        current_row <- current_row + 1
                   
                        if(current_row > nrow(res_mat)){
                            res_mat <- rbind(res_mat, row_tmp)
                            col_par_names <- c(col_par_names, par_)
                        } else{
                            res_mat[current_row,1:(length(par_names)+1)] <- row_tmp
                            col_par_names[current_row] <- par_
                        }
                      count_par <- count_par + 1
                    }
                    ## Adding parameters in the negative direction
                    
                    if(verbose){
                        message(paste("Profiling", par_, "in the negative direction."))
                    }
                    ## Reordering
                    count_par <- 1

                    res_mat[row_next_par:current_row, ] <- res_mat[current_row:row_next_par,]
                    
                    initial_par <- initial_par_bkp

                    z <- 0
                    while((z > -cutoff) && count_par < maxpts){
                        numer <- res_mat[current_row, 1 + col_num_par] - res_mat[current_row - 1, 1 + col_num_par]
                        denom <- res_mat[current_row, 1] - res_mat[current_row - 1, 1]
                        if(numer == 0){
                                numer <- minstep * sign(denom)
                        }                        
                        max_step <- abs(numer * maxmult)                        
                        if(denom == 0){
                                step <- max_step  
                        } else{
                            step <- delta * numer / denom
                        }
                        current_par <- res_mat[current_row, 1 + col_num_par]
                        if(step < 0){
                            step <- minstep
                        }

                        if(abs(step) > max_step){
                            step <- sign(step) * max_step
                        }
                        new_par <- current_par - step

                        dev_list <- withCallingHandlers(tryCatch(dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls), error = function(e){return(NA)}), 
                                    warning = function(w){invokeRestart("muffleWarning")})

                        possible_methods <- c("Nelder-Mead", "L-BFGS-B", "BFGS", "CG")
                        possible_methods <- setdiff(possible_methods, optim_method)
                        new_method <- optim_method

                        while(all(is.na(dev_list)) && (length(possible_methods) >= 1)){
                            if(verbose){
                                message(paste(new_method, "failed for one profile. Trying", possible_methods[1]))
                            }                            
                            new_method <- possible_methods[1]
                            dev_list <- withCallingHandlers(tryCatch(dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls), error = function(e){return(NA)}), 
                                warning = function(w){invokeRestart("muffleWarning")})
                            possible_methods <- setdiff(possible_methods, new_method)
                        }

                        if(all(is.na(dev_list))){
                            if(verbose){
                                message("The optimizers failed, we will try again with a new step size.")
                            }
                                # Trying a very small step
                                new_par <- current_par + minstep

                                dev_list <- withCallingHandlers(tryCatch(dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls), error = function(e){return(NA)}), 
                                        warning = function(w){invokeRestart("muffleWarning")})

                                possible_methods <- c("Nelder-Mead", "L-BFGS-B", "BFGS", "CG")
                                possible_methods <- setdiff(possible_methods, optim_method)
                                new_method <- optim_method
                                while(all(is.na(dev_list)) && (length(possible_methods) >= 1)){
                                    if(verbose){
                                        message(paste(new_method, "failed for one profile. Trying", possible_methods[1]))
                                    }
                                    new_method <- possible_methods[1]
                                    dev_list <- withCallingHandlers(tryCatch(dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls), error = function(e){return(NA)}), 
                                    warning = function(w){invokeRestart("muffleWarning")})
                                    possible_methods <- setdiff(possible_methods, new_method)
                            }
                        }                              

                        if(all(is.na(dev_list))){
                            stop("The profiling failed. Try changing the controls parameters.")
                        } 

                        initial_par <- dev_list[["par"]]
                        z <- dev_list[["score"]]
                        if(col_num_par == 1){
                            row_par <- c(new_par, initial_par)
                        } else if (col_num_par == n_par_full){
                            row_par <- c(initial_par, new_par)
                        } else{
                            row_par <- c(initial_par[1:(col_num_par-1)], new_par, initial_par[col_num_par:(n_par_full-1)])
                        }

                        dev <- dev_list[["deviance"]]
                        if(dev < 0){
                            warning("Larger likelihood found. Updating. The larger likelihood will be saved as an attribute.")
                            mle_par <- row_par
                            max_loglik <- dev_list[["new_lik"]]
                        }                        
                        row_tmp <- c(z, row_par)
                        current_row <- current_row + 1
                   
                        if(current_row > nrow(res_mat)){
                            res_mat <- rbind(res_mat, row_tmp)
                            col_par_names <- c(col_par_names, par_)
                        } else{
                            res_mat[current_row,1:(length(par_names)+1)] <- row_tmp
                            col_par_names[current_row] <- par_
                        }
                      count_par <- count_par + 1
                    }

                    row_next_par <- current_row + 1
                }
                colnames(res_mat) <- c(".zeta", par_names)
                res_mat <- as.data.frame(res_mat)
                res_mat[[".par"]] <- col_par_names
                ord_idx <- order(res_mat[[".par"]], res_mat[[".zeta"]])
                res_mat <- res_mat[ord_idx,]
                rownames(res_mat) <- 1:nrow(res_mat)
                attr(res_mat, "max_loglik") <- max_loglik
                attr(res_mat, "mle") <- mle_par

                # Rescaling back
                if(fitted$latent_model$type %in% c("graphLaplacian", "WhittleMatern")){
                    res_mat[,2] <- exp(res_mat[,2])
                    if(is.null(fitted$stationary)){
                        res_mat[,3:(length(fitted$coeff$random_effects)+1)] <- exp(res_mat[,3:(length(fitted$coeff$random_effects)+1)])
                        if(fitted$latent_model$type %in% c("graphLaplacian", "WhittleMatern")){
                            res_mat[,3] <- 1/res_mat[,3]
                        }
                    } else if(fitted$stationary){
                        if(fitted$estimate_nu){
                            res_mat[,3] <- exp(res_mat[,3]) + 0.5
                            res_mat[,4:(length(fitted$coeff$random_effects)+1)] <- exp( res_mat[,4:(length(fitted$coeff$random_effects)+1)])

                        } else{
                            res_mat[,3:(length(fitted$coeff$random_effects)+1)] <- exp( res_mat[,3:(length(fitted$coeff$random_effects)+1)])
                        }                    
                    } else{
                        if(fitted$estimate_nu){
                            res_mat[,3] <- exp(res_mat[,3]) + 0.5
                        }
                    }
                } else{
                    res_mat[,2:(length(fitted$coeff$random_effects)+1)] <- exp(res_mat[,2:(length(fitted$coeff$random_effects)+1)])
                }
            return(res_mat)
        }
