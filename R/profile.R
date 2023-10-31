    # initial_par, the vector of the MLE of the remaining parameters, to 
    #           use as starting values for the optimization
    # fixed_par the value at the fixed parameter will be fixed
    # coord_fixed_par coordinate of the fixed parameter
    # loglik_fun the log likelihood function
    # max_loglik the maximum value of the log likelihood
    # returns the signed square root of the deviance and the profile mle of the remaining parameters
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
        dev <- max(0, dev)
        return(list(score = sign(fixed_par - base_par) * sqrt(dev), par = res$par, deviance = dev, new_lik = lik1))
    }




# ### Estrategia para random effects

# step_size -> comeca com 1.01 na primeira vez

# na etapa seguinte: fix_{n+1} - fix_n -> numerador 
#                    zeta_{n+1} - zeta_{n} -> denominador

#                    absstep * num/denom

# absstep <- cutoff.delta (ou seja 1/8) *  sqrt(qchisq(1 - alphamax, DFlogLik)) # usar o mesmo alphamax (isso eh delta) -> absstep eh delta




# zeta = sgn(fix_{n+1} - fix_ORIGINAL) * sqrt(dev_diff) # Dev diff eh o deviance

# nptot -> DF loglik

# cutoff <- sqrt(qchisq(1 - alphamax, nptot)) -> Esse eh o cutoff verdadeiro



# cap step at max step:  abs(maxmult * num) # num eh o numerador

# maxmult = 10 eh um parametro

# maxpts diz o tamanho inicial da matrix


# minstep (se step for negativo)

### Estrategia para fixed effects -> mesma


# which_par can be a character vector with the names of the parameters, 
# if which par is numeric it must obey the order c(measurement_error, random effects parameters, fixed_effects parameters)

# if which_par contains ".fixed", this stands for all fixed parameters
# and ".random" contains stands for all the random effects.

profile.graph_lme <- function(fitted, which_par = NULL, alphamax = 0.01, maxpts = 100,
            delta = NULL, delta_cutoff = 1/8, maxmult = 10, minstep = 1e-6, optim_method = "L-BFGS-B", parallel = FALSE, n_cores = parallel::detectCores()-1,optim_controls = list()){

                max_loglik <- logLik(fitted)
                df_model <- attr(logLik(fitted),"df")
                par_names <- c("std. dev", names(fitted$coeff$random_effects), names(fitted$coeff$fixed_effects))

                n_par_full <- length(par_names)

                #### REMEMBER TO TAKE EXPONENTIAL IN THE END!!!!!!!
                if(is.null(fitted$stationary)){
                    mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$coeff$random_effects), fitted$coeff$fixed_effects)
                    if(fitted$latent_model$type %in% c("graphLaplacian", "WhittleMatern")){
                        mle_par[2] <- 1/mle_par[2]
                    }
                } else if(fitted$stationary){
                    if(fitted$estimate_nu){
                        mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$coeff$random_effects[1]-0.5), log(fitted$coeff$random_effects[2:length(fitted$coeff$random_effects)]), fitted$coeff$fixed_effects)
                    } else{
                        mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$coeff$random_effects), fitted$coeff$fixed_effects)
                    }                    
                } else{
                    if(fitted$estimate_nu){
                        mle_par <- c(log(fitted$coeff$measurement_error), log(fitted$coeff$random_effects[1]-0.5), fitted$coeff$random_effects[2:length(fitted$coeff$random_effects)], fitted$coeff$fixed_effects)
                    } else{
                        mle_par <- c(log(fitted$coeff$measurement_error), fitted$coeff$random_effects, fitted$coeff$fixed_effects)
                    }
                }


                loglikfun <- fitted$lik_fun

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
                         which_par <- c(which_par, names(fitted$coeff$random_effects))                       
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
                for(par_ in which_par){
                    
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
                        print(res_mat)
                        if(current_row == row_next_par){
                            if(base_par == 0){
                                new_par <- 0.001
                            } else {
                                new_par <- base_par + sign(base_par) * 0.01 * base_par
                            }
                        } else{
                            numer <- res_mat[current_row, 1 + col_num_par] - res_mat[current_row - 1, 1 + col_num_par]
                            denom <- res_mat[current_row, 1] - res_mat[current_row - 1, 1]
                            step <- delta * numer / denom
                            current_par <- res_mat[current_row, 1 + col_num_par]
                            if(step < 0){
                                step <- minstep
                            }
                            max_step <- abs(numer * maxmult)
                            if(abs(step) > max_step){
                                step <- max_step
                            }
                            new_par <- current_par + step
                        }

                        dev_list <-  dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls)
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

                    ## Reordering
                    count_par <- 1

                    res_mat[row_next_par:current_row, ] <- res_mat[current_row:row_next_par,]
                    
                    initial_par <- initial_par_bkp

                    z <- 0
                    while((z > -cutoff) && count_par < maxpts){
                        numer <- res_mat[current_row, 1 + col_num_par] - res_mat[current_row - 1, 1 + col_num_par]
                        denom <- res_mat[current_row, 1] - res_mat[current_row - 1, 1]
                        if(denom == 0){
                            step <- minstep
                        } else{
                            step <- delta * numer / denom
                        }
                        current_par <- res_mat[current_row, 1 + col_num_par]
                        if(step < 0){
                            step <- minstep
                        }
                        max_step <- abs(numer * maxmult)
                        if(abs(step) > max_step){
                            step <- max_step
                        }
                        new_par <- current_par - step

                        dev_list <-  dev_fun_score(initial_par = initial_par, fixed_par = new_par, base_par = base_par, coord_fixed_par = col_num_par, loglik_fun = loglikfun, max_loglik = max_loglik, optim_method = optim_method, optim_controls = optim_controls)
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
            }
