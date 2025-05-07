# ABC rejection algorithm
library(EasyABC)


set.seed(1)
toy_model <- function(x){ 2 * x + 5 + rnorm(1,0,0.1) }
toy_prior <- list(c("unif",0,10))
sum_stat_obs <- 6.5
tolerance <- 0.2
n <- 1000 

ABC_rej <- ABC_rejection(model=toy_model, prior=toy_prior, 
                         nb_simul=n, 
                         summary_stat_target=sum_stat_obs, 
                         tol=tolerance)

mean(ABC_rej$param)
boxplot(ABC_rej$param, 
        main="Estimation result of ABC Rejection algorithm",
        ylab="Parameter Value",
        col="lightblue")

true_param <- (sum_stat_obs - 5)/2
abline(h=true_param, col="red", lty=2)
legend("topright", legend="Theoretical True Value", col="red", lty=2)

cat("Number of accepted particles:", length(ABC_rej$param), "\n")
cat("Acceptance rate:", length(ABC_rej$param)/n, "\n")
cat("Mean of posterior:", mean(ABC_rej$param), "\n")
cat("Theoretical true value:", true_param, "\n")








# ABC-PMC algorithm
par(mfrow = c(1, 1))
set.seed(1)
toy_model <- function(x){ 2 * x + 5 + rnorm(1,0,0.1) }
toy_prior <- list(c("unif",0,10))
sum_stat_obs <- 6.5
tolerance <- c(1.5, 0.5,0.1)

ABC_Beaumont <- ABC_sequential(method="Beaumont", 
                               model=toy_model, 
                               prior=toy_prior,
                               nb_simul=1000, 
                               summary_stat_target=sum_stat_obs, 
                               tolerance_tab=tolerance,
                               verbose=TRUE) 
mean(ABC_Beaumont$param)
sd(ABC_Beaumont$param)
n_iterations <- length(ABC_Beaumont$intermediary)
cat("Number of iterations:", n_iterations, "\n")
param_values <- list()
stat_values <- list()


for (i in 1:n_iterations) {
  if (!is.null(ABC_Beaumont$intermediary[[i]]$posterior)) {
    posterior <- ABC_Beaumont$intermediary[[i]]$posterior
    param_values[[i]] <- posterior[, 2]
    stat_values[[i]] <- posterior[, 3]
  } else {
    cat("Structure of iteration", i, ":\n")
    print(str(ABC_Beaumont$intermediary[[i]]))
  }
}

par(mfrow=c(1,1))

boxplot(param_values, 
        main="Estimation result of ABC-PMC",
        xlab="Iteration", 
        ylab="Parameter Value",
        col="lightblue")









# Adaptive ABC-PMC
set.seed(1)
toy_model <- function(x){ 2 * x + 5 + rnorm(1,0,0.1) }
toy_prior <- list(c("unif",0,10))
sum_stat_obs <- 6.5
pacc <- 0.4
ABC_Lenormand <- ABC_sequential(method="Lenormand", model=toy_model, prior=toy_prior,
                                nb_simul=1000, summary_stat_target=sum_stat_obs, 
                                p_acc_min=pacc, verbose=TRUE)


mean(ABC_Lenormand$param)
n_iterations <- length(ABC_Lenormand$intermediary)

param_values <- list()
stat_values <- list()

for (i in 1:n_iterations) {

  posterior <- ABC_Lenormand$intermediary[[i]]$posterior
  

  param_values[[i]] <- posterior[, 2]
  
  stat_values[[i]] <- posterior[, 3]
}


par(mfrow=c(1,1))  # Set 2x1 plot layout, two plots vertically

boxplot(param_values, 
        main="Estimation result of adaptive ABC-PMC",
        xlab="Iteration", 
        ylab="Parameter Value",
        col="lightblue")
abline(h=(sum_stat_obs - 5)/2, col="red", lty=2)  # Theoretical true parameter value
legend("topleft", legend="Theoretical True Value", col="red", lty=2)

