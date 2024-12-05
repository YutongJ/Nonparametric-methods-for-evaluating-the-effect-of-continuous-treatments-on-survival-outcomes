#-----------------------
# simulation data 
#-----------------------
#' Data Generate in Simulation
#'
#' Creates a dataset including Y, Delta, S, Trt, W1 and W2.
#'
#' @param seed seed used in generating dataset.
#' @param n sample size in a study. Default is 300.
#' @param n_a the number of treatment groups in a study. Default is 5
#' @param contA whether if the treatment is continuous. Default is FALSE
#' @param effect Treatment effect as A increases from -1 to 1. Three available options: "null", "monotone", "concave".
#'
#' @return A data.frame consists of Y, A, W1 and W2.
#'
#' @examples
#' dat <- data_generate(123)
#'
#' @export

# seed=123;n=1000
data_generate <- function(seed, 
                          n = 1000, 
                          n_a = NA, 
                          contA = FALSE, 
                          effect){
  # inverse cumulative distribution function (ICDF)
  inverse_F <- function(n, beta) {
    # inverse cumulative distribution 
    out <- (1/beta) * log(exp((log(exp(beta) + 1) - log(exp(-beta) + 1)) * u +
                                log(exp(-beta) + 1)) - 1)
    #if beta == 0, generate from a uniform dist'n
    out[beta == 0] <- runif(sum(beta == 0), -1, 1) 
    return(out)
  }
  
  # set random seed
  set.seed(seed)
  
  # baseline characteristics
  w1 <- runif(n, 1,2)
  w2 <- rbinom(n, 1, 0.5)
  
  # f(w)
  f_w <- -3 + 0.3*w1 + 1.1*w2
  
  # treatment (discrete and purely random)
  if (contA){
    if (effect == "null"){
      a <- runif(n, min = -1, max = 1)
      
      # g(a)
      g_a <- function(a){a}
      
      # outer function for generating event T
      gamma_t <- function(f, g){exp(0.2*f)}
      # observed event time
      t <- rexp(n, rate = gamma_t(f_w, g_a(a)))*10
      # hist(t)
      
      # outer function for generating event T
      gamma_c <- function(f, g){exp(-0.2 + 0.4*f)}
      # observed right-censoring time
      c <- rexp(n, rate = gamma_c(f_w, g_a(a)))*9
      # hist(c)
    }
    if (effect == "monotone"){
      # inverse transform sampling
      u <- runif(n)
      p <- 5*(plogis(-1 + w1 - w2) - 0.5) # plogis - expit
      # treatment
      a <- inverse_F(n, beta = p)
      
      # g(a)
      g_a <- function(a){a}
      
      # outer function for generating event T
      gamma_t <- function(f, g){exp(0.6*f - 0.75*g)}
      # observed event time
      t <- rexp(n, rate = gamma_t(f_w, g_a(a)))*10*0.35
      # hist(t)
      
      # outer function for generating event T
      gamma_c <- function(f, g){exp(-1.2 + 0.4*f - 0.5*g)}
      # observed right-censoring time
      c <- rexp(n, rate = gamma_c(f_w, g_a(a)))*9*0.35
      # hist(c)
    }
    if (effect == "concave"){
      # inverse transform sampling
      u <- runif(n)
      p <- 5*(plogis(-1 + w1 - w2) - 0.5) # plogis - expit
      # treatment
      a <- inverse_F(n, beta = p)
      
      # g(a)
      g_a <- function(a){-2*a^2 + 1.2}
      
      # outer function for generating event T
      gamma_t <- function(f, g){exp(0.6*f - 0.75*g)}
      # observed event time
      t <- rexp(n, rate = gamma_t(f_w, g_a(a)))*10*0.35
      # hist(t)
      
      # outer function for generating event T
      gamma_c <- function(f, g){exp(-1.2 + 0.4*f - 0.5*g)}
      # observed right-censoring time
      c <- rexp(n, rate = gamma_c(f_w, g_a(a)))*9*0.35
      # hist(c)
    }
  }else{
    a <- sample(seq(n_a)-1, size = n, replace = TRUE)
  }
  
  
  # make outcome integer
  c <- ceiling(c)
  t <- ceiling(t)
  # hist(c)
  c[c > 35] <- 35
  
  # observed time-to-event
  y <- pmin(t, c)
  ftype <- as.numeric(y == t)
  # hist(y, breaks = 20)
  # range(t);range(c);table(ftype)
  
  # full data
  dat <- data.frame(id = 1:n,
                    Y = y,
                    T = t,
                    C = c,
                    delta = ftype,
                    A = a,
                    W1 = w1,
                    W2 = w2)
  return(dat)
}


# dat <- data_generate(123, n = 500, n_a = 5)
# hist(dat$Y)



# # for developing purpose
# source(paste0(dir_fun, "DGP_complex_ceiling.R"))
# dat0 <- data_generate(seed = 1, n = 5000, contA = TRUE, effect = "null")
# dat1 <- data_generate(seed = 1, n = 5000, contA = TRUE, effect = "monotone")
# dat2 <- data_generate(seed = 1, n = 5000, contA = TRUE, effect = "concave")
# 
# t0 <- 25
# 
# par(mfrow=c(1,3))
#
# hist(dat0$A)
# hist(dat1$A)
# hist(dat2$A)
#
# plot(x = dat0$A, y = dat0$Y)
# plot(x = dat1$A, y = dat1$Y)
# plot(x = dat2$A, y = dat2$Y)
# 
# test <- ksmooth(x = dat0$A, y = ifelse(dat0$T >= t0, 1, 0))
# plot(test$x, test$y, ylim = c(0, 1), type = "l")
# test <- ksmooth(x = dat1$A, y = ifelse(dat1$T >= t0, 1, 0))
# plot(test$x, test$y, ylim = c(0, 1), type = "l")
# test <- ksmooth(x = dat2$A, y = ifelse(dat2$T >= t0, 1, 0))
# plot(test$x, test$y, ylim = c(0, 1), type = "l")
# 
# 
# plot(x = dat0$A, y = ifelse(dat0$Y >= t0, 1, 0))
# plot(x = dat1$A, y = ifelse(dat1$Y >= t0, 1, 0))
# plot(x = dat2$A, y = ifelse(dat2$Y >= t0, 1, 0))


