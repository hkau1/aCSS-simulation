library(dplyr)
library(truncnorm)
library(ggplot2) 

mu_0 = 0
sigma_1 = 1
sigma_2 = 2
sigma = 1
tol = 1e-6
M = 10
n0 = 50
n1 = 50
N = 1000

signal = seq(from = 0, to = 1, by = 0.2)
nsignal = length(signal)
contam = seq(from = 0, to = 5, by = 1)
n_contam = length(contam)

solve_quadratic <- function(a, b, c) {
  if (a == 0) {
    if (b == 0) {
      stop("Invalid equation: a and b cannot both be zero.")
    }
    root <- -c / b
    return (root)
  }

  discriminant <- b^2 - 4 * a * c
  
  if (discriminant > 0) {
    root <- (-b + sqrt(discriminant)) / (2 * a)
    return (root)
  } else {
    stop("Invalid equation: discriminant not positive.")
  }
}

hessian <- function(X1, X2, thetahat){
  mu = thetahat[1]
  sig2_0 = thetahat[2]
  sig2_1 = thetahat[3]
  m0 = length(X1)
  m1 = length(X2)
  mat = matrix(c(m0/sig2_0+m1/sig2_1,
           1/sig2_0^2*sum(X1-mu),
           1/sig2_1^2*sum(X2-mu),
           1/sig2_0^2*sum(X1-mu),
           -m0/2/sig2_0^2+1/sig2_0^3*sum((X2-mu)^2),
           0,
           1/sig2_1^2*sum(X2-mu),
           0,
           -m1/2/sig2_1^2+1/sig2_1^3*sum((X2-mu)^2)),3,3)
  return (det(mat))
}

gradient = function(X1, X2, thetahat){
  mu = thetahat[1]
  sig2_0 = thetahat[2]
  sig2_1 = thetahat[3]
  m0 = length(X1)
  m1 = length(X2)
  c(-1/sig2_0*sum(X1-mu)
    -1/sig2_1*sum(X2-mu),
    m0/2/sig2_0-1/2/sig2_0^2*sum((X1-mu)^2),
    m1/2/sig2_1-1/2/sig2_1^2*sum((X2-mu)^2))
}

SSOSP <- function(X1_copies, X2_copies, thetahat, trim_value,  sigma, tol = 1e-6){
  mu = thetahat[1]
  SSOSP = 0
  w_new = -gradient(X1_copies, X2_copies, thetahat)/sigma
  thetah_new = theta_hat(X1_copies, X2_copies,w_new,sigma, rep(0,2))$theta_hat
  if((max(abs(X1_copies - mu)) < trim_value[1] - tol) && (max(abs(X2_copies - mu)) < trim_value[2] - tol) && (hessian(X1_copies, X2_copies, thetahat) > 0) && (sum(((thetah_new)-thetahat)^2)<=tol)){
    SSOSP = hessian(X1_copies, X2_copies, thetahat)
  }
  return (SSOSP * exp(-3 * sum(gradient(X1_copies, X2_copies, thetahat)^2)/2*sigma^2))
}

pval_compute <- function(thetahat, trim_value, X1, X2, sigma, N = 1000){
  T0 = mean(X2) - mean(X1)
  SSOSP_0 = SSOSP(X1, X2, thetahat, trim_value, sigma)
  Tcopies = SSOSPcopies = rep(0, N)
  if (SSOSP_0 == 0){
    return (1)
  }
  else{  for(i in 1:N){
    X1_copies = rtruncnorm(length(X1), a = thetahat[1] - trim_value[1], b = thetahat[1] + trim_value[1], mean =  thetahat[1], sd = sqrt(thetahat[2]))
    X1_copies = thetahat[1] + sqrt(sum((X1-thetahat[1])^2)/sum((X1_copies-thetahat[1])^2))*(X1_copies-thetahat[1])
    X2_copies = rtruncnorm(length(X2), a =thetahat[1] - trim_value[2], b = thetahat[1] + trim_value[2], mean =  X2, sd = thetahat[3])
    X2_copies = thetahat[1] + sqrt(sum((X2-thetahat[1])^2)/sum((X2_copies-thetahat[1])^2))*(X2_copies-thetahat[1])
    SSOSPcopies[i] = SSOSP(X1_copies, X2_copies, thetahat, trim_value, sigma)
    Tcopies[i] = (mean(X2_copies) - mean(X1_copies)>=T0)
  }
  return ((sum(Tcopies * SSOSPcopies)+SSOSP_0)/(sum(SSOSPcopies)+SSOSP_0))
  }
}

theta_hat <- function(X1, X2, W, sigma, s, max_iter = 500, tol = 1e-6) {
  mu <- (mean(X1) + mean(X2))/2
  sigma_1 <- mean((X1 - mu)^2)
  sigma_2 <- mean((X2 - mu)^2)
  log_likelihood <- -Inf
  
  for (iter in 1:max_iter) {
    trimmed_X1 <- X1
    trimmed_X2 <- X2
    
    if(s[1] > 0){
      weight_1 <- abs(X1 - mu)
      indice_1 <- order(weight_1, decreasing = TRUE)
      keep_indice_1 <- indice_1[-(1:s[1])]
      trimmed_X1 <- X1[keep_indice_1]
    }
    
    if(s[2] > 0){
      weight_2 <- abs(X2 - mu)
      indice_2 <- order(weight_2, decreasing = TRUE)
      keep_indice_2 <- indice_2[-(1:s[2])]
      trimmed_X2 <- X2[keep_indice_2]
    }
    
    len_1 <- length(trimmed_X1)
    len_2 <- length(trimmed_X2)
    new_mu <- (sum(trimmed_X1)/sigma_1 + sum(trimmed_X2)/sigma_2 - sigma * W[1])/(len_1/sigma_1 + len_2/sigma_2)
    new_sigma_1 <- solve_quadratic(2 * sigma * abs(W[2]), len_1, -sum((trimmed_X1 - new_mu)^2))
    new_sigma_2 <- solve_quadratic(2 * sigma * abs(W[3]), len_2, -sum((trimmed_X2 - new_mu)^2))
    new_log_likelihood <- sum(log(dnorm(trimmed_X1, mean = new_mu, sd = sqrt(new_sigma_1)))) + sum(log(dnorm(trimmed_X2, mean = new_mu, sd = sqrt(new_sigma_2))))
    if (abs(new_mu - mu) < tol) {
      # cat("Converged in", iter, "iterations.\n")
      break
    }
    mu <- new_mu
    sigma_1 <- new_sigma_1
    sigma_2 <- new_sigma_2
    log_likelihood <- new_log_likelihood
  }
  trim_value = rep(Inf, 2)
  if(s[1]>0){
    trim_value[1] = weight_1[indice_1[s[1]]]
  }
  if(s[2]>0){
    trim_value[2] = weight_2[indice_2[s[2]]]
  }
  return (list(theta_hat = c(mu, sigma_1, sigma_2), trimmed_X1 = trimmed_X1, trimmed_X2 = trimmed_X2, trim_value = trim_value))
}

score_test <- function(X1,X2){
  thetahat = theta_hat(X1, X2, rep(0,3), 0, s = rep(0,2))$theta_hat
  teststat = (mean(X2)-mean(X1))^2/(thetahat[2]/length(X1)+thetahat[3]/length(X2))
  pval = 1-pchisq(teststat,1)
}

pval_MLE = pval_MTLE = pval_MTLE_oracle = pval_Ttest = pval_Ttest_oracle = pval_score = pval_score_oracle = rep(0, M * nsignal * n_contam)
dim(pval_MLE) = dim(pval_MTLE) = dim(pval_MTLE_oracle) = dim(pval_Ttest) = dim(pval_Ttest_oracle) = dim(pval_score) = dim(pval_score_oracle) = c(M, nsignal, n_contam)

starttime = Sys.time()
for (i in 1:M){
  if (i%%1 == 0){
    print(i)
    print(Sys.time() - starttime)
    starttime = Sys.time()
  }
  for(j in 1:nsignal){
    for(k in 1:n_contam){
      set.seed(1 + 10 * k + 100 * j + 1000 * i)
      s = c(contam[k],0)
      X1 = c(rnorm(n0 - s[1], mean = mu_0, sd = sqrt(sigma_1)), 3 + abs(rt(s[1], 1)))
      X2 = c(rnorm(n1 - s[2], mean = signal[j], sd = sqrt(sigma_2)), 3 + rt(s[2],1))
      W = rnorm(3, sd = sqrt(1/3))
      
      # MLE part
      res_MLE = theta_hat(X1, X2, W, sigma, rep(0,2))
      pval_MLE[i,j,k] = pval_compute(res_MLE$theta_hat, res_MLE$trim_value, res_MLE$trimmed_X1, res_MLE$trimmed_X2, sigma)
      
      # MTLE part
      res_MTLE_oracle = theta_hat(X1, X2, W, sigma,s)
      pval_MTLE_oracle[i,j,k] = pval_compute(res_MTLE_oracle$theta_hat, res_MTLE_oracle$trim_value, res_MTLE_oracle$trimmed_X1, res_MTLE_oracle$trimmed_X2, sigma)
      
      res_MTLE = theta_hat(X1, X2, W, sigma, c(5, 0))
      pval_MTLE[i,j,k] = pval_compute(res_MTLE$theta_hat, res_MTLE$trim_value, res_MTLE$trimmed_X1, res_MTLE$trimmed_X2, sigma)
      
      pval_Ttest[i,j,k] = t.test(X1, X2, alternative = "less")$p.value
      pval_Ttest_oracle[i,j,k] = t.test(X1[1:(n0 - s[1])], X2[1:(n1 - s[2])], alternative = "less")$p.value
      
      pval_score[i,j,k]=score_test(X1, X2)
      pval_score_oracle[i,j,k] = score_test(X1[1:(n0 - s[1])], X2[1:(n1 - s[2])])
    }
  }
}


alpha = 0.1
power_MLE = apply(pval_MLE, c(2,3), function(x) sum(as.numeric(x<alpha))/ M)
power_MTLE = apply(pval_MTLE, c(2,3), function(x) sum(as.numeric(x<alpha))/ M)
power_MTLE_oracle = apply(pval_MTLE_oracle, c(2,3), function(x) sum(as.numeric(x<alpha))/ M)
power_Ttest = apply(pval_Ttest, c(2,3), function(x) sum(as.numeric(x<alpha))/ M)
power_Ttest_oracle = apply(pval_Ttest_oracle, c(2,3), function(x) sum(as.numeric(x<alpha))/ M)
power_score = apply(pval_score, c(2,3), function(x) sum(as.numeric(x<alpha))/ M)
power_score_oracle = apply(pval_score_oracle, c(2,3), function(x) sum(as.numeric(x<alpha))/ M)

print("power_MLE")
print(power_MLE)
print("power_MTLE")
print(power_MTLE)
print("power_MTLE_oracle")
print(power_MTLE_oracle)
print("power_Ttest")
print(power_Ttest)
print("power_Ttest_oracle")
print(power_Ttest_oracle)

for (i in 1:n_contam){
  plot_data <- data.frame(
    signal = rep(signal, times = 5),
    power = c(unlist(power_MLE[,i]), unlist(power_MTLE[,i]), unlist(power_MTLE_oracle[,i]), unlist(power_Ttest[,i]),  unlist(power_Ttest_oracle[,i])),
    method =rep(c("MLE", "MTLE", "MTLE_oracle", "t_test",  "t_test_oracle"), each = 6)) 
  plot_data$method <- factor(plot_data$method, levels = c("MLE", "MTLE", "MTLE_oracle", "t_test",  "t_test_oracle"))
  plot_data$delta <- sqrt(plot_data$power * (1 - plot_data$power) / M) 
  p <- ggplot(plot_data, aes(x = signal, y = power, linetype = method, color = method)) +
    geom_line(aes( color = method), size = 0.5) +
    geom_ribbon(aes(ymin = power - delta, ymax = power +  delta, fill = method), alpha = 0.3, color = NA) +
    scale_color_manual(values=c("#D8C033","#B9A2A3","#225F6D","#B0946F",
                                "#D8C033","#B9A2A3","#225F6D")) +
    labs(x = expression(beta[0]), y = "power") +
    xlim(0, 1) +
    ylim(0, 1) +
    geom_hline(yintercept = 0.1, linetype = "dotted", color = "red", size = 0.5) +
    theme_bw() +
    theme(
      legend.text = element_text(size = 14),  # Font size for legend text
      legend.title = element_text(size = 16),  # Font size for legend title
      axis.title.x = element_text(size = 16),  # Font size for x-axis label
      axis.title.y = element_text(size = 16),  # Font size for y-axis label
      axis.text.x = element_text(size = 14),  # Font size for x-axis tick labels
      axis.text.y = element_text(size = 14)  # Font size for y-axis tick labels
    )
  print(p) 
  myfile=paste("number of contam = ",contam[i],".pdf",sep="")
  ggsave( p, filename =myfile, height = 5, width = 7, units = "in")
}






