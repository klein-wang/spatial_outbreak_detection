library(MCMCpack)
library(MASS)
read_folder = 'data_maryland/'
trial = 2

# load data and parameters
load(paste0(read_folder,"data_maryland.RData"))  # data,x,y,pop,d_inv,param_name,startvalue

alpha = startvalue[1]
beta = startvalue[2]
alpha.min = log(0.5) # try (log(0),log(2))
beta.min = log(0.5) # try (log(0),log(2))
gamma = startvalue[3]
Delta = startvalue[4]
a_theta = c(gamma,Delta,1-gamma-Delta)

iterations = 10000 #2000, running time 1hr55mins
burnIn = 500 #200


# local_prob function 
local_prob <- function(x_t,d_inv,theta){
  # x_t: a vector of Xs indicating an outbreak or not for each location at time t
  # d_inv: an inverse distance metrix among all locations within the region
  # j: the index of the location
  # theta: a vector of gamma, Delta and 1-gamma-Delta
  gamma = theta[1]
  Delta = theta[2]
  m <- dim(d_inv)[1] # number of locations 
  prob_t <- rep(0,m) # probability for each location at time t
  
  if (length(x_t) != 1) {
    for (j in 1:m){
      x_t = x_t[-j] # removing x_j
      d_j = d_inv[j,-j] # the jth row excluding the jth location, giving a vector where j!=k
      weighted_d_j = d_j/sum(d_j) # compute the weight to each entry
      
      prob_t[j] <- gamma+Delta*sum(weighted_d_j*x_t)
    }
  }
  else {prob_t = gamma+Delta*x_t} # if x_t is single value
  return(prob_t)
}



############### Construct Posterior Distribution ################
# log density function of the Dirichlet distribution
lddirichlet <- function(x,a) {
  # x: current value of theta, a: a_theta
  wh<-which(a!=0) # allow a for having 0s for some a_i
  return(lgamma(sum(a[wh]))+sum((a[wh]-1)*log(x[wh])-lgamma(a[wh])))
}


# Log Prior density
prior <- function(param,x,sum=TRUE){
  # x: a matrix of x_it, outbreak indicator in location i at time t
  alpha = param[1]
  beta = param[2]
  gamma = param[3]
  Delta = param[4]
  theta = c(gamma,Delta,1-gamma-Delta)
  
  alpha_prior = dgamma(alpha-alpha.min, shape = 6, log = T) # shape size determined the support range of alpha
  beta_prior = dgamma(beta-beta.min, shape=6, log = T) # shape size determined the support range of beta
  
  x <- as.matrix(x)
  x_prior <- matrix(data = NA, nrow = 1, ncol = dim(x)[2])  
  for (t in 1:dim(x)[1]){
    p_t = local_prob(x[t,],d_inv,theta) # compute localized probability at time t
    x_t_prior <- dbinom(x[t,],1,p_t,log = T)  # a vector of density for each x_t
    x_prior <- rbind(x_prior,x_t_prior)
  }
  x_prior <- x_prior[-1,] # remove the 1st empty row
  
  gamma_Delta_prior = lddirichlet(theta,a_theta)
  
  if (sum==T){res = alpha_prior+beta_prior+sum(x_prior)+gamma_Delta_prior}
  else {res = alpha_prior+beta_prior+x_prior+gamma_Delta_prior}
  return(res)
}


# Log Likelihood
# likelihood of y given x 
likelihood <- function(param,x){
  alpha = param[1]
  beta = param[2]
  pred_per_person = exp(alpha + beta*x) # lambda per person
  t = dim(y)[1]
  m = dim(y)[2]
  
  if(length(x)==1){
    pred = matrix(0,t,m)
    for (j in 1:m){
      pred[,j] = pred_per_person*pop[j]
    }
  } else {
    pred = matrix(0,t,m)
    for (j in 1:m){
      pred[,j] = pred_per_person[,j]*pop[j]
    }
  }
  
  singlelikelihoods = dpois(y, lambda = pred, log = TRUE) # predicted lambda, y comes from above
  return(singlelikelihoods) # a matrix of likelihoods 
}

# likelihood of x given p(theta)
likelihood_x <- function(x,theta){ 
  likelihood_x <- matrix(data = NA, nrow = 1, ncol = dim(x)[2]) 
  for (i in 1:dim(x)[1]){
    res = dbinom(x[i,], 1, local_prob(x[i,],d_inv,theta), log = T)
    likelihood_x <- rbind(likelihood_x,res)
  }
  likelihood_x <- likelihood_x[-1,] # remove the 1st empty row
  return(sum(likelihood_x)) # likelihood for x given theta
}


# Posterior 
posterior <- function(param,x,sum=TRUE){
  if (sum==T){res = sum(likelihood(param,x)) + prior(param,x,T)}
  else {res = likelihood(param,x) + prior(param,x,F)}
  return (res) # matrix + matrix
}


### Metropolis Hasting algorithm ###

load("trials_mcmc/maryland/maryland_ln(0.52)_ln(4.36).RData")
sigma <- cov(chain[-1,c(1,2)])
Sigma <- (2.38^2/dim(sigma)[1])*sigma
proposalfunction <- function(param){
  # param: alpha,beta,gamma,Delta
  tmp = mvrnorm(n=1,mu=param[c(1,2)],Sigma=Sigma) # proposal distribution of alpha,beta
  param[1] = tmp[1]
  param[2] = tmp[2]
  # print(paste(param[c(1,2)]))
  return(param)
} 



run_metropolis_MCMC <- function(startvalue,x,iterations,burnIn){
  # create chain to store alpha,beta,gamma,Delta
  chain = matrix(NA,iterations+1,length(startvalue)) # each row: c(alpha,beta,gamma,Delta)
  chain[1,] = startvalue # 1st row 
  
  # create chain_x to store x
  t <- dim(x)[1] # total time
  m <- dim(x)[2] # total locations
  chain_x <- array(NA, c(iterations+1,t,m)) # create a 3 dimensional array
  chain_x[1,,] <- x # startvalue for x at time 0
  a <- 1 # scale parameter in Dirichlet distribution Dir(theta+a*theta_)
  prob_alpha_beta = 0
  prob_theta = 0
  
  # iterations
  for (i in 1:iterations){
    if(i %% 100==0){print(paste('iteration',i))}
    
    # MH for alpha, beta
    proposal = proposalfunction(chain[i,]) # propose new alpha, beta
    if (proposal[2]>beta.min & proposal[1]>alpha.min){ # reject if it's smaller than min
      probab = exp(posterior(proposal,x,T) - posterior(chain[i,],x,T)) # ratio
      prob_alpha_beta[i] = probab  # record probability ratio
      if (is.na(probab) || runif(1) >= probab){
        chain[i+1,] = chain[i,] # reject
      } else {
        chain[i+1,] = proposal # accept alpha,beta
      }
    }
    else{
      prob_alpha_beta[i] = 0
      chain[i+1,] = chain[i,] # reject
    }
    
    # Dirichlet random walk
    theta <- c(chain[i,3],chain[i,4],1-chain[i,3]-chain[i,4]) # current theta
    theta_ = rdirichlet(1,a*theta+a_theta) # proposal for new theta
    prob1 <- likelihood_x(x,theta_) - likelihood_x(x, theta) # ratio of likelihoods
    prob2 <- lddirichlet(theta_, a_theta)-lddirichlet(theta, a_theta) # ratio of proposal density
    prob3 <- lddirichlet(theta,a_theta+a*theta_) - lddirichlet(theta_,a_theta+a*theta) # ratio of proposal density
    full_prob <- exp(prob1+prob2+prob3)
    prob_theta[i] = full_prob  # record probability ratio
    if (is.na(full_prob) || runif(1) >= full_prob){
      a = a+1
    } else {
      chain[i+1,3] = theta_[1] # accept gamma
      chain[i+1,4] = theta_[2] # accept Delta
      a = max(0,a-3) # adaptive approach on choosing scaling parameter a, to get close to 25% ratio
    }
    
    # update x, Gibbs sampler
    param <- chain[i+1,]
    # normalized_c <- 1/(exp(posterior(param,1,F))+exp(posterior(param,0,F)))
    log_p1 <- posterior(param,1,F)
    log_p0 <- posterior(param,0,F) # R = exp(log_p1)/(exp(log_p1)+exp(log_p0)), i.e. p1/(p1+p0)
    M <- pmax(log_p0,log_p1)
    R <- exp(log_p1-M)/(exp(log_p1-M)+exp(log_p0-M)) # re-write R, using log-sum-exp trick
    
    for (j in 1:m){
      chain_x[i+1,,j] = rbinom(t,1,R[,j]) # vector of proposed x for one location
    }
  }
  
  # return acceptances
  acceptance_alpha_beta = 1-mean(duplicated(chain[-(1:burnIn),c(1:2)]))
  acceptance_theta = 1-mean(duplicated(chain[-(1:burnIn),c(3:4)]))
  
  # return(chain)
  setClass(Class="res",
           representation(
             chain="matrix",
             chain_x="array",
             prob_alpha_beta="numeric",
             prob_theta="numeric",
             acceptance_alpha_beta="numeric",
             acceptance_theta="numeric"
           )
  )
  return(new('res',
             chain = chain,
             chain_x = chain_x,
             prob_alpha_beta = prob_alpha_beta,
             prob_theta = prob_theta,
             acceptance_alpha_beta = acceptance_alpha_beta,
             acceptance_theta = acceptance_theta))
}

# trace plots & histgrams
# plot the histogram of the parameter distributions
plot_hist <- function(chain,burnIn,param_index,param_name,param_true){
  hist(chain[-(1:burnIn),param_index],nclass=30, main=paste("Posterior of",param_name[param_index]), xlab="True value = red line" )
  abline(v = mean(chain[-(1:burnIn),param_index]),col='blue')
}
# plot the MCMC chains for the parameter
plot_chain <- function(chain,burnIn,param_index,param_name,param_true){
  plot(chain[-(1:burnIn),param_index], type = "l", xlab="True value = red line" , main = paste("Chain values of",param_name[param_index]))
  abline(h = mean(chain[-(1:burnIn),param_index]),col='blue')
}
# visualizing x
plot_x <- function(chain_x,x){
  # chain_x: simulated x from MCMC
  # x: true x in dataset
  
  s = apply(chain_x[-(1:burnIn),,],2:3,mean) # average x over iterations (row:time;col:location)
  s = as.vector(s)
  quantile.no = 400
  segments <- quantile(s,probs = seq(0,1,1/quantile.no))
  point <- unname(segments[c(1,quantile.no-1)])
  plot(s, xlab="Week Number,Location" , ylab="Outbreak probability", main = "Probability of having an outbreak")
  abline(h = point[1], col="red")
  abline(h = point[2], col='blue')
  # compare with real x
  wh <- which(x==1)
  for (j in wh){ abline(v = j, col='green') }
  
  # ROC curve
  df <- data.frame(s*100,x)
  names(df) <- c('chain_x','x')
  library('pROC')
  roc(df$x,df$chain_x,plot=TRUE, print.thres=TRUE, print.auc=TRUE)
}



find_beta <- function(beta_){
  startvalue[2] = beta_ # update the startvalue of beta
  start_time <- Sys.time()
  res <- run_metropolis_MCMC(startvalue,x,iterations,burnIn)
  end_time <- Sys.time()
  run_time <- end_time-start_time
  print(paste('Running time of MCMC:',round(run_time,2),'hrs.'))
  
  chain <- res@chain
  chain_x <- res@chain_x
  prob_alpha_beta <- res@prob_alpha_beta
  prob_theta <- res@prob_theta
  acceptance_alpha_beta <- res@acceptance_alpha_beta
  acceptance_theta <- res@acceptance_theta
  print(paste('Acceptance rate of alpha, beta:',round(100*acceptance_alpha_beta,4),'%.'))
  print(paste('Acceptance rate of theta:',round(100*acceptance_theta,4),'%.'))
  save(chain,chain_x,prob_alpha_beta,prob_theta,acceptance_alpha_beta,
       acceptance_theta,file = paste0("trials_mcmc/maryland/maryland_ln(",round(exp(alpha),2),
                                      ")_ln(",round(exp(beta_),2),")_trial",trial,"_.RData"))
  
  param_estimate = rep(0,length(param_name))
  for (i in 1:length(param_name)){
    param_estimate[i] <- mean(chain[-(1:burnIn),i])
  }
  print(paste('We have the starting values of parameters as:',list(round(startvalue,2))))
  print(paste('We have the estimated values of parameters as:',list(round(param_estimate,2))))
  
  pdf(paste0("trials_mcmc/maryland/maryland_ln(",round(exp(alpha),2),
             ")_ln(",round(exp(beta_),2),")_trial",trial,"_.pdf"))
  par(mfrow = c(3,4))
  for (i in 1:length(param_name)){
    plot_hist(chain,burnIn = burnIn, param_index = i, param_name = param_name, param_true = param_true)
  }
  for (i in 1:length(param_name)){
    plot_chain(chain,burnIn = burnIn, param_index = i, param_name = param_name, param_true = param_true)
  }
  
  # visualizing x
  plot_x(chain_x,x) # return average x per week, location
  
  plot(0:1,0:1,t="n") # make the plot transparent using axes
  legend("topright",c('0.25% quantile','99.75% quantile','real outbreak week'),col=c('red','blue','green'),lty=1)
  
  dev.off()
}


######################################
# run different betas
# betas <- c(log(1.5),log(6),log(10),log(20))
# for (i in betas){ find_beta(i) }
find_beta(beta)
# plot_x(chain_x,x)
