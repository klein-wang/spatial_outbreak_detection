read_folder = 'data_generate/'

# load data and parameters
load(paste0(read_folder,'spatial.RData')) # param_true,param_name,y,lambda,x,p,m,d_inv
alpha = param_true[1]
beta = param_true[2]
alpha.min = log(0.5) # try (log(0),log(2))
beta.min = log(0.5) # try (log(0),log(2))
a = 7
b = 52-a
x = as.vector(x)
y = as.vector(y)
p = mean(p) # > gamma, which only accounts for isolated prop

param_true = c(alpha,beta,p)
param_name = c('alpha','beta','p')
startvalue = param_true # initial value for alpha, beta and p (we pick the true value)
iterations = 30000
burnIn = 5000

# Log Likelihood
likelihood <- function(param,x){
  alpha = param[1]
  beta = param[2]
  pred = exp(alpha + beta*x)
  # each week has a lambda based on x_i
  singlelikelihoods = dpois(y, lambda = pred, log = TRUE) # an array of likelihoods for each week
  return(singlelikelihoods) 
}

# Log Prior density
prior <- function(param,x,sum=TRUE){
  alpha = param[1]
  beta = param[2]
  p = param[3]
  
  alpha_prior = dgamma(alpha-alpha.min, shape=6, log = T) # shape size determined the support range of alpha
  beta_prior = dgamma(beta-beta.min, shape=6, log = T) 
  x_prior = dbinom(x,1,p,log = T)  # a vector of density for each x
  p_prior = dbeta(p,a,b, log = T) # need to define a,b, in line8,9
  
  if (sum==T){
    res = alpha_prior+beta_prior+sum(x_prior)+p_prior # gives a value
  } else { res = alpha_prior+beta_prior+x_prior+p_prior} # gives a vector
  return(res)
}

# Log Posterior 
posterior <- function(param,x,sum=TRUE){
  if (sum==T){
    res = sum(likelihood(param,x)) + prior(param,x,T) # value + value, gives a value
  } else { res = likelihood(param,x) + prior(param,x,F)} # vector + vector, gives a vector
  return(res)
}



### Metropolis Hasting algorithm ###

proposalfunction <- function(param){
  param[1] = rnorm(1,mean = param[1], sd= .2) # propose alpha
  param[2] = rnorm(1,mean = param[2], sd= .1) # propose beta
  return(param)
}

run_metropolis_MCMC <- function(startvalue,x,iterations,burnIn){
  # create chain to store alpha,beta,p
  param_no = length(startvalue) # number of parameters
  chain = array(dim = c(iterations+1,param_no))
  chain[1,] = startvalue # 1st row 
  
  # create chain_x to store x
  length_x <- length(x) # size of x in an iteration, 52 here
  chain_x <- array(dim = c(iterations+1,length_x))
  chain_x[1,] <- x # startvalue for x
  
  p_index = 3
  
  # iterations
  for (i in 1:iterations){
    
    # MH for alpha, beta
    proposal = proposalfunction(chain[i,]) # proposal for alpha, beta
    if (proposal[2]>beta.min & proposal[1]>alpha.min){ # reject beta if it's smaller than beta.min
      probab = exp(posterior(proposal,x,T) - posterior(chain[i,],x,T)) # ratio
      if (runif(1) < probab){ 
        chain[i+1,1] = proposal[1] # accept alpha
        chain[i+1,2] = proposal[2] # accept beta
      }else{
        chain[i+1,] = chain[i,] # reject
      }
    }else{ 
      chain[i+1,] = chain[i,] 
    } # reject 
    
    # update p
    n_case <- sum(chain_x[i,]) # sum of x in ith iteration
    p <- rbeta(1,a+n_case,b-n_case+length_x)
    chain[i+1,p_index] = p # update p here using Gibbs Simpler
    param <- chain[i+1,] # set as current proposal
    
    # update x
    log_p1 <- posterior(param,1,F)
    log_p0 <- posterior(param,0,F) # R = exp(log_p1)/(exp(log_p1)+exp(log_p0)), i.e. p1/(p1+p0)
    M <- pmax(log_p0,log_p1)
    R <- exp(log_p1-M)/(exp(log_p1-M)+exp(log_p0-M)) # re-write R, using log-sum-exp trick
    
    x <- rbinom(length_x,1,R) # vector of proposed x, R = exp(posterior(param,1))*normalized_c
    chain_x[i+1,] = x # update proposal_x     
  }
  
  # return acceptance for alpha,beta
  acceptance = 1-mean(duplicated(chain[-(1:burnIn),c(1,2)]))
  
  # return(chain)
  setClass(Class="res",
           representation(
             chain="matrix",
             chain_x="matrix",
             acceptance="numeric"
           )
  )
  return(new('res',
             chain = chain,
             chain_x = chain_x,
             acceptance = acceptance))
}



# trace plots & histgrams
# plot the histogram of the parameter distributions
plot_hist <- function(chain,burnIn,param_index,param_name,param_true){
  hist(chain[-(1:burnIn),param_index],nclass=30, main=paste("Posterior of",param_name[param_index]), xlab="True value = red line" )
  abline(v = mean(chain[-(1:burnIn),param_index]),col='blue')
  abline(v = param_true[param_index], col="red" )
}
# plot the MCMC chains for the parameter
plot_chain <- function(chain,burnIn,param_index,param_name,param_true){
  plot(chain[-(1:burnIn),param_index], type = "l", xlab="True value = red line" , main = paste("Chain values of",param_name[param_index]))
  abline(h = param_true[param_index], col="red" )
  abline(h = mean(chain[-(1:burnIn),param_index]),col='blue')
}
# visualizing x
plot_x <- function(chain_x,x){
  # chain_x: simulated x from MCMC
  # x: true x in dataset
  
  s = apply(chain_x[-(1:burnIn),],2,sum) # accumulating column sum for each week of x
  quantile.no = 400
  segments <- quantile(s/(iterations-burnIn),probs = seq(0,1,1/quantile.no))
  point <- unname(segments[c(1,quantile.no-1)])
  plot(s/(iterations-burnIn), xlab="Week Number" , ylab="Outbreak probability",main = "Probability of having an outbreak")
  abline(h = point[1], col="red")
  abline(h = point[2], col='blue')
  # compare with real x
  wh <- which(x==1)
  for (j in wh){ abline(v = j, col='green') }
  
  # ROC curve
  df <- data.frame(s/(iterations-burnIn)*100,x)
  names(df) <- c('chain_x','x')
  library('pROC')
  roc(df$x,df$chain_x,plot=TRUE, print.thres=TRUE, print.auc=TRUE)
}




############### beta ################
find_beta <- function(beta_){
  startvalue[2] = beta_
  start_time <- Sys.time()
  res <- run_metropolis_MCMC(startvalue,x,iterations,burnIn)
  end_time <- Sys.time()
  run_time <- end_time-start_time
  print(paste('Running time of MCMC:',round(run_time,2),'seconds.'))
  
  chain <- res@chain
  chain_x <- res@chain_x
  acceptance <- res@acceptance
  print(paste('Acceptance rate of MCMC:',round(100*acceptance,4),'%.'))
  
  param_estimate = rep(0,length(param_name))
  for (i in 1:length(param_name)){
    param_estimate[i] <- mean(chain[-(1:burnIn),i])
  }
  print(paste('We have the starting values of parameters as:',list(round(startvalue,2))))
  print(paste('We have the estimated values of parameters as:',list(round(param_estimate,2))))
  
  pdf(paste0('trials_mcmc/simple_with_spatial/simple_with_spatial_ln(',exp(alpha),')_ln(',exp(beta_),").pdf"))
  par(mfrow = c(3,3))
  for (i in 1:length(param_name)){
    plot_hist(chain,burnIn = burnIn, param_index = i, param_name = param_name, param_true = param_true)
  }
  for (i in 1:length(param_name)){
    plot_chain(chain,burnIn = burnIn, param_index = i, param_name = param_name, param_true = param_true)
  }
  
  # visualizing x
  plot_x(chain_x,x) # return cumulative sum s
  
  plot(0:1,0:1,t="n") # make the plot transparent using axes
  legend("topright",c('0.25% quantile','99.75% quantile','real outbreak week'),col=c('red','blue','green'),lty=1)
  
  dev.off()
  # return(chain_x)
  save(chain,chain_x,acceptance,file = paste0("trials_mcmc/simple_with_spatial/simple_with_spatial_ln(",exp(alpha),
                                      ")_ln(",exp(beta_),").RData"))
}

# store R data

######################################
# run different betas
betas <- c(log(1.5),log(6),log(10),log(20),log(40),log(80),log(200))
# for (i in betas){ find_beta(i) }
find_beta(beta)
# apply(chain_x[-(1:burnIn),],2,sum) # cumulative sum of x 
# plot_x(chain_x,x)
paste("Results have been saved at 'trials_mcmc/simple_with_spatial_beta_log...'")
