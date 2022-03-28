# set parameters
alpha = log(4) # = ln(3)
beta = log(10) # = ln(6)
gamma = 7/52
max_p = 0.8 # maximum p when all neighboring locations are having an outbreak
Delta = max_p - gamma
a_theta <- c(gamma,Delta,1-gamma-Delta)
n = 10  # number of weeks


# load distance information
library("readxl")
d <- read_excel('london/inner_london.xlsx', sheet=2) # distance (in miles)
d <- d[-1] #drop the 1st index column
d_inv <- 1/d # inverse distance metric
names(d_inv) <- c('Her','Cam','Isl','Hac','Ham','Ken','Wes','Cit','Tow','New','Wan','Lam','Sou','Lew')
m <- dim(d)[1] # number of boroughs


# define probability function at each time t
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



#################### Generate data for Spatially dependent model #####################
# Gibbs sampling
set.seed(2022)
x_0 <- rep(0,14) # looking at 14 boroughs in inner London area at time 0
x_0[c(8,13)] = c(1,1) # City of London-8, Southwark-13
x_i <- x_0 
p_i <- local_prob(x_i,d_inv,a_theta)

# running time: 2mins
for (k in 1:300){ # let p converge to have neighboring effect after 300 updates
  for (j in 1:m){
    x_i[j] <- rbinom(1,size = 1,prob = p_i[j]) #p in location j at time i
    p_i <- local_prob(x_i,d_inv,a_theta) # update p_i
  }
}

### Generate data
x <- x_i
p <- p_i

for (i in 1:n){
  for (j in 1:m){
    x_i[j] <- rbinom(1,size = 1,prob = p_i[j]) #p in location j at time i
    p_i <- local_prob(x_i,d_inv,a_theta) # update p_i
  }
  x <- rbind(x,x_i)
  # update p_i using x_i
  p_i <- local_prob(x_i,d_inv,a_theta)
  p <- rbind(p,p_i)
}

x <- x[-1,] # a matrix of x_it with time and location, remove the 1st empty column
p <- p[-1,] # a matrix of p_it with time and location, remove the 1st empty column



# set parameters
lambda = exp(alpha + beta*x) # excluding the noise, returns a matrix
y <- matrix(data = NA, nrow = 1, ncol = m)  
for (i in 1:n){
  y_i <- rpois(1,lambda = lambda[i,1]) # sampled observations for the 1st location at time i
  for (j in 2:m){
    y_i[j] = rpois(1,lambda = lambda[i,j]) # sampled observations for one location
  }
  y <- rbind(y,y_i)
}

y <- y[-1,] # a matrix of y_it with time and location, remove the 1st empty row

param_true <- c(alpha,beta,gamma,Delta)
param_name <- c('alpha','beta','gamma','Delta')
save(param_true,param_name,y,lambda,x,p,m,d_inv, file = "spatial.RData")
