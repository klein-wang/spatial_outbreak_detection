# Simulate parameters p,x,y
a = 7
b = 52-a # "a" weeks in a year would has an outbreak (Assumption: a/(a+b))
# p <- mean(rbeta(100, shape1 = a, shape2 = b)) # probability of having an disease (once every year) 
p = a/(a+b)
alpha = log(2)
beta = log(4)

# generate data 
set.seed(2022)
n = 52  # number of weeks
x <- rbinom(n,size = 1,prob = p)  # an array of 0 or 1 in 52 weeks
lambda = exp(alpha + beta*x) # excluding the noise
y = rpois(n,lambda = lambda) # construct y (observations)
model <- glm(y~x,family=poisson(link = 'log'))
# summary(model)


param_true = c(alpha,beta,p) # true value for alpha, beta and p
param_name = c('alpha','beta','p')

save(param_true,param_name,y,lambda,x,a,b, file = "simple.RData")
save(param_true,param_name,y,lambda,x,a,b, file = paste0("simple_ln(",exp(alpha),")_ln(",exp(beta),").RData"))
