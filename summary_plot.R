model_choice = c("simple","spatial","simple_with_spatial")
model = model_choice[3] # choosing the model

folder = paste0("trials_mcmc/",model,"/",model,"_")
exp_alpha = c(1,2)
exp_beta = c(2,3,4)
plot_trace = F

if (model=="spatial"){
  burnIn = 500 # 500 for spatial
} else if (model=="simple"){
  burnIn = 5000 # 5000 for simple
} else if (model=="simple_with_spatial"){
  burnIn = 5000 # 5000 for simple
}

# visualizing x
plot_x_sum <- function(chain_x,model="simple"){
  # chain_x: simulated x from MCMC, 1st row being the true x in dataset
  
  if (model=="simple" | model=="simple_with_spatial"){
    x = chain_x[1,]
    chain_x = chain_x[-1,]
    s = apply(chain_x[-(1:burnIn),],2,mean) # accumulating column sum for each week of x
  } else if (model=="spatial"){
    x = chain_x[1,,]
    chain_x = chain_x[-1,,]
    s = apply(chain_x[-(1:burnIn),,],2:3,mean) # average x over iterations (row:time;col:location)
    s = as.vector(s)
  } 
  
  quantile.no = 400
  segments <- quantile(s,probs = seq(0,1,1/quantile.no))
  point <- unname(segments[c(1,quantile.no-1)])
  plot(s, xlab="Week Number,Location" , ylab="Outbreak probability")
  abline(h = point[1], col="red")
  abline(h = point[2], col='blue')
  # compare with real x
  wh <- which(x==1)
  for (k in wh){ abline(v = k, col='green') }
  
  #plot(0:1,0:1,t="n") # make the plot transparent using axes
  #legend("topright",c('0.25% quantile','99.75% quantile','real outbreak week'),col=c('red','blue','green'),lty=1)
}

plot_roc <- function(chain_x,model="simple"){
  # ROC curve
  if (model=="simple" | model=="simple_with_spatial"){
    x = chain_x[1,]
    chain_x = chain_x[-1,]
    s = apply(chain_x[-(1:burnIn),],2,mean) # accumulating column sum for each week of x
  } else if (model=="spatial"){
    x <- as.vector(chain_x[1,,])
    chain_x = chain_x[-1,,]
    s = apply(chain_x[-(1:burnIn),,],2:3,mean) # average x over iterations (row:time;col:location)
    s = as.vector(s)
  }
  
  df <- data.frame(s*100,x)
  names(df) <- c('chain_x','x')
  library('pROC')
  roc(df$x,df$chain_x,plot=TRUE, print.thres=TRUE, print.auc=TRUE)
}


par(mfrow = c(length(exp_alpha),length(exp_beta)))
for (i in exp_alpha){
  for (j in exp_beta){
    filename = paste0(folder,"ln(",i,")_ln(",j,").RData")
    load(filename)
    #plot_roc(chain_x,model=model)
    plot_x_sum(chain_x,model=model)
    title(paste0('alpha = ln',i,', beta = ln',j))
  }
}


#1plot for spatial_ln(4)_ln(2).RData
load(paste0(folder,"ln(4)_ln(2).RData"))
#plot_roc(chain_x,model=model)
plot_x_sum(chain_x,model=model)
title(paste0('alpha = ln',4,', beta = ln',2))



########### Plot Trace Plot #############
plot_chain <- function(chain,burnIn,param_index,param_name,param_true){
  plot(chain[-(1:burnIn),param_index], type = "l", xlab="True(red) & Estimated(blue)" , ylab="Value", main = paste("Chain values of",param_name[param_index]))
  abline(h = param_true[param_index], col="red" )
  abline(h = mean(chain[-(1:burnIn),param_index]),col='blue')
}

if (plot_trace==T){
  load('data_generate/spatial.RData')
  load(paste0(folder,"ln(1)_ln(4).RData"))
  par(mfrow = c(2,2))
  for (i in 1:length(param_name)){
    plot_chain(chain,burnIn = burnIn, param_index = i, param_name = param_name, param_true = param_true)
  }
}
