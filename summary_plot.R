model_choice = c("simple","spatial")
model = model_choice[2]

if (model=="spatial"){
  folder = "trials_mcmc/spatial/spatial_"
  exp_alpha = c(1)
  exp_beta = c(4)
  burnIn = 30 # 500 for spatial
} else if (model=="simple"){
  folder = "trials_mcmc/simple/simple_"
  exp_alpha = c(1,2)
  exp_beta = c(2,3,4)
  burnIn = 5000 # 5000 for simple
}

# visualizing x
plot_x_sum <- function(chain_x,x,model="simple"){
  # chain_x: simulated x from MCMC
  # x: true x in dataset
  if (model=="simple"){
    s = apply(chain_x[-(1:burnIn),],2,mean) # accumulating column sum for each week of x
  } else if (model=="spatial"){
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
  title(paste0('alpha = ln',i,', beta = ln',j))
}

plot_roc <- function(chain_x,x,model="simple"){
  # ROC curve
  if (model=="simple"){
    s = apply(chain_x[-(1:burnIn),],2,mean) # accumulating column sum for each week of x
  } else if (model=="spatial"){
    s = apply(chain_x[-(1:burnIn),,],2:3,mean) # average x over iterations (row:time;col:location)
    s = as.vector(s)
    x <- as.vector(x)
  }
  
  df <- data.frame(s*100,x)
  names(df) <- c('chain_x','x')
  library('pROC')
  roc(df$x,df$chain_x,plot=TRUE, print.thres=TRUE, print.auc=TRUE)
  title(paste0('alpha = ln',i,', beta = ln',j))
}


par(mfrow = c(length(exp_alpha),length(exp_beta)))
for (i in exp_alpha){
  for (j in exp_beta){
    filename = paste0(folder,"ln(",i,")_ln(",j,").RData")
    load(filename)
    plot_roc(chain_x,x,model=model)
    plot_x_sum(chain_x,x,model=model)
  }
}