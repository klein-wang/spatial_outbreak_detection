read_folder <- 'data_maryland/'

# data source: https://catalog.data.gov/dataset/md-covid-19-cases-by-county
raw <- read.csv(paste0(read_folder,'MD_COVID-19_-_Cases_by_County.csv'),header = T)
raw[is.na(raw)] <- 0 # replace NA with 0
raw <- raw[,-dim(raw)[2]] # remove last column

# population by 24 counties: https://worldpopulationreview.com/us-counties/states/md
# dynamic population: https://data.bts.gov/Research-and-Statistics/Trips-by-Distance-Maryland-Counties/k4y4-vf3n
# land area: http://www.mgs.md.gov/geology/areas_and_lengths.html
pop <- read.csv(paste0(read_folder,'MD_Population_by_County.csv'),header = T)
pop <- as.vector(pop[,2]) # vectorize by counties
pop <- round(pop/1000,2) # in thousands(k) 

# build dataset with time and locations
dense = 7 # condense data by every n days
t_include <- c(623:715) # from 2021/11/27(ID:623) to 2022/2/27(ID:715)
nrow <- dim(raw)[1]
ncol = dim(raw)[2]
data <- as.data.frame(matrix(, nrow,ncol))
colnames(data) = colnames(raw)
data[1,] = raw[1,] # keep 1st day data
data[,1:2] = raw[,1:2] #keep first two columns

for (i in 3:ncol){
  for (j in 2:nrow){ # start with 2nd row
    if (raw[j,i] > raw[j-1,i]){
      data[j,i] = raw[j,i] - raw[j-1,i]
    } else {
      data[j,i] = 0
    }
  }
}
data <- data[t_include,]
data_full <- data # used to calculate segments
nrow = as.integer(dim(data)[1]/dense) # condense into weekly data # update nrow after choosing the time
wh <- 0
for(i in 1: nrow){wh[i] = (i-1)*dense+1}
data <- data[wh,]

########### x #############
data_x = data
quantile.no = 100
for (i in 3:ncol){
  segments <- quantile(data_full[,i],probs = seq(0,1,1/quantile.no))
  for (j in 1:nrow){
    if (data_x[j,i] > segments["86%"]){ # around 1-7/52, 7/52 being the assumed outbreak prob
      data_x[j,i] = 1
    } else {data_x[j,i] = 0}
  }
}
x <- data_x[,-(1:2)]
########### rate,lambda ##########
y <- data[,-(1:2)]
rate = apply(y,2,mean) # estimated lambda for each location, using MLE
lambda = rate/pop # lambda by counties

########### distance matrix ##########
# distance matrix
d <- read.csv(paste0(read_folder,'MD_Distance_by_County.csv'),header = T)
d_inv <- 1/d[,-1]
d_inv[3,4] = 0 # Baltimore to Baltimore city, fix the Inf
d_inv[4,3] = 0 # Baltimore to Baltimore city, fix the Inf


########### Find proper startvalue #############
param_name <- c('alpha','beta','gamma','Delta')

y_per_person <- y
for (j in 1:dim(y)[2]){ # each column
  y_per_person[,j] = y_per_person[,j]/pop[j]
}
lambda_outbreak <- as.matrix(y_per_person)[which(x==1)]
lambda_no_outbreak <- as.matrix(y_per_person)[which(x==0)]

max_p = 0.8
alpha <- log(mean(lambda_no_outbreak)) # estimate on alpha
beta <- log(mean(lambda_outbreak)) - alpha # estimate on beta
gamma <- sum(x)/(dim(x)[1]*dim(x)[2]) # estimate on independent probability p
Delta <- max_p - gamma
startvalue <- c(alpha,beta,gamma,Delta)

x <- as.matrix(x)
y <- as.matrix(y)

save(data,x,y,lambda,pop,d_inv,param_name,startvalue,file = paste0(read_folder,"data_maryland.RData"))
