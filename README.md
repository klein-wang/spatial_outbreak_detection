# Spatial Outbreak Detection
Evaluate spatio-temporal effect in predicting possible disease outbreaks.


## Project Information

4 steps:
- Generate simple data and implement MCMC using a spatially independent model
- Generate spatial data and implement MCMC using both spatially independent and dependent models 
- Evaluate the spatio-temporal effect in terms of MCMC estimation
- Implement MCMC on Maryland Covid data using spatial dependent model 

## Code Structure
```
├── data_generate # generate data for both models
│   ├── generate_data_simple.R # for simple independent model
│   ├── generate_data_spatial.R # for spatially dependent model
│   ├── london # geographical information of inner London area
│   ├── simple.RData # data for simple independent model
│   └── spatial.RData # data for spatially dependent model
├── data_maryland
│   ├── data_maryland.RData
│   ├── MD_COVID-19_-_Cases_by_County.csv
│   ├── MD_Distance_by_County.csv
│   ├── MD_Population_by_County.csv
│   └── Trips_by_Distance_-_Maryland_Counties.csv
├── data_process.R # data process of Maryland raw data
├── mcmc_maryland.R
├── mcmc_simple_with_spatial.R # implement mcmc on spatial data using simple model
├── mcmc_simple.R
├── mcmc_spatial.R 
├── diagnosis.Rmd # notebook for generated data
├── maryland.Rmd # notebook for maryland data
├── trials_mcmc # store mcmc trials
│   ├── simple # mcmc trials of simple (independent) model
│   ├── simple_with_spatial
│   ├── spatial
│   └── maryland # mcmc trials of spatial model on Maryland data
├── summary_plot.R # generate summary plots using results from mcmc_trails
├── plots # store summary plots
└── README.md
```

## MCMC Implementation ##

| Model | Parameters | Iterations | Burn-In |
| ----- | ---------- | ---------- | ------- |
| simple  | &alpha;, &beta;, p | 30000 | 5000 |
| spatial | &alpha;, &beta;, &gamma;, &Delta; | 3000 | 500 |


### Maryland Covid Data

**Model parameters**

- &alpha; = log(0.52)
- &beta; = log(4.36)
- &gamma; = 0.22
- &Delta; = 0.58


## Evaluation

| Trail | &alpha; | &beta; | &gamma; | ROC (spatial) | ROC (simple) |
| ----- | ------- | ------ | ------- | ------------- | ------------ |
|       | ln(1)   | ln(2)  | 7/52    | 0.701 |  0.668 |
|       | ln(1)   | ln(3)  | 7/52    | 0.822 |  0.809 |
|       | ln(1)   | ln(4)  | 7/52    | 0.923 |  0.922 |
|       | ln(2)   | ln(2)  | 7/52    | 0.797 |  0.797 |
|       | ln(2)   | ln(3)  | 7/52    | 0.910 |  0.913 |
|       | ln(4)   | ln(2)  | 7/52    | 0.864 |  0.861 |


