### Risk, Return, and Equilibrium: Empirical Tests
### Replication of Fama-MacBeth (1973)

### Author: Marcos Diniz
### April-2022

## Load libraries --------------------------------------------------------------------
rm(list=ls(all=TRUE))

library(dplyr)
library(tidyverse)
library(PerformanceAnalytics)
library(quantmod)
library(TTR)
library(zoo)
library(sandwich)
library(lmtest)


## Import data and calculate returns -------------------------------------------------

## We have three important datasets in FM estimation: portfolio betas, 
## portfolio standard errors, portfolio returns.
## Our goal is to create this from individual securities data for five portfolios
## In FM paper, portfolios are defined based on individual betas for a period and 
## the portfolio information for next period (port.betas, port.returns, port.se)
## are used make the regressions.

## Data used in this replication can be downloaded at: encurtador.com.br/lAPR0

prices <- read_delim("./Data/data.csv", 
                     delim = ";", 
                     escape_double = FALSE, 
                     col_types = cols(ym = col_skip()), 
                     locale = locale(grouping_mark = " "), 
                     trim_ws = TRUE)

glimpse(prices)

## Basic descriptive functions
str(prices)
head(prices)
summary(prices)
dim(prices)

## Transform prices into returns, omit the first row
## Declare first the prices to be a time series object
prices <- ts(data=prices, frequency=12, start=c(1969, 12))
returns <- Return.calculate(prices)
returns <- na.omit(returns)

## Keep track of rows (time)
returns.t <- nrow(returns)

## Risk-free rate: read straight from FRED database and transform into monthly returns for our time period
getSymbols("TB3MS", src="FRED")
rf <- TB3MS[paste("1970-02-01", "2014-12-01", sep="/")]
rf <- (1+(rf/100))^(1/12)-1
rfts <- ts(data=rf, frequency=12, start=c(1970, 1))

## Finally calculate the market return factor
rmrf <- returns[,1]-rfts




### Calculating individual betas and standard errors for each period -----------------
## Beta time period
mo <- 60

## Update Betas and SEs yearly
update <- 12

## Keep track of time
beta.t <- returns.t-mo

## Create beta matrix
betas <- returns
betas[!is.na(betas)] <- NA

## Rolling betas
for (i in 1:ncol(betas)){
  for (j in 1:(nrow(betas)-mo)){
    betas[j+mo-1, i] <- coef(lm(returns[j:(j+mo-1),i] ~ returns[j:(j+mo-1), 1]))[2]
  }
}

## Remove NAs
betas <- na.omit(betas)

## Betas updated yearly matrix
betas.y <- betas

# Betas updated yearly
for (i in seq(from=1, to=nrow(betas.y), by=update)) {
  for (j in 1:update-1) {
    betas.y[i+j,] <- betas.y[i,]
  }
}

## Rolling SEs
se <- returns
se[!is.na(se)] <- NA

## Estimate SEs
for (i in 1:ncol(se)){
  for (j in 1:(nrow(se)-mo)){
    
    # standard deviation of residuals from estimation of betas
    se[j+mo-1, i] <- sd(summary(lm(returns[j:(j+mo-1),i] ~ returns[j:(j+mo-1), 1]))$residuals)
  }
}

## Omit NAs
se <- na.omit(se)

## SEs updated yearly
se.y <- se

## Update SEs yearly
for (i in seq(from=1, to=nrow(se.y), by=update)) {
  for (j in 1:update-1) {
    se.y[i+j,] <- se.y[i,]
  }
}




### Portfolio structure based on betas ---------------------------------------------------

## Portfolios are created based on individual betas and fixed for
## five years when they are recomputed

## Number of portfolios individual securities will be divided
n.port <- 5

## Define portfolio-betas
betas.p <- betas

## Copy the beta of portfolio mo until the reforming of portfolio
for (i in seq(from=1, to=nrow(betas), by=mo)) {
  for (j in 1:mo-1) {
    betas.p[i+j,] <- betas[i,]
  }
}

## Apply row-wise rank - higher beta, higher rank
## Classify the betas into five portfolios
ranks <- t(apply(betas.p, 1, ntile, n.port))

## Assign the assets into definied portfolios
beta.port <- list()

for (s in 1:n.port){
  
  # create five beta portfolios of logical values
  beta.port[[paste0('beta.p', s)]] <- t(apply(ranks, 1, function(x){x == s}))
  
  }


### Portfolio betas  -------------------------------------------------------------
## Estimated betas
beta.est <- list()

## Assign individual betas to portfolios
for(s in 1:n.port){
  
  beta.est[[paste0('beta.p', s)]] <- betas.y * beta.port[[s]]
  
  # Average individual betas within portfolios
  beta.est[[s]] <- as.matrix(rowSums(beta.est[[s]])/rowSums(beta.est[[s]] != 0))
  
}


### Portfolio betas squared  -------------------------------------------------------------
## Estimated betas
beta.sq.est <- list()

## Assign individual betas to portfolios
for(s in 1:n.port){
  
  beta.sq.est[[paste0('beta.sq.p', s)]] <- (betas.y^2) * beta.port[[s]]
  
  # Average individual betas within portfolios
  beta.sq.est[[s]] <- as.matrix(rowSums(beta.sq.est[[s]])/rowSums(beta.sq.est[[s]] != 0))
  
}

### Portfolio returns -------------------------------------------------------

## Calculate returns for the five portfolios
## Remove sixty first observations to get same size as ranks data frame
ret <- returns[-c(1:mo),] 
ret <- ts(data=ret, frequency=12, start=c(1975, 1))

## Returns of portfolios
ret.port <- list()

for (i in 1:n.port){
  
  # Multiply logical portfolio by ret portfolio
  ret.port[[paste0('ret.p', i)]] <- beta.port[[i]] * ret
  
  # Average individual returns within portfolios
  ret.port[[i]] <- as.matrix(rowSums(ret.port[[i]])/rowSums(ret.port[[i]] != 0))
  
}

### Portfolio standard errors -----------------------------------------------

## Calculate SEs for each portfolio
se.est <- list()
for(s in 1:n.port){
  
  # Assign individual SEs to portfolios
  se.est[[paste0('se.p', s)]] <- se.y * beta.port[[s]]
  
  # Average individual betas within portfolios
  se.est[[s]] <- as.matrix(rowSums(se.est[[s]])/rowSums(se.est[[s]] != 0))
  
}



### Dataset with estimations ------------------------------------------------

## Now we have the portfolio data (port. returns, port.betas,
## port. betas squared and port.SEs )

## Join all data
estimates <- cbind( rmrf[-c(1:mo),], do.call('cbind', c(ret.port, beta.est, beta.sq.est, se.est)))

## Replace colnames
names <- c( 'rmrf', names(c(ret.port, beta.est, beta.sq.est, se.est)) )
colnames(estimates) <- names

## As a timeseries
estimates <- ts(estimates, start = c(1974, 12), frequency = 12)

# Plot data timeseries
plot(estimates[,2:6], main = 'Portfolio Returns')
plot(estimates[,7:11], main = 'Portfolio Betas') # betas
plot(estimates[,12:16], main = 'Portfolio Squared Betas') # betas squared
plot(estimates[,17:21], main = 'Portfolio Standard Errors') # SEs



#### Run the regressions -----------------------------------------------------

# Estimates as data frame
estimates <- as.data.frame(estimates)

# Formulas to be used in linear regressions
fm <- matrix(nrow = n.port, ncol = 1)

for (i in 1:n.port){ 
  fm[i,] <- paste0('ret.p',i,' ~ beta.p', i, ' + beta.sq.p', i, ' + se.p', i)
}

## Matrix for allocation of different periods estimation
est.list <- list()

## Estimation of coefficients for each period
est.coef <- matrix(nrow = n.port, ncol = (ncol(estimates)-1)/n.port)
rownames(est.coef) <- c(paste0('p.', 1:5))
colnames(est.coef) <- c('intercept', 'coef.beta', 'coef.beta.sq', 'coef.se')

## Run a regression for each portfolio
for (j in 1:8){ # split data in 8 estimation periods
for (i in 1:n.port){ # number of portfolios
  
  # run regressions for each estimation period 
  # Roll the panel with a interval of 60 (mo) rows
  est.coef[i,] <- coef( lm(formula = fm[i,],
                           data = estimates[ ((j*mo)-mo+1):(j*mo), ]) )
  
  est.list[[paste0('est.per', j)]] <- est.coef

}
}







