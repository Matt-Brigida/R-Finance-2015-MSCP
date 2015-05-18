# server.R

library(EIAdata)
library(quantmod)

shinyServer(function(input, output) {

  key <- source("./markov_switching_coin_eq/key")$value  

  dataInput <- reactive({

    ticker1 <- switch(input$symb1,
		      "1" = "PET.RWTC.M",
		      "2" = "PET.RBRTE.M",
		      "3" = "NG.RNGWHHD.M",
		      "4" = "PET.EER_EPD2F_PF4_Y35NY_DPG.M",
		      "5" = "PET.EER_EPMRU_PF4_RGC_DPG.M",
		      "6" = "PET.EER_EPMRR_PF4_Y05LA_DPG.M",
		      "7" = "PET.EER_EPJK_PF4_RGC_DPG.M",
		      "8" = "PET.EER_EPD2DXL0_PF4_RGC_DPG.M")

    ticker2 <- switch(input$symb2,
		      "1" = "PET.RWTC.M",
		      "2" = "PET.RBRTE.M",
		      "3" = "NG.RNGWHHD.M",
		      "4" = "PET.EER_EPD2F_PF4_Y35NY_DPG.M",
		      "5" = "PET.EER_EPMRU_PF4_RGC_DPG.M",
		      "6" = "PET.EER_EPMRR_PF4_Y05LA_DPG.M",
		      "7" = "PET.EER_EPJK_PF4_RGC_DPG.M",
		      "8" = "PET.EER_EPD2DXL0_PF4_RGC_DPG.M")

    lhs <- getEIA(ticker1, key = key)
  
    rhs  <- getEIA(ticker2, key = key) 

    lhs <- lhs[paste(input$dates[1], "/", input$dates[2], sep = "")]
    rhs <- rhs[paste(input$dates[1], "/", input$dates[2], sep = "")]

    data <- merge.xts(lhs, rhs, join = "inner")
  })

  output$plot <- renderPlot({   
    data <- dataInput()
p.0 <- c(input$init, (1 - input$init))

## Create transition matrix

P <- matrix(c(p.0[1], 1-p.0[1], 1- p.0[2], p.0[2]), nrow=2, ncol=2)

## Define log oil and log ng ##

lnoil <- log(as.vector(data[,2]))
lnng <- log(as.vector(data[,1]))

##################################
## Log-Likelihood Function     ###
##################################

lik <- function(theta, lnoil, lnng){

alpha1 <- theta[1]
alpha2 <- theta[2]
alpha3 <- theta[3]
alpha4 <- theta[4]
alpha5 <- theta[5]
alpha6 <- theta[6]
p11 <- 1 / (1 + exp(-theta[7]))
p22 <- 1 / (1 + exp(-theta[8]))

dist.1 <- 0
dist.1 <- (1/(alpha5*sqrt(2*pi)))*exp((-(lnng-alpha1-alpha3*lnoil)^2)/(2*alpha5^2))


dist.2 <- 0
dist.2 <- (1/(alpha6*sqrt(2*pi)))*exp((-(lnng-alpha2-alpha4*lnoil)^2)/(2*alpha6^2))


dist <- cbind(dist.1, dist.2)

o.v <- c(1,1)

P <- matrix(c(p11, 1-p11, 1- p22, p22), nrow=2, ncol=2)

xi.a <- rep(0,2*length(lnoil))
xi.a <- matrix(xi.a, nrow=(length(lnoil)),ncol=2)
xi.b <- rep(0,2*length(lnoil))
xi.b <- matrix(xi.a, nrow=(length(lnoil)),ncol=2)
model.lik <- rep(0, length(lnoil))

xi.a[1,] <- (c(p11,p22)*dist[1,])/(o.v%*%(c(p11,p22)*dist[1,]))

for (i in 1:(length(lnoil)-1)){
  xi.b[i+1,] <- P%*%xi.a[i,]
  xi.a[i+1,] <- (xi.b[i+1,]*dist[i+1,])/(o.v%*%(xi.b[i+1,]*dist[i+1,]))
  model.lik[i+1] <- o.v%*%(xi.b[i+1,]*dist[i+1,])
}

logl <- sum(log(model.lik[2:length(model.lik)]))
return(-logl)
} 

##################################
## Max Log-Likelihood Function ###
##################################

###################
## Unconstrained ##
###################

theta.start <- c(-.05, .01, .2, .4, .1, .2, input$init, (1 - input$init))
max.lik.optim <- optim(theta.start, lik, lnoil=lnoil, lnng=lnng, hessian=F)

#OI <- solve(max.lik.optim$hessian)
#se <- sqrt(diag(OI))
#t <- max.lik.optim$par/se
#pval <- 2*(1-pt(abs(t), nrow(data)-2))
#results <- cbind(max.lik.optim$par, se, t, pval)
#colnames(results) <- c("Coef", "se", "t-stat", "pval")
#rownames(results) <- c("beta_01", "beta_02", "beta_11", "beta_12", "sigma2_1", "sigma2_2","p11", "p22")
# print(results, digits=4)


######################################################################
## Run filter to get state probabilities from estiamted parameters  ##
######################################################################

## Vector of Model Parameters, beta 0, beta 1, sigma ###

alpha.hat <- c(max.lik.optim$par[1], max.lik.optim$par[2], max.lik.optim$par[3], max.lik.optim$par[4], max.lik.optim$par[5], max.lik.optim$par[6])

p.0.hat <- c(( 1 / (1 + exp(-max.lik.optim$par[7]))), ( 1 / (1 + exp(-max.lik.optim$par[8]))))

## Create transition matrix

P.hat <- matrix(c(p.0.hat[1], 1-p.0.hat[1], 1- p.0.hat[2], p.0.hat[2]), nrow=2, ncol=2)


### Create the distributions for both states ###

## State 1

dist.1.hat <- 0
dist.1.hat <- (1/(alpha.hat[5]*sqrt(2*pi)))*exp((-(lnng-alpha.hat[1]-alpha.hat[3]*lnoil)^2)/(2*alpha.hat[5]^2))


## State 2

dist.2.hat <- 0
dist.2.hat <- (1/(alpha.hat[6]*sqrt(2*pi)))*exp((-(lnng-alpha.hat[2]-alpha.hat[4]*lnoil)^2)/(2*alpha.hat[6]^2))


## Create 2x1 distribution vectors at time t
## this actually creates a Tx2 matrix
## so refer to the t vector by dist[t,]

dist.hat <- cbind(dist.1.hat, dist.2.hat)

## Create 2x1 one vector

o.v <- c(1,1)

############
## Filter ##
############

xi.a.hat <- rep(0,2*length(lnoil))
xi.a.hat <- matrix(xi.a.hat, nrow=(length(lnoil)),ncol=2)
xi.b.hat <- rep(0,2*length(lnoil))
xi.b.hat <- matrix(xi.b.hat, nrow=(length(lnoil)),ncol=2)
model.lik.hat <- rep(0, length(lnoil))

xi.a.hat[1,] <- (p.0.hat*dist.hat[1,])/(o.v%*%(p.0.hat*dist.hat[1,]))

for (i in 1:(length(lnoil)-1)){
  xi.b.hat[i+1,] <- P.hat%*%xi.a.hat[i,]
  xi.a.hat[i+1,] <- (xi.b.hat[i+1,]*dist.hat[i+1,])/(o.v%*%(xi.b.hat[i+1,]*dist.hat[i+1,]))
	model.lik.hat[i+1] <- o.v%*%(xi.b.hat[i+1,]*dist.hat[i+1,])
}

# plot(xi.a.hat[,2], type='l', xlab='Month', ylab='Probability', main='Filtered Probability of Being in State 2')

## Log-Likelihood Value
# max.log.lik.value <- sum(log(model.lik.hat[2:length(model.lik.hat)]))
# print(max.log.lik.value)
#
Filtered_Probability_of_State_1 <- as.xts(xi.a.hat[,2], order.by = index(data))    
       chartSeries(Filtered_Probability_of_State_1, theme = chartTheme("white"), 
                    type = "line", TA = NULL)
  })
})
