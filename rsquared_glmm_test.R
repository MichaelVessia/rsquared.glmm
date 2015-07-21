library(MASS)
library(ggplot2)
#Generate a mock dataset
set.seed(4)
data <- data.frame(y.nbinom=rnegbin(100,5,2))
data <- cbind(data, 
              data.frame(fixed1=c(runif(50, 0, 5)),
                         fixed2=c(runif(50, 10, 50)),
                         rand1=c(rnorm(50, 5, 10)),
                         rand2=c(rpois(50, 10))) )



#Generate models using mock dataset
library(lme4)
mod5 <- glmer(y.nbinom ~ fixed1*fixed2 + (1|rand2/rand1), family="negative.binomial", data)
mod5.sqrt<- update(mod5, family=negative.binomial(link = "sqrt"))
#Get values for all kinds of models
(lme4.models <- rsquared.glmm(list(mod5, mod5.sqrt)))


#Uniform distribution
plot(data$fixed1, xlim=c(0, 100), ylim=c(0, 100))

#Negative Binomial distribution
plot(data$y.nbinom, xlim=c(0, 100), ylim=c(0, 100))
hist(data$y.nbinom)

plot(data$rand1, xlim=c(0,100), ylim=c(0,100))

x   <- seq(-15,15,length=1000)
y <- dnorm(x, mean=0, sd=3)
#plot(x,y, type="l", lwd=1)

ggplot(data=data, aes(x,y)) + geom_point()
