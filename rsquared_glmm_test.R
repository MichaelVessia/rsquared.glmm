#Generate a mock dataset
set.seed(4)
data <- data.frame(y.nbinom=rnegbin(100,5,2))
data <- cbind(data, 
              data.frame(fixed1=data$y+c(runif(50, 0, 5),runif(50, 10, 50)),
                         fixed2=c("Treatment1", "Treatment2"),
                         rand1=LETTERS[1:2],
                         rand2=rep(LETTERS[23:26],each=25)) )



#Generate models using mock dataset
library(lme4)
mod5 <- glmer(y.nbinom ~ fixed1*fixed2 + (1|rand2/rand1), family="negative.binomial", data)
mod5.sqrt<- update(mod5, family=negative.binomial(link = "sqrt"))
#Get values for all kinds of models
(lme4.models <- rsquared.glmm(list(mod5, mod5.sqrt)))
