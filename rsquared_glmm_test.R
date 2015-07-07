# #Creates a negative binomial dataset using the sample Quine data
# quine.nb <- glmer(Days ~ 1|Sex, family=negative.binomial(theta = 2, link = "log"),  data = quine)
# 
# #Calculate the rsquared value
# rsquared.glmm(quine.nb)


set.seed(4)
data <- data.frame(y=rnorm(100, 5, 10), y.binom=rbinom(100, 1, 0.5), y.poisson=rpois(100, 5), y.nbinom=rnegbin(100,5,2))
data <- cbind(data, 
              data.frame(fixed1=data$y+c(runif(50, 0, 5),runif(50, 10, 50)),
                         fixed2=c("Treatment1", "Treatment2"),
                         rand1=LETTERS[1:2],
                         rand2=rep(LETTERS[23:26],each=25)) )



library(lme4)
#Linear model
mod0 <- lm(y ~ fixed1, data)
#Linear mixed effects model
mod1 <- lmer(y ~ fixed1 + (1|rand2/rand1), data)
rsquared.glmm(mod1)
mod1.1 <- lmer(y ~ fixed1 + (fixed1|rand2/rand1), data)
mod2 <- lmer(y ~ fixed1 + fixed2 + (1|rand2/rand1), data)
rsquared.glmm(list(mod0, mod1, mod1.1, mod2))
#Generalized linear mixed effects model (binomial)
mod3 <- glmer(y.binom ~ fixed1*fixed2 + (1|rand2/rand1), family="binomial", data)
mod3.prob <- update(mod3, family = binomial(link = "probit"))
rsquared.glmm(list(mod3, mod3.prob))
#Generalized linear mixed effects model (poisson)
mod4 <- glmer(y.poisson ~ fixed1*fixed2 + (1|rand2/rand1), family="poisson", data)
mod4.sqrt <- update(mod4, family = poisson(link = "sqrt"))
rsquared.glmm(list(mod4, mod4.sqrt))
mod5 <- glmer(y.nbinom ~ fixed1*fixed2 + (1|rand2/rand1), family="negative.binomial", data)
rsquared.glmm(mod5)
#Get values for all kinds of models
(lme4.models <- rsquared.glmm(list(mod0, mod1, mod1.1, mod2, mod3, mod3.prob, mod4, mod4.sqrt, mod5)))
