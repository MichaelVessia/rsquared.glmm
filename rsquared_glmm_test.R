#Creates a negative binomial dataset using the sample Quine data
quine.nb <- glmer(Days ~ 1|Sex, family=negative.binomial(theta = 2, link = "log"),  data = quine)

#Calculate the rsquared value
rsquared.glmm(quine.nb)
