#' R-squared and pseudo-rsquared for a list of (generalized) linear (mixed) models
#'
#' This function calls the generic \code{\link{r.squared}} function for each of the
#' models in the list and rbinds the outputs into one data frame
#'
#' @param a single model or a list of fitted (generalized) linear (mixed) model objects
#' @return a dataframe with one row per model, and "Class",
#'         "Family", "Marginal", "Conditional", "AIC", and "BIC" columns  # BIC added JPB
rsquared.glmm <- function(modlist) {
  if( class(modlist) != "list" ) modlist = list(modlist) else modlist
  # Iterate over each model in the list
  do.call(rbind, lapply(modlist, r.squared))
}

#' R-squared and pseudo-rsquared for (generalized) linear (mixed) models
#'
#' This generic function calculates the r squared and pseudo r-squared for
#' a variety of(generalized) linear (mixed) model fits.
#' Currently implemented for \code{\link{lm}}, \code{\link{lmerTest::merMod}},
#' and \code{\link{nlme::lme}} objects.
#' Implementing methods usually call \code{\link{.rsquared.glmm}}
#'
#' @param mdl a fitted (generalized) linear (mixed) model object
#' @return Implementing methods usually return a dataframe with "Class",
#'         "Family", "Marginal", "Conditional", "AIC" and "BIC" columns  # BIC added JPB
r.squared <- function(mdl){
  UseMethod("r.squared")
}

#' Marginal r-squared for lm objects
#'
#' This method uses r.squared from \code{\link{summary}} as the marginal.
#' Contrary to other \code{\link{r.squared}} methods, 
#' this one doesn't call \code{\link{.rsquared.glmm}}
#'
#' @param mdl an lm object (usually fit using \code{\link{lm}},
#' @return a dataframe with with "Class" = "lm", "Family" = "gaussian",
#'        "Marginal" = unadjusted r-squared, "Conditional" = NA, "AIC" and "BIC" columns #BIC added JPB
r.squared.lm <- function(mdl){
  data.frame(Class=class(mdl), Family="gaussian", Link="identity",
             Marginal=summary(mdl)$r.squared,
             Conditional=NA, AIC=AIC(mdl), BIC=BIC(mdl))   # BIC added JPB
}

#' Marginal and conditional r-squared for merMod objects
#'
#' This method extracts the variance for fixed and random effects, residuals,
#' and the fixed effects for the null model (in the case of Poisson family),
#' and calls \code{\link{.rsquared.glmm}}
#' call with: source("C:\\Users\\Jimb\\Documents\\Birds\\rsquared_glmm.R")
#'
#' @param mdl an merMod model (usually fit using \code{\link{lme4::lmer}},
#'        \code{\link{lme4::glmer}}, \code{\link{lmerTest::lmer}},
#'        \code{\link{blme::blmer}}, \code{\link{blme::bglmer}}, etc)
r.squared.merMod <- function(mdl){
  # Get variance of fixed effects by multiplying coefficients by design matrix
#print("Here 1")
#print(paste(lme4::fixef(mdl)))
#print("Here 1a")
#print(paste("Length t(X)=",length(t(mdl@pp$X))))
#print("Here 1b")
  VarF <- var(as.vector(lme4::fixef(mdl) %*% t(mdl@pp$X)))
#print(paste("VarF =",VarF,"Length =",length(VarF)))
  # Get variance of random effects by extracting variance components
  # Omit random effects at the observation level, variance is factored in later
  DispNames<-gsub("([//(|//)])","",names(ranef(mdl))) # Added JPB - Include intercept?
  MargNames<-names(fixef(mdl)[(1:(length(names(fixef(mdl)))))])
#  MargNames<-names(fixef(mdl)[(2:(length(names(fixef(mdl)))))])
#print("Here 1c")
VarFixList<-(lme4::fixef(mdl))
for(i in names(lme4::fixef(mdl))) {
           VarFixList[i]<-var(skim.glmm2@pp$X[,i]*lme4::fixef(skim.glmm2)[i])
            }
#  VarFixList<-sapply(     # Note \\)|:\\(| added to RegEx string J.P.B - replaced with DispNames
#      VarCorr(mdl)[!sapply(unique(unlist(strsplit(MargNames,":|/"))), function(l) length(unique(mdl@frame[,l])) == nrow(mdl@frame))],
#      function(Sigma) {
#        X <- model.matrix(mdl)
#        Z <- X[,rownames(Sigma)]
#        print(sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X))   #For debugging  
#        sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } 
#        )
#print(paste("Here 1d",VarFixList))
  VarFixList<-data.frame(VarFixList)   # Added for partial R2 - JPB
  namesList<-unique(unlist(strsplit(MargNames,":|/")))
  namesList2<-as.list(unique(unlist(strsplit(MargNames,":|/"))))
#  print(paste("namesList is class",class(namesList)))
#  print(paste("namesList is length",length(namesList)))
#  print(paste("namesList2 is class",class(namesList2)))
#  print(paste("namesList2 is length",length(namesList2)))
#  print(paste("splitting namesList makes it class",class(strsplit(namesList," "))))
#  print(paste("c(\"logNDIN\",\"scCTern\",\"Time\") is class",class(c("logNDIN","scCTern","Time"))))
#  print(c(namesList))
#  print(paste(t(namesList)))
  listNums<-c(1:length(namesList2))
  #print(c(listNums))
   VarFixList<-t(VarFixList)
  for(i in 1:length(namesList2)) {
      #print(paste(i,":",listNums[i],"<-",namesList2[i]))
      names(VarFixList)[i]<-namesList2[i]
      }
#  names(VarFixList)<-c(t(namesList))
#   print(paste(VarFixList))
#   print(paste(t(VarFixList)))
#print(paste("Here 1c, length VarFixList=",length(VarFixList),"variables:",paste(names(VarFixList))))
  VarRand <- sum(
    sapply(     # Note \\)|:\\(| added to RegEx string J.P.B - replaced with DispNames
      VarCorr(mdl)[!sapply(unique(unlist(strsplit(DispNames,":|/"))), function(l) length(unique(mdl@frame[,l])) == nrow(mdl@frame))],
      function(Sigma) {
        X <- model.matrix(mdl)
        Z <- X[,rownames(Sigma)]
        #print(sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X))   #For debugging  
        sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } 
        ) )
#print(paste("Here 2, VarRand=",VarRand))
  # Get the dispersion variance     # Note )&( removed before sending to RegEx string J.P.B
#print(paste("Here 3a, VarCorr=",VarCorr(mdl)))
                           # Added split character blank ('\ :') to strsplit in next line; is this == or != ????? JPB
#        DispNames<-gsub("([//(|//)])","",names(ranef(mdl))) # Added JPB & moved up
#print(paste("Here 3a, VarCorr=",paste(DispNames)))
               # Requiring the # unique Random Variables == total table makes no sense!
               # DispNames with any parr removed added JPB
 VarDisp <- unlist(VarCorr(mdl)[sapply(unique(unlist(strsplit(DispNames,":|/"))), function(l) length(unique(mdl@frame[,l])) == nrow(mdl@frame))])
            # the else part of this statement should not be needed JPB
  if(is.null(VarDisp)) VarDisp = 0 else VarDisp = VarDisp
#print(paste("Here 3b, VarDisp=",VarDisp))
  if(inherits(mdl, "lmerMod")){
    # Get residual variance
    VarResid <- attr(lme4::VarCorr(mdl), "sc")^2
    # Get ML model AIC & BIC
    mdl.aic <- AIC(update(mdl, REML=F))
    mdl.bic <- BIC(update(mdl, REML=F))
    # Model family for lmer is gaussian
    family <- "gaussian"
    # Model link for lmer is identity
    link <- "identity"
  }
  else if(inherits(mdl, "glmerMod")){
    # Get the model summary
    mdl.summ <- summary(mdl)
    # Get the model's family, link, AIC, and BIC # BIC added JPB
    family <- mdl.summ$family
    link <- mdl.summ$link
    mdl.aic <- AIC(mdl)
    mdl.bic <- BIC(mdl)
    # Pseudo-r-squared for poisson also requires the fixed effects of the null model
    if(family=="poisson") {
      # Get random effects names to generate null model
      rand.formula <- reformulate(sapply(findbars(formula(mdl)),
                                         function(x) paste0("(", deparse(x), ")")),
                                  response=".")
      # Generate null model (intercept and random effects only, no fixed effects)
      null.mdl <- update(mdl, rand.formula)
      # Get the fixed effects of the null model
      null.fixef <- as.numeric(lme4::fixef(null.mdl))
    }
   # New Code Here - James P. Browne Dec. 7, 2014
   # for lgmerMod class using the Negative Binomial family from MASS 
   else if(length(i <- grep("Negative Binomial", family))) {
     cat("'Neg-Bin' appears at least once in\n\t", family, "\n")
     i # 2 and 4
#     txt[i]
      # Get random effects names to generate null model
      rand.formula <- reformulate(sapply(findbars(formula(mdl)),
                                         function(x) paste0("(", deparse(x), ")")),
                                  response=".")
      # Generate null model (intercept and random effects only, no fixed effects)
      null.mdl <- update(mdl, rand.formula)
      # Get the fixed effects of the null model
      null.fixef <- as.numeric(lme4::fixef(null.mdl))
   }
#print("Here 4")
  }
  # Call the internal function to do the pseudo r-squared calculations
  .rsquared.glmm(VarF, VarRand, VarResid, VarDisp, family = family, link = link,
                 mdl.aic = mdl.aic,
                 mdl.bic = mdl.bic,
                 mdl.class = class(mdl),
                 null.fixef = null.fixef,
                 VarFix.list = VarFixList) # added JPB for partial R2 calculation
}

#' Marginal and conditional r-squared for lme objects
#'
#' This method extracts the variance for fixed and random effects,
#' as well as residuals, and calls \code{\link{.rsquared.glmm}}
#'
#' @param mdl an lme model (usually fit using \code{\link{nlme::lme}})
r.squared.lme <- function(mdl){
  # Get design matrix of fixed effects from model
  Fmat <- model.matrix(eval(mdl$call$fixed)[-2], mdl$data)
  # Get variance of fixed effects by multiplying coefficients by design matrix
  VarF <- var(as.vector(nlme::fixef(mdl) %*% t(Fmat)))
  # First, extract variance-covariance matrix of random effects
  Sigma.list = VarCorr(mdl)[!grepl(" =",rownames(VarCorr(mdl))) & rownames(VarCorr(mdl)) != "Residual", colnames(VarCorr(mdl))=="Variance", drop=F]
  corr.list = as.numeric(VarCorr(mdl)[!grepl(" =",rownames(VarCorr(mdl))) & rownames(VarCorr(mdl)) != "Residual" & rownames(VarCorr(mdl)) != "(Intercept)",colnames(VarCorr(mdl))=="Corr",drop=F])
  Sigma.list2 = split(as.numeric(Sigma.list), cumsum(rownames(Sigma.list) == "(Intercept)"), drop=F)
  Sigma.list2 = lapply(1:length(Sigma.list2), function(i) { 
    mat = matrix(prod(Sigma.list2[[i]])*abs(corr.list[i]), ncol=length(Sigma.list2[[i]]), nrow=length(Sigma.list2[[i]]))
    diag(mat) = Sigma.list2[[i]]
    colnames(mat) = rownames(Sigma.list)[1:sum(cumsum(rownames(Sigma.list) == "(Intercept)") == 1)]
    rownames(mat) = colnames(mat)
    return(mat) } )
  # Calculate variance of random effects
  VarRand = sum(
    sapply(
      Sigma.list2,
      function(Sigma) {
        Z <- Fmat[,colnames(Sigma),drop=F]
        sum(diag(Z %*% Sigma %*% t(Z)))/nrow(Fmat) } ) )
#print(paste("Here A, VarRand=",VarRand))  
  # Get residual variance
  VarResid <- as.numeric(nlme::VarCorr(mdl)[rownames(nlme::VarCorr(mdl))=="Residual", 1])
#print(paste("Here B, VarResid=",VarResid))
  # Call the internal function to do the pseudo r-squared calculations
  .rsquared.glmm(VarF, VarRand, VarResid, VarDisp, family = "gaussian", link = "identity",
                 mdl.aic = AIC(update(mdl, method="ML")),
                 mdl.bic = BIC(update(mdl, method="ML")), # Added JPB
                 mdl.class = class(mdl))
}

#' Marginal and conditional r-squared for glmm given fixed and random variances
#'
#' This function is based on Nakagawa and Schielzeth (2013). It returns the marginal
#' and conditional r-squared, as well as the AIC for each glmm.
#' Users should call the higher-level generic "r.squared", or implement a method for the
#' corresponding class to get varF, varRand and the family from the specific object
#'
#' @param varF Variance of fixed effects
#' @param varRand Variance of random effects
#' @param varResid Residual variance. Only necessary for "gaussian" family
#' @param family family of the glmm (currently works with gaussian, binomial and poisson)
#' @param link model link function. Working links are: gaussian: "identity" (default);
#'        binomial: "logit" (default), "probit"; poisson: "log" (default), "sqrt"
#' @param mdl.aic The model's AIC
#' @param mdl.bic The model's BIC   # BIC added JPB
#' @param mdl.class The name of the model's class
#' @param null.fixef Numeric vector containing the fixed effects of the null model.
#'        Only necessary for "poisson" family
#' @return A data frame with "Class", "Family", "Marginal", "Conditional", "AUC", and "BIC" columns
#     VarFix.list added for partial R2 - JPB
.rsquared.glmm <- function(varF, varRand, varResid = NULL, varDisp = NULL, family, link,
                           mdl.aic, mdl.bic, mdl.class, null.fixef = NULL, VarFix.list = NULL){
  if(family == "gaussian"){
    # Only works with identity link
    if(link != "identity")
      family_link.stop(family, link)
    # Calculate marginal R-squared (fixed effects/total variance)
    Rm <- varF/(varF+varRand+varResid)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+varResid)
  }
  else if(family == "binomial"){
    # Get the distribution-specific variance
    if(link == "logit")
      varDist <- (pi^2)/3
    else if(link == "probit")
      varDist <- 1
    else
      family_link.stop(family, link)
    # Calculate marginal R-squared
    Rm <- varF/(varF+varRand+varDist+varDisp)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+varDist+varDisp)
  }
  else if(family == "poisson"){
    # Get the distribution-specific variance
    if(link == "log")
      varDist <- log(1+1/exp(null.fixef))
    else if(link == "sqrt")
      varDist <- 0.25
    else
      family_link.stop(family, link)
    # Calculate marginal R-squared
    Rm <- varF/(varF+varRand+varDist+varDisp)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+varDist+varDisp)
  }
  # Added working code Jim Browne based on poisson
   else if(length(i <- grep("Negative Binomial", family))) {
#print(paste("Here 5, Link=",link,", null model fixed effect=",null.fixef))
    if(link == "log")
#      if(is.na(null.fixef) | ????? # Add error mesage if no null.fixef
      varDist <- log(1+1/exp(null.fixef))
    else if(link == "sqrt")  # Why is sqrt varDist = 0.25 ?
      varDist <- 0.25
#    else if(link == "identity")  # What should this be??? Only for Gausian!
#      varDist <- 1.0
    else
      family_link.stop(family, link)
#print(paste("Here 5a, varDist=",varDist,", varF=",varF,", varRand=",varRand,", varDisp=",varDisp))

    # Calculate marginal R-squared
    Rm <- varF/(varF+varRand+varDist+varDisp)
    # Calculate conditional R-squared (fixed effects+random effects/total variance)
    Rc <- (varF+varRand)/(varF+varRand+varDist+varDisp)
#print(paste("Here 6, Rm=",Rm,"Rc=",Rc))
#print("Here 6a calculate partial R2") # Calculation and printing of Partial R2 added JPB
#print(paste(names(VarFix.list)))
#print(paste(VarFix.list))
  for(i in names(VarFix.list)) {
       Rp<-VarFix.list[i]/(varF+varRand+varDist+varDisp)
       print(paste("Partial R2 ",i," = ",Rp,", ",(Rp*100/Rm),"% of marginal"),sep="")
        }

   }
  else
    family_link.stop(family, link)
  # Bind R^2s into a matrix and return with AIC values
  #  BIC added, JPB
  data.frame(Class=mdl.class, Family = family, Link = link,
             Marginal=Rm, Conditional=Rc, AIC=mdl.aic, BIC=mdl.bic)
}

#' stop execution if unable to calculate variance for a given family and link
family_link.stop <- function(family, link){
  stop(paste("Don't know how to calculate variance for",
             family, "family and", link, "link."))
}

HedgesG.merMod<-function(mdl) {
    Coef.mdl<-summary(mdl)$coefficients
    rows<-length(Coef.mdl[,1])
    fix_t<-Coef.mdl[2:rows,3]
    df<-1 # Get improvement
    n1<-length(residuals(skim.glmm2))   # Placeholder - correct this!!!
    n2<-length(residuals(skim.glmm2))   # Placeholder - correct this!!!
    r.part<-sqrt(fix_t^2/(fix_t^2+df))
    names(r.part)<-names(fix_t)
    print(r.part)
    HeG<-(r.part/sqrt(1-r.part^2))/sqrt(df*(n1+n2)/n1*n2)
    print(HeG)
   }

# Code from Christopher Moore on September 7, 2010
# http://blog.lib.umn.edu/moor0554/canoemoore/
# NOTE: nHe claims that multiple regression doesn't work, 
# but points to a list question about supressing residuals?
# corections from Laurie Samuels, Vanderbilt University
# The p.value.LRT[i] line now needs to be "p.value.LRT[i]  ????
#Steve Brady | October 20, 2010 2:49 PM | Reply
#This works well for singular random effects.
#Can the function be modified to handle nested random effects? At the moment, when I try it with nested random effects, I receive the following error:
#Error in names(data.ranef) 'names' attribute [6] must be the same length as the vector [2]

p.values.lmer <- function(x) {
  summary.model <- summary(x)
  data.lmer <- data.frame(model.matrix(x))
  names(data.lmer) <- names(fixef(x))
  names(data.lmer) <- gsub(pattern=":", x=names(data.lmer), replacement=".", fixed=T)
  names(data.lmer) <- ifelse(names(data.lmer)=="(Intercept)", "Intercept", names(data.lmer))
  string.call <- strsplit(x=as.character(x@call), split=" + (", fixed=T)
  var.dep <- unlist(strsplit(x=unlist(string.call)[2], " ~ ", fixed=T))[1]
  vars.fixef <- names(data.lmer)
  formula.ranef <- paste("+ (", string.call[[2]][-1], sep="")
  formula.ranef <- paste(formula.ranef, collapse=" ")
  formula.full <- as.formula(paste(var.dep, "~ -1 +", paste(vars.fixef, collapse=" + "), 
                  formula.ranef))
  data.ranef <- data.frame(x@frame[, 
                which(names(x@frame) %in% names(ranef(x)))])
  names(data.ranef) <- names(ranef(x))
  data.lmer <- data.frame(x@frame[, 1], data.lmer, data.ranef)
  names(data.lmer)[1] <- var.dep
  out.full <- lmer(formula.full, data=data.lmer, REML=F)
  p.value.LRT <- vector(length=length(vars.fixef))
  for(i in 1:length(vars.fixef)) {
    formula.reduced <- as.formula(paste(var.dep, "~ -1 +", paste(vars.fixef[-i], 
                       collapse=" + "), formula.ranef))
    out.reduced <- lmer(formula.reduced, data=data.lmer, REML=F)
    print(paste("Reduced by:", vars.fixef[i]))
    print(out.LRT <- data.frame(anova(out.full, out.reduced)))
    p.value.LRT[i] <- round(out.LRT[2, 7], 3)
  }
  summary.model$coefficients <- cbind(summary.model$coefficients, p.value.LRT)
  summary.model$methTitle <- c("\n", summary.model$methTitle, 
                           "\n(p-values from comparing nested models fit by maximum likelihood)")
  print(summary.model)
}

plotEffects.merMod <- function(fit1) {
    require(ggplot2)
    randoms<-ranef(fit1, postVar = TRUE)
    qq <- attr(ranef(fit1, postVar = TRUE)[[1]], "postVar")
    rand.interc<-randoms$Batch
    df<-data.frame(Intercepts=randoms$Batch[,1],
              sd.interc=2*sqrt(qq[,,1:length(qq)]),
              lev.names=rownames(rand.interc))
    df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)])
    p <- ggplot(df,aes(lev.names,Intercepts,shape=lev.names))

    #Added horizontal line at y=0, error bars to points and points with size two
    p <- p + geom_hline(yintercept=0) +geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0,color="black") + geom_point(aes(size=2)) 

    #Removed legends and with scale_shape_manual point shapes set to 1 and 16
    p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=c(1,1,1,16,16,16))

    #Changed appearance of plot (black and white theme) and x and y axis labels
    p <- p + theme_bw() + xlab("Levels") + ylab("")

    #Final adjustments of plot
    p <- p + theme(axis.text.x=element_text(size=rel(1.2)),
                   axis.title.x=element_text(size=rel(1.3)),
                   axis.text.y=element_text(size=rel(1.2)),
                   panel.grid.minor=element_blank(),
                   panel.grid.major.x=element_blank())

    #To put levels on y axis you just need to use coord_flip()
    p <- p+ coord_flip()
    print(p)
  }
## re = object of class ranef.mer
ggCaterpillar <- function(re, QQ=TRUE, likeDotplot=TRUE) {
    require(ggplot2)
    f <- function(x) {
        pv   <- attr(x, "postVar")
        cols <- 1:(dim(pv)[1])
        se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
        ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
        pDf  <- data.frame(y=unlist(x)[ord],
                           ci=1.96*se[ord],
                           nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                           ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                           ind=gl(ncol(x), nrow(x), labels=names(x)))

        if(QQ) {  ## normal QQ-plot
            p <- ggplot(pDf, aes(nQQ, y))
            p <- p + facet_wrap(~ ind, scales="free")
            p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
        } else {  ## caterpillar dotplot
            p <- ggplot(pDf, aes(ID, y)) + coord_flip()
            if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
                p <- p + facet_wrap(~ ind)
            } else {           ## different scales for random effects
                p <- p + facet_grid(ind ~ ., scales="free_y")
            }
            p <- p + xlab("Levels") + ylab("Random effects")
        }

        p <- p + theme(legend.position="none")
        p <- p + geom_hline(yintercept=0)
        p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
        p <- p + geom_point(aes(size=1.2), colour="blue") 
        return(p)
    }

    lapply(re, f)
}
##################### find overdispersion
# From DRAFT r-sig-mixed-models FAQ
# http://glmm.wikidot.com/faq
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
#################################
## re = object of class ranef.mer
ggCaterpillar <- function(re, QQ=TRUE, likeDotplot=TRUE) {
    require(ggplot2)
    f <- function(x) {
        pv   <- attr(x, "postVar")
        cols <- 1:(dim(pv)[1])
        se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
        ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
        pDf  <- data.frame(y=unlist(x)[ord],
                           ci=1.96*se[ord],
                           nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                           ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                           ind=gl(ncol(x), nrow(x), labels=names(x)))

        if(QQ) {  ## normal QQ-plot
            p <- ggplot(pDf, aes(nQQ, y))
            p <- p + facet_wrap(~ ind, scales="free")
            p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
        } else {  ## caterpillar dotplot
            p <- ggplot(pDf, aes(ID, y)) + coord_flip()
            if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
                p <- p + facet_wrap(~ ind)
            } else {           ## different scales for random effects
                p <- p + facet_grid(ind ~ ., scales="free_y")
            }
            p <- p + xlab("Levels") + ylab("Random effects")
        }

        p <- p + theme(legend.position="none")
        p <- p + geom_hline(yintercept=0)
        p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
        p <- p + geom_point(aes(size=1.2), colour="blue") 
        return(p)
    }

    lapply(re, f)
}
#determinant(summary(skim.glmm2)$vcov)
#names(summary(skim.glmm2))
# [1] "methTitle"    "objClass"     "devcomp"      "isLmer"       "useScale"    
# [6] "logLik"       "family"       "link"         "ngrps"        "coefficients"
#[11] "sigma"        "vcov"         "varcor"       "AICtab"       "call"        
#[16] "residuals"   

# as.function(skim.glmm2)
#function (pars) 
#{
#    resp$setOffset(baseOffset)
#    resp$updateMu(lp0)
#    pp$setTheta(as.double(pars[dpars]))
#    spars <- as.numeric(pars[-dpars])
#    offset <- if (length(spars) == 0) 
#        baseOffset
#    else baseOffset + pp$X %*% spars
#    resp$setOffset(offset)
#    p <- pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, 
#        verbose)
#    resp$updateWts()
#    p
#}
#<environment: 0x23d3662c>


# Modified theta.m** functions from MASS for use with {lme4}
# Use residual degrees freedom estimate from overdispersion function (above)
# from: http://glmm.wikidot.com/faq

df.residual.glmerMod<-function(y) {
       vpars <- function(m) {
          nrow(m)*(nrow(m)+1)/2
           }
        mu <- fitted(y)
        # Next 2 lines from: DRAFT r-sig-mixed-models FAQ
        model.df <- sum(sapply(VarCorr(y),vpars))+length(fixef(y))
        rdf <- nrow(model.frame(y))-model.df
        rdf
}

theta.glmer.mj<-function(mod) { #,theta.seed=NULL) {
#     res<-residuals(mod)
#     fit<-fitted(mod)
      mnres<-mean(residuals(mod))
      mdres<-median(residuals(mod))
      theta.seed<-as.numeric(t(data.frame(strsplit(gsub("\\)|\\(",",",summary(mod)$family),",")))[2])
##      print(paste("old theda =",theta.seed,", mean res =",mnres,", median res =",mdres))
#      theta.new<-theta.seed * mean(res)/median(res)
#      theta.new<-theta.seed * 0.25 * mean(residuals(mod))/median(residuals(mod))
      theta.new <- theta.seed+(mdres-mnres)/(5*mnres)
#      theta.new2 <- theta.seed+(mdres-mnres)/(5*mdres)
#      c(theta.new,theta.new2)
     theta.new
     }

Loop.fit<-function(mod, lim=0.0025) {
     newTheta<-theta.glmer.mj(mod)
     print(paste("New theta =",newTheta))
     oldTheta<-as.numeric(t(data.frame(strsplit(gsub("\\)|\\(",",",summary(mod)$family),",")))[2])
     while(abs(as.numeric(oldTheta)-as.numeric(newTheta))>lim) {
#        print(paste("old theta =",oldTheta,", New theta=",newTheta))
         newCallParts<-gsub(oldTheta,newTheta,summary(mod)$call)
         newCall<-paste("glmer(",newCallParts[2],", data =",newCallParts[3],", family =",newCallParts[4],")")
        print(newCall)
         mod<-eval(parse(text=newCall))
         oldTheta<-newTheta
         newTheta<-theta.glmer.mj(mod)
      }
     newTheta
     } 

theta.glmer.md <-
    function(y, mu=NULL, dfr=NULL, weights=NULL, limit = 20, eps = .Machine$double.eps^0.25)
{
    #require(pbkrtest) #KRmodcomp() also fails with "glmerMod"
#    print("Here Td1")
    if(inherits(y, "glmerMod")) {
       vpars <- function(m) {
          nrow(m)*(nrow(m)+1)/2
           }
#       print("Here Td2")
        mu <- fitted(y)
#       print("Here Td3")
        # Next 2 lines from: DRAFT r-sig-mixed-models FAQ
        model.df <- sum(sapply(VarCorr(y),vpars))+length(fixef(y))
        rdf <- nrow(model.frame(y))-model.df
        dfr<-rdf
        #dfr <- y$df.residual
#        print("Here Td4")
##        y <- if(is.null(y$y)) mu + residuals(y) else y$y
          y <- mu + residuals(y)
#       print("Here Td5")
    }
    if(missing(weights)) weights <- rep(1, length(y))
    n <- sum(weights)
#   print(paste("length y=",length(y),", length mu=",length(mu)))
    t0 <- n/sum(weights*(y/mu - 1)^2)
    a <- 2 * sum(weights*y * log(pmax(1, y)/mu)) - dfr
    it <- 0
    del <- 1
    while((it <- it + 1) < limit && abs(del) > eps) {
        t0 <- abs(t0)
        tmp <- log((y + t0)/(mu + t0))
        top <- a - 2 * sum(weights*(y + t0) * tmp)
        bot <- 2 * sum(weights*((y - mu)/(mu + t0) - tmp))
        del <- top/bot
        t0 <- t0 - del
    }
    if(t0 < 0) {
        t0 <- 0
        warning("estimate truncated at zero")
        attr(t0, "warn") <- gettext("estimate truncated at zero")
    }
    t0
}

theta.glmer.ml <-
    function(y, mu, n = sum(weights), weights, limit = 10,
             eps = .Machine$double.eps^0.25,
             trace = FALSE)
{
    score <- function(n, th, mu, y, w)
        sum(w*(digamma(th + y) - digamma(th) + log(th) +
               1 - log(th + mu) - (y + th)/(mu + th)))
    info <- function(n, th, mu, y, w)
        sum(w*( - trigamma(th + y) + trigamma(th) - 1/th +
               2/(mu + th) - (y + th)/(mu + th)^2))
    if(inherits(y, "glmerMod")) {
        mu <- fitted(y)
#        y <- if(is.null(y$y)) mu + residuals(y) else y$y
         y <- mu + residuals(y)
    }
    if(missing(weights)) weights <- rep(1, length(y))
    t0 <- n/sum(weights*(y/mu - 1)^2)
    it <- 0
    del <- 1
    if(trace) message(sprintf("theta.ml: iter %d 'theta = %f'",
                              it, signif(t0)), domain = NA)
    while((it <- it + 1) < limit && abs(del) > eps) {
        t0 <- abs(t0)
        del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, mu, y, weights))
        t0 <- t0 + del
        if(trace) message("theta.ml: iter", it," theta =", signif(t0))
    }
    if(t0 < 0) {
        t0 <- 0
        warning("estimate truncated at zero")
        attr(t0, "warn") <- gettext("estimate truncated at zero")
    }
    if(it == limit) {
        warning("iteration limit reached")
        attr(t0, "warn") <- gettext("iteration limit reached")
    }
    attr(t0, "SE") <- sqrt(1/i)
    t0
}

theta.glmer.mm <- function(y, mu, dfr, weights, limit = 10,
                     eps = .Machine$double.eps^0.25)
{
    if(inherits(y, "glmerMod")) {
       vpars <- function(m) {
          nrow(m)*(nrow(m)+1)/2
           }
        mu <- fitted(y)
        # Next 2 lines from: DRAFT r-sig-mixed-models FAQ
        model.df <- sum(sapply(VarCorr(y),vpars))+length(fixef(y))
        rdf <- nrow(model.frame(y))-model.df
        dfr<-rdf
        #dfr <- y$df.residual
        #y <- if(is.null(y$y)) mu + residuals(y) else y$y
        y <- mu + residuals(y)
    }
    if(missing(weights)) weights <- rep(1, length(y))
    n <- sum(weights)
    t0 <- n/sum(weights*(y/mu - 1)^2)
    it <- 0
    del <- 1
        print(paste(" mu=",mu," dfr=",dfr))
    while((it <- it + 1) < limit && abs(del) > eps) {
        t0 <- abs(t0)
        print(paste(" t0=",t0))
        del <- (sum(weights*((y - mu)^2/(mu + mu^2/t0))) - dfr)/
            sum(weights*(y - mu)^2/(mu + t0)^2)
        print(paste(" del=",del))
        t0 <- t0 - del
    }
    if(t0 < 0) {
        t0 <- 0
        warning("estimate truncated at zero")
        attr(t0, "warn") <- gettext("estimate truncated at zero")
    }
    t0
}



# Get printout of this Fox page
#http://rstudio-pubs-static.s3.amazonaws.com/10858_19492a31e76245e292aaaae235d3efa8.html

