#### Last updated: 12 December 2019; 13 Jan 2020

#### JLH 13 Jan 2020: Make sure that we have removed any negative feeding rates estimates in the calculation below ####

rm(list = ls())

library(magrittr)
library(pomp)

## load data
#setwd("~/Dropbox/Alaina/MaxLikelihood")
data <- read.csv("means_4MaxLikelihood_AJV2.csv", header = T, sep = ",")
names(data)
## 'data' is a dataframe with columns:
## 'CloneID' host genotype
## 'S0' number of susceptible hosts (does not vary)
## 'P0' parasite/pathogen spore dose
## 'Tube' vial number
## 'plate_treatment' an identifier for each clone-by-spore treatment
## 'A0' Chla reading in a sample WITHOUT an animal: control reading - note, each genotyope-spore combination has a unique control reading
## 'mm' gives the size of the animal (convert to surface area: L^2)
## 'A1' Chla reading in a sample WITH an animal
## 'I1' infection outcome ('S1' and 'E1' are meaningless and can be ignored)

## remove one row that is problematic because it is missing data
data <- data[-which(is.na(data),arr.ind=TRUE)[1,'row'],]


#### All of this code needs to loop through multiple genotypes, and needs to match each individual's feeding rate with that animal's length, infection status (at the end of the experiment), spore exposure and plate-specific control (A0) value
#### Bolker's book has this cool trick which will help with some of this.

## 'params' is a named vector of parameters f, sd, and u
## For the purposes of the fitting, these parameters should be
## able to take any values from -Inf to Inf, so they must be
## transformed using the log (for f and sd) or logit (for u)
## It is useful to sweep across a large number of parameter
## combinations before attempting any ML fitting. The code below
## sets up 'nseq' initial parameter estimates spread across
## a hyperdimensional volume bounded by lower and upper estimates
## of each parameter.
## CEC: I removed the A01-A04 "parameters" because these are not
## actually estimated - they do not appear in any of the calculations
## in nll2 below. Rather, these are the plate-specific initial algae
## concentrations.....####  JLH 13 Jan 2020, they aren't estimated but I think they are needed to calcuate the actual feeding rate?? ####
## CEC: I also realized that, with the current model formulation, the
## genotype-specific feeding parameters (f11-f14) need to be transformed
## in order to vary between -Inf and Inf

#### provide initial values and bounds for the MLE optimizer ####
sobolDesign(lower = c(u1 = log(0.003/(1-0.003)),## genotype-specific per-spore susceptibility
                      u2 = log(0.003/(1-0.003)),
                      u3 = log(0.003/(1-0.003)),
                      u4 = log(0.003/(1-0.003)),
                      f11 = log(0.0001), ## genotype-specific feeding rate
                      f12 = log(0.0001),
                      f13 = log(0.0001),
                      f14 = log(0.0001),
                      a = -0.01, ## per-spore affect on feeding rate
                      sd = log(0.1)), ## observation error in feeding
            upper = c(u1 = log(0.3/(1-0.3)),
                      u2 = log(0.3/(1-0.3)),
                      u3 = log(0.3/(1-0.3)),
                      u4 = log(0.3/(1-0.3)),
                      f11 = log(0.1),
                      f12 = log(0.1),
                      f13 = log(0.1),
                      f14 = log(0.1),
                      a = 0.01,
                      sd = log(0.1)),
            nseq =200) -> pars

## compute the negative log-likelihood
## 'params' is a named vector of parameter values
## 'data' is the full dataset
## 'mod' is a character string giving the name of the feeding model to fit
nll2 <- function(params, data, mod) {
    T <- 7/24 ## length of feeding/exposure in days

    ## For each CloneID value, compute the NLL of its set of parameter values and add this to the existing total
    nll <- 0 ## for storage of the final NLL
    for (id in unique(data$CloneID)) {
        ## transform genotype-specific feeding rate back to the natural scale and choose the fhat value for this genotype
        fhat <- c(exp(params["f11"]),
                  exp(params["f12"]),
                  exp(params["f13"]),
                  exp(params["f14"]))[id] %>% unname
        ## per-spore effect on feeding can already vary between -Inf (a<0 implies a decrease in feeding rate with spores) and Inf (a>0 implies an increase in feeding rate with spores)
        a  <- params["a"] %>% unname
        ## transform observation error to the natural scale
        sd <- exp(params["sd"]) %>% unname
        ## transform genotype-specific susceptibility back to the natural scale and choose the u value for this genotype
        u  <- c(exp(params["u1"])/(1 + exp(params["u1"])),
                exp(params["u2"])/(1 + exp(params["u2"])),
                exp(params["u3"])/(1 + exp(params["u3"])),
                exp(params["u4"])/(1 + exp(params["u4"])))[id] %>% unname

        ## extract the data for this genotype
        data2 <- subset(data, CloneID == id)

        ## compute the size- and spore-dependent feeding rate
        ## CEC: if you wanted to use 'mod' to consider different possible feeding models, here is how you would do it:
        ## switch(mod,
        ##        mod1 = fhat*(data2$mm^2)*exp(a*data2$P0),
        ##        mod2 = fhat*(data2$mm^2)*(1+a*data2$P0),
        ##        mod3 = ...) -> f_pred
        ## where the function variable 'mod' is a character string that must
        ## match the model names in the switch function ('mod1', 'mod2', etc.)
        ## You can specify as many different feeding models as you want
        ## I'll just use the two hypothetical models above (one with an
        ## exponential change in feeding with spores, and one with a linear
        ## change) as an example:
        ## JLH 13 Jan 2020: I really thought this should be negative exponential function... shoot, I don't remember why you changed this??
        ## CEC 15 Jan 2020: the parameter a can take on any value from -Inf to Inf, so in practice, it likely will be a negative exponential (because a will most likely be a negative number), but I didn't see any reason to force a to be a negative number.
        f_pred <- switch(mod,
                         mod1=fhat*(data2$mm^2)*exp(a*data2$P0),
                         mod2=fhat*(data2$mm^2)*(1+a*data2$P0))

        ## if any of the feeding rates are negative, this can generate a negative prob of infection,  breaking dbinom and generating errors
        ## these are bad parameter sets, so simply assign them very high -logLik values. Note: you shouldn't assign a log-likelihood of 0, as although the likelihood must be a positive number, the log-likelihood can be negative if the likelihood is small, which it often will be when the fitting algorithm gets started. Thus, assigning a value of 0 might actually be saying that this parameter set isn't terrible whereas assigning it a value of Inf says that this parameter set is the worst.
        if (any (f_pred < 0))
            nll <- nll + Inf

        ## Otherwise, estimate the amount of algae that is remaining
        else {

            ## Compute the expected amount of food remaining, assuming that food in the vial changes according to the ODE dA/dt = -f_pred*A, which has the solution A(t) = A(0)*exp(-f_pred*t). Then, after T time units, A(T) = A(0)*exp(-f_pred*T). Taking logarithms, you have log(A(T)) = log(A(0)) - f_pred*T
            exp.A1 <- log(data2$A0) - f_pred*T

            ## compute the expected probability of infection, given feeding and per-spore susceptibility
            lambda <- 1 - exp(-u*f_pred*data2$P0)
            ## Compute the negative log-liklihood of observing log(data2$A1) algae remaining, under the expectation exp.A1, and the negative log-likelihood of observing the binary infection outcome data2$I1, under the expected probability of infection lambda; sum these two NLLs, and add that sum to the total NLL across all genotypes
            nll <- -sum(dnorm(log(data2$A1),
                              mean = exp.A1,
                              sd = sd,
                              log = TRUE) %>% sum,
                        dbinom(data2$I1,
                               size = 1,
                               prob = lambda,
                               log = TRUE) %>% sum) + nll
        }
    }
    return(nll)
}


## calculate the likelihood of all the parameter sets using nll2
## CEC: to do this we will do the following
## the function apply(pars,1,function(p)...) says:
## 1. Take the information in the dataframe 'pars' and
## 2. For each row of that dataframe (specify performing an option on the row by setting the second value passed to apply equal to 1)
## 3. Pass the values in that row to the function 'nll2'. Thus 'p' in function(p) is the current row in pars. You have to do the unlist(p) because pars is a dataframe, whereas nll2 is expected 'pars' to be a named numeric vector. (You can see this by running class(pars[1,]) versus class(unlist(pars[1,]))

## Conceptually, what apply is doing is doing the following for loop:
## ll <- vector(length=nrow(pars))
## for (i in 1:nrow(pars)) {
##      p <- unlist(pars[i,])
##      ll[i] <- nll2(p, data, mod='mod1')
## }
apply(pars,
      1,
      function(p) nll2(unlist(p), data, mod ='mod1')) -> ll #### JLH: 13 Jan 2020: This should stay data?
##JLH the error message says object f not found....
## CEC: This is because there was a spot in nll2 that required an object named 'f' but no such variable existed. Whenever you see this error it means that there's something being called inside your function that R doesn't have a value for.

## from tutorial:
## optimize only from the other starting points
## use Nelder-Mead with fixed RNG seed
# fit <- optim(
#   par=c(log(2), log(1), log(0.9/(1-0.9))),
#   est=c("Beta","mu_I","rho"),
#   fn=neg.ll,
#   method="Nelder-Mead",
#   control=list(maxit=400,trace=0)
# )

## remove all parameter sets whose -logLik was Inf
pars2 <- pars[which(ll!=Inf),]
## CEC: for each of the parameter guesses that produced a finite NLL,
## start a Nelder-Mead optimizer to find the parameter values that minimize
## the NLL. Again, you need to use unlist to turn the parameter dataframe into
## a numeric vector.
## create a list to store the output from each Nelder-Mead optimization
fits  <- vector(mode = 'list', length = nrow(pars2))
for(i in 1:nrow(pars2)) {
    print(i)
    optim(par = unlist(pars2[i,]), ##JLH what's this loop doing? It's going through each row of the "useable" NLL estimates and unlisting them?
          fn = nll2,
          data = data,
          mod='mod1',
          control=list(maxit = 5000)
          ) -> fits[[i]] ## then what's this doing?
}

## CEC: before proceeding, there's something else we need to do here, which is check for convergence. For each output from optim stored in fits, we need to look at the convergence value. If convergence is anything other than 0, it means that the optimizer had not yet converged. For example, if fits[[1]]$convergence = 1, it means that the optimizer hit the max. number of function iterations and gave up without converging. When maxit was set to 200 previously, I was getting this convergence value for every parameter set. When I increased maxit to 5000 (which obviously really slows everything down), I started to see convergence. Of course, this greatly, greatly increased how long the for loop takes to run!

## The following line tells which of the fits actually converged. Those are the only ones to pay attention to.
good.ests <- which(unlist(lapply(fits, function(f) f$convergence))==0)
fits2 <- fits[good.ests]

## Extract the parameter estimates from fits2
## We'll do this using the function lapply, since fits2 is a list
## What lapply(fits2, function(f) f$pars) is doing is essentially
## the following for loop:
## f <- vector(mode='list', length=length(fits))
## for (i in 1:length(fits2))
##      f[[i]] <- fits2[[i]]$pars
## What I want to do with those parameter estimates is put them into a data.frame where each row is the a set of optimized parameter values
## I am doing this by taking the list f and unlisting it, which turns it from a list into a vector with a length equal to (the number of estimated parameters) * (the length of fits). I want to just put this into a matrix with a number of columns that is equal to the number of parameters that I am estimating. The 'matrix' function will essentially take this really long vector and "wrap" it into the matrix by taking the first (number of parameters) values and putting them into the first row, then taking values (number of parameters)+1 to 2*(number of parameters) and putting them into the second row, and so on. I then turn the matrix into a dataframe.## I could do this with the extra line
## parsets <- as.data.frame(matrix(unlist(f), ncol=ncol(pars), byrow=T))
## I've just put all of that into the following more concise form

lapply(fits2, function(f) f$par) %>% ## JLH should this be f or f_hat?? CEC: it should be f, not f_hat!
    unlist %>%
    matrix(., ncol = ncol(pars), byrow = T) %>% ## JLH why 5 columns? can we make this more general?
    ## CEC: Yes, you are right. ncol should be equal to the number of parameters you are trying to estimate, which is the same as the number of columns of pars
    as.data.frame -> parsets
colnames(parsets) <- names(fits[[1]]$par)

library(plyr)
## transform the parameter estimates back to the natural scale
mutate(parsets,
       u1 = exp(u1)/(1 + exp(u1)),
       u2 = exp(u2)/(1 + exp(u2)),
       u3 = exp(u3)/(1 + exp(u3)),
       u4 = exp(u4)/(1 + exp(u4)),
       f11 = exp(f11),
       f12 = exp(f12),
       f13 = exp(f13),
       f14 = exp(f14),
       sd = exp(sd)) -> parsets
## One last thing to do: pull out the NLL for each parameter set. Each initial guess converges to a different final spot, and these may differ a lot in NLL. We're only interested in parameter sets with high likelihood.
parsets$nll <- lapply(fits2, function(f) f$value) %>% unlist
range(parsets$nll) ## you can see that there is a *huge* range here

## pull out just the best parameter set
bestpars <- unlist(parsets[which.min(parsets$nll),])

## Using this best-fitting parameter set, let's compare the predicted and observed algae levels.
## storage for the predicted algae measurements
pred.A1 <- vector(length=nrow(data))
for (i in 1:nrow(data)) {
    ## the current genotype
    g <- data$CloneID[i]
    ## pull the genotype-specific parameter estimate for fhat
    fhat <- unname(bestpars[paste0("f1",g)])
    ## pull the estimate of a
    a <- unname(bestpars["a"])
    ## pull the observed length
    l <- data$mm[i]
    ## pull the observed initial algae measurement
    A0 <- data$A0[i]
    ## pull the observed final algae measurement
    A1 <- data$A1[i]
    ## pull the observed spore load
    P0 <- data$P0[i]
    ## compute the predicted algae
    pred.A1[i] <- log(A0) - (fhat*l^2*exp(a*P0))*(7/24)
}

## plot the observed A1 against the predicted A1
plot(log(data$A1), pred.A1, col=data$CloneID, pch=21, bg=data$CloneID)
abline(a=0,b=1)

## We can also compute the probability of infection for each genotype at each spore level and compare that against the observed
## Of course, the probability of infection depends on length, and every individual was a slightly different length at the time of exposure, but we'll just do it for a fixed length, just to see.
## Here is the feeding rate for the mean length individual when there are 150 spores around
f.g1.p150 <- bestpars['f11']*mean(data$mm)^2*exp(bestpars['a']*150)
## Here is what that translates into in terms of infection probability
lambda.g1.p150 <- 1-exp(-bestpars['u1']*f.g1.p150*150)
## Note that we must multiply by 150 because there are no "150 values" - nll2 doesn't estimate different values of the parameters as spore load changes.
##** JLH  why multiplying for 150 here? why not just select...the 150 values?

## Do the same thing for the 300 spore dose
f.g1.p300 <- bestpars['f11']*mean(data$mm)^2*exp(bestpars['a']*300)
lambda.g1.p300 <- 1-exp(-bestpars['u1']*f.g1.p300*300)

## Repeat for the other three genotypes
f.g2.p150 <- bestpars['f12']*mean(data$mm)^2*exp(bestpars['a']*150)
lambda.g2.p150 <- 1-exp(-bestpars['u2']*f.g2.p150*150)
f.g2.p300 <- bestpars['f12']*mean(data$mm)^2*exp(bestpars['a']*300)
lambda.g2.p300 <- 1-exp(-bestpars['u2']*f.g2.p300*300)

f.g3.p150 <- bestpars['f13']*mean(data$mm)^2*exp(bestpars['a']*150)
lambda.g3.p150 <- 1-exp(-bestpars['u3']*f.g3.p150*150)
f.g3.p300 <- bestpars['f13']*mean(data$mm)^2*exp(bestpars['a']*300)
lambda.g3.p300 <- 1-exp(-bestpars['u3']*f.g3.p300*300)

f.g4.p150 <- bestpars['f14']*mean(data$mm)^2*exp(bestpars['a']*150)
lambda.g4.p150 <- 1-exp(-bestpars['u4']*f.g4.p150*150)
f.g4.p300 <- bestpars['f14']*mean(data$mm)^2*exp(bestpars['a']*300)
lambda.g4.p300 <- 1-exp(-bestpars['u4']*f.g4.p300*300)

## Observed infection prevalences for each genotype, at each spore load
prev.g1.p150 <- sum(subset(data, CloneID==1 & P0==150)$I1)/nrow(subset(data, CloneID==1 & P0==150))
prev.g1.p300 <- sum(subset(data, CloneID==1 & P0==300)$I1)/nrow(subset(data, CloneID==1 & P0==300))
prev.g2.p150 <- sum(subset(data, CloneID==2 & P0==150)$I1)/nrow(subset(data, CloneID==2 & P0==150))
prev.g2.p300 <- sum(subset(data, CloneID==2 & P0==300)$I1)/nrow(subset(data, CloneID==2 & P0==300))
prev.g3.p150 <- sum(subset(data, CloneID==3 & P0==150)$I1)/nrow(subset(data, CloneID==3 & P0==150))
prev.g3.p300 <- sum(subset(data, CloneID==3 & P0==300)$I1)/nrow(subset(data, CloneID==3 & P0==300))
prev.g4.p150 <- sum(subset(data, CloneID==4 & P0==150)$I1)/nrow(subset(data, CloneID==4 & P0==150))
prev.g4.p300 <- sum(subset(data, CloneID==4 & P0==300)$I1)/nrow(subset(data, CloneID==4 & P0==300))

## Plot observed against predicted
plot(c(prev.g1.p150,
       prev.g1.p300,
       prev.g2.p150,
       prev.g2.p300,
       prev.g3.p150,
       prev.g3.p300,
       prev.g4.p150,
       prev.g4.p300),
     c(lambda.g1.p150,
       lambda.g1.p300,
       lambda.g2.p150,
       lambda.g2.p300,
       lambda.g3.p150,
       lambda.g3.p300,
       lambda.g4.p150,
       lambda.g4.p300),
     pch=21,
     cex=2,
     bg=c(1,1,2,2,3,3,4,4),
     xlab="Obseved prevalence",
     ylab="predicted prevalence")
abline(0,1)


# https://kingaa.github.io/short-course/pfilter/pfilter.html
## JLH from the tutorial
pomp(flu,toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     paramnames=c("Beta","mu_I","mu_R1","rho")) -> flu ## JLH Clay, what does this -> do? It just assigns all this to 'flu' ...but it's at the end instead of the front? <-
## then run nll like above
mle <- flu
coef(mle,c("Beta","mu_I","rho"),transform=TRUE) <- fit$par
coef(mle)

lls <- replicate(n=5,logLik(pfilter(mle,Np=20000)))
ll <- logmeanexp(lls,se=TRUE); ll
simulate(mle,nsim=10,as.data.frame=TRUE,include.data=TRUE) -> sims





## JLH Modification for our code based on the tutorial
pomp(whatgoeshere?nll2?, what goes here?
     and what goes here?,
     paramnames=c("f_pred", "u", "lambda")) -> whatgoeshere?nll2?
## then run nll like above
mle <- whatgoeshere?nll2?
coef(mle,c("f_pred", "u", "lambda"), transform = TRUE) <- fit$parsets ##JLH - is this right? parsets?
coef(mle)

lls <- replicate( n = 5, logLik(pfilter(mle, Np = 20000))) ##JLH why n = 5?
ll <- logmeanexp(lls, se = TRUE); ll
simulate(mle,nsim = 10,as.data.frame = TRUE, include.data = TRUE) -> sims
