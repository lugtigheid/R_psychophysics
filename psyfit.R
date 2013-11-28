d <- read.table("data.out")
names(d) <- c("stimval", "resp", "n");

# this is a function to fit the cumulative normal
boot.normfit <- function(x,p)  pnorm(x,mean=p[1],sd=p[2])

# this computes the log likelihood
# s: the sample we're fitting 
# p: a two-element vector where p[1] = mu, and p[2] = sigma
# d: the data we're fitting these data on
boot.loglikeli <- function(p,d,s) -sum(log(dbinom(s,d$n, boot.normfit(d$stimval,p))))

# this is the bootstrap function
boot.fun <- function (data, nboot, cuts = c(.025, 0.159, .5, 0.841 ,.975 ) ) {
  
  # resample these data from a binomial distribution
  samples <- matrix(rbinom(nrow(d)*nboot, d$n, d$resp/d$n), ncol=nboot)
  
  # we use the sapply function to replicate the optimisation
  sim <- sapply(seq_len(nboot),function(n) {
    
    # this is the optimisation result for this run
    r <- optim(c(1,1), boot.loglikeli, d=data, s=samples[,n])
    
    # return only the parameters
    return(r$par)
  })
  
  # calculate the percentiles (CI's) for mu and sigma
  mu = quantile(sim[1,], cuts)
  sg = quantile(sim[2,], cuts)
  
  # return the object
  return( list(mu=mu, sg=sg, sim=sim) );
}

# here we actually do the bootstrap
results <- boot.fun(d,999)

# print(data.frame(mu=mu,sigma=sg))

par(mfrow=c(1,2)) 

# let's plot this fucker
plot(d$stimval, d$resp/d$n, ylim=c(0,1), pch=16,
     xlab="Retinal size ratio", ylab="P[ref]")

# create the x-values along which we plot the function
x <- seq(min(d$stimval), max(d$stimval), len=100)

# calculate the y-values
pred.nrm <- pnorm(x, mean=results$mu[3], sd=results$sg[3])
lines(x,pred.nrm)

# plot the scatterplot of simulated points
plot(results$sim[1,],results$sim[2,],pch=16, col='grey', xlab="mu", ylab="sigma")
abline(h = results$sg[3], v = results$mu[3], col = "black", lty=3, lw=1)
# abline(h = results$sg[3], v = results$mu[3], col = "red", lty=1, lw=2)

# done

# ----------------------------------------------------------------

# this also computes the log likelihood but in a different manner - answer is the same
# loglikeli1 <- function(p,d) {
#   pr <- pnorm(d$stimval, mean=p[1], sd=p[2])
#   -sum(d$nYes * log(pr) + d$nNo * log(1-pr))
# }
