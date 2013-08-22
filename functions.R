##' commandArgs parsing
##' 
##' return a named list of command line arguments
##'
##' Usage:
##' call the R script thus
##'   ./myfile.R --args --myarg=something
##' or
##'   R CMD BATCH --args --myarg=something myfile.R
##'
##' Then in R do
##'   myargs <- getArgs()
##' and myargs will be a named list
##' > str(myargs)
##' List of 2
##' $ file : chr "myfile.R"
##' $ myarg: chr "something"
##'
##' @title getArgs
##' @param verbose print verbage to screen 
##' @param defaults a named list of defaults, optional
##' @param numeric names of arguments that should be converted from character to numeric, optional
##' @return a named list
##' @author Chris Wallace
getArgs <- function(verbose=FALSE, defaults=NULL, numeric=NULL) {
  myargs <- gsub("^--","",commandArgs(TRUE))
  setopts <- !grepl("=",myargs)
  if(any(setopts))
    myargs[setopts] <- paste(myargs[setopts],"=notset",sep="")
  myargs.list <- strsplit(myargs,"=")
  myargs <- lapply(myargs.list,"[[",2 )
  names(myargs) <- lapply(myargs.list, "[[", 1)
  
  ## logicals
  if(any(setopts))
    myargs[setopts] <- TRUE
  
  ## defaults
  if(!is.null(defaults)) {
    defs.needed <- setdiff(names(defaults), names(myargs))
    if(length(defs.needed)) {
      myargs[ defs.needed ] <- defaults[ defs.needed ]
    }
  }
  
  ## numerics
  if(!is.null(numeric)) {
    numeric <- intersect(numeric, names(myargs))
    if(length(numeric))
      myargs[numeric] <- lapply(myargs[numeric], as.numeric)
  }
  
  ## verbage
  if(verbose) {
    cat("read",length(myargs),"named args:\n")
    print(myargs)
  }
  myargs
}

library(mlogitBMA)
calc.q1 <- function(q0,p,N0,N1) {
  X2 <- qchisq(p,df=1,lower=FALSE)
  ## find q1
  f <- function(q1) {
    q <- (N0*q0 + N1*q1)/(N0 + N1)
    qq <- q*(1-q)
    abs(X2 - (2*N0*(q0-q)^2/qq + 2*N1*(q1-q)^2/qq))
  }
  optimise(f,interval=c(q0,0.5))$minimum
}
Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}
approx.bf <- function(p,f,type, N, s, suffix) {
  if(type=="quant") {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  } else {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF # I want ln scale to compare in log natural scale with LR diff
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ret <- data.frame(V,z,r,lABF)
  colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)
}
logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}
calc.pp.multi <- function(df1) {
  f1 <- as.formula("Y ~ 1 | X")
  binmod<-mlogit2logit(f1,data=df1,choices=0:2,base.choice = 1)$data
  
  index<-grep("z",colnames(binmod))
  binX<-binmod[,index]
                                        #remove z_1 (we do not want it in our model matrix)
  binX<-binX[,-1]
                                        #extract the new reponse
  binY<-binmod[,"Y.star"]
  models<-cbind(c(0,1,0,1),c(0,0,1,1),rep(1,4))
                                        #run glib
  mods1 <- glib(x=binX, y=binY, error="binomial", link="logit",models=models)
                                        #log the bf with a flat prior
  B10=mods1$bf$twologB10[,2]/2
  B10 <- B10 - B10[1]
  
  pp <- log(priors) + B10 - log(priors[1])
  exp(pp - logsum(pp))
}
calc.pp.indep <- function(df1,q0,N0,N1,N2) {
  p1 <- anova(glm(Y~X,data=df1,subset=Y!=2,family="binomial"),test="Chisq")$P[2]
  df2 <- subset(df1,Y!=1)
  df2$Y <- df2$Y/2
  p2 <- anova(glm(Y~X,data=df2,family="binomial"),test="Chisq")$P[2]
  
  bf1 <- approx.bf(p1,q0,"cc",N1+N0,N1/(N1+N0),"1")[1,4]
  bf2 <- approx.bf(p2,q0,"cc",N2+N0,N2/(N2+N0),"2")[1,4]
  B10 <- c(1,bf1,bf2,bf1+bf2)
  pp <- log(priors) + B10 - log(priors[1])
  pp.indep <- c(obs.p1=p1,obs.p2=p2,exp(pp - logsum(pp)))
}
cprob <- function(pp) {
  pp["p12"]/sum(pp[c("p1","p12")])
}

gfq <- function(q) {
  c((1-q)^2,2*q*(1-q),q^2)
}
indep <- function(N0,N1,N2,q0,p1,p2) {  
  bf1 <- approx.bf(p1, q0, "cc", N0+N1, N1/(N0+N1), "1")[1,4]
  bf2 <- approx.bf(p2, q0, "cc", N0+N2, N2/(N0+N2), "2")[1,4]
  B10 <- c(1,bf1,bf2,bf1+bf2)
  pp <- log(priors) + B10 - log(priors[1])
  exp(pp - logsum(pp))
}

runsim <- function(N0,N1,N2,q0,p1,p2,B=100) {
  
  ## DERIVED VARIABLES
  q1 <- calc.q1(q0=q0,p=p1,N0=N0,N1=N1)
  q2 <- calc.q1(q0=q0,p=p2,N0=N0,N1=N2)
  
    ## GENERATE DATA
    df1 <- rbind(data.frame(Y=0,X=rep(0:2,times=round(N0*gfq(q0)))),
                 data.frame(Y=1,X=rep(0:2,times=round(N1*gfq(q1)))),
                 data.frame(Y=2,X=rep(0:2,times=round(N2*gfq(q2)))))
    
    ## POSTERIOR PROBS
    pp.multi <- calc.pp.multi(df1)
    pp.indep <- calc.pp.indep(df1,q0,N0,N1,N2)
    
    ## RETURN
    c(N0=N0,N1=N1,N2=N2,q0=q0,p1=p1,p2=p2,
      pp.multi=pp.multi,pp.indep=pp.indep)
  
}
