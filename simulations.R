#!/usr/bin/Rscript

source("functions.R")


## PRIORS
priors <- c(H1=1e-4,H2=1e-4,H12=1e-5)
priors <- c(H0=1-sum(priors),priors)

## INPUT VARIABLES
## N0 shared controls
## N1 cases disease 1
## N2 cases disease2
## q0 MAF in controls
## p1 p value in trait 1
## p2 p value in trait 2
svals <- c("N0","N1","N2","q0","p1","p2")
args <- getArgs(defaults=list(N0=3000,N1=2000,N2=2000,q0=0.2,p1=5e-8,p2=1e-4))
args <- args[!duplicated(names(args))] ## trick
if(is.null(args$file))
  args$file <- paste(unlist(args),collapse="-")

args[svals] <- lapply(args[svals], function(a) as.numeric(strsplit(as.character(a),",")[[1]]))

## DO CALCS
B <- prod(sapply(args[svals],length))
results <- vector("list",B)
i <- 1
for(N0 in args$N0) {
  for(N1 in args$N1) {
    for(N2 in args$N2) {
      for(q0 in args$q0) {
        for(p1 in args$p1) {
          for(p2 in args$p2) {
            results[[i]] <- runsim(N0=N0,N1=N1,N2=N2,q0=q0,p1=p1,p2=p2)
            i <- i+1
          }
        }
      }
    }
  }
}
results <- do.call("rbind",results)

## SAVE
save(results,file=sprintf("results/sim-%s.RData",args$file))

