#!/usr/bin/Rscript
source("functions.R")

## NOVEL REGIONS
regions <- read.table("novel_regions.csv",header=TRUE,sep=",",as.is=TRUE)


## PRIORS
priors <- c(H1=1e-4,H2=1e-4,H12=1e-5)
priors <- c(H0=1-sum(priors),priors)

## INPUT VARIABLES
N<-list( #case,control
    'Celiac'=c(12041,12228), ## http://www.immunobase.org/page/Overview/display/study_id/GDXHsS00009
    'PBC'=c(2861,8514), ## http://www.immunobase.org/page/Overview/display/study_id/GDXHsS00020
    'RA'=c(11475,15870), ## http://www.immunobase.org/page/Overview/display/study_id/GDXHsS00019
    'T1D'=c(6670,9416), 
    'JIA'=c(2816,13056) ## curation not verified http://www.immunobase.org/page/Overview/display/study_id/GDXHsS00023
)

results <- vector("list",nrow(regions))
for(i in 1:nrow(regions)) {
  d <- regions[i,"disease"]
  N0 <- min(N[[d]][2],N[["T1D"]][2])
  N1 <- N[[d]][1]
  N2 <- N[["T1D"]][1]
  p1 <- regions[i,"p.val"]
  p2 <- regions[i,"t1d_cooper"]
  q0 <- regions[i,"MAF"]  
  results[[i]] <- runsim(N0=N0,N1=N1,N2=N2,q0=q0,p1=p1,p2=p2)
}
results <- do.call("rbind",results)
results <- as.data.frame(results)
results$cp.multi <- results$pp.multi.H12/(results$pp.multi.H1 + results$pp.multi.H12)

results <- cbind(regions,results)
results[,c("index.snp","disease","q0","p1","p2","N0","N1","N2","cp.multi")]
