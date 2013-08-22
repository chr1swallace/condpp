library(ggplot2)
library(gridExtra)
sfiles <- list.files("results",pattern="^sim",full=TRUE)
ifiles <- list.files("results",pattern="^indep",full=TRUE)

read.sfile <- function(f) {
  cat(f,"\n")
  (load(f)) # results
  pp <- as.data.frame(results)
  pp$N <- paste(pp$N0,pp$N1,pp$N2,sep="-")
  pp$cp.multi <- pp$pp.multi.H12/(pp$pp.multi.H1+pp$pp.multi.H12)
  pp$cp.indep <- pp$pp.indep.H12/(pp$pp.indep.H1+pp$pp.indep.H12)
  pp$p <- paste(pp$p1,pp$p2,sep="-")
  cbind(as.data.frame(pp[c("N","N0","N1","N2","q0","p","p1","p2","cp.multi","cp.indep")]),file=gsub("results/sim-|.RData","",f))
}
myplot2 <- function(df,f) {
  df <- subset(sdata,file==f)
  keys <- c("F","q0","logp1","logp2")
  v <- apply(df[,keys],2,function(x) length(unique(x)))
  v <- keys[v>1]
  if(length(v)!=1)
    stop("expect exactly one column to vary from",keys,"\n")
  ggplot(df, aes_string(x=v,y="cp",col="analysis",group="analysis")) +
         geom_point() + geom_line() + ggtitle(f) + ylim(0,1) + theme(legend.position="none")
}

## VARY SPECIFIC PARAMETERS
sfiles <- list.files("results",pattern="wtccc|ichip",full=TRUE)
sdata <- do.call("rbind", c(lapply(sfiles, read.sfile)))
multi <- sdata[,-grep("indep",colnames(sdata))]
  indep <- sdata[,-grep("multi",colnames(sdata))]
  multi$analysis <- "multinomial"
  indep$analysis <- "independent"
  colnames(indep) <- sub(".indep","",colnames(indep))
  colnames(multi) <- sub(".multi","",colnames(multi))
  sdata <- rbind(indep,multi)

sdata$logp1 <- -log10(sdata$p1)
sdata$logp2 <- -log10(sdata$p2)
sdata$F <- sdata$N0/(sdata$N0 + sdata$N1+sdata$N2)
sdata$F2 <- sdata$N0/sdata$N2

levels(sdata$file) <-  sub("maf","q0",levels(sdata$file))
levels(sdata$file) <-  sub("ichip","ImmunoChip",levels(sdata$file))
levels(sdata$file) <-  sub("wtccc","WTCCC",levels(sdata$file))

plots <- lapply(levels(sdata$file)[c(1,7) + rep(0:5,each=2)], function(f) myplot2(sdata,f))

theme_set(theme_bw())
do.call("grid.arrange",c(plots,ncol=2))
pdf("scenarios.pdf",height=12,width=10)
do.call("grid.arrange",c(plots,ncol=2))
dev.off()

