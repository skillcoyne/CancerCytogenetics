
rm(list=ls())

setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")

roll<-function(probs)
  {
  rand = runif(1,0,1)
  match = probs[probs >= rand]
  row = match[1]
  return(which(probs == row))
  }

set.probs<-function(probs, df)
  {
  for(i in length(probs):1)
    df[i,'p']<-sum(probs[i:1])
  return(df)
  }

norm.ks<-function(t)
  {
  ks = ks.test(t, pnorm, mean(t), sd(t))
  return(ks)
  }

pois.ks<-function(t)
  {
  ks = ks.test(t, ppois, mean(t))
  return(ks)
  }
  
probs.in.ranges<-function(ranges, counts, decimals=5, total_counts)
  {
  rframe = as.data.frame(matrix( rep(0, length(ranges)), ncol=1, nrow=length(ranges), dimnames=list(sapply( ranges, function(s) paste(s[1], s[2], sep="-") ), c("count")) ))
  for (n in 1:length(ranges))
    {
    x = unlist(ranges[n])
    rframe[ paste(x[1], x[2], sep="-"), 1] =  sum(counts[ which( as.numeric(names(counts)) >= x[1] & as.numeric(names(counts)) <= x[2]  ) ])
    }
  rframe$prob = round(rframe[,'count']/total_counts, decimals)
  rframe$prob = round(rframe[,'prob']/sum(rframe[,'prob']), decimals)
  
  if (sum(rframe$prob) < 1)
    rframe$prob[rframe$prob == min(rframe$prob)] = signif(rframe$prob[rframe$prob == min(rframe$prob)] + (1-sum(rframe$prob)), decimals)
  if (sum(rframe$prob) > 1)
    rframe$prob[rframe$prob == min(rframe$prob)] = signif(rframe$prob[rframe$prob == min(rframe$prob)] - (sum(rframe$prob)-1), decimals)
  
  rframe = rframe[order(-rframe$prob),]
  return(rframe)
  }




datadir = "~/Data/sky-cgh/output"
setwd(datadir)
outdir = "~/Analysis/Database/cancer"

abrkt = read.table("current/abr_per_kt.txt", header=T, sep="\t", row.names=1)
total_karyotypes = nrow(abrkt)

norm.ks(abrkt[,'aneuploidy.count'])
pois.ks(abrkt[,'aneuploidy.count'])

bpkt = read.table("current/bp_per_kt.txt", sep="\t")
colnames(bpkt) = c('bps','kt.count') 
# Don't much care about those that have no breakpoints
bpkt = bpkt[bpkt[,'bps'] > 0, ]
bpkt = bpkt[order(bpkt[,1]),]

plot((bpkt$kt.count), type='o', col='blue', xlab="Number of bps", ylab="(Karyotype Count)")

# neither fits very well - so use a simple probability
norm.ks(bpkt[,2])  
pois.ks(bpkt[,2])  

bpkt$prob = round( bpkt[,2]/total_karyotypes, 5)
bpkt = bpkt[ order(bpkt$prob), ]

# group them for more continuous probabilities
bp_tbl = bpkt[,2]
names(bp_tbl) = bpkt[,1]
bp_probs = probs.in.ranges( list(c(1,5), c(6,10), c(11,20), c(21,100)), bp_tbl, 5, total_karyotypes )

## Aberrations and aneuploidy per karyotype
abrs = read.table("current/abr_per_kt.txt", header=T)

pdy_probs = probs.in.ranges(list(c(0,3), c(4,10), c(11,20), c(21,35)), table(abrs[,2]), 5, total_karyotypes)
abr_probs = probs.in.ranges(list(c(0,2), c(3,7), c(8,14), c(15,20), c(21,55)), table(abrs[,3]), 5, total_karyotypes)

bp_probs$type = "breakpoint"
pdy_probs$type = "aneuploidy" 
abr_probs$type = "aberration"

#write("## General karyotype probabilities given for counts in ranges", file=paste(outdir, "karyotype-probs.txt", sep="/"), app=F)
#for (pf in list(bp_probs, pdy_probs, abr_probs))
#  write.table(cbind(rownames(pf), pf[,c('prob', 'type')]), row.name=F,col.names=F, quote=F, sep="\t", app=T, file=paste(outdir, "karyotype-probs.txt", sep="/"))


