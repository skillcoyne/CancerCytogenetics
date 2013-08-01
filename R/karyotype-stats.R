
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
  

total_karyotypes = 100240
datadir = "~/Data/sky-cgh/output"
setwd(datadir)
outdir = "~/Analysis/Database/cancer"

abrkt = read.table("current/abr_per_kt.txt", header=T, sep="\t", row.names=1)
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
bp = as.data.frame(matrix(c(0,0,0,0), ncol=1, nrow=4, dimnames=list(c("1-5","5-10", "10-20", "20-100"), c("count"))))

for(r in 1:nrow(bpkt))
  {
  if (bpkt[r,'bps'] <= 5) bprow = "1-5"
  if (bpkt[r,'bps'] > 5 & bpkt[r, 'bps'] <= 10) bprow = "5-10"
  if (bpkt[r,'bps'] > 10 & bpkt[r, 'bps'] <= 20) bprow = "10-20"
  if (bpkt[r,'bps'] > 20 & bpkt[r, 'bps'] <= 100) bprow = "20-100"
  
  bp[bprow, "count"] = bp[bprow, "count"] + bpkt[r,"kt.count"]
  }
bp$prob = round(bp$count/total_karyotypes, 5)

bp$prob = round(bp$prob/sum(bp$prob), 5)

if (sum(bp$prob) < 1)
  bp$prob[bp$prob == min(bp$prob)] = signif(bp$prob[bp$prob == min(bp$prob)] + (1-sum(bp$prob)), 5)
if (sum(bp$prob) > 1)
  bp$prob[bp$prob == min(bp$prob)] = signif(bp$prob[bp$prob == min(bp$prob)] - (sum(bp$prob)-1), 5)

bp = bp[order(-bp$prob),]
#bp$bps = rownames(bp)


write.table(cbind(rownames(bp), bp$prob), row.name=F,col.names=c('bps','probs'), quote=F, sep="\t", filename=)




chrkt = read.table("current/chr_per_kt.txt", sep="\t")
colnames(chrkt) = c('chrs', 'kt.count')
chrkt = chrkt[order(chrkt[,1]),]

plot((chrkt$kt.count), type='o', col='blue', xlab="Number of chrs", ylab="(Karyotype Count)")

norm.ks(chrkt[,2]) # just about normal
pois.ks(chrkt[,2]) # again, pois fits better

pdy = read.table("current/ploidy_per_kt.txt", sep="\t")
colnames(pdy) = c('ploidy','kt.count')
pdy = pdy[order(pdy$ploidy),]

plot((pdy$kt.count), type='o', col='blue', xlab="Total Ploidy", ylab="Karyotype Count")

norm.ks(pdy[,2])
pois.ks(pdy[,2])


#probability_list = pnorm(pdy[,2],mean(pdy[,2]),sd(pdy[,2]))
#x = ppois(pdy[,2], mean(pdy[,2]))

#x = rpois(1000, mean(pdy[,2]))
#x = rnorm(1000, mean(pdy[,2]), sd(pdy[,2]))


#


