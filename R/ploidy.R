rm(list=ls())

simple.prob<-function(d)
  {
  d = d[ order( d$chromosome, decreasing=T ), ]
  
  d$gain.chance = round(d$gain/(d$gain + d$loss), 3)
  d$loss.chance = 1-d$gain.chance
  
  d$overall.chance = round( d$karyotypes/total_karyotypes, 3 )
  return(d[order(d$overall.chance, decreasing=T)])
  }

setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")


datadir = "~/Data/sky-cgh/output/current"
setwd(datadir)

total_karyotypes = 100240
## Ploidy changes ##
pdy = read.table("noleuk-ploidy.txt", header=T, sep="\t")

lpdy = read.table("leuk-ploidy.txt", header=T, sep="\t")

all = pdy
all$gain = all$gain + lpdy$gain
all$loss = all$loss + lpdy$loss
all = simple.prob(all)

# non leukemia
pdy = pdy[ order( pdy$chromosome, decreasing=T ), ]

pdy$gain.chance = round(pdy$gain/(pdy$gain + pdy$loss), 3)
pdy$loss.chance = 1-pdy$gain.chance

pdy$overall.chance = round( pdy$karyotypes/total_karyotypes, 3 )



pdygain = pdy[pdy$'gain',]
pdygain = pdygain[order(-pdygain$count),]

pdyloss = pdy[pdy$class == 'loss',]
pdyloss = pdyloss[order(-pdyloss$count),]

# leukemia
lpdygain = lpdy[lpdy$class == 'gain',]
lpdygain = lpdygain[order(-lpdygain$count),]

lpdyloss = lpdy[lpdy$class == 'loss',]
lpdyloss = lpdyloss[order(-lpdyloss$count),]

# merged
#gains = merge(pdygain, lpdygain, by.x=c('class','chromosome'), by.y=c('class','chromosome'))
#gains = gains[order(-gains$count.x, -gains$count.y), ]

#loss = merge(pdyloss, lpdyloss, by.x=c('class','chromosome'), by.y=c('class','chromosome'))
#loss = loss[order(-loss$count.x, -loss$count.y), ]

#total_gain_prob = gains[,c(1,2)] 
#total_gain_prob$prob = round((gains$count.x+gains$count.y)/total_karyotypes, 4)
#total_gain_prob = total_gain_prob[order(-total_gain_prob$prob),]


#total_loss_prob = loss[,c(1,2)] 
#total_loss_prob$prob = round((loss$count.x+loss$count.y)/total_karyotypes, 4)
#total_loss_prob = total_loss_prob[order(-total_loss_prob$prob),]



par(mfrow=c(2,1))
plot(total_gain_prob$prob, type='h', col='blue', xaxt='n', main="Chr Gain", ylab="Prob. to gain", xlab="chr")
lines(gains$count.x/total_karyotypes, type='p', col='green', xaxt='n')
lines(gains$count.y/total_karyotypes, type='p', col='red')
axis(1, at=1:nrow(gains), label=gains$chromosome)


plot( total_loss_prob$prob, type='h', col='blue', xaxt='n', main="Chr Loss", ylab="Prob. to lose", xlab="chr")
lines(loss$count.x/total_karyotypes, type='p', col='green', xaxt='n')
lines(loss$count.y/total_karyotypes, type='p', col='red')
axis(1, at=1:nrow(loss), label=loss$chromosome)


probs = merge(total_gain_prob[,c(2,3)], total_loss_prob[,c(2,3)], by.x='chromosome', by.y='chromosome', suffixes=c('.gain', '.loss'))
write.table(probs, quote=F, sep="\t", file="aneuploidy-probs.txt")


