rm(list=ls())
setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")


datadir = "~/Data/sky-cgh/output/current"
setwd(datadir)

total_karyotypes = 100240
## Ploidy changes ##
pdy = read.table("noleuk-ploidy.txt", header=T, sep="\t")

lpdy = read.table("leuk-ploidy.txt", header=T, sep="\t")

# non leukemia
pdygain = pdy[pdy$class == 'gain',]
pdygain = pdygain[order(-pdygain$count),]


gains = merge(pdygain, lpdygain, by.x=c('class','chromosome'), by.y=c('class','chromosome'))
gains = gains[order(-gains$count.x, -gains$count.y), ]



pdyloss = pdy[pdy$class == 'loss',]
pdyloss = pdyloss[order(-pdyloss$count),]


loss = merge(pdyloss, lpdyloss, by.x=c('class','chromosome'), by.y=c('class','chromosome'))
loss = loss[order(-loss$count.x, -loss$count.y), ]

par(mfrow=c(2,1))
plot( (gains$count.x+gains$count.y)/total_karyotypes, type='h', col='blue', xaxt='n', main="Chr Gain", ylab="Prob. to gain", xlab="chr")
lines(gains$count.x/total_karyotypes, type='p', col='green', xaxt='n')
lines(gains$count.y/total_karyotypes, type='p', col='red')
axis(1, at=1:nrow(gains), label=gains$chromosome)


plot( (loss$count.x+loss$count.y)/total_karyotypes, type='h', col='blue', xaxt='n', main="Chr Loss", ylab="Prob. to lose", xlab="chr")
lines(loss$count.x/total_karyotypes, type='p', col='green', xaxt='n')
lines(loss$count.y/total_karyotypes, type='p', col='red')
axis(1, at=1:nrow(loss), label=loss$chromosome)



# leukemia
lpdygain = lpdy[lpdy$class == 'gain',]
lpdygain = lpdygain[order(-lpdygain$count),]

lpdyloss = lpdy[lpdy$class == 'loss',]
lpdyloss = lpdyloss[order(-lpdyloss$count),]




