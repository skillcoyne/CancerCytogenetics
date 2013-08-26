rm(list=ls())

## This sets up the basic probability of each chromosome either being deleted, or gained. ##
adjust.to.one<-function(p, r=5)
  {
  adjusted = round(p/sum(p), r) 
  
  if (sum(adjusted) > 1)
    adjusted[ which(adjusted == min(adjusted)) ] = adjusted[ which(adjusted == min(adjusted)) ] - (sum(adjusted) - 1)
  
  if (sum(adjusted) < 1)
    adjusted[ which(adjusted == min(adjusted)) ] = adjusted[ which(adjusted == min(adjusted)) ] + (1 - sum(adjusted))
  
  return(adjusted)
  }

simple.prob<-function(d, total)
  {
  d = d[ order( d$chromosome, decreasing=T ), ]
  
  d$gain.prob = d[,'gain']/d[,'karyotypes']
  d$loss.prob = d[,'loss']/d[,'karyotypes']
  
  for (r in 1:nrow(d))
    d[r,c('gain.prob','loss.prob')] = adjust.to.one(d[r,c('gain.prob', 'loss.prob')])
  
  overall =  d$karyotypes/total
  d$prob = adjust.to.one(overall/sum(overall))
  
  return(d[ order(-d$prob), ])
  }

datadir = "~/Data/sky-cgh/output/current"
setwd(datadir)

total_karyotypes = 100240
## Ploidy changes ##
pdy = read.table("noleuk-ploidy.txt", header=T, sep="\t")
lpdy = read.table("leuk-ploidy.txt", header=T, sep="\t")

all = pdy
all$gain = all$gain + lpdy$gain
all$loss = all$loss + lpdy$loss
all$karyotypes = all$karyotypes + lpdy$karyotypes


chrinfo = read.table("~/Data/sky-cgh/genomic_info/chromosome_gene_info_2012.txt", header=T, sep="\t")
chrinfo = chrinfo[ -which(chrinfo[,'Chromosome'] == 'mtDNA'), ]

ploidy_info = merge(all, chrinfo[,c(1,3,5)], by.x=c('chromosome'), by.y=c('Chromosome'))

# length? nope
cor.test(ploidy_info[,'karyotypes'], ploidy_info[, 'Base.pairs'])
cor.test(ploidy_info[,'loss'], ploidy_info[, 'Base.pairs'])

# genes? nope
cor.test(ploidy_info[,'karyotypes'], ploidy_info[, 'Confirmed.proteins']) 
cor.test(ploidy_info[,'loss'], ploidy_info[, 'Confirmed.proteins']) 

# gain *might* be related to length or number of genes
cor.test(ploidy_info[,'gain'], ploidy_info[, 'Base.pairs']) 
cor.test(ploidy_info[,'gain'], ploidy_info[, 'Confirmed.proteins']) 

## if there's no length correlation, is there still a gene correlation? - nope
adj = ploidy_info[,'gain']/ploidy_info[,'Base.pairs']
cor.test(adj, ploidy_info[,'Base.pairs']) 
# is there still a gene correlation?
cor.test(adj, ploidy_info[,'Confirmed.proteins']) 


## So back to a very simple probability
all = simple.prob(all, total_karyotypes)



# non leukemia
pdy = simple.prob(pdy)
# leukemia
lpdy = simple.prob(lpdy)

plot(all[,'overall.prob'], xaxt='n')
axis(1, at=1:nrow(all), label=all$chromosome)

setwd("~/Analysis/Database/cancer")
write.table( all[,c('chromosome','prob', 'gain.prob', 'loss.prob')], quote=F, sep="\t", row.name=F, file="aneuploidy-probs.txt" )

write.table( all[,c('chromosome','prob', 'gain.prob', 'loss.prob')], quote=F, sep="\t", row.name=F )



