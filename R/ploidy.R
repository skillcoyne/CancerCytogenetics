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

simple.prob<-function(d)
  {
  d = d[ order( d$chromosome, decreasing=T ), ]
  
  d$gain.prob = signif( adjust.to.one(d[,'gain']/d[,'karyotypes']), 5)
  d$loss.prob = signif( adjust.to.one(d[,'loss']/d[,'karyotypes']), 5)
  
  overall =  d$karyotypes/total_karyotypes
  d$overall.prob = adjust.to.one(overall/sum(overall))
  
  return(d[ order(-d$overall.prob), ])
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
all = simple.prob(all)



# non leukemia
pdy = simple.prob(pdy)
# leukemia
lpdy = simple.prob(lpdy)

plot(all[,'overall.prob'], xaxt='n')
axis(1, at=1:nrow(all), label=all$chromosome)

setwd("~/Analysis/Database/cancer")
write.table( all[,c('chromosome','gain.prob', 'loss.prob')], file="aneuploidy-probs.txt", quote=F, sep="\t", row.name=F )

write.table( all[,c('chromosome','gain.prob', 'loss.prob')], quote=F, sep="\t", row.name=F )



