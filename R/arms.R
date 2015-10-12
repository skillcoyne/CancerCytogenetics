# Script looks at breakpoints within arms. Centromeres are known to be unstable so looking at chromosomal instability 
# in just the arms could provide a more accurate or interesting view on overall instability.

rm(list=ls())

setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")

datadir = "~/Data/sky-cgh/output"
outdir = "~/Analysis/Database/cancer"
setwd(datadir)

leuk = T

# Load files
bp = load.breakpoints("current/noleuk-breakpoints.txt")
if (leuk)
  bp = load.breakpoints(c("current/noleuk-breakpoints.txt", "current/leuk-breakpoints.txt"))


nrow(bp)
chrinfo = loadChromosomeInfo("../genomic_info/chromosome_gene_info_2012.txt")

bandinfo = read.table("../genomic_info/band_genes.txt",  header=T, sep="\t")
bandinfo$length = bandinfo$end - bandinfo$start
## Most of these are driven by breakages at the centromeres.  What happens when we look just at the arms
# how the hell did I decide this?


cmrows = grep("(p|q)(11)",bp$band)
arms = grep("(p|q)(11)",bp$band,invert=T)


## CENTROMERES ## ignore X,Y
centromeres = bp[ cmrows, ]
centinfo = bandinfo[ grep("(p|q)(10|11|12)",bandinfo$band),]

# drop centromeric bands
bp = bp[ -cmrows, ]
bandinfo = bandinfo[ grep("(p|q)(10|11|12)",bandinfo$band,invert=T),]

## Centromeres are removed, so now how long are the arms
# ignoring X and Y
chrs = c(1:22)#, 'X', 'Y')
arm_counts=as.data.frame(matrix(nrow=length(chrs),ncol=3,dimnames=list(chrs,c('count','length','genes'))))
for(i in chrs)
  arm_counts[i,] = c(sum(bp[bp$chr == i,'total.karyotypes']), sum(bandinfo[bandinfo$chr == i, 'length']), sum(bandinfo[bandinfo$chr == i, 'gene.count']))

cent_counts=as.data.frame(matrix(nrow=length(chrs),ncol=3,dimnames=list(chrs,c('count','length','genes'))))
for(i in chrs)
  cent_counts[i,] = c(sum(centromeres[centromeres$chr == i,'total.karyotypes']), sum(centinfo[centinfo$chr == i, 'length']), sum(centinfo[centinfo$chr == i, 'gene.count']))



## does it matter how far the band is from the centromere?
distance_cor = matrix(nrow=length(chromosomes), ncol=2, dimnames=list(chromosomes, c('p','q')))
for (chr in chrs)
  {
  scores = bp[bp$chr == chr,]
  scores = scores[order(scores$band),]
  scores = na.omit(scores)
  # get the arms in the correct order p -> q
  parm = grep("p", scores$band )

  q = scores[ grep("q", scores$band),]
  q$c.distance = seq(1:nrow(q))
  #q$c.distance = sample(c(1:100), nrow(q))
  
  p = scores[grep("p", scores$band ),]
  if (nrow(p) > 0)
    {
    p$c.distance = seq(1:nrow(p))
    #p$c.distance = sample(c(1:100), nrow(p))
    p = p[order(p$band, decreasing=T),]
    }
  scores = rbind(p,q)
  scores
  tryCatch({
    test = cor.test(q$total.karyotypes, q$c.distance, method="pearson")
    distance_cor[chr, 'q'] = round(test$estimate, 4)
    
    test = cor.test(p$total.karyotypes, p$c.distance, method="pearson")
    distance_cor[chr, 'p'] = round(test$estimate, 4)
  }, error = function(e) {
    print(e)
  })
  
  }  

# CENT  chr1 and 22 lower the correlation between length and count
inc=c(2:21)
inc=c(1:22)

plot(cent_counts$count, cent_counts$length, type='n')
text(cent_counts$count, cent_counts$length, labels=rownames(cent_counts)[inc])
# 22 appears to be an outlier, but doesn't seem to drive the correlation
cor.test(cent_counts$count, cent_counts$length) # low correlation but it is there

cor.test(cent_counts$count, cent_counts$length) # low correlation but it is there
cor.test(cent_counts$count, cent_counts$genes) # low correlation but it is there

# without 1 and 22, 17 is the outlier
length_adj=0.9
cent_adj_scores = cent_counts$count/(cent_counts$length^length_adj)
cor.test(cent_adj_scores[inc], cent_counts$length[inc])
plot(cent_adj_scores[inc], type='n')
text(cent_adj_scores[inc], labels=rownames(cent_counts)[inc])


# ARMS highly correlated still, so perhaps not driven by the centromere?
cor.test(arm_counts$count,arm_counts$length) 
cor.test(arm_counts$count,arm_counts$genes) 
length_adj=0.7
arm_adj_scores = arm_counts[,'count']/(arm_counts[,'length']^length_adj)
cor.test(arm_adj_scores, arm_counts[,'length'])

plot(arm_adj_scores, arm_counts[,'length'], type="n", main="Chromosome Arms vs Calculated base pairs", xlab="Breakpoint counts", ylab="Base pair length for arms")
text(arm_adj_scores, arm_counts[,'length'], labels=rownames(arm_counts))


## So centromeres don't have protein coding genes so if we compare (the way we did with whole chromosomes) to information encoded in total on just the arms...
#and now we can (sorta) say that chromosome instability related directly to the amount of information encoded on that chromosome
cor.test(arm_adj_scores, arms[,'genes'])
cor.test(cent_adj_scores, arms[,'genes']) # no correlation at all for centromeres


# Again look at distribution - it's normal, not sure that I'd expect that to change really
# Can't be compared to the one from bp-analysis though as that score is calculated based on information content rather than length.  
ks.test(arm_adj_scores, pnorm, mean(arm_adj_scores), sd(arm_adj_scores))
arm_probability_list=pnorm(arm_adj_scores,mean(arm_adj_scores),sd(arm_adj_scores))

plot(function(x) dnorm(x,mean(arm_adj_scores),sd(arm_adj_scores)), min(arm_adj_scores)-sd(arm_adj_scores)*1,max(arm_adj_scores)+sd(arm_adj_scores)*1,
     ylab="Density",xlab="Chromosome Arms Instability Score")
xpos=vector("numeric",2)
ypos=vector("numeric",2)
ypos[1]=0
density_pos=dnorm(arm_adj_scores,mean(arm_adj_scores),sd(arm_adj_scores))
for(i in 1:length(arm_adj_scores))
	{
	xpos[1]= arm_adj_scores[i]
	xpos[2]= arm_adj_scores[i]
	ypos[2]=density_pos[i]
	lines(xpos,ypos)
	text(xpos[2],ypos[2],labels=i, pos=3)
	}

instability_score=arm_probability_list


plot(instability_score,arm_counts$genes/(arm_counts$length^length_adj), main="Chromosome Instability, Arms", xlab="Chromosome instability score",
     ylab="Length normalized protein count", type="n")
text(instability_score,arm_counts$genes/(arm_counts$length^length_adj),labels=rownames(arm_counts),col='blue')

plot(instability_score, chrinfo[c,"Confirmed.proteins"]/(chrinfo[c,"Base.pairs"]^length_adjust), 
     main="Chromosome Instability",
     sub=corr_str,xlab="Chromosome instability score",ylab="Length normalised protein count", type="n")
text(instability_score, chrinfo[c,"Confirmed.proteins"]/(chrinfo[c,"Base.pairs"]^length_adjust),labels=names(instability_score), col='blue')


hc = hclust(dist(instability_score))
groups = cutree(hc, k=3)
plot(hc, main="Chromosome Instability", cex=2)
rect.hclust(hc, k=3, border=c("red", "blue", "green"))

plot(hc, col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071",
     col.axis = "#F38630", lwd = 3, lty = 3, sub = '', hang = -1, axes = FALSE)
rect.hclust(hc, k=3, border=c("red", "blue", "green"))

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")

A2Rplot(hc, k = 3, boxes = FALSE, col.up = "gray50", col.down = c("red", "blue", "purple"), main="Chromosome Instability")




#filename = "~/Analysis/Database/cancer/arm_chr_instability_prob.txt"
#write("# Normal distribution, probability score per chromosome. Each score is independent of the other chromosomes", file=filename, app=F)
#write.table( arm_probability_list/sum(arm_probability_list), quote=F, col.names=F, sep="\t", app=T, file=filename)

ks.test(cent_adj_scores, pnorm, mean(cent_adj_scores), sd(cent_adj_scores)) # on the other hand, centromere bps are not normally distributed

write.table(bp[1:10,c('name','total.karyotypes')], row.names=F, quote=F,sep="\t")



