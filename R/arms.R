# Script looks at breakpoints within arms. Centromeres are known to be unstable so looking at chromosomal instability 
# in just the arms could provide a more accurate or interesting view on overall instability.

rm(list=ls())

setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")

datadir = "~/Data/sky-cgh/output"
outdir = "~/Analysis/Database/cancer"
setwd(datadir)

total_karyotypes = 100240
leuk = FALSE

# Load files
bp = read.table("current/noleuk-breakpoints.txt", sep="\t", header=T)

if (leuk)
  {
  bp2 = read.table("current/leuk-breakpoints.txt", sep="\t", header=T)
  cols = c('chr','band','start','end', 'total.karyotypes')
  merged = merge(bp[,cols], bp2[,cols], by=cols[1:4])

  merged$total.karyotypes = merged[,'total.karyotypes.x'] + merged[,'total.karyotypes.y']
  merged$total.karyotypes.x = NULL
  merged$total.karyotypes.y = NULL
  bp = merged
  }

nrow(bp)
chrinfo = loadChromosomeInfo("../genomic_info/chromosome_gene_info_2012.txt")

bp = bp[bp$total.karyotypes > 5,]

bandinfo = read.table(paste(datadir,"band_genes.txt", sep="/"), header=T, sep="\t")

## Most of these are driven by breakages at the centromeres.  What happens when we look just at the arms
## CENTROMERES ## ignore X,Y
cmrows = grep("(11|12)", bp$band)
centromeres = bp[ cmrows, ]

# drop centromeric bands
bp = bp[ -cmrows, ]

cmrows = grep("(11|12)", bandinfo$band)
bandinfo = bandinfo[ -cmrows, ]
bandinfo$length = bandinfo$end - bandinfo$start

## Centromeres are removed, so now how long are the arms
chromosomes = c(1:22)
cols = c('bps', 'genes', 'length')
arms = matrix(nrow=length(chromosomes), ncol=length(cols), dimnames=list(chromosomes, cols))
distance_cor = matrix(nrow=length(chromosomes), ncol=2, dimnames=list(chromosomes, c('p','q')))
for (chr in chromosomes)
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
  
  arms[chr, 'bps'] = sum(scores$total.karyotypes)
  
  # ---- #
  bands = bandinfo[bandinfo$chr == chr,]
  bands = bands[order(bands$band),]
  qarm = grep("q", bands$band)
  arms[chr, 'length'] = sum(bands[ qarm, 'length' ]) + sum(bands[ -qarm, 'length' ])
  ## just a sanity check
  if (arms[chr, 'length'] > chrinfo[chr, 'Base.pairs'])
    stop( paste(chr, "arm lengths incorrect") )
  arms[chr, 'genes'] = sum(bands$gene.count)
  }  

## aaaand, this suggests no correlation? So would this mean most of the instability is due to the centromeres and length so not very good example
arm_adj_scores = arms[,'bps']/arms[,'length']
cor.test(arm_adj_scores, arms[,'length'])

plot(arm_adj_scores, arms[,'length'], type="n", main="Chromosome Arms vs Calculated base pairs", xlab="Breakpoint counts", ylab="Base pair counts for arms")
text(arm_adj_scores, arms[,'length'], labels=names(arm_adj_scores))

adj_factor = 0.7
## So centromeres don't have protein coding genes so if we compare (the way we did with whole chromosomes) to information encoded in total on just the arms...
arm_adj_scores = arms[,'bps']/(arms[,'length']^adj_factor)
cor.test(arm_adj_scores, arms[,'length'])
#and now we can (sorta) say that chromosome instability related directly to the amount of information encoded on that chromosome
cor.test(arm_adj_scores, arms[,'genes'])


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
cor.test(instability_score,arms[,'genes'])

ins_fit = hclust(dist(instability_score))
groups = cutree(ins_fit, k=3)
plot(ins_fit, main="Chromosome Arm Instability")
rect.hclust(ins_fit, k=3, border=c("red", "blue", "green"))

#filename = "~/Analysis/Database/cancer/arm_chr_instability_prob.txt"
#write("# Normal distribution, probability score per chromosome. Each score is independent of the other chromosomes", file=filename, app=F)
#write.table( arm_probability_list/sum(arm_probability_list), quote=F, col.names=F, sep="\t", app=T, file=filename)

