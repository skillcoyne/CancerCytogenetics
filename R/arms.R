# Script looks at breakpoints within arms. Centromeres are known to be unstable so looking at chromosomal instability 
# in just the arms could provide a more accurate or interesting view on overall instability.
resetwd()
source("R/lib/load_files.R")
source("R/lib/wd.R")

setDataDirectory(date = '09012013')

# Load files
bp = loadBreakpoints("breakpoints.txt")
nrow(bp)
chrinfo = loadChromosomeInfo("../../genomic_info/chromosome_gene_info_2012.txt")

# This just gets me a list of bands per chromosome, probably useful, ignoring X/Y
knownbands = read.table("../../genomic_info/bands_by_chr.txt", header=T, sep="\t")
knownbands = dropSexChr(knownbands)


## sample the most frequent leukemia/lymphoma cases 
sample_leukemia_bps = FALSE
if (sample_leukemia_bps)
  {
  leuks = c('Acute myeloid leukemia', 'Acute lymphoblastic leukemia', "Non-hodgkin's lymphoma", 'Chronic myelogenous leukemia', 'Chronic lymphocytic leukemia')
  bp = sampleCancers(bp, "Cancer", leuks)
  }

bpfreq = table(bp$Breakpoint)



## Most of these are driven by breakages at the centromeres.  What happens when we look just at the arms
## CENTROMERES ## ignore X,Y
centromeres = vector(mode="character")
for (i in 1:22)
 	{
	b = 11
  centromeres = c(centromeres, paste(i, "q", b, sep=""))
  centromeres = c(centromeres, paste(i, "p", b, sep=""))
	}

# Centromere
bpcent=vector("numeric", 22)
#bpcent=vector("numeric",length(centromeres))
names(bpcent)=1:22
for (i in 1:22)  
	{
  p = paste(i,"p", 11, sep="")
  q = paste(i,"q", 11, sep="")
  bpcent[i] = sum(bp$Breakpoint == p) + sum(bp$Breakpoint == q)
	}

`%nin%`=Negate(`%in%`) # cute way to create my own 'not' operator
# drop centromeric bands
knownbands = knownbands[ knownbands$Band %nin% centromeres, ]

# Sum the breakpoints seen per band
for (band in knownbands$Band)
  {
  knownbands[knownbands$Band == band, 'Counts'] = nrow(bp[bp$Breakpoint == band,])
  }
sorted = knownbands[order(-knownbands$Count), ]
write.table(sorted, quote=F, sep="\t", file="arm-bp-freq.txt")


# OK so ultimately all that does is let me pick out the bands that aren't centromeric.  I can't find exact numbers on this.  
# So far what I've read is the estimates of centromere size in base pairs is about 1% 
chrinfo$Base.pairs.arms = floor(chrinfo$Base.pairs-(chrinfo$Base.pairs*.01))
# so then...
arm_counts=vector("numeric", 22)
names(arm_counts)=1:22
for(i in 1:22)
	{
	arm_counts[i]=sum(knownbands[knownbands$Chr == i, 'Counts'])
	}

## aaaand, this suggests no correlation? So would this mean most of the instability is due to the centromeres and length so not very good example
arm_adj_scores=arm_counts/(chrinfo[1:22,"Base.pairs.arms"])
cor.test(arm_adj_scores[1:22],chrinfo[1:22,"Base.pairs.arms"])
dev.new()
plot(arm_adj_scores[1:22],chrinfo[1:22,"Base.pairs.arms"], type="n", main="Chromosome Arms vs Calculated base pairs", xlab="Breakpoint counts", ylab="Base pair counts for arms")
text(arm_adj_scores[1:22],chrinfo[1:22,"Base.pairs.arms"],labels=names(arm_adj_scores))


## So centromeres don't have protein coding genes so if we compare (the way we did with whole chromosomes) to information encoded in total on just the arms...
adjusted_scores=arm_counts/(chrinfo[1:22,"Base.pairs.arms"]^(.7))
cor.test(adjusted_scores[1:22],chrinfo[1:22,"Base.pairs.arms"])
#and now we can (sorta) say that chromosome instability related directly to the amount of information encoded on that chromosome
cor.test(adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"])


# Again look at distribution - it's normal, not sure that I'd expect that to change really
# Can't be compared to the one from bp-analysis though as that score is calculated based on information content rather than length.  NEED TO FIND THAT
ks.test(arm_adj_scores,pnorm,mean(arm_adj_scores),sd(arm_adj_scores))
arm_probability_list=pnorm(arm_adj_scores,mean(arm_adj_scores),sd(arm_adj_scores))
dev.new()
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
cor.test(instability_score,chrinfo[1:22,"Confirmed.proteins"])

ins_fit = hclust(dist(instability_score))
groups = cutree(ins_fit, k=3)
plot(ins_fit, main="Chromosome Arm Instability")
rect.hclust(ins_fit, k=3, border=c("red", "blue", "green"))



