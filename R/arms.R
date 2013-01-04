# Script looks at breakpoints within arms. Centromeres are known to be unstable so looking at chromosomal instability 
# in just the arms could provide a more accurate or interesting view on overall instability.
#


# This just gets me a list of bands per chromosome, probably useful
knownbands=read.table("/Users/sarah.killcoyne/workspace/KaryotypeAnalysis/resources/bands_by_chr.txt", header=T, sep="\t")
## Ignoring X/Y
knownbands=knownbands[which(knownbands$Chr!="X" & knownbands$Chr!="Y"),]

# Actual breakpoints
basedir="/Users/sarah.killcoyne/Data/sky-cgh"
bp = read.table(paste(basedir, "/output/26112012/","breakpoints.txt", sep=""), sep="\t", comment="#", header=T)
## Lets ignore subbands (11.1) and just group them by the major band designation (11)
bp$Breakpoint = sub("\\.[0-9]+", "", bp$Breakpoint)
# Not using sex chromosomes generally, 
bp=bp[which(bp$Chr!="X" & bp$Chr!="Y"),]
# but here they are if we need them
bp_sex=bp[which(bp$Chr=="X" | bp$Chr=="Y" )]

# 
chrinfo = read.table(paste(basedir, "/chromosome_gene_info_2012.txt", sep=""), sep="\t", row.names=1, header=T)
# don't need the mtDNA row
chrinfo = chrinfo[ -(nrow(chrinfo)), ]

# ignoring X and Y
chromosome_counts=vector("numeric",22)
names(chromosome_counts)=1:22
for(i in 1:22)
	{
	chromosome_counts[i]=sum(bp$Chr==i)
	}

#adjust the length in a non-linear manner
adjusted_scores=chromosome_counts/(chrinfo[1:22,"Base.pairs"]^(.7))
#chromosome instability score (not as a function of length)
cor.test(adjusted_scores[1:22],chrinfo[1:22,"Base.pairs"])  
# instability is (sorta) related to the proteins encoded by that chromosome (rather than length)
cor.test(adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"])
#plot(adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"], type="n")
#text(adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"],labels=names(adjusted_scores))
probability_list=pnorm(adjusted_scores,mean(adjusted_scores),sd(adjusted_scores))

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
bpcent=vector("numeric",length(centromeres))
names(bpcent)=centromeres
for (i in centromeres)  
	{
	bpcent[i] = sum(bp$Breakpoint == i)
	}
	

`%nin%`=Negate(`%in%`) # cute way to create my own 'not' operator
# Arm
knownbands$Counts=0
bparm = bp[which(bp$Breakpoint %nin% centromeres),]  
bparm_f=sort(table(bparm$Breakpoint))
for (name in names(bparm_f))
	{
	knownbands[knownbands$Band == name, 'Counts'] = sum(bparm$Breakpoint == name)
	}


# OK so ultimately all that does is let me pick out the bands that aren't centromeric.  I can't find exact numbers on this.  So far what I've read is the estimates of centromere size in base pairs is about 1%
chrinfo$Base.pairs.arms = floor(chrinfo$Base.pairs-(chrinfo$Base.pairs*.01))
# so then...
arm_counts=vector("numeric", 22)
names(arm_counts)=1:22
for(i in 1:22)
	{
	arm_counts[i]=sum(knownbands[knownbands$Chr == i, 'Counts'])
	}
## aaaand, this suggests no correlation? So would this mean most of the instability is due to the centromeres?
arm_adj_scores=arm_counts/(chrinfo[1:22,"Base.pairs.arms"])
cor.test(arm_adj_scores[1:22],chrinfo[1:22,"Base.pairs.arms"])
dev.new()
plot(arm_adj_scores[1:22],chrinfo[1:22,"Base.pairs.arms"], type="n", main="Chromosome Arms vs Calculated base pairs", xlab="Breakpoint counts", ylab="Base pair counts for arms")
text(arm_adj_scores[1:22],chrinfo[1:22,"Base.pairs.arms"],labels=names(arm_adj_scores))


	
### TODO should really find another file with information regarding what is encoded for by band 

# Again look at distribution - it's normal, not sure that I'd expect that to change really
# Can't be compared to the one from bp-analysis though as that score is calculated based on information content rather than length.  NEED TO FIND THAT
ks.test(arm_adj_scores,pnorm,mean(arm_adj_scores),sd(arm_adj_scores))
arm_probability_list=pnorm(arm_adj_scores,mean(arm_adj_scores),sd(arm_adj_scores))
dev.new()
plot(function(x) dnorm(x,mean(arm_adj_scores),sd(arm_adj_scores)), min(arm_adj_scores)-sd(arm_adj_scores)*1,max(arm_adj_scores)+sd(arm_adj_scores)*1,ylab="Density",xlab="Instability Score")
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
	text(xpos[2],ypos[2],labels=i)
	}







