rm(list=ls())
setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")

datadir = "~/Data/sky-cgh/output/29052013"
setwd(datadir)

`%nin%`=Negate(`%in%`) 
chrinfo = loadChromosomeInfo()  

bp = read.table("breakpoints.txt", sep="\t", header=T)


plots = TRUE

# drop anyting that shows up in < 5 samples
#bp = bp[ bp$total.samples >= 5, ]
nrow(bp)

total_major_bps = nrow(bp)
total_count_bps = sum(bp$total.samples)

# turn it into a table
bp_freq = vector(mode="numeric", length=nrow(bp))
names(bp_freq) = paste(bp$chr, bp$band, sep="")
bp_freq[1:length(bp_freq)] = bp$total.samples

# ignoring X and Y
chromosome_counts=vector("numeric",24)

names(chromosome_counts)= c(1:22, 'X', 'Y')
for(i in c(1:22, 'X', 'Y'))
	{
  curr = bp[bp$chr == i,]
  chromosome_counts[i] = sum(curr$patients+curr$cell.lines)
	}
c=1:24

# Obviously correlates to length
cor.test(chromosome_counts[c],chrinfo[c,"Base.pairs"])

#-- so we adjust the length in a non-linear manner (with a bit of mucking arround found ^.7 gave nearly no correlation) --#
adjusted_scores=chromosome_counts/(chrinfo[c,"Base.pairs"]^(.7))
#and now we don;t have a correlation with length- this is our chromosome instability score
cor.test(adjusted_scores[c],chrinfo[c,"Base.pairs"])
#and now we can (sorta) say that chromosome instability related directly to the amount of information encoded on that chromosome
cor.test(adjusted_scores[c],chrinfo[c,"Confirmed.proteins"])


if (plots)
  {
  #although this is being drived by chromosome "1", which means we will have to think about how to normalise the score a bit more
  plot(adjusted_scores[c],chrinfo[c,"Confirmed.proteins"], type="n", xlab="Chromosome Counts adjusted for length", ylab="Protein Counts")
  text(adjusted_scores[c],chrinfo[c,"Confirmed.proteins"],labels=names(adjusted_scores))
  }

#anyway now as we wanna be clever we can map it to a probability distribution -first we check to see if adjusted_scores are normal
ks.test(adjusted_scores,pnorm,mean(adjusted_scores),sd(adjusted_scores))

## -- Probabilities are probably the most useful -- ##
#they are normal so we generate the proabilities of these scores occuring if they follow the normal dist
probability_list=pnorm(adjusted_scores,mean(adjusted_scores),sd(adjusted_scores))

if (plots)
  {
  #and now a (vaguely) pretty pic just for ewe
  plot(function(x) dnorm(x,mean(adjusted_scores),sd(adjusted_scores)), min(adjusted_scores)-sd(adjusted_scores)*1,max(adjusted_scores)+sd(adjusted_scores)*1,
     ylab="Density",xlab="Chromosome Instability Score")
  xpos=vector("numeric",2)
  ypos=vector("numeric",2)
  ypos[1]=0
  density_pos=dnorm(adjusted_scores,mean(adjusted_scores),sd(adjusted_scores))
  for(i in 1:length(adjusted_scores))
	  {
	  xpos[1]=adjusted_scores[i]
	  xpos[2]=adjusted_scores[i]
	  ypos[2]=density_pos[i]
	  lines(xpos,ypos)
	  text(xpos[2],ypos[2],labels=i,pos=3)
	  }
  }

#this just gives us the list of how unusual the score is (stable or instable) -not that useful really
one_sided_probability_list=((0.5-probability_list)^2)^.5 *2 
names(chromosome_counts)=names(one_sided_probability_list)
sorted_one_sided_probability_list=sort(one_sided_probability_list)
plot(1:length(sorted_one_sided_probability_list),sorted_one_sided_probability_list,xlab="Chromosome Number",ylab="Unusual Chromosome stability/instability", type="n",  xaxt = "n")
text(1:length(sorted_one_sided_probability_list),sorted_one_sided_probability_list,labels=names(sorted_one_sided_probability_list))

#so the question is how does our instability score relate to other chromosome function (and does mapping to a dist give us anything)
#this is a bit dubious as we aren't adjusting for length for the protein count
instability_score=probability_list

cor.test(instability_score,chrinfo[c,"Confirmed.proteins"])
cor.test(instability_score,chrinfo[c,"Total.Prot.RNA"])

# if we don;t adjust for length we get ok results- fitting to dist doesn't give us much more apart from a pretty pic and soemthing that sounds clever...
# should be noted that the mapped to dist and normal scores are basically the same (very highly correlated as the scores are basically normally distributed) 
# but working with probs might be easier
if (plots)
  {
  plot(instability_score,chrinfo[c,"Confirmed.proteins"], main="Cancer instability and known protein counts", 
     sub="pearson cor=.31 (pval=0.16)",xlab="Chromosome instability score",ylab="Protein count", type="n")
  text(instability_score,chrinfo[c,"Confirmed.proteins"],labels=names(instability_score))
  #plot(instability_score,chrinfo[c,"Total.Prot.RNA"], main="Cancer instability and summed protein/RNA counts",sub="pearson cor=.31 (pval=0.16)",xlab="Chromosome instability score",ylab="Protein count", type="n")
  #text(instability_score,chrinfo[c,"Total.Prot.RNA"],labels=names(instability_score))
  }

# we can adjust for length
ins_corr = cor.test(instability_score,chrinfo[c,"Confirmed.proteins"]/chrinfo[c,"Base.pairs"]^.7)
corr_str = paste("pearson cor=", round(ins_corr$estimate, 2), "  (pval=", round(ins_corr$p.value, 2), ")", sep="")
if (plots)
  {
  #instability_score = sort(instability_score)
  plot(instability_score, chrinfo[c,"Confirmed.proteins"]/chrinfo[c,"Base.pairs"]^.7, 
      main="Cancer instability and known protein counts",
      sub=corr_str,xlab="Chromosome instability score",ylab="Length normalised protein count", type="p", col='blue')
  text(instability_score, chrinfo[c,"Confirmed.proteins"]/chrinfo[c,"Base.pairs"]^.7,labels=names(instability_score),pos=3)
  }

#but now 19 is an outlier- but if we ignore it we get slightly more significant results (reasonably high instability with very high protein encoding)
exclude_nineteen<-(names(instability_score) != "19")
cor.test(instability_score[names(instability_score[exclude_nineteen])],
  chrinfo[names(instability_score[exclude_nineteen]),"Confirmed.proteins"]/chrinfo[names(instability_score[exclude_nineteen]),"Base.pairs"]^.7)
cor.test(adjusted_scores[names(instability_score[exclude_nineteen])],
  chrinfo[names(instability_score[exclude_nineteen]),"Confirmed.proteins"]/chrinfo[names(instability_score[exclude_nineteen]),"Base.pairs"]^.7)

#so- what next...
dev.new()

plot(hclust(dist(chromosome_counts)))
plot(hclust(dist(adjusted_scores)))

ins_fit = hclust(dist(instability_score))
groups = cutree(ins_fit, k=3)
plot(ins_fit, main="Chromosome Instability")
rect.hclust(ins_fit, k=3, border=c("red", "blue", "green"))

for (i in 1:3)
  {
  print(paste("Group", i))
  print(groups[groups == i])  
  }

filename = "~/Analysis/Database/cancer/chr_instability_prob.txt"
write("# Normal distribution, probability score per chromosome. Each score is independent of the other chromosomes", file=filename)
write(paste("# Stats for scores used to generated scores - Mean:", round(mean(adjusted_scores), 3), "SD:", round(  sd(adjusted_scores), 3)), app=T, file=filename)
write.table( round(instability_score, 3), quote=F, col.names=F, sep="\t", app=T, file=filename)

