
bp = read.table("breakpoints.txt", sep="\t", comment="#", header=T)
chrinfo = read.table("../../chromosome_gene_info_2012.txt", sep="\t", row.names=1, header=T)

# don't need the mtDNA row
chrinfo = chrinfo[ -(nrow(chrinfo)), ]

bpfreq = table(bp$Breakpoint)
sorted = sort(bpfreq)

# ignoring X and Y
chromosome_counts=vector("numeric",22)
names(chromosome_counts)=1:22
for(i in 1:22)
	{
	sum(bp$Chr==i)->chromosome_counts[i]
	}

dev.new()
par(mfrow=c(3,1))
#shock/horror number of events correlates with length of chromosome
cor.test(chromosome_counts,chrinfo[1:22, "Base.pairs"])
plot(chromosome_counts, chrinfo[1:22,"Base.pairs"], type="n")
text(chromosome_counts, chrinfo[1:22,"Base.pairs"],labels=row.names(chrinfo)[1:22])

#and as length correlates with number of proteins this also correlates- but the score is higher
cor.test(chromosome_counts,chrinfo[1:22,"Confirmed.proteins"])
plot(chromosome_counts, chrinfo[1:22,"Confirmed.proteins"], type="n")
text(chromosome_counts, chrinfo[1:22,"Confirmed.proteins"],labels=row.names(chrinfo)[1:22])

#again length correlates with number of psudeogenes- although this is a very high correlation and so is interesting
cor.test(chromosome_counts,chrinfo[1:22,"Pseudogenes"])
plot(chromosome_counts, chrinfo[1:22,"Pseudogenes"], type="n")
text(chromosome_counts, chrinfo[1:22,"Pseudogenes"],labels=row.names(chrinfo)[1:22])

#problem is that this strong correlation is largely being driven by length of chromosome
cor.test(chromosome_counts,chrinfo[1:22,"Variations"])

#so we adjust the length in a non-linear manner (with a bit of mucking arround found ^.7 gave nearly no correlation)
adjusted_scores=chromosome_counts/(chrinfo[1:22,"Base.pairs"]^(.7))
#and now we don;t have a correlation with length- this is our chromosome instability score
cor.test(adjusted_scores[1:22],chrinfo[1:22,"Base.pairs"])
#and now we can (sorta) say that chromosome instability related directly to the amount of information encoded on that chromosome
cor.test(adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"])


#although this is being drived by chromosome "1", which means we will have to think about how to normalise the score a bit more
dev.new()
plot(adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"], type="n")
text(adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"],labels=names(adjusted_scores))

#anyway now as we wanna be clever we can map it to a probability distribution -first we check to see if adjusted_scores are normal
ks.test(adjusted_scores,pnorm,mean(adjusted_scores),sd(adjusted_scores))

#they are normal so we generate the proabilities of these scores occuring if they follow the normal dist
probability_list=pnorm(adjusted_scores,mean(adjusted_scores),sd(adjusted_scores))
dev.new()
#and now a (vaguely) pretty pic just for ewe
plot(function(x) dnorm(x,mean(adjusted_scores),sd(adjusted_scores)), min(adjusted_scores)-sd(adjusted_scores)*1,max(adjusted_scores)+sd(adjusted_scores)*1,ylab="Density",xlab="Instability Score")
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
	text(xpos[2],ypos[2],labels=i)
	}


#this just gives us the list of how unusual the score is (stable or instable) -not that useful really
one_sided_probability_list=((0.5-probability_list)^2)^.5 *2 
names(chromosome_counts)=names(one_sided_probability_list)
sorted_one_sided_probability_list=sort(one_sided_probability_list)

dev.new()


plot(1:length(sorted_one_sided_probability_list),sorted_one_sided_probability_list,xlab="Chromosome Number",ylab="Unusual Chromosome stability/instability", type="n",  xaxt = "n")
text(1:length(sorted_one_sided_probability_list),sorted_one_sided_probability_list,labels=names(sorted_one_sided_probability_list))

#so the question is how does our instability score relate to other chromosome function (and does mapping to a dist give us anything)
#this is a bit dubious as we aren't adjusting for length for the protein count
instability_score=probability_list
names(instability_score)=names(chromosome_counts)

cor.test(instability_score,chrinfo[1:22,"Confirmed.proteins"])

#if we don;t adjust for length we get ok results- fitting to dist doesn't give us much more apart from a pretty pic and soemthing that sounds clever...
#should be noted that the mapped to dist and normal scores are basically the same (very highly correlated as the scores are basically normally distributed) but working with probs might be easier

dev.new()
par(mfrow=c(2,1))
plot(instability_score,chrinfo[1:22,"Confirmed.proteins"], main="Cancer instability and known protein counts", sub="pearson cor=.31 (pval=0.16)",xlab="Chromosome instability score",ylab="Protein count", type="n")
text(instability_score,chrinfo[1:22,"Confirmed.proteins"],labels=names(instability_score))

#we can adjust for length
cor.test(instability_score,chrinfo[1:22,"Confirmed.proteins"]/chrinfo[1:22,"Base.pairs"]^.7)
plot(instability_score, chrinfo[1:22,"Confirmed.proteins"]/chrinfo[1:22,"Base.pairs"]^.7, 
main="Cancer instability and known protein counts",
sub="pearson cor=.31 (pval=0.15)",xlab="Chromosome instability score",ylab="Length normalised protein count", type="n")
text(instability_score, chrinfo[1:22,"Confirmed.proteins"]/chrinfo[1:22,"Base.pairs"]^.7,labels=names(instability_score))

#but now 19 is an outlier- but if we ignore it we get more significant results (reasonably high instability with very high protein encoding)
exclude_nineteen<-(names(instability_score) != "19")
cor.test(instability_score[names(instability_score[exclude_nineteen])],
chrinfo[names(instability_score[exclude_nineteen]),"Confirmed.proteins"]/chrinfo[names(instability_score[exclude_nineteen]),"Base.pairs"]^.7)
cor.test(adjusted_scores[names(instability_score[exclude_nineteen])],
chrinfo[names(instability_score[exclude_nineteen]),"Confirmed.proteins"]/chrinfo[names(instability_score[exclude_nineteen]),"Base.pairs"]^.7)

#so- what next...
dev.new()
plot(hclust(dist(chromosome_counts)))
plot(hclust(dist(adjusted_scores)))
plot(hclust(dist(instability_score)))

