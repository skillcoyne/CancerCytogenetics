## This script does an analysis of the instability of each chromosome. The current result is a discrete distribution of probabilities ##


rm(list=ls())
setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")

datadir = "~/Data/sky-cgh/output"
setwd(datadir)

`%nin%`=Negate(`%in%`) 
chrinfo = loadChromosomeInfo("../genomic_info/chromosome_gene_info_2012.txt")  
chrinfo$total.genes = rowSums(chrinfo[,c('Confirmed.proteins','Pseudogenes','miRNA','rRNA','snRNA','snoRNA','Misc.ncRNA')])
chrinfo$gene.density.mb = chrinfo[,'total.genes']/(chrinfo[,'Base.pairs']/1E6)

plot(chrinfo$Base.pairs/1E6, chrinfo$gene.density.mb, type='p', pch=19, col='blue', ylab='Gene Density per Mb', xlab='Chr length in Mb', main='Gene Density')
text(chrinfo$Base.pairs/1E6, chrinfo$gene.density.mb, labels=rownames(chrinfo), pos=4)

plot(chrinfo[,'Base.pairs'],chrinfo[,'Confirmed.proteins'], type='n')
text(chrinfo[,'Base.pairs'],chrinfo[,'Confirmed.proteins'], labels=rownames(chrinfo))

leuk = T
plots = T

# Load files
bp = load.breakpoints("current/noleuk-breakpoints.txt")
if (leuk)
  bp = load.breakpoints(c("current/noleuk-breakpoints.txt", "current/leuk-breakpoints.txt"))


total_major_bps = nrow(bp)

# ignoring X and Y
c = c(1:22)#, 'X', 'Y')
chromosome_counts=vector("numeric",length(c))
names(chromosome_counts)=c
for(i in names(chromosome_counts))
  chromosome_counts[i] = sum(bp[bp$chr == i,'total.karyotypes'])

# Obviously correlates to length
cor.test(chromosome_counts,chrinfo[c,"Base.pairs"])
# and to the number of protein coding 
cor.test(chromosome_counts,chrinfo[c,"Confirmed.proteins"])
cor.test(chromosome_counts,chrinfo[c,'total.genes'])

# But not to the gene density
cor.test(chromosome_counts,chrinfo[c,"gene.density.mb"])


length_adjust = 0.7
#-- so we adjust the length in a non-linear manner (with a bit of mucking arround found ^.7 gave nearly no correlation) --#
adjusted_scores=chromosome_counts/(chrinfo[c,"Base.pairs"]^length_adjust)

# but adjusting with gene density doesn't get rid of the correlation so...
#adjusted_scores = chromosome_counts/chrinfo[c,'gene.density.mb']

#and now we don;t have a correlation with length- this is our chromosome instability score
cor.test(adjusted_scores[c],chrinfo[c,"Base.pairs"])
#and now we can (sorta) say that chromosome instability related directly to the amount of information encoded on that chromosome
cor.test(adjusted_scores[c],chrinfo[c,"total.genes"])
cor.test(adjusted_scores[c],chrinfo[c,"Confirmed.proteins"])


if (plots)
  {
#  par(mfrow=c(2,1))
  #although this is being drived by chromosome "1", which means we will have to think about how to normalise the score a bit more
  plot(adjusted_scores[c],chrinfo[c,"Confirmed.proteins"], type="n", xlab="Chromosome Counts adjusted for length", ylab="Protein Counts")
  text(adjusted_scores[c],chrinfo[c,"Confirmed.proteins"],labels=names(adjusted_scores))
  
#  plot(adjusted_scores[c],chrinfo[c,"total.genes"], type="n", xlab="Chromosome Counts adjusted for length", ylab="Gene Counts")
#  text(adjusted_scores[c],chrinfo[c,"total.genes"],labels=names(adjusted_scores))
  
  
  }

#anyway now as we wanna be clever we can map it to a probability distribution -first we check to see if adjusted_scores are normal
ks.test(adjusted_scores,pnorm,mean(adjusted_scores),sd(adjusted_scores))

## Normal distribution
if (plots)
  {
  #and now a (vaguely) pretty pic just for ewe
  plot(function(x) dnorm(x,mean(adjusted_scores),sd(adjusted_scores)), min(adjusted_scores)-sd(adjusted_scores)*1,max(adjusted_scores)+sd(adjusted_scores)*1,
     ylab="Density",xlab="Chromosome Instability Score", col='blue',lwd=2)
  xpos=vector("numeric",2)
  ypos=vector("numeric",2)
  ypos[1]=0
  density_pos=dnorm(adjusted_scores,mean(adjusted_scores),sd(adjusted_scores))
  for(i in 1:length(adjusted_scores))
	  {
	  xpos[1]=adjusted_scores[i]
	  xpos[2]=adjusted_scores[i]
	  ypos[2]=density_pos[i]
	  lines(xpos,ypos,col='blue',lwd=2)
    #if (i %% 2 == 0) ypos[2] = ypos[2]+1
	  text(xpos[2],ypos[2],labels=names(adjusted_scores[i]),pos=3)
	  }
  }

## -- Probabilities are probably the most useful -- ##
#they are normal so we generate the proabilities of these scores occuring if they follow the normal dist
probability_list = pnorm(adjusted_scores,mean(adjusted_scores),sd(adjusted_scores))

#this just gives us the list of how unusual the score is (stable or instable) -not that useful really
one_sided_probability_list=((0.5-probability_list)^2)^.5 *2 
names(chromosome_counts)=names(one_sided_probability_list)
sorted_one_sided_probability_list=sort(one_sided_probability_list)
plot(1:length(sorted_one_sided_probability_list),sorted_one_sided_probability_list,xlab="Chromosome Number",ylab="Unusual Chromosome stability/instability", type="n", xaxt="n")
text(1:length(sorted_one_sided_probability_list),sorted_one_sided_probability_list,labels=names(sorted_one_sided_probability_list))

#so the question is how does our instability score relate to other chromosome function (and does mapping to a dist give us anything)
#this is a bit dubious as we aren't adjusting for length for the protein count
instability_score=probability_list


ins_corr = cor.test(instability_score[c], chrinfo[c,"Confirmed.proteins"])
corr_str = paste("pearson cor=", round(ins_corr$estimate, 2), "  (pval=", round(ins_corr$p.value, 2), ")", sep="")


# if we don;t adjust for length we get ok results- fitting to dist doesn't give us much more apart from a pretty pic and soemthing that sounds clever...
# should be noted that the mapped to dist and normal scores are basically the same (very highly correlated as the scores are basically normally distributed) 
# but working with probs might be easier
if (plots)
  {
  plot(instability_score,chrinfo[c,"Confirmed.proteins"], main="Chromosome Instability", 
     sub=corr_str,xlab="Chromosome instability score",ylab="Protein count", type="n")
  text(instability_score,chrinfo[c,"Confirmed.proteins"],labels=names(instability_score),col='blue')
  }


# we can adjust for length
ins_corr = cor.test(instability_score,chrinfo[c,"Confirmed.proteins"]/(chrinfo[c,"Base.pairs"]^length_adjust))
corr_str = paste("pearson cor=", round(ins_corr$estimate, 2), "  (pval=", round(ins_corr$p.value, 2), ")", sep="")
if (plots)
  {
  #instability_score = sort(instability_score)
  plot(instability_score, chrinfo[c,"Confirmed.proteins"]/(chrinfo[c,"Base.pairs"]^length_adjust), 
      main="Chromosome Instability", col='blue',
      sub=corr_str,xlab="Chromosome instability score",ylab="Length normalised protein count", type="p", pch=19)
  text(instability_score, chrinfo[c,"Confirmed.proteins"]/(chrinfo[c,"Base.pairs"]^length_adjust),labels=names(instability_score),  pos=2)
  }


plot(inst_all, chrinfo[c,"Confirmed.proteins"]/(chrinfo[c,"Base.pairs"]^length_adjust), 
     main="Chromosome Instability", xlab="Chromosome instability score",ylab="Length normalised protein count", col='blue', type="p", pch=19)
text(inst_all, chrinfo[c,"Confirmed.proteins"]/(chrinfo[c,"Base.pairs"]^length_adjust),labels=names(inst_all), col='blue', pos=4)

points(inst_noleuk, chrinfo[c,"Confirmed.proteins"]/(chrinfo[c,"Base.pairs"]^length_adjust), col='red', type="p", pch=19)
text(inst_noleuk, chrinfo[c,"Confirmed.proteins"]/(chrinfo[c,"Base.pairs"]^length_adjust),labels=names(inst_noleuk), col='red',pos=2)
legend('topright', fill=c('blue','red'), legend=c("All data","No leukemia"), border=F)


plot(sort(inst_all-inst_noleuk), main="Instability Difference",type='p',col='purple',pch=19,ylab="All-No Leukemia")
text(sort(inst_all-inst_noleuk), labels=names(sort(inst_all-inst_noleuk)), pos=4)

#but now 19 is an outlier- but if we ignore it we get slightly more significant results (reasonably high instability with very high protein encoding)
enin<-which(names(instability_score) %nin% c("19"))
cor.test(instability_score, chrinfo[c,"Confirmed.proteins"]/chrinfo[c,"Base.pairs"]^length_adjust)
cor.test(instability_score[enin], chrinfo[enin,"Confirmed.proteins"]/chrinfo[enin,"Base.pairs"]^length_adjust)
cor.test(adjusted_scores[enin], chrinfo[enin,"Confirmed.proteins"]/chrinfo[enin,"Base.pairs"]^length_adjust)

#so- what next...
dev.off()

hc = hclust(dist(instability_score))
groups = cutree(hc, k=3)
plot(hc, main="Chromosome Instability", cex=2)
rect.hclust(hc, k=3, border=c("red", "blue", "green"))

plot(hc, col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071",
     col.axis = "#F38630", lwd = 3, lty = 3, sub = '', hang = -1, axes = FALSE)
rect.hclust(hc, k=3, border=c("red", "blue", "green"))

source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")

A2Rplot(hc, k = 3, boxes = FALSE, col.up = "gray50", col.down = c("red", "blue", "purple"), main="Chromosome Instability")


for (i in 1:3)
  {
  print(paste("Group", i))
  print(groups[groups == i])  
  print(instability_score[groups == i])
  }


write.table(bp[order(-bp$total.karyotypes),][1:10, c('chr','band','total.karyotypes')], quote=F,sep="\t", row.names=F)

#filename = "~/Analysis/Database/cancer/chr_instability_prob.txt"
#write("# Normal distribution, probability score per chromosome. Each score is independent of the other chromosomes", file=filename, app=F)
#write.table( round(probability_list/sum(probability_list), 5), quote=F, col.names=F, sep="\t", app=T, file=filename)

