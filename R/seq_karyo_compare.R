#source("R/lib/load_files.R")
#source("R/lib/wd.R")

count_chr<-function(df, colname = 'Chr')
  {
  # ignoring X and Y, and the odd chromosomes that are likely just bad file formatting
  chromosome_counts=vector("numeric",22)
  names(chromosome_counts)=1:22
  
  df = df[df[[colname]] %in% names(chromosome_counts),]
  freq = table(df[[colname]])
  for (i in 1:22)
    {
    chromosome_counts[i] = freq[[as.character(i)]]
    }
  return(chromosome_counts)
  }


setwd("~/Data/TCGA")
`%nin%`=Negate(`%in%`) 
if (vars != NULL)
  {
  vars = read.table("all_variants.txt", header=T, sep="\t") 
  ## Not looking at SNPs may help
  vars = vars[vars$VarType %nin% c('SNP'),]
  }

chrinfo = loadChromosomeInfo("../sky-cgh/genomic_info/chromosome_gene_info_2012.txt")  

seq_chr_count = count_chr(vars)
# looks like a clear correlation with length
plot(seq_chr_count[1:22],type="n", xlab="Chromosome", ylab="Variant count")
text(seq_chr_count[1:22],labels=names(seq_chr_count))
# yep, correlated with length
cor.test(seq_chr_count[1:22],chrinfo[1:22,"Base.pairs"])
#-- so we adjust the length in a non-linear manner (with a bit of mucking arround found ^.7 gave nearly no correlation) --#
seq_adjusted_scores=seq_chr_count/(chrinfo[1:22,"Base.pairs"]^(.7))
cor.test(seq_adjusted_scores[1:22],chrinfo[1:22,"Base.pairs"])

#and now we don;t have a correlation with length- this is our chromosome instability score
cor.test(seq_adjusted_scores[1:22],chrinfo[1:22,"Base.pairs"])

# and now we can (sorta) say that chromosome instability related directly to the amount of information encoded on that chromosome
# except that 1 is still a very obvious outlier...
plot(seq_adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"], type="n", xlab="Chromosome Counts adjusted for length", ylab="Protein Counts")
text(seq_adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"],labels=names(seq_adjusted_scores))

#anyway now as we wanna be clever we can map it to a probability distribution -first we check to see if adjusted_scores are normal
ks.test(seq_adjusted_scores,pnorm,mean(seq_adjusted_scores),sd(seq_adjusted_scores))

## -- Probabilities are probably the most useful -- ##
# they are normal so we generate the proabilities of these scores occuring if they follow the normal dist
probability_list=pnorm(seq_adjusted_scores,mean(seq_adjusted_scores),sd(seq_adjusted_scores))

# and now a (vaguely) pretty pic which is fairly different from the karyotype one...
plot(function(x) dnorm(x,mean(seq_adjusted_scores),sd(seq_adjusted_scores)), min(seq_adjusted_scores)-sd(seq_adjusted_scores)*1,max(seq_adjusted_scores)+sd(seq_adjusted_scores)*1,
     ylab="Density",xlab="Chromosome Instability Score")
xpos=vector("numeric",2)
ypos=vector("numeric",2)
ypos[1]=0
density_pos=dnorm(seq_adjusted_scores,mean(seq_adjusted_scores),sd(seq_adjusted_scores))
for(i in 1:length(seq_adjusted_scores))
  {
  xpos[1]=seq_adjusted_scores[i]
  xpos[2]=seq_adjusted_scores[i]
  ypos[2]=density_pos[i]
  lines(xpos,ypos)
  text(xpos[2],ypos[2],labels=i,pos=3)
  }

instability_score=probability_list
names(instability_score)=names(seq_chr_count)

# and 1 is still out there, but the correlations are still ok
cor.test(instability_score,chrinfo[1:22,"Confirmed.proteins"])
plot(instability_score,chrinfo[1:22,"Confirmed.proteins"], main="Cancer instability and known protein counts", 
     sub="pearson cor=.31 (pval=0.16)",xlab="Chromosome instability score",ylab="Protein count", type="n")
text(instability_score,chrinfo[1:22,"Confirmed.proteins"],labels=names(instability_score))

#we can adjust for length
ins_corr = cor.test(instability_score,chrinfo[1:22,"Confirmed.proteins"]/chrinfo[1:22,"Base.pairs"]^.7)
corr_str = paste("pearson cor=", round(ins_corr$estimate, 2), "  (pval=", round(ins_corr$p.value, 2), ")", sep="")

# and now it's 19...which also happens in the karyotypes.  However, 17 is on the high end now which is good
plot(instability_score, chrinfo[1:22,"Confirmed.proteins"]/chrinfo[1:22,"Base.pairs"]^.7, 
       main="Cancer instability and known protein counts",
       sub=corr_str,xlab="Chromosome instability score",ylab="Length normalised protein count", type="n")
text(instability_score, chrinfo[1:22,"Confirmed.proteins"]/chrinfo[1:22,"Base.pairs"]^.7,labels=names(instability_score))

#but now 19 is an outlier- but if we ignore it we get more significant results (reasonably high instability with very high protein encoding)
exclude_nineteen<-(names(instability_score) != "19")
cor.test(instability_score[names(instability_score[exclude_nineteen])],
         chrinfo[names(instability_score[exclude_nineteen]),"Confirmed.proteins"]/chrinfo[names(instability_score[exclude_nineteen]),"Base.pairs"]^.7)
cor.test(seq_adjusted_scores[names(instability_score[exclude_nineteen])],
         chrinfo[names(instability_score[exclude_nineteen]),"Confirmed.proteins"]/chrinfo[names(instability_score[exclude_nineteen]),"Base.pairs"]^.7)


plot(hclust(dist(seq_chr_count)))  # fairly different from the karyotypes
plot(hclust(dist(seq_adjusted_scores)))  # 17/19 are the highest cluster

# The clusters are different, but share some similarities.  In cluster 1: 11,17,19 are shared by cluster 1 in the karyotype analysis
ins_fit = hclust(dist(instability_score))
groups = cutree(ins_fit, k=3)
plot(ins_fit, main="Chromosome Instability")
rect.hclust(ins_fit, k=3, border=c("red", "blue", "green"))



