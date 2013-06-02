resetwd()
source("R/lib/load_files.R")
source("R/lib/wd.R")

chromosome_count<-function(bp)
	{
	# ignoring X and Y
	chromosome_counts=vector("numeric",22)
	names(chromosome_counts)=1:22
	for(i in 1:22)
		{
		sum(bp$Chr==i)->chromosome_counts[i]
		}
	return(chromosome_counts)
	}

adjust_for_length<-function(chromosome_counts, chrinfo)
	{
	#-- so we adjust the length in a non-linear manner (with a bit of mucking arround found ^.7 gave nearly no correlation) --#
	adjusted_scores=chromosome_counts/(chrinfo[1:22,"Base.pairs"]^(.7))
	#and now we don;t have a correlation with length- this is our chromosome instability score
	ct = cor.test(adjusted_scores[1:22],chrinfo[1:22,"Base.pairs"])
	message( paste("Base.pairs\n\t", ct$method, "\n\tCorrelation:", round(ct$estimate, 2), " P-value:", round(ct$p.value, 2), sep="")  )
	#and now we can (sorta) say that chromosome instability related directly to the amount of information encoded on that chromosome
	ct = cor.test(adjusted_scores[1:22],chrinfo[1:22,"Confirmed.proteins"])
	message( paste("Confirmed.proteins\n\t",ct$method, "\n\tCorrelation:", round(ct$estimate, 2), " P-value:", round(ct$p.value, 2), sep="")  )
	return (adjusted_scores)
	}
	
instability<-function()
	{
	instability_score=probability_list
names(instability_score)=names(chromosome_counts)

cor.test(instability_score,chrinfo[1:22,"Confirmed.proteins"])
cor.test(instability_score,chrinfo[1:22,"Total.Prot.RNA"])

	}


norm_dist_plot<-function(adjusted_scores, label)
	{
	#and now a (vaguely) pretty pic just for ewe
  	plot(function(x) dnorm(x,mean(adjusted_scores),sd(adjusted_scores)), min(adjusted_scores)-sd(adjusted_scores)*1,max(adjusted_scores)+sd(adjusted_scores)*1,
     ylab="Density",xlab=paste("Chromosome Instability Score", label, sep=" "))
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

instability_plot<-function(instability_score, chrinfo, label)
	{
	#ct = cor.test(instability_score,chrinfo[1:22,"Confirmed.proteins"])
	#corr_str = paste("pearson cor=", round(ct$estimate, 2), "  (pval=", round(ct$p.value, 2), ")", sep="")
	#xlabel = paste("Chromosome instability score", label, sep=" ")
	#plot(instability_score,chrinfo[1:22,"Confirmed.proteins"], main="Cancer instability and known protein counts", 
    #sub=corr_str,xlab=xlabel,ylab="Protein count", type="n")
  	#text(instability_score,chrinfo[1:22,"Confirmed.proteins"],labels=names(instability_score))
	ct = cor.test(instability_score,chrinfo[1:22,"Confirmed.proteins"]/chrinfo[1:22,"Base.pairs"]^.7)
	corr_str = paste("pearson cor=", round(ct$estimate, 2), "  (pval=", round(ct$p.value, 2), ")", sep="")
  	plot(instability_score, chrinfo[1:22,"Confirmed.proteins"]/chrinfo[1:22,"Base.pairs"]^.7, 
      main=paste("Cancer Instability", label, sep=" "),
      	sub=corr_str,xlab="Chromosome instability score",ylab="Length normalised protein count", type="n")
  	text(instability_score, chrinfo[1:22,"Confirmed.proteins"]/chrinfo[1:22,"Base.pairs"]^.7,labels=names(instability_score))
	}

cluster_plot<-function(instability_score, label)
	{
	ins_fit = hclust(dist(instability_score))
	groups = cutree(ins_fit, k=3)
	plot(ins_fit, main=paste("Chromosome Instability", label, sep=" "))
	rect.hclust(ins_fit, k=3, border=c("red", "blue", "green"))
	}

chrinfo = loadChromosomeInfo()  
`%nin%`=Negate(`%in%`) 


setDataDirectory("09012013")
bp_09 = loadBreakpoints("breakpoints.txt")

setDataDirectory()
bp_19 = loadBreakpoints("breakpoints.txt")

nrow(bp_09)
nrow(bp_19)

## sample the most frequent leukemia/lymphoma cases 
sample_leukemia_bps = TRUE
if (sample_leukemia_bps)
  {
  leuks = c('Acute myeloid leukemia', 'Acute lymphoblastic leukemia', "Non-hodgkin's lymphoma", 'Chronic myelogenous leukemia', 'Chronic lymphocytic leukemia')
  bp_09 = sampleCancers(bp_09, "Cancer", leuks)
  bp_19 = sampleCancers(bp_19, "Cancer", leuks)
  }






# ignoring X and Y
chr_cnt_09=chromosome_count(bp_09)
chr_cnt_19=chromosome_count(bp_19)

adj_09=adjust_for_length(chr_cnt_09, chrinfo)
adj_19=adjust_for_length(chr_cnt_19, chrinfo)


ks.test(adj_09,pnorm,mean(adj_09),sd(adj_09))
ks.test(adj_19,pnorm,mean(adj_19),sd(adj_19))

dev.new()
par(mfrow=c(1,2))
norm_dist_plot(adj_09, "09")
norm_dist_plot(adj_19, "19")


prob_09=pnorm(adj_09,mean(adj_09),sd(adj_09))
prob_19=pnorm(adj_19,mean(adj_19),sd(adj_19))

instability_score_09=prob_09
names(instability_score_09)=names(chr_cnt_09)

instability_score_19=prob_19
names(instability_score_19)=names(chr_cnt_19)

dev.new()
par(mfrow=c(1,2))
instability_plot(instability_score_09, chrinfo, "09")
instability_plot(instability_score_19, chrinfo, "19")

## The clusters are different!
dev.new()
par(mfrow=c(1,2))
cluster_plot(instability_score_09, "09")
cluster_plot(instability_score_19, "19")

dev.new()
par(mfrow=c(1,2))
plot(hclust(dist(chr_cnt_09)))
plot(hclust(dist(chr_cnt_19)))

dev.new()
par(mfrow=c(1,2))
plot(hclust(dist(adj_09)))
plot(hclust(dist(adj_19)))

## TODO:  Leukemia skews the results. May be worth doing on a per-cancer basis but the cancers should be clustered by type somehow

chr_cnts = cbind(chr_cnt_09, chr_cnt_19)



