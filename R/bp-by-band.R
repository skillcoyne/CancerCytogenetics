## This pretty much just supports the initial by chromosome analysis.  Which means I need to score the bands by the number of genes perhaps?

rm(list=ls())

instabilityPlot<-function(x, y, names, corr="")
  {
  plot(x, y, main="Band instability and known gene counts", sub=corr, xlab="Band instability score", ylab="Length normalised gene count", type="n")
  text(x, y, labels=names)
  }


setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")

datadir = "~/Data/sky-cgh/output/current"
setwd(datadir)

`%nin%`=Negate(`%in%`) 
chrinfo = loadChromosomeInfo()  
bandinfo = read.table("~/Data/sky-cgh/output/band_genes.txt", header=T, sep="\t")

df = read.table("breakpoints.txt", sep="\t", header=T)

chromosomes = c(1:22, 'X')

bi = data.frame()
for (chr in chromosomes)
  {
  bp = df[df$chr == chr,]
  bp$band = clearSubbands(bp$band)

  uniq_bands = unique(bp$band)
  tmp = data.frame()
  for (i in 1:length(uniq_bands))
    {
    b = uniq_bands[i]
    curr = bp[bp$band == b,]
    
    tmp[i,'chr'] = chr
    tmp[i,'band'] = b
    tmp[i,'bp.length'] = curr[1,'end'] - curr[1,'start']
    tmp[i,'total.abr'] = sum(curr$total.aberrations)
    tmp[i,'total.breaks'] = sum(curr[,'total.breaks'])
    tmp[i,'total.recombination'] = sum(curr[,'total.recombination'])
    }
  bi = rbind(bi, tmp)
  }

break_info = merge(bi, bandinfo, by.x=c('chr', 'band'), by.y=c('chr','band'))

# All bands with breaks but NO genes are in the centromere regions...pull them out see if correlations improve
# Note that not all centromeres lack genes apparently
all_centromeres = grep("(11|12)", break_info$band)
centromeres = which(break_info$gene.count <= 0)
break_info = break_info[ -centromeres,] 

test_cols = c('length', 'gene.count', 'length.normalized', 'gene.normalized', 'ks')

bk_tests = as.data.frame( matrix( nrow=length(chromosomes), ncol=length(test_cols)) )
colnames(bk_tests) = test_cols
rownames(bk_tests) = chromosomes

rc_tests = as.data.frame( matrix( nrow=length(chromosomes), ncol=length(test_cols)) )
colnames(rc_tests) = test_cols
rownames(rc_tests) = chromosomes

for (chr in chromosomes)
  {
  scores = break_info[break_info$chr == chr,]
  scores = scores[order(scores$band),]

  # get the arms in the correct order p -> q
  parm = grep("p", scores$band )
  p = scores[parm,]
  p = p[order(p$band, decreasing=T),]
  scores[parm,] = p

  # do the number of breaks correlate to the band length? 
  ct = cor.test(scores$bp.length, scores$total.breaks)
  bk_tests[chr, 1] = round(ct$estimate, 4)
  
  ct = cor.test(scores$bp.length, scores$total.recombination)
  rc_tests[chr, 1] = round(ct$estimate, 4)
  
  # so what about gene counts?
  ct = cor.test(scores$gene.count, scores$total.breaks)
  bk_tests[chr, 2] = round(ct$estimate, 4)

  ct = cor.test(scores$gene.count, scores$total.recombination)
  rc_tests[chr, 2] = round(ct$estimate, 4)
  
  # Normalize the total breaks/recombinations for length 
  breaks_adj = scores$total.breaks/(scores$bp.length^.7)
  recomb_adj = scores$total.recombination/(scores$bp.length^.7)
  
  ct = cor.test(breaks_adj, scores$bp.length)
  bk_tests[chr, 3] = round(ct$estimate, 4)
  
  ct = cor.test(recomb_adj, scores$bp.length)
  rc_tests[chr, 3] = round(ct$estimate, 4)
  
  ct = cor.test(breaks_adj, scores$gene.count/scores$bp.length^.7)
  bk_tests[chr, 4] = round(ct$p.value, 4)
  
  ct = cor.test(recomb_adj, scores$gene.count/scores$bp.length^.7)
  rc_tests[chr, 4] = round(ct$p.value, 4)
  
  # are the adjusted scores are normal
  ct = ks.test(breaks_adj, pnorm, mean(breaks_adj), sd(breaks_adj))
  bk_tests[chr, 5] = round(ct$p.value, 4)
  
  ct = ks.test(recomb_adj, pnorm, mean(recomb_adj), sd(recomb_adj))
  rc_tests[chr, 5] = round(ct$p.value, 4)
  }

chrinfo = loadChromosomeInfo()  

# Not exactly the same as the chromosome level instability analysis, but pretty close
par(mfrow=c(2,3))
for (i in 1:ncol(bk_tests))
  {
  label = "cor estimate"
  if (i == 4) title = "p.value"
  
  plot(bk_tests[,i], chrinfo[1:23, 'Confirmed.proteins']/chrinfo[1:23,'Base.pairs']^.7,  col='blue', type='p', main=colnames(bk_tests)[i], 
       ylab="Length normalized protein count", xlab=label, sub="blue=break, red=recombination")
  text(bk_tests[,i], chrinfo[1:23, 'Confirmed.proteins']/chrinfo[1:23,'Base.pairs']^.7 ,labels=rownames(bk_tests),pos=3)
  lines(rc_tests[,i], chrinfo[1:23, 'Confirmed.proteins']/chrinfo[1:23,'Base.pairs']^.7,  col='red', type='p')
  text(rc_tests[,i], chrinfo[1:23, 'Confirmed.proteins']/chrinfo[1:23,'Base.pairs']^.7 ,labels=rownames(rc_tests),pos=3)
  
  }
