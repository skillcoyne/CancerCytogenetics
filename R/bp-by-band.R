## This pretty much just supports the initial by chromosome analysis.  
# Which means I need to score the bands by the number of genes perhaps?

## Three classes of breakpoints for the purposes of insilico generation
# 1. centromeres
# 2. Top breakpoints 
# 3. The rest based on correlation with gene count.
## Note that chromosomes 19,20,22 don't have enough observations due to the fact that most/all of their breakpoints occur within centromeres

rm(list=ls())

instabilityPlot<-function(x, y, names, corr="")
  {
  plot(x, y, main="Band instability and known gene counts", sub=corr, xlab="Band instability score", ylab="Length normalised gene count", type="n")
  text(x, y, labels=names)
  }

corTests<-function(bpinfo)
  {
  chromosomes = c(1:22, 'X')
  
  test_cols = c('length', 'gene.count', 'ks.gene', 'length.normalized', 'gene.normalized', 'ks.normalized')
  
  bk_tests = as.data.frame( matrix( nrow=length(chromosomes), ncol=length(test_cols)) )
  colnames(bk_tests) = test_cols
  rownames(bk_tests) = chromosomes
  
  rc_tests = as.data.frame( matrix( nrow=length(chromosomes), ncol=length(test_cols)) )
  colnames(rc_tests) = test_cols
  rownames(rc_tests) = chromosomes
  
  for (chr in chromosomes)
    {
    scores = bpinfo[bpinfo$chr == chr,]
    scores = scores[order(scores$band),]
    scores = na.omit(scores)
    # get the arms in the correct order p -> q
    parm = grep("p", scores$band )
    p = scores[parm,]
    p = p[order(p$band, decreasing=T),]
    scores[parm,] = p
    
    tryCatch({
      # do the number of breaks correlate to the band length? 
      ct = cor.test(scores$bp.length, scores$total.breaks)
      bk_tests[chr, 'length'] = round(ct$estimate, 4)
    
      ct = cor.test(scores$bp.length, scores$total.recombination)
      rc_tests[chr, 'length'] = round(ct$estimate, 4)
    
      # so what about gene counts?
      ct = cor.test(scores$gene.count, scores$total.breaks)
      bk_tests[chr, 'gene.count'] = round(ct$estimate, 4)
    
      ct = cor.test(scores$gene.count, scores$total.recombination)
      rc_tests[chr, 'gene.count'] = round(ct$estimate, 4)
    
      ct = ks.test(scores$total.breaks, pnorm, mean(scores$total.breaks), sd(scores$total.breaks))
      bk_tests[chr, 'ks.gene'] = round(ct$p.value, 4)
      
      ct = ks.test(scores$total.recombination, pnorm, mean(scores$total.recombination), sd(scores$total.recombination))
      rc_tests[chr, 'ks.gene'] = round(ct$p.value, 4)

      # Normalize the total breaks/recombinations for length 
      breaks_adj = scores$total.breaks/(scores$bp.length^.7)
      recomb_adj = scores$total.recombination/(scores$bp.length^.7)
    
      ct = cor.test(breaks_adj, scores$bp.length)
      bk_tests[chr, 'length.normalized'] = round(ct$estimate, 4)
    
      ct = cor.test(recomb_adj, scores$bp.length)
      rc_tests[chr, 'length.normalized'] = round(ct$estimate, 4)
    
      ct = cor.test(breaks_adj, scores$gene.count/scores$bp.length^.7)
      bk_tests[chr, 'gene.normalized'] = round(ct$p.value, 4)
    
      ct = cor.test(recomb_adj, scores$gene.count/scores$bp.length^.7)
      rc_tests[chr, 'gene.normalized'] = round(ct$p.value, 4)
    
      # are the adjusted scores are normal
      ct = ks.test(breaks_adj, pnorm, mean(breaks_adj), sd(breaks_adj))
      bk_tests[chr, 'ks.normalized'] = round(ct$p.value, 4)
    
      ct = ks.test(recomb_adj, pnorm, mean(recomb_adj), sd(recomb_adj))
      rc_tests[chr, 'ks.normalized'] = round(ct$p.value, 4)
      }, error = function(e) {
        print(e)
      })
    }
  return(list("bk" = bk_tests, "rc" = rc_tests))
  }

plotCor<-function(bk_tests, rc_tests)
  {
  chrinfo = loadChromosomeInfo("../genomic_info/chromosome_gene_info_2012.txt")  
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
  }

setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")

datadir = "/Volumes/Spark/Data/sky-cgh/output"
setwd(datadir)

`%nin%`=Negate(`%in%`) 
bandinfo = read.table(paste(datadir,"band_genes.txt", sep="/"), header=T, sep="\t")

df = read.table("current/breakpoints.txt", sep="\t", header=T)


## To simplify, the bands with < 5 aberrations are all subbands.  We're only dealing with major bands and these
# don't make a difference in the analysis
df = df[-which(df$total.aberrations <= 5),]


break_info = merge(df, bandinfo[,c(1,2,5)], by.x=c('chr', 'band'), by.y=c('chr','band'))
break_info$bp.length = break_info$end - break_info$start

# All bands with breaks but NO genes are in the centromere regions...pull them out see if correlations improve
# Note that not all centromeres lack genes apparently

## CLASS 1 ##
cmrows = grep("(11|12)", break_info$band)
centromeres = break_info[ cmrows, ]
centromeres$chance.to.break = round(centromeres[,'total.aberrations']/sum(centromeres[,'total.aberrations']), 4)

# Not using this for anything, it mostly just confirms the chromosome instability analysis
break_info = break_info[ -cmrows,] 
allcors = corTests(break_info)

# So lets see the range of breakpoints -- AGAIN leukemia bias is a problem here
break_info = break_info[ order(break_info$total.aberrations, decreasing=T),]

sd(break_info$total.aberrations)
summary(break_info$total.aberrations)

## CLASS 2 ##
topBP = break_info[ break_info$total.aberrations >= mean(break_info$total.aberrations), ]
topBP$chance.to.break = round(topBP[,'total.aberrations']/sum(topBP[,'total.aberrations']), 4)

## CLASS 3 ##
remainder = break_info[ break_info$total.aberrations < mean(break_info$total.aberrations), ]

remCor = corTests(remainder)
plotCor(remCor$bk, remCor$rc)

remainder$chance.to.break = round(remainder[,'total.aberrations']/sum(remainder[,'total.aberrations']), 4)



