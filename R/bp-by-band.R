## This pretty much just supports the initial by chromosome analysis.  
# Which means I need to score the bands by the number of genes perhaps?

## Three classes of breakpoints for the purposes of insilico generation
# 1. centromeres
# 2. Top breakpoints 
# 3. The rest based on correlation with gene count.
## Note that chromosomes 19,20,22 don't have enough observations due to the fact that most/all of their breakpoints occur within centromeres

rm(list=ls())

cor.tests<-function(scores, col="total.karyotypes")
  {
  test_cols = c('length', 'gene.count', 'ks.gene', 'length.normalized', 'gene.normalized', 'ks.normalized')
  bk_tests = as.data.frame( matrix( nrow=1, ncol=length(test_cols)) )
  
  bk_tests = vector("numeric", length(test_cols))
  names(bk_tests) = test_cols
  
  adj_scores = as.data.frame(matrix(ncol=3, nrow=nrow(scores), dimnames=list(c(1:nrow(scores)), c('chr','band','scores'))))
  adj_scores$chr = scores$chr
  adj_scores$band = scores$band
  
  tryCatch({
    # do the number of breaks correlate to the band length? - yep
    ct = cor.test(scores[,'bp.length'], scores[[col]])
    bk_tests['length'] = round(ct$estimate, 4)
  
    # so what about gene counts? - yep
    ct = cor.test(scores[,'gene.count'], scores[[col]])
    bk_tests['gene.count'] = round(ct$estimate, 4)
  
    ct = ks.test(scores[, col], pnorm, mean(scores[[col]]), sd(scores[[col]]))
    bk_tests['ks.gene'] = round(ct$p.value, 4)
    }, error = function(e) {
      print(e)
    })
    
    length_adj = 0.7
    # Normalize the total for length 
    adj_scores[,'scores'] = scores[[col]]/(scores[,'bp.length']^length_adj)
    breaks_adj = adj_scores[,'scores']

  tryCatch({
    ct = cor.test(breaks_adj, scores[,'bp.length'])
    bk_tests['length.normalized'] = round(ct$estimate, 4)
  
    # length correlation normalized
    ct = cor.test(breaks_adj, scores[,'gene.count']/(scores[,'bp.length']^length_adj))
    bk_tests['gene.normalized'] = round(ct$p.value, 4)
  
    # are the adjusted scores are normal
    ct = ks.test(breaks_adj, pnorm, mean(breaks_adj), sd(breaks_adj))
    bk_tests['ks.normalized'] = round(ct$p.value, 4)
    }, error = function(e) {
      print(e)
  })
      
  return(list("tests" = bk_tests, "scores" = adj_scores))
  }

bp.cor.tests<-function(bpinfo, col)
  {
  chromosomes = c(1:22, 'X')
  
  test_cols = c('length', 'gene.count', 'ks.gene', 'length.normalized', 'gene.normalized', 'ks.normalized')
  bk_tests = as.data.frame( matrix( nrow=length(chromosomes), ncol=length(test_cols)) )
  
  colnames(bk_tests) = test_cols
  rownames(bk_tests) = chromosomes
  
  adj_scores = list()
  
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
      ct = cor.tests(scores, col)
      
      bk_tests[chr, ] = ct$tests
      as = ct$scores
      as[,c('band','scores')]
      adj_scores[[chr]] = as$scores
      names(adj_scores[[chr]]) = as$band
      
      }, error = function(e) {
        print(e)
      })
    }
  return(list("tests" = bk_tests, "scores" = adj_scores))  
  }

plot.cor<-function(tests)
  {
  chrinfo = loadChromosomeInfo("../genomic_info/chromosome_gene_info_2012.txt")  
  par(mfrow=c(2,3))
  for (i in 1:ncol(tests))
    {
    label = "cor estimate"
    if (i == 4) title = "p.value"
    plot(tests[,i], chrinfo[1:23, 'Confirmed.proteins']/chrinfo[1:23,'Base.pairs']^.7,  col='blue', type='p', main=colnames(tests)[i], 
         ylab="Length normalized protein count", xlab=label)
    text(tests[,i], chrinfo[1:23, 'Confirmed.proteins']/chrinfo[1:23,'Base.pairs']^.7 ,labels=rownames(tests), pos=3)
    }
  }

plot.norm<-function(adjusted_scores, title)
  {
  plot(function(x) dnorm(x,mean(adjusted_scores),sd(adjusted_scores)), min(adjusted_scores)-sd(adjusted_scores)*1,max(adjusted_scores)+sd(adjusted_scores)*1,
       ylab="Density",xlab="Instability Score", main=title)
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
    text(xpos[2],ypos[2],labels=names(adjusted_scores[i]),pos=3)
    }
  }

setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")

datadir = "~/Data/sky-cgh/output"
outdir = "~/Analysis/Database/cancer"
setwd(datadir)

total_karyotypes = 100240

bandinfo = read.table("~/Data/sky-cgh/genomic_info/band_genes.txt", header=T, sep="\t")

df = read.table("current/noleuk-breakpoints.txt", sep="\t", header=T)
#df = read.table("current/leuk-breakpoints.txt", sep="\t", header=T)  # sampled leuk bp's to cut down on bias

## To simplify, the bands with < 5 aberrations are all subbands.  We're only dealing with major bands and these
# don't make a difference in the analysis
df = df[-which(df[,'total.karyotypes'] <= 5),]

break_info = merge(df, bandinfo[,c(1,2,5)], by.x=c('chr', 'band'), by.y=c('chr','band'))
break_info$bp.length = break_info$end - break_info$start

# All bands with breaks but NO genes are in the centromere regions...pull them out see if correlations improve
# Note that not all centromeres lack genes apparently

# Not using this for anything, it mostly just confirms the chromosome instability analysis
allcors = bp.cor.tests(break_info, "total.karyotypes")
plot.cor(allcors$tests)

## CLASS 1 ##
cmrows = grep("(11|12)", break_info$band)
centromeres = break_info[ cmrows, ]
centromeres = centromeres[order(centromeres$chr, decreasing=T),]

# centromeres show a correlation with length, not surprised
cor.test(scores[,'bp.length'], centromeres[,'total.karyotypes'])

length_adj = 0.6
# Normalize the total for length  and the correlation drops
cent_adj = centromeres[,'total.karyotypes']/(centromeres[,'bp.length']^length_adj)
cor.test(cent_adj, scores[,'bp.length'])

# not really normal, so adjusted scores maybe best I can do -- er, this is actually the same thing I use in all the others too...
ks.test(cent_adj, pnorm, mean(cent_adj), sd(cent_adj))

## 22q11 is an outlier, and now it's normal.  
drop22 = cent_adj[ cent_adj < max (cent_adj) ]
ks.test(drop22, pnorm, mean(drop22), sd(drop22))

centromeres$bp.prob = round(cent_adj/sum(cent_adj), 5)

break_info = break_info[ -cmrows,] 
# If you compare this to the "all" correlations above the gene correlations get better (for obvious reasons)
nocm = bp.cor.tests(break_info, "total.karyotypes")
plot.cor(nocm$tests)
chr_adj_scores = nocm$scores

# scores per breakpoint within each chromosome
for (i in 1:length(chr_adj_scores))
  {
  chr = names(chr_adj_scores[i])
  probability_list = pnorm(chr_adj_scores[[chr]], mean(chr_adj_scores[[chr]]), sd(chr_adj_scores[[chr]]))
  pf = as.data.frame( round(probability_list/sum(probability_list), 5) )  # prob adjusted to 0-1
  pf$band = row.names(pf)
  
  per_chr = merge(break_info[break_info$chr == chr,], pf, by=c('band'))
  names(per_chr)[length(per_chr)] = 'per.chr.prob'
  if (exists("temp_bk"))  temp_bk = rbind(temp_bk, per_chr) 
  else  temp_bk = per_chr 
  }
break_info = temp_bk
rm(temp_bk)


# scores per breakpoint across entire dataset
arms = cor.tests(break_info, "total.karyotypes")
bp_adj_scores = arms$scores
probability_list = pnorm(bp_adj_scores[,'scores'], mean(bp_adj_scores[,'scores']), sd(bp_adj_scores[,'scores']))
bp_adj_scores[,'bp.prob'] = round(probability_list/sum(probability_list), 5)

# these two are outliers and known leukemia bps
#adjusted_scores = bp_adj_scores$scores
#names(adjusted_scores) = paste(bp_adj_scores$chr, bp_adj_scores$band, sep="")
#adjusted_scores = adjusted_scores[ -which(names(adjusted_scores) == "14q32" | names(adjusted_scores) == "9q34")]
#plot.norm(adjusted_scores, "Breakpoint instability by gene count")

break_info = merge(break_info, bp_adj_scores[,c('chr','band','bp.prob')], by=c('chr','band'))


# So lets see the range of breakpoints 
break_info = break_info[ order(break_info[,'total.karyotypes'], decreasing=T),]

sd(break_info[,'total.karyotypes'])
summary(break_info[,'total.karyotypes'])

## CLASS 2 ##
#topBP = break_info[ break_info[,'total.karyotypes'] >= mean(break_info[,'total.karyotypes']), ]

## CLASS 3 ##
#remainder = break_info[ break_info[,'total.karyotypes'] < mean(break_info[,'total.karyotypes']), ]

cols = c('chr','band','start','end','bp.prob')
write.table(centromeres[,cols], row.name=F, quote=F, sep="\t", file=paste(outdir, "centromeres-probs.txt", sep="/"))

#cols = c('chr','band','start','end','per.chr.prob')
#write.table(topBP[,cols], row.name=F, quote=F, sep="\t", file=paste(outdir, "class2-topbp.txt", sep="/"))
#write.table(remainder[,cols], row.name=F, quote=F, sep="\t", file=paste(outdir, "class3-remainder.txt", sep="/"))


write.table(break_info[,c('chr','band','bp.prob', 'per.chr.prob')], quote=F, row.name=F, sep="\t", file=paste(outdir, "all-bp-prob.txt", sep="/"))


