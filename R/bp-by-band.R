## This pretty much just supports the initial by chromosome analysis.  
# Which means I need to score the bands by the number of genes perhaps?

## Three classes of breakpoints for the purposes of insilico generation
# 1. centromeres
# 2. Top breakpoints 
# 3. The rest based on correlation with gene count.
## Note that chromosomes 19,20,22 don't have enough observations due to the fact that most/all of their breakpoints occur within centromeres

rm(list=ls())
setwd("~/workspace/CancerCytogenetics/R")
source("lib/load_files.R")
source("lib/wd.R")


cor.adj<-function(scores, col="total.karyotypes")
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

bp.cor.tests<-function(bpinfo, col, chromosomes = c(1:22, 'X'))
  {
  #chromosomes = c(1:22, 'X')
  
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
      ct = cor.adj(scores, col)
      
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
  chrinfo = read.table("../genomic_info/chromosome_gene_info_2012.txt", header=T, sep="\t")  
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

adjust.to.one<-function(p, r=5)
  {
  adjusted = round(p/sum(p), r) 
  
  if (sum(adjusted) > 1)
    adjusted[ which(adjusted == min(adjusted)) ] = adjusted[ which(adjusted == min(adjusted)) ] - (sum(adjusted) - 1)
  
  if (sum(adjusted) < 1)
    adjusted[ which(adjusted == min(adjusted)) ] = adjusted[ which(adjusted == min(adjusted)) ] + (1 - sum(adjusted))
  
  return(adjusted)
  }



datadir = "~/Data/sky-cgh/output"
outdir = "~/Analysis/Database/cancer"
setwd(datadir)

total_karyotypes = 100240

bandinfo = read.table("~/Data/sky-cgh/genomic_info/band_genes.txt", header=T, sep="\t")
bandinfo$bp.length = bandinfo$end - bandinfo$start
bandinfo$chr = as.character(bandinfo$chr)
bandinfo$band = as.character(bandinfo$band)


correlations = matrix(nrow=3,ncol=3,dimnames=list(c('cor.all','cor.nophil','adj.genes'),c('all','centromeres','arms')))

leuk = T

# Load files
df = load.breakpoints("current/noleuk-breakpoints.txt")
if (leuk)
  df = load.breakpoints(c("current/noleuk-breakpoints.txt", "current/leuk-breakpoints.txt"))

all = merge(df,bandinfo,by=c('chr','band'))
all$start=NULL
all$end=NULL

# All bands with breaks but NO genes are in the centromere regions...pull them out see if correlations improve
# Note that not all centromeres lack genes apparently

# Not using this for anything, it mostly just confirms the chromosome instability analysis
allcors = bp.cor.tests(all, "total.karyotypes")
plot.cor(allcors$tests)
dev.off()

ct = cor.test(all$total.karyotypes, all$bp.length)
correlations['cor.all','all'] = ct$estimate


plot(all$total.karyotypes, all$gene.count,type='n', main="Per Band", ylab="Gene", xlab="# Breakpoints")
text(all$total.karyotypes, all$gene.count,labels=paste(all$chr,all$band,sep=""))

noPhil = which(paste(all$chr,all$band,sep="") != '22q11' & paste(all$chr,all$band,sep="") != "9q34")

ct = cor.test(all$total.karyotypes[noPhil], all$bp.length[noPhil])
correlations['cor.nophil','all'] = ct$estimate

all_adj = all$total.karyotypes/(all$bp.length^0.7)
cor.test(all_adj, all$bp.length)
ks.test(all_adj, pnorm, mean(all_adj), sd(all_adj))
plot.norm(all_adj, "all")
ct = cor.test(all_adj, all$gene.count)
correlations['adj.genes','all'] = ct$estimate



## CLASS 1 ##
cmrows = grep("(p|q)11", all$band)
centromeres = all[ cmrows, ]
centromeres = centromeres[order(centromeres$chr, decreasing=T),]

dev.off()
# centromeres show a correlation with length, not surprised -- 22q11 skews this a little but not significantly
ct = cor.test(centromeres[,'total.karyotypes'],centromeres[,'bp.length'])
correlations['cor.all','centromeres'] = ct$estimate

par(mfrow=c(2,1))
plot(centromeres[,'total.karyotypes'],centromeres[,'gene.count'], type='n', main="Per Centromere", ylab="Gene count", xlab="# Breakpoints")
text(centromeres[,'total.karyotypes'],centromeres[,'gene.count'], labels=unlist(apply(centromeres[,c('chr','band')], 1, paste, collapse="")))

plot(centromeres[,'total.karyotypes'],centromeres[,'bp.length'], type='n', main="Per Centromere", ylab="Length", xlab="# Breakpoints")
text(centromeres[,'total.karyotypes'],centromeres[,'bp.length'], labels=unlist(apply(centromeres[,c('chr','band')], 1, paste, collapse="")))

noPhil = which(centromeres$name != '22q11')
ct = cor.test(centromeres[noPhil,'total.karyotypes'],centromeres[noPhil,'bp.length'])
correlations['cor.nophil','centromeres'] = ct$estimate

length_adj = 0.4
# Normalize for length  and the correlation drops

cent_adj = centromeres[,'total.karyotypes']/(centromeres[,'bp.length']^length_adj)
cor.test(cent_adj, centromeres[,'bp.length'])
cor.test(cent_adj[noPhil], centromeres[noPhil,'bp.length'])

ct = cor.test(cent_adj, centromeres[,'gene.count'])
correlations['adj.genes','centromeres'] = ct$estimate



# not enough observations per chromosome to test normality per chromosome
# not really normal, so adjusted scores may be best I can do 
ks.test(cent_adj, pnorm, mean(cent_adj), sd(cent_adj))

## 22q11 is an outlier, and it's better but still not normal  
drop22 = cent_adj[ cent_adj < max (cent_adj) ]
ks.test(drop22, pnorm, mean(drop22), sd(drop22))

# get probabilities from pnorm then adjust to 1 
cent_adj = unlist(bp.cor.tests(centromeres,'total.karyotypes', c(1:22,'X','Y'))$scores)
names(cent_adj) =  sub("\\.","", names(cent_adj))
centromeres[match(names(cent_adj), centromeres$name), 'bp.prob'] =  adjust.to.one(cent_adj, 5)

## Arms
arms = all[-cmrows,]
arms = arms[order(arms$chr,decreasing=T),]
ct = cor.test(arms$total.karyotypes,arms$bp.length)
correlations['cor.all','arms'] = ct$estimate


plot(arms$total.karyotypes,arms$bp.length, type='n')
text(arms$total.karyotypes,arms$bp.length, labels=unlist(apply(arms[,c('chr','band')], 1, paste, collapse="")))

noPhil = which(arms$name != '9q34')
ct = cor.test(arms$total.karyotypes[noPhil],arms$bp.length[noPhil])
correlations['cor.nophil','arms'] = ct$estimate

length_adj = 0.6
arm_adj = arms[,'total.karyotypes']/(arms[,'bp.length']^length_adj)
cor.test(arm_adj, arms[,'bp.length'])

ct = cor.test(arm_adj, arms[,'gene.count'])
correlations['adj.genes','arms'] = ct$estimate

dev.off()
colors=c('firebrick3','blue','green3')
bp = barplot((correlations), beside=T, col=colors, border=NA, ylim=c(0,0.45), main="Karyotype Count Correlations", ylab="Pearson's correlation", cex=2, cex.axis=1.5)
legend('topleft', legend=c('All(v)Length','-22q11,9q34(v)Length', '(v)Length Adj. Gene'), fill=colors, border=NA)
for (i in 1:nrow(bp))
  text(bp[i,], round(correlations[i,]+.01,2), labels=round(correlations[i,], 2), cex=1.5)

## So if I look at the chromosomal instability as a whole without centromeres?
chrs = c(1:22,'X','Y')
arm_chr = as.data.frame(matrix(nrow=length(chrs), ncol=3, dimnames=list(chrs,c('total.karyotypes', 'gene.count','length'))))
arm_chr[,'total.karyotypes'] = unlist(lapply(chrs, function(x)  sum(arms[arms$chr == x, 'total.karyotypes'])))
arm_chr[,'gene.count'] = unlist(lapply(chrs, function(x)  sum(arms[arms$chr == x, 'gene.count'])))
arm_chr[,'length'] = unlist(lapply(chrs, function(x)  sum(arms[arms$chr == x, 'bp.length'])))


cor.test(arm_chr$total.karyotypes,arm_chr$length)
cor.test(arm_chr$total.karyotypes,arm_chr$gene.count)


# If you compare this to the "all" correlations above the gene correlations get better (for obvious reasons)
arm.cor =  bp.cor.tests(arms, "total.karyotypes")
plot.cor(arm.cor$tests)
arm_adj_scores = arm.cor$scores

chrs=c(1:22)
length_adj=0.7
cor.test(arm_chr$total.karyotypes[chrs], arm_chr$length[chrs])
arm_adj_scores = arm_chr$total.karyotypes[chrs]/(arm_chr$length[chrs]^length_adj)
cor.test(arm_adj_scores, arm_chr$length[chrs])
cor.test(arm_adj_scores, arm_chr$gene.count[chrs])

instability_score = pnorm(arm_adj_scores,mean(arm_adj_scores),sd(arm_adj_scores))
names(instability_score) = rownames(arm_chr)[chrs]

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





# scores per breakpoint within each chromosome arms
for (i in 1:length(arm_adj_scores))
  {
  chr = names(arm_adj_scores[i])
  probability_list = pnorm(arm_adj_scores[[chr]], mean(arm_adj_scores[[chr]]), sd(arm_adj_scores[[chr]]))
  pf = as.data.frame( round(probability_list/sum(probability_list), 5) )  # prob adjusted to 0-1
  pf$band = row.names(pf)
  
  per_chr = merge(all[all$chr == chr,], pf, by=c('band'))
  names(per_chr)[length(per_chr)] = 'per.chr.prob'
  if (exists("temp_bk"))  temp_bk = rbind(temp_bk, per_chr) else  temp_bk = per_chr 
  }
arm_break_info = temp_bk
rm(temp_bk)


# scores per arm breakpoint across entire dataset
arms = cor.adj(arm_break_info, "total.karyotypes")
arm_bp_adj_scores = arms$scores
probability_list = pnorm(arm_bp_adj_scores[,'scores'], mean(arm_bp_adj_scores[,'scores']), sd(arm_bp_adj_scores[,'scores']))
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

cols = c('chr','band','start','end','bp.prob')
write.table(centromeres[,cols], row.name=F, quote=F, sep="\t", file=paste(outdir, "centromeres-probs.txt", sep="/"))

## Note...the scores don't change if the centromeres are't removed from the list so for simplicity in the selection process centromeres will be put into this list
# What I really need for the database though is breakpoint probabilities per chromosomes
all_scores = bp.cor.tests(all, "total.karyotypes")$scores
# scores per breakpoint within each chromosome
for (i in 1:length(all_scores))
  {
  chr = names(all_scores[i])
  probability_list = pnorm(all_scores[[chr]], mean(all_scores[[chr]]), sd(all_scores[[chr]]))
  pf = as.data.frame( adjust.to.one(probability_list/sum(probability_list), 5) )  # prob adjusted to 0-1
  pf$band = row.names(pf)
  
  per_chr = merge(break_info[break_info$chr == chr,], pf, by=c('band'))
  names(per_chr)[length(per_chr)] = 'per.chr.prob'
  if (exists("temp_bk"))  temp_bk = rbind(temp_bk, per_chr) 
  else  temp_bk = per_chr 
  }
all = temp_bk
rm(temp_bk)




#write.table(all[,c('chr','band','bp.prob', 'per.chr.prob')], quote=F, row.name=F, sep="\t", file=paste(outdir, "all-bp-prob.txt", sep="/"))


chr_ins = read.table(paste(outdir, "chr_instability_prob.txt", sep="/"), header=F)
colnames(chr_ins) = c("chr", "prob")

chr_ins = merge(chr_ins, carm_probs, by.x="chr", by.y="row.names")


#write.table(chr_ins, quote=F, row.name=F, sep="\t", file=paste(outdir, "chr_instability_prob.txt", sep="/"))

