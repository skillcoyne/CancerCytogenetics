source("R/lib/load_files.R")
source("R/lib/wd.R")

setDataDirectory(date = '09012013')

# Load files
bp = loadBreakpoints("breakpoints.txt")

sample_leukemia_bps = FALSE
if (sample_leukemia_bps)
  {
  leuks = c('Acute myeloid leukemia', 'Acute lymphoblastic leukemia', "Non-hodgkin's lymphoma", 'Chronic myelogenous leukemia', 'Chronic lymphocytic leukemia')
  bp = sampleCancers(bp, "Cancer", leuks)
  }

# File lists only protein coding genes. Ignore subbands (11.1) and just group them by the major band designation (11)
gene_loc = read.table("../../genomic_info/Hs_gene_band.txt", header=T, sep="\t")
gene_loc$Band = clearSubbands(gene_loc$Band)

# Reorder
bp = bp[ order(bp$Chr, bp$Breakpoint), ]
gene_loc = gene_loc[ order(gene_loc$Chromosome, gene_loc$Band), ]
gene_loc = gene_loc[ gene_loc$Band %in% bp$Breakpoint,] # select genes where breakpoints occur, drops only about 115 entries

breaks = table(bp$Breakpoint)

genes_per_band = vector("numeric", length(breaks))
names(genes_per_band) = names(breaks)
for (band in names(genes_per_band))
  {
  genes_per_band[band] = nrow(gene_loc[gene_loc$Band == band,])  
  }


## Clustering the instability scores had given us 3 major clusters need to know what they share
# !@! This isn't actually useful.  I could just go look up all the genes on each chromosome.  What could be more useful is per band  !@!
c1 = c(22,9,17,1,14,11,21)
c2 = c(4,10,2,13,5,20,15,16)
c3 = c(6,18,8,12,7,3,19)
clusters = list('1' = c1, '2' = c2, '3' = c3)

file_names = c("instability_cluster_1.txt", "instability_cluster_2.txt", "instability_cluster_3.txt")
for (i in 1:length(clusters))
  {
  cls = clusters[[i]]
  curr_gene = gene_loc[ gene_loc$Chromosome %in% cls, ]
  cancers = bp[ bp$Chr %in% cls,]
  curr_gene = curr_gene[ order(curr_gene$Chromosome, curr_gene$Band), ]
  #write.table(curr_gene, file = file_names[i], row.names=F, quote=F, sep="\t")
  }







