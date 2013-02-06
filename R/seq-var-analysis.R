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
vars = read.table("all_variants.txt", header=T, sep="\t") 
## Not looking at SNPs may help
vars = vars[vars$VarType %nin% c('SNP'),]

vars$Length = vars$EndPosition - vars$StartPosition



ins = vars[vars$VarType == 'INS',]
del = vars[vars$VarType == 'DEL',]

table(del$Length)
table(ins$Length)


