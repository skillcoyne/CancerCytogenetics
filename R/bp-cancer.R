source("R/lib/load_files.R")
source("R/lib/wd.R")

setDataDirectory()
bp = loadBreakpoints("breakpoints.txt")

cancers = table(bp$Cancer)
cancers = rownames(cancers)

cancer_counts=vector("numeric",length(cancers))
names(cancer_counts) = cancers
for(i in 1:length(cancers))
  {
  c = cancers[i]
  cancer_counts[i] = sum(bp$Cancer==c)
  }
#cancer_counts = log(cancer_counts)


# This probably is not useful since I have classified the cancers in the dataset.
# I would need to drop all the of the cambridge samples as well since they are cell line
# only with no clinical information
sorted_cc = sort(cancer_counts)
plot(sorted_cc)
axis(1, at=c(1:length(cancer_counts)), labels=names(cancer_counts))

# not sure this tells me anything actually, particularly since I didn't log how many samples 
# of each cancer are represented in the data...
cor.test(cancer_counts,c(1:length(cancers)))

# Perhaps more interested to look at which events correlate to which cancers?

# drop low frequence breakpoints first as there are a fair number of low frequency 
bp_freq = log(table(bp$Breakpoint))
# take all those above one standard dev, this drops about 200 breakpoints
bp_freq = bp_freq[ which(bp_freq > mean(bp_freq)) ] 
nrow(bp)
bp = bp[bp$Breakpoint %in% names(bp_freq),]  
nrow(bp) 

# > 1std drops about 600 observations from 151k
# > mean drops about 1k

c_vs_bp = table(bp$Cancer, bp$Breakpoint)
c_vs_bp = t(c_vs_bp) # cancer names as the columns

for (cancer in colnames(c_vs_bp))
  {
#  cancer_table = c_vs_bp[c_vs_bp[,cancer] > 0, cancer]
#  plot(cancer_table)
#  print(cancer_table)
  print (cancer)
  print(sum(c_vs_bp[,cancer])) 
  break
  }



