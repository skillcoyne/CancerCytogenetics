bp = read.table("breakpoints.txt", sep="\t", comment="#", header=T)

## Lets ignore subbands (11.1) and just group them by the major band designation (11)
bp$Breakpoint = sub("\\.[0-9]+", "", bp$Breakpoint)


bpfreq = table(bp$Breakpoint)
sorted = sort(bpfreq)

cancers = table(bp$Cancer)
cancers = rownames(cancers)

cancer_counts=vector("numeric",length(cancers))
names(cancer_counts) = cancers
for(i in 1:length(cancers))
  {
  c = cancers[i]
  cancer_counts[i] = sum(bp$Cancer==c)
  }
cancer_counts = log(cancer_counts)

# Not super useful
plot(hclust(dist(cancer_counts)))


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

# drop low frequence breakpoints first (just for computational purposes)
bp_freq = log(table(bp$Breakpoint))
bph = bp[bp$Breakpoint %in% names(bp_freq),]  
nrow(bph)



c_vs_bp = table(bp$Cancer, bp$Breakpoint)
c_vs_bp = t(c_vs_bp) # cancer names as the columns

for (cancer in colnames(c_vs_bp))
  {
  print(c_vs_bp[,cancer])
  print (cancer)
  
  cancer_table = c_vs_bp[,cancer]
  #cancer_table = cancer_table[cancer_table > 0]
  print(cancer_table)
  break
  }


