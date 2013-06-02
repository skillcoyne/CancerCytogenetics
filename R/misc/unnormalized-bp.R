source("R/lib/load_files.R")
source("R/lib/wd.R")

setDataDirectory(date = NA)
`%nin%`=Negate(`%in%`) 

bp = loadBreakpoints("breakpoints.txt")
chrinfo = loadChromosomeInfo()  


bpfreq = table(bp$Breakpoint)
sorted = sort(bpfreq)

# ignoring X and Y
chromosome_counts=vector("numeric",22)
names(chromosome_counts)=1:22
for(i in 1:22)
  {
  sum(bp$Chr==i)->chromosome_counts[i]
  }


#shock/horror number of events correlates with length of chromosome
cor.test(chromosome_counts,chrinfo[1:22, "Base.pairs"])
plot(chromosome_counts, chrinfo[1:22,"Base.pairs"], type="n")
text(chromosome_counts, chrinfo[1:22,"Base.pairs"],labels=row.names(chrinfo)[1:22])

#and as length correlates with number of proteins this also correlates- but the score is higher
cor.test(chromosome_counts,chrinfo[1:22,"Confirmed.proteins"])
plot(chromosome_counts, chrinfo[1:22,"Confirmed.proteins"], type="n")
text(chromosome_counts, chrinfo[1:22,"Confirmed.proteins"],labels=row.names(chrinfo)[1:22])

#again length correlates with number of psudeogenes- although this is a very high correlation and so is interesting
cor.test(chromosome_counts,chrinfo[1:22,"Pseudogenes"])
plot(chromosome_counts, chrinfo[1:22,"Pseudogenes"], type="n")
text(chromosome_counts, chrinfo[1:22,"Pseudogenes"],labels=row.names(chrinfo)[1:22])

# just to show that the correlation is still to length
cor.test(chromosome_counts,chrinfo[1:22,"Total.Prot.RNA"])
plot(chromosome_counts, chrinfo[1:22,"Total.Prot.RNA"], type="n")
text(chromosome_counts, chrinfo[1:22,"Total.Prot.RNA"],labels=row.names(chrinfo)[1:22])


#problem is that this strong correlation is largely being driven by length of chromosome
cor.test(chromosome_counts,chrinfo[1:22,"Variations"])

