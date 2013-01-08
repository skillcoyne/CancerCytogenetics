source("R/lib/load_files.R")
source("R/lib/wd.R")

setDataDirectory(date = NA)

bp = loadBreakpoints("breakpoints.txt")
chrinfo = loadChromosomeInfo()  

# frequency tables
bpfreq = table(bp$Breakpoint)
lowfreq = bpfreq[bpfreq < mean(bpfreq)]
highfreq = bpfreq[bpfreq > mean(bpfreq)]

#write.table(sort(highfreq), file="/Users/sarah.killcoyne/Data/sky-cgh/analysis/high-freq-bp.txt", quote=F, col.names=F, sep="\t")

dev.new()
par(mfrow=c(3,1))
plot(bpfreq, main="Breakpoint log frequency", ylab="BP frequency", xlab="Breakpoints")
# TODO would be niceish to get these labeled
plot(lowfreq, main=paste("Frequencies below the mean (", mean(bpfreq), ")", sep=""), type="h")
plot(highfreq, main=paste("Frequencies above the mean (", mean(bpfreq), ")", sep=""), type="h")

## Plot by chromosome  
ymax = max(log(bpfreq))
dev.new()
par(mfrow=c(4,6))
for (i in c(1:22, c("X", "Y")))
  {
  nd = subset(bp, bp$Chr == i)
  freq = table(as.vector(nd$Breakpoint))
  title = paste("Chr", i, sep=" ")
  plot(log(freq), ylim=c(0, ymax), main=title, xlab="breakpoints", ylab="log(freq)")
  }


## CENTROMERES ## ignoring X,Y
centromeres = vector(mode="character")
for (i in 1:22)
	{
	b = 11 ## p11/q11 are the centromeric bands, 10/12 are adjacent 
    centromeres = c(centromeres, paste(i, "q", b, sep=""))
    centromeres = c(centromeres, paste(i, "p", b, sep=""))
	}

# Centromere
bpcent=vector("numeric",length(centromeres))
names(bpcent)=centromeres
for (i in centromeres)  
	{
	bpcent[i] = sum(bp$Breakpoint == i)
	}
#write.table(sort(bpcent), file=""/Users/sarah.killcoyne/Data/sky-cgh/analysis/centromere-frequencies.txt", quote=F, sep="\t", col.names=F)

`%nin%`=Negate(`%in%`) # cute way to create my own 'not' operator
# Arm
bparm = bp[which(bp$Breakpoint %nin% centromeres),]  
sum(bpcent)/sum(bpfreq)
nrow(bparm)/nrow(bp)

dev.new()
sorted_cent_freq=sort(bpcent)
sorted_cent_freq=log(sorted_cent_freq)
plot(sorted_cent_freq[1:length(sorted_cent_freq)], main="Centromere log frequency", ylab="log(centromeric frequency)", xlab="Chr centromeres", pch=19, xaxt="n")
text(sorted_cent_freq[1:length(sorted_cent_freq)],labels=names(sorted_cent_freq), pos=3)

