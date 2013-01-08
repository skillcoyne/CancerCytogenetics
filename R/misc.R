source("R/lib/load_files.R")
source("R/lib/wd.R")

setDataDirectory(date = NA)

bp = loadBreakpoints("breakpoints.txt")
chrinfo = loadChromosomeInfo()  

bpfreq = table(bp$Breakpoint)
sorted = sort(bpfreq)

ymax = max(log(bpfreq))

## Plot by chromosome  
dev.new()
par(mfrow=c(4,6))
for (i in c(1:22, c("X", "Y")))
  {
  nd = subset(bp, bp$Chr == i)
  freq = table(as.vector(nd$Breakpoint))
#pn = qnorm(log(freq))
  title = paste("Chr", i, sep=" ")
#  plot(pn, ylim=c(0, 1), main=title, xlab="breakpoints", ylab="log(freq)")
  plot(log(freq), ylim=c(0, ymax), main=title, xlab="breakpoints", ylab="log(freq)")
  }


## CENTROMERES ##

centromeres = vector(mode="character")
#for (i in c(1:22, c("X", "Y")))
for (i in 1:22)
  {
#  for (b in c(10:11))  ## p11/q11 are the centromeric bands, 10/12 are adjacent 
#    {
	b = 11
    centromeres = c(centromeres, paste(i, "q", b, sep=""))
    centromeres = c(centromeres, paste(i, "p", b, sep=""))
#    }
  }

# Centromere
bpcent=vector("numeric",length(centromeres))
names(bpcent)=centromeres
for (i in centromeres)  
	{
	bpcent[i] = sum(bp$Breakpoint == i)
	}

`%nin%`=Negate(`%in%`) # cute way to create my own 'not' operator
# Arm
bparm = bp[which(bp$Breakpoint %nin% centromeres),]  
sum(bpcent)/sum(bpfreq)
nrow(bparm)/nrow(bp)

total_percentage_centromeric=round(length(centromeres)/nrow(bpfreq), digits=2)*100

dev.new()
sorted_cent_freq=sort(bpcent)
sorted_cent_freq=log(sorted_cent_freq)
plot(sorted_cent_freq[1:length(sorted_cent_freq)], main="Centromere log frequency", ylab="log(centromeric frequency)", xlab="Chr centromeres", pch=19)
text(sorted_cent_freq[1:length(sorted_cent_freq)],labels=names(sorted_cent_freq), pos=3)

