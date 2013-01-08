source("R/lib/load_files.R")
source("R/lib/wd.R")

setDataDirectory(date = NA)

fg = loadFragments("fragments.txt")

freq=table(fg$Chr)

# filter these out
for (i in 1:22)
	{
	print(i)
	print(nrow(fg$Chr==i))
	# p arm
	nt=fg[fg$Chr==i,]
	nt=nt[nt$Start!=ter[i, "p"] & nt$End!="pter",] 
	nt=nt[fg$Start!="pter" & fg$End!=ter[i, "p"],] 
#	fg=fg[which(fg$Chr==i & fg$Start!=ter[i, "p"] & fg$End!="pter"),] 
#	fg=fg[which(fg$Chr==i & fg$Start!="pter" & fg$End!=ter[i, "p"]),] 
	fg$Chr
	print(nrow(fg$Chr==i))
	}
